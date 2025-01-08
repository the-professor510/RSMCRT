module kernels
!! Contains the main program and scattering loop. Calls all other routine to setup, run and break down the simulation.

    implicit none
    
    private
    public :: weight_scatter, pathlength_scatter, pathlength_scatter2, test_kernel

contains
!###############################################################################
!                   KERNELS
    subroutine weight_scatter(input_file)

        !Shared data
        use iarray
        use constants, only : wp, CHANCE, THRESHOLD

        !subroutines
        use detectors,     only : dect_array
        use detector_mod,  only : hit_t
        use historyStack,  only : history_stack_t
        use inttau2,       only : tauint2, update_voxels
        use photonMod,     only : photon
        use piecewiseMod
        use random,        only : ran2, init_rng
        use sdfs,          only : sdf
        use sim_state_mod, only : state
        use utils,         only : pbar
        use vec4_class,    only : vec4
        use vector_class,  only : vector
        
        !external deps
        use tev_mod, only : tevipc
        use tomlf,   only : toml_table
#ifdef _OPENMP
        use omp_lib
#endif
        character(len=*), intent(in)   :: input_file
        
        integer                        :: numproc, id, j, i
        type(history_stack_t)          :: history
        type(pbar)                     :: bar
        type(photon)                   :: packet
        type(toml_table)               :: dict
        real(kind=wp),     allocatable :: distances(:), image(:,:,:)
        type(hit_t)                    :: hpoint
        type(vector)                   :: dir
        type(dect_array),  allocatable :: dects(:)
        type(sdf),         allocatable :: array(:)
        real(kind=wp)                  :: nscatt, start, weight_absorb
        type(tevipc)                   :: tev
        integer                        :: celli, cellj, cellk
        type(spectrum_t)               :: spectrum

        call setup(input_file, tev, dects, array, packet, spectrum, dict, distances, image, nscatt, start)

#ifdef _OPENMP
        !is state%seed private, i dont think so...
        !$omp parallel default(none) shared(dict, array, numproc, start, state, bar, jmean, emission, tev, dects, spectrum)&
        !$omp& private(id, distances, image, dir, hpoint, history, weight_absorb, cellk, cellj, celli) &
        !$omp& reduction(+:nscatt) firstprivate(packet)
        numproc = omp_get_num_threads()
        id = omp_get_thread_num()
        if(numproc > state%nphotons .and. id == 0)print*,"Warning, simulation may be underministic due to low photon count!"
        if(state%trackHistory)history = history_stack_t(state%historyFilename, id)
#elif MPI
    !nothing
#else
        numproc = 1
        id = 0
        if(state%trackHistory)history = history_stack_t(state%historyFilename, id)
#endif
        if(id == 0)print("(a,I3.1,a)"),'Photons now running on', numproc,' cores.'

        ! set seed for rnd generator. id to change seed for each process
        call init_rng(state%iseed, fwd=.true.)

        bar = pbar(state%nphotons/ 10)
        !$OMP BARRIER
        !$OMP do
        !loop over photons 
        do j = 1, state%nphotons
            if(mod(j, 10) == 0)call bar%progress()

            ! Release photon from point source
            call packet%emit(spectrum, dict)
            packet%step = 0
            packet%id = id
            distances = 0._wp
            do i = 1, size(distances)
                distances(i) = array(i)%evaluate(packet%pos)
                !if(distances(i) > 0._wp)distances(i)=-999.0_wp
            end do
            packet%layer=maxloc(distances,dim=1, mask=(distances<=0._wp))
            if(state%trackHistory)call history%push(vec4(packet%pos, packet%step))
            ! Find scattering location
            call tauint2(state%grid, packet, array, dects, history)

            do while(.not. packet%tflag)
                if(state%trackHistory)call history%push(vec4(packet%pos, packet%step))
                
                weight_absorb = packet%weight * (1._wp - array(packet%layer)%getAlbedo())

                packet%weight = packet%weight - weight_absorb
                call update_voxels(state%grid, &
                packet%pos + vector(state%grid%xmax, state%grid%ymax, state%grid%zmax), celli, cellj, cellk)

                if(celli < 1)then
                    packet%tflag = .true.
                    exit
                end if
                if(cellj < 1)then
                    packet%tflag = .true.
                    exit
                end if
                if(cellk < 1)then
                    packet%tflag = .true.
                    exit
                end if
                !$omp atomic
                jmean(celli,cellj,cellk) = jmean(celli,cellj,cellk) + weight_absorb
                call packet%scatter(array(packet%layer)%gethgg(), array(packet%layer)%getg2(), dects)
                if(packet%weight < THRESHOLD)then
                    if(ran2() < CHANCE)then
                        packet%weight = packet%weight / CHANCE
                    else
                        packet%tflag = .true.
                        exit
                    end if
                end if

                ! !Find next scattering location
                call tauint2(state%grid, packet, array, dects, history)
            end do

            dir = vector(packet%nxp, packet%nyp, packet%nzp)
            hpoint = hit_t(packet%pos, dir, packet%weight, packet%layer, packet%weight)
            do i = 1, size(dects)
                call dects(i)%p%record_hit(hpoint, history)
            end do

            if(id == 0 .and. mod(j,1000) == 0)then
                if(state%tev)then
!$omp critical
                    image = reshape(jmean(:,100:100,:), [state%grid%nxg,state%grid%nzg,1])
                    call tev%update_image(state%experiment, real(image(:,:,1:1)), ["I"], 0, 0, .false., .false.)

                    image = reshape(jmean(100:100,:,:), [state%grid%nyg,state%grid%nzg,1])
                    call tev%update_image(state%experiment, real(image(:,:,1:1)), ["J"], 0, 0, .false., .false.)

                    image = reshape(jmean(:,:,100:100), [state%grid%nxg,state%grid%nyg,1])
                    call tev%update_image(state%experiment, real(image(:,:,1:1)), ["K"], 0, 0, .false., .false.)
!$omp end critical
                end if
            end if
        end do

#ifdef _OPENMP
!$OMP end  do
!$OMP end parallel
#endif
    call finalise(dict, dects, nscatt, start, history)
    end subroutine weight_scatter


    !Full weight reduction
    subroutine pathlength_scatter(input_file)

        !Shared data
        use iarray
        use constants, only : wp

        !subroutines
        use detector_mod,  only : hit_t
        use detectors,     only : dect_array
        use historyStack,  only : history_stack_t
        use inttau2,       only : tauint2
        use photonMod,     only : photon
        use piecewiseMod
        use random,        only : ran2, init_rng, seq
        use sdfs,          only : sdf
        use sim_state_mod, only : state
        use utils,         only : pbar
        use vec4_class,    only : vec4
        use vector_class,  only : vector
        use writer_mod,    only : checkpoint

        !external deps
        use tev_mod, only : tevipc
        use tomlf,   only : toml_table
#ifdef _OPENMP
        use omp_lib
#endif
        character(len=*), intent(in) :: input_file
        
        integer                       :: numproc, id, j, i
        type(history_stack_t)         :: history
        type(pbar)                    :: bar
        type(photon)                  :: packet
        type(toml_table)              :: dict
        real(kind=wp),    allocatable :: distances(:), image(:,:,:)
        type(hit_t)                   :: hpoint
        type(vector)                  :: dir
        type(dect_array), allocatable :: dects(:)
        type(sdf),        allocatable :: array(:)
        real(kind=wp)                 :: ran, nscatt, start, weight_absorb
        type(tevipc)                  :: tev
        type(seq)                     :: seqs(2)
        type(spectrum_t)              :: spectrum
        real :: tic, toc

        integer :: nphotons_run,pos
        character(len=128) :: line
        character(len=:), allocatable :: checkpt_input_file

        if(state%loadckpt)then
            call setup(input_file, tev, dects, array, packet, spectrum, dict, distances, image, nscatt, start, .false.)
            open(newunit=j,file=state%ckptfile, access="stream", form="formatted")
            read(j,"(a)")line
            pos = scan(line, "=")
            checkpt_input_file = trim(line(pos+1:))

            read(j,"(a)")line
            pos = scan(line, "=")
            read(line(pos+1:),*) nphotons_run

            inquire(j,pos=pos)
            close(j)

            open(newunit=j,file=state%ckptfile, access="stream", form="unformatted")
            read(j,pos=pos)jmean
            close(j)

            call setup(checkpt_input_file, tev, dects, array, packet, spectrum, dict, distances, image, nscatt, start, .true.)
            state%iseed=state%iseed*101
            state%nphotons = state%nphotons - nphotons_run
        else
            call setup(input_file, tev, dects, array, packet, spectrum, dict, distances, image, nscatt, start, .true.)
        end if

#ifdef _OPENMP
        tic=omp_get_wtime()
!$omp parallel default(none)& 
        !$omp& shared(dict, array, numproc, start, bar, jmean, emission, absorb, input_file, phasor, tev, dects, spectrum)& 
        !$omp& private(ran, id, distances, image, dir, hpoint, history, seqs)& 
        !$omp& reduction(+:nscatt) firstprivate(state, packet)
        numproc = omp_get_num_threads()
        id = omp_get_thread_num()
        if(numproc > state%nphotons .and. id == 0)print*,"Warning, simulation may be underministic due to low photon count!"
        if(state%trackHistory)history = history_stack_t(state%historyFilename, id)
#elif MPI
    !nothing
#else
        call cpu_time(tic)
        numproc = 1
        id = 0
        if(state%trackHistory)history = history_stack_t(state%historyFilename, id)
#endif
        if(id == 0)print("(a,I3.1,a)"),'Photons now running on', numproc,' cores.'
        state%iseed = state%iseed + id
        ! set seed for rnd generator. id to change seed for each process
        call init_rng(state%iseed, fwd=.true.)
        seqs = [seq((id+1)*(state%nphotons/numproc), 2),&
                seq((id+1)*(state%nphotons/numproc), 3)]

        bar = pbar(state%nphotons/ 10)

        !$OMP BARRIER
        !$OMP do
        !loop over photons
        do j = 1, state%nphotons
            if(mod(j, 10) == 0)call bar%progress()
            if(mod(j, state%ckptfreq) == 0 .and. id==0)call checkpoint(input_file, state%ckptfile, j, .true.)

            ! Release photon from point source
            call packet%emit(spectrum, dict, seqs)

            do while (packet%xcell < 1 .or. packet%xcell > state%grid%nxg .or. &
                        packet%ycell < 1 .or. packet%ycell > state%grid%nyg .or. &
                        packet%zcell < 1 .or. packet%zcell > state%grid%nzg)
                call packet%emit(spectrum, dict, seqs)
            end do

            if(state%render_source)call recordEmission(packet)
            packet%step = 0
            packet%id = id
            distances = 0._wp
            do i = 1, size(distances)
                distances(i) = array(i)%evaluate(packet%pos)
                !if(distances(i) > 0._wp)distances(i)=-999.0_wp
            end do
            packet%layer=maxloc(distances,dim=1, mask=(distances<=0._wp))
            
            if(state%trackHistory)call history%push(vec4(packet%pos, packet%step))
            ! Find scattering location
            call tauint2(state%grid, packet, array, dects, history)

            do while(.not. packet%tflag)
                if(state%trackHistory)call history%push(vec4(packet%pos, packet%step))
                ran = ran2()

                if(ran < array(packet%layer)%getAlbedo()) then !interacts with tissue
                    call packet%scatter(array(packet%layer)%gethgg(), &
                                        array(packet%layer)%getg2())
                    nscatt = nscatt + 1
                    packet%step = packet%step + 1
                else
                    packet%tflag = .true.
                    call recordWeight(packet, 1.0_wp)
                    exit
                end if
                ! Find next scattering location
                call tauint2(state%grid, packet, array, dects, history)
            end do

            if(id == 0 .and. mod(j,1000) == 0)then
                if(state%tev)then
!$omp critical
                    image = reshape(jmean(:,100:100,:), [state%grid%nxg,state%grid%nzg,1])
                    call tev%update_image(state%experiment, real(image(:,:,1:1)), ["I"], 0, 0, .false., .false.)

                    image = reshape(phasor(100:100,:,:), [state%grid%nyg,state%grid%nzg,1])
                    call tev%update_image(state%experiment, real(image(:,:,1:1)), ["J"], 0, 0, .false., .false.)

                    image = reshape(phasor(:,:,100:100), [state%grid%nxg,state%grid%nyg,1])
                    call tev%update_image(state%experiment, real(image(:,:,1:1)), ["K"], 0, 0, .false., .false.)
!$omp end critical
                end if
            end if
        end do

#ifdef _OPENMP
!$OMP end  do
!$OMP end parallel
toc=omp_get_wtime()
#else
    call cpu_time(toc)
#endif
    print*,"Photons/s: ",(state%nphotons / (toc - tic))

    call finalise(dict, dects, nscatt, start, history)
    end subroutine pathlength_scatter


    !Partial weight reduction with survival biasing as a variance reduction technique
    subroutine pathlength_scatter2(input_file)

        !Shared data
        use iarray
        use constants, only : wp, CHANCE, THRESHOLD

        !subroutines
        use detector_mod,  only : hit_t
        use detectors,     only : dect_array
        use historyStack,  only : history_stack_t
        use inttau2,       only : tauint2
        use photonMod,     only : photon
        use piecewiseMod
        use random,        only : ran2, init_rng, seq
        use sdfs,          only : sdf
        use sim_state_mod, only : state
        use utils,         only : pbar
        use vec4_class,    only : vec4
        use vector_class,  only : vector
        use writer_mod,    only : checkpoint

        !external deps
        use tev_mod, only : tevipc
        use tomlf,   only : toml_table
#ifdef _OPENMP
        use omp_lib
#endif
        character(len=*), intent(in) :: input_file
        
        integer                       :: numproc, id, j, i
        type(history_stack_t)         :: history
        type(pbar)                    :: bar
        type(photon)                  :: packet
        type(toml_table)              :: dict
        real(kind=wp),    allocatable :: distances(:), image(:,:,:)
        type(hit_t)                   :: hpoint
        type(vector)                  :: dir
        type(dect_array), allocatable :: dects(:)
        type(sdf),        allocatable :: array(:)
        real(kind=wp)                 :: ran, nscatt, start, weight_absorb
        type(tevipc)                  :: tev
        type(seq)                     :: seqs(2)
        type(spectrum_t)              :: spectrum
        real :: tic, toc

        integer :: nphotons_run,pos
        character(len=128) :: line
        character(len=:), allocatable :: checkpt_input_file

        if(state%loadckpt)then
            call setup(input_file, tev, dects, array, packet, spectrum, dict, distances, image, nscatt, start, .false.)
            open(newunit=j,file=state%ckptfile, access="stream", form="formatted")
            read(j,"(a)")line
            pos = scan(line, "=")
            checkpt_input_file = trim(line(pos+1:))

            read(j,"(a)")line
            pos = scan(line, "=")
            read(line(pos+1:),*) nphotons_run

            inquire(j,pos=pos)
            close(j)

            open(newunit=j,file=state%ckptfile, access="stream", form="unformatted")
            read(j,pos=pos)jmean
            close(j)

            call setup(checkpt_input_file, tev, dects, array, packet, spectrum, dict, distances, image, nscatt, start, .true.)
            state%iseed=state%iseed*101
            state%nphotons = state%nphotons - nphotons_run
        else
            call setup(input_file, tev, dects, array, packet, spectrum, dict, distances, image, nscatt, start, .true.)
        end if

#ifdef _OPENMP
        tic=omp_get_wtime()
!$omp parallel default(none)& 
        !$omp& shared(dict, array, numproc, start, bar, jmean, emission, absorb, input_file, phasor, tev, dects, spectrum)&
        !$omp& private(ran, id, distances, image, dir, hpoint, history, seqs, weight_absorb)& 
        !$omp& reduction(+:nscatt) firstprivate(state, packet)
        numproc = omp_get_num_threads()
        id = omp_get_thread_num()
        if(numproc > state%nphotons .and. id == 0)print*,"Warning, simulation may be underministic due to low photon count!"
        if(state%trackHistory)history = history_stack_t(state%historyFilename, id)
#elif MPI
    !nothing
#else
        call cpu_time(tic)
        numproc = 1
        id = 0
        if(state%trackHistory)history = history_stack_t(state%historyFilename, id)
#endif
        if(id == 0)print("(a,I3.1,a)"),'Photons now running on', numproc,' cores.'
        state%iseed = state%iseed + id
        ! set seed for rnd generator. id to change seed for each process
        call init_rng(state%iseed, fwd=.true.)
        seqs = [seq((id+1)*(state%nphotons/numproc), 2),&
                seq((id+1)*(state%nphotons/numproc), 3)]

        bar = pbar(state%nphotons/ 10)

        !$OMP BARRIER
        !$OMP do
        !loop over photons
        do j = 1, state%nphotons
            if(mod(j, 10) == 0)call bar%progress()
            if(mod(j, state%ckptfreq) == 0 .and. id==0)call checkpoint(input_file, state%ckptfile, j, .true.)

            ! Release photon from point source
            call packet%emit(spectrum, dict, seqs)

            do while (packet%xcell < 1 .or. packet%xcell > state%grid%nxg .or. &
                        packet%ycell < 1 .or. packet%ycell > state%grid%nyg .or. &
                        packet%zcell < 1 .or. packet%zcell > state%grid%nzg)
                call packet%emit(spectrum, dict, seqs)
            end do

            if(state%render_source)call recordEmission(packet)
            packet%step = 0
            packet%id = id
            distances = 0._wp
            do i = 1, size(distances)
                distances(i) = array(i)%evaluate(packet%pos)
                !if(distances(i) > 0._wp)distances(i)=-999.0_wp
            end do
            packet%layer=maxloc(distances,dim=1, mask=(distances<=0._wp))
            
            if(state%trackHistory)call history%push(vec4(packet%pos, packet%step))
            ! Find scattering location
            call tauint2(state%grid, packet, array, dects, history)

            do while(.not. packet%tflag)
                if(state%trackHistory)call history%push(vec4(packet%pos, packet%step))
                ran = ran2()

                !Reduce the packet weight
                weight_absorb = packet%weight * (1._wp - array(packet%layer)%getAlbedo())
                packet%weight = packet%weight - weight_absorb
            
                call recordWeight(packet, weight_absorb)

                ! is the packet weight below a threshold
                if(packet%weight < THRESHOLD)then
                    !yes, then put through roulette
                    if(ran < CHANCE)then
                        ! survive, continue emission with higher weight
                        packet%weight = packet%weight / CHANCE
                    else
                        !doesn't survive, don't re-emit
                        packet%tflag = .true.
                        exit
                    end if
                end if

                ! scatter the particle
                call packet%scatter(array(packet%layer)%gethgg(), array(packet%layer)%getg2())
                nscatt = nscatt + 1
                packet%step = packet%step + 1

                ! Find next scattering location
                call tauint2(state%grid, packet, array, dects, history)
            end do

            ! Used to test if the detectors were working, they are
            !if(packet%weight > THRESHOLD)then
            !    !the packet has left the scene, not been absorbed
            !    !record the weight as either diffuse transmission or diffuse reflection
            !    !print*, packet%xcell, packet%ycell, packet%zcell 
            !    !print*, packet%pos%x, packet%pos%y, packet%pos%z
            !
            !    if(packet%xcell == -1) then
            !        if(packet%pos%x >= 0.0_wp) then
            !            packet%xcell = state%grid%nxg
            !        else if(packet%pos%x < 0.0_wp) then
            !            packet%xcell = 1
            !        end if
            !    end if
            !    if (packet%ycell == -1) then
            !        if(packet%pos%y >= 0.0_wp) then
            !            packet%ycell = state%grid%nyg
            !        else if(packet%pos%y < 0.0_wp) then
            !            packet%ycell = 1
            !        end if
            !    end if
            !    if (packet%zcell == -1) then
            !        if(packet%pos%z >= 0.0_wp) then
            !            packet%zcell = state%grid%nzg
            !        else if(packet%pos%z < 0.0_wp) then
            !            packet%zcell = 1
            !        end if
            !    end if
            !    !print*, packet%xcell, packet%ycell, packet%zcell
            !    !print*, packet%pos%x, packet%pos%y, packet%pos%z
            !    call recordWeight(packet, packet%weight)
            !end if

            if(id == 0 .and. mod(j,1000) == 0)then
                if(state%tev)then
!$omp critical
                    image = reshape(jmean(:,100:100,:), [state%grid%nxg,state%grid%nzg,1])
                    call tev%update_image(state%experiment, real(image(:,:,1:1)), ["I"], 0, 0, .false., .false.)

                    image = reshape(phasor(100:100,:,:), [state%grid%nyg,state%grid%nzg,1])
                    call tev%update_image(state%experiment, real(image(:,:,1:1)), ["J"], 0, 0, .false., .false.)

                    image = reshape(phasor(:,:,100:100), [state%grid%nxg,state%grid%nyg,1])
                    call tev%update_image(state%experiment, real(image(:,:,1:1)), ["K"], 0, 0, .false., .false.)
!$omp end critical
                end if
            end if
        end do

#ifdef _OPENMP
!$OMP end  do
!$OMP end parallel
toc=omp_get_wtime()
#else
    call cpu_time(toc)
#endif
    print*,"Photons/s: ",(state%nphotons / (toc - tic))

    call finalise(dict, dects, nscatt, start, history)
    end subroutine pathlength_scatter2

    subroutine test_kernel(input_file, end_early)

        !Shared data
        use iarray
        use constants, only : wp

        !subroutines
        use detectors,     only : dect_array
        use historyStack,  only : history_stack_t
        use inttau2,       only : tauint2
        use photonMod,     only : photon
        use piecewiseMod
        use random,        only : ran2, init_rng
        use sdfs,          only : sdf
        use sim_state_mod, only : state
        use utils,         only : pbar
        use vector_class,  only : vector
        
        !external deps
        use tev_mod, only : tevipc
        use tomlf,   only : toml_table
#ifdef _OPENMP
        use omp_lib
#endif
        character(len=*), intent(in) :: input_file
        
        integer :: numproc, id, j, i
        type(history_stack_t) :: history
        ! type(pbar)        :: bar
        type(photon)      :: packet
        type(toml_table)  :: dict
        real(kind=wp), allocatable :: distances(:), image(:,:,:)
        type(dect_array),  allocatable :: dects(:)

        type(sdf),   allocatable :: array(:)
        real(kind=wp) :: ran, nscatt, start
        type(tevipc)      :: tev
        type(vector)  :: pos(4), pos2(4)
        logical, intent(in) :: end_early
        type(spectrum_t) :: spectrum

        pos = vector(0.0_wp, 0.0_wp, 0.0_wp)
        pos2 = vector(0.0_wp, 0.0_wp, 0.0_wp)
        call setup(input_file, tev, dects, array, packet, spectrum, dict, distances, image, nscatt, start)

        numproc = 1
        id = 0

        if(id == 0)print("(a,I3.1,a)"),'Photons now running on', numproc,' cores.'

        ! set seed for rnd generator. id to change seed for each process
        call init_rng(state%iseed, fwd=.true.)

        ! bar = pbar(state%nphotons/ 10)
        !loop over photons 
        do j = 1, state%nphotons
            ! if(mod(j, 10) == 0)call bar%progress()

            ! Release photon from point source
            call packet%emit(spectrum, dict)
            packet%step = 0
            packet%id = id
            distances = 0._wp
            do i = 1, size(distances)
                distances(i) = array(i)%evaluate(packet%pos)
                !if(distances(i) > 0._wp)distances(i)=-999.0_wp
            end do
            packet%layer=maxloc(distances,dim=1, mask=(distances<=0._wp))

            ! Find scattering location
            call tauint2(state%grid, packet, array, dects, history)


            do while(.not. packet%tflag)
                ran = ran2()
                if(ran < array(packet%layer)%getalbedo())then!interacts with tissue
                    call packet%scatter(array(packet%layer)%gethgg(), &
                                        array(packet%layer)%getg2())
                    nscatt = nscatt + 1
                    packet%step = packet%step + 1
                    if(packet%step == 1)then
                        pos(1) = pos(1) + packet%pos
                        pos2(1) = pos2(1) + packet%pos**2
                    elseif(packet%step == 2)then
                        pos(2) = pos(2) + packet%pos
                        pos2(2) = pos2(2) + packet%pos**2
                    elseif(packet%step == 3)then
                        pos(3) = pos(3) + packet%pos
                        pos2(3) = pos2(3) + packet%pos**2
                    elseif(packet%step == 4)then
                        pos(4) = pos(4) + packet%pos
                        pos2(4) = pos2(4) + packet%pos**2
                    else
                        if(end_early)packet%tflag = .true.
                    end if
                else
                    packet%tflag = .true.
                    exit
                end if
                ! !Find next scattering location
                call tauint2(state%grid, packet, array, dects, history)
            end do
        end do

    open(newunit=j,file="positions.dat")
    do i = 1, 4
        write(j,*)10.*pos(i)%x/state%nphotons,10.*pos(i)%y/state%nphotons,10.*pos(i)%z/state%nphotons
    end do
    do i = 1,4
        write(j,*)100.*pos2(i)%x/state%nphotons,100.*pos2(i)%y/state%nphotons,100.*pos2(i)%z/state%nphotons
    end do
    close(j)
    call finalise(dict, dects, nscatt, start, history)
    end subroutine test_kernel

    subroutine recordEmission(packet)
        !! record emission using path length estimators. Uses voxel grid
        use photonMod
        use iarray,     only: phasor, jmean, emission, absorb
        use constants , only : sp
        
        !> packet stores the photon related variables
        type(photon),    intent(IN) :: packet
        
        integer       :: celli, cellj, cellk
        celli = packet%xcell
        cellj = packet%ycell
        cellk = packet%zcell

!$omp atomic
        emission(celli,cellj,cellk) = emission(celli,cellj,cellk) + real(1.0, kind=sp)
    end subroutine recordEmission

    subroutine recordWeight(packet, weightAbsorbed)
        !! record weight absorbed
        use photonMod
        use iarray,     only: phasor, jmean, emission, absorb
        use constants , only : wp

        !> packet stores the photon related variables
        type(photon),    intent(IN) :: packet
        !> weight absorbed at this point in space
        real(kind=wp),   intent(IN) :: weightAbsorbed
        
        integer       :: celli, cellj, cellk
        celli = packet%xcell
        cellj = packet%ycell
        cellk = packet%zcell

!$omp atomic
        absorb(celli, cellj, cellk) = absorb(celli, cellj, cellk) + weightAbsorbed
    end subroutine recordWeight


!####################################################################################################
!                           Setup and break down helper routines
    subroutine setup(input_file, tev, dects, array, packet, spectrum, dict, distances, image, nscatt, start, display)
        !! setup simulation by reading in setting file, and setup variables to be used.
        
        !shared data
        use iarray
        use constants, only : wp
        
        !subroutines
        use detectors,     only : dect_array
        use parse_mod,     only : parse_params
        use photonMod,     only : photon
        use random,        only : init_rng
        use piecewiseMod

        use sdfs,          only : sdf, render
        use sim_state_mod, only : state
        use setupMod,      only : setup_simulation, directory
        use utils,         only : get_time, print_time, str
        use vector_class,  only : vector
        ! !external deps
        use tev_mod, only : tevipc, tev_init
        use tomlf,   only : toml_table, toml_error
        
        !> Filename for toml settings to be used
        character(*),                  intent(in)  :: input_file
        !> array of SDF objects that create the geometry
        type(sdf),        allocatable, intent(out) :: array(:)
        !> array of photon detectors
        type(dect_array), allocatable, intent(out) :: dects(:)
        !> toml table of meta-data to be written to output files.
        type(toml_table),              intent(out) :: dict
        !> handle for communicating with TEV
        type(tevipc),                  intent(out) :: tev
        !> photon that is to be simulated
        type(photon),                  intent(out) :: packet
        real(kind=wp),    allocatable, intent(out) :: distances(:), image(:,:,:)
        real(kind=wp),                 intent(out) :: nscatt, start
        type(spectrum_t),              intent(out) :: spectrum
        !> flag to display simulation init settings
        logical, optional,             intent(in) :: display
        
        ! mpi/mp variables
        integer       :: id
        real(kind=wp) :: chance, threshold
        type(toml_error), allocatable :: error
        logical :: disp

        if(present(display))then
            disp = display
        else
            disp = .true.
        end if

        chance = 1._wp/10._wp
        threshold = 1e-6_wp
        
        call directory()

        ! Read in the toml
        dict = toml_table()
        call parse_params("res/"//trim(input_file), packet, dects, spectrum, dict, error)
        if(allocated(error))then
            print*,error%message
            stop 1
        end if
        allocate(image(state%grid%nxg,state%grid%nzg,1))
        
        if(disp)call display_settings(state, input_file, packet, "Pathlength")

        if(state%tev)then
            !init TEV link
            tev = tevipc()
            call tev%close_image(state%experiment)
            call tev%create_image(state%experiment, state%grid%nxg, state%grid%nzg, ["I", "J", "K"], .true.)
        end if

        nscatt = 0._wp
        !call init_rng(state%iseed, fwd=.true.)
        
        call setup_simulation(array, dict)
        ! render geometry to voxel format for debugging
        if(state%render_geom)then
            print*,"Rendering geometry to file"
            call render(array, state)
        end if
        
        allocate(distances(size(array)))
        
        start = get_time()
        id = 0            
        
        if(id == 0)then
           print*,'# of photons to run',state%nphotons
        end if
end subroutine setup

subroutine finalise(dict, dects, nscatt, start, history)
    !! Routine writes out simulation data, deallocates arrays and prints total runtime
    use constants,     only : wp, fileplace
    use detectors,     only : dect_array
    use historyStack,  only : history_stack_t
    use iarray,        only : phasor, phasorGLOBAL, jmean, jmeanGLOBAL, absorb, absorbGLOBAL, emission, emissionGLOBAL
    use sim_state_mod, only : state
    use setupMod,      only : dealloc_array
    use writer_mod,    only : normalise_fluence, write_data, write_detected_photons
    
    use utils, only : get_time, print_time, str
    use tomlf, only : toml_table, set_value

    !> Total number of scattered photon packets
    real(kind=wp),         intent(in)    :: nscatt
    !> Start time of simulation. Used to calculate total runtime.
    real(kind=wp),         intent(in)    :: start
    !> Detector array
    type(dect_array),      intent(in)    :: dects(:)
    !> Photon histyor object
    type(history_stack_t), intent(in)    :: history
    !> Dictionary of metadata
    type(toml_table),      intent(inout) :: dict

    integer       :: id, numproc, i
    real(kind=wp) :: nscattGLOBAL, time_taken

    id = 0
    numproc = 1

#ifdef MPI
    ! collate fluence from all processes
    call mpi_reduce(jmean, jmeanGLOBAL, size(jmean),MPI_DOUBLE_PRECISION, MPI_SUM,0,MPI_COMM_WORLD)
    call mpi_reduce(absorb, absorbGLOBAL, size(absorb),MPI_DOUBLE_PRECISION, MPI_SUM,0,MPI_COMM_WORLD)
    call mpi_reduce(phasor, phasorGLOBAL, size(phasor),MPI_DOUBLE_COMPLEX, MPI_SUM,0,MPI_COMM_WORLD)
    call mpi_reduce(emission, emissionGLOBAL, size(emission),MPI_DOUBLE_PRECISION, MPI_SUM,0,MPI_COMM_WORLD)
    call mpi_reduce(nscatt,nscattGLOBAL,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD)
#else
    jmeanGLOBAL = jmean
    absorbGLOBAL = absorb
    phasorGLOBAL = phasor
    emissionGLOBAL = emission
    nscattGLOBAL = nscatt
#endif

    if(id == 0)then
#ifdef _OPENMP
        print*,'Average # of scatters per photon:',nscattGLOBAL/(state%nphotons)
#else
        print*,'Average # of scatters per photon:',nscattGLOBAL/(state%nphotons*numproc)
        ! for testing purposes
        open(newunit=i,file="nscatt.dat")
        write(i,*)nscattGLOBAL/(state%nphotons)
        close(i)
#endif
        !write out files
        !create dict to store metadata and nrrd hdr info
        call set_value(dict, "grid_data", "fluence map")
        call set_value(dict, "real_size", str(state%grid%xmax,7)//" "//str(state%grid%ymax,7)//" "//str(state%grid%zmax,7))
        call set_value(dict, "nphotons", state%nphotons)
        call set_value(dict, "source", state%source)
        call set_value(dict, "experiment", state%experiment)

#ifdef pathlength
        call normalise_fluence(state%grid, jmeanGLOBAL, state%nphotons)
        call write_data(jmeanGLOBAL, trim(fileplace)//"jmean/"//state%outfile, state, dict)
#endif

        call normalise_fluence(state%grid, emissionGLOBAL, state%nphotons)
        call write_data(emissionGLOBAL, trim(fileplace)//"emission/"//state%rendersourcefile, state, dict)

        call write_data(absorbGLOBAL, trim(fileplace)//"absorb/"//"absorb.nrrd", state, dict)
        ! if(state%absorb)call write_data(absorbGLOBAL, trim(fileplace)//"deposit/"//state%outfile_absorb, state, dict)
        !INTENSITY
        ! call write_data(abs(phasorGLOBAL)**2, trim(fileplace)//"phasor/"//state%outfile, state, dict)    
    end if
    
    !write out detected photons
    if(size(dects) > 0)then
        call write_detected_photons(dects)
        block
            logical :: mask(size(dects))
            do i = 1, size(dects)
                mask(i) = dects(i)%p%trackHistory
            end do
            if(state%trackHistory)call history%finish()
        end block
    end if

    time_taken = get_time() - start
    call print_time(time_taken, 4)
#ifdef MPI
    call MPI_Finalize()
#endif
    call dealloc_array()
end subroutine finalise

subroutine display_settings(state, input_file, packet, kernel_type)
    !! Displays the settings used in the current simulation run

    use sim_state_mod, only : settings_t
    use photonMod,     only : photon
    use utils,         only : str

    !> Simulation state
    type(settings_t), intent(IN) :: state
    !> Input filenname
    character(*),     intent(IN) :: input_file
    !> Kernel type to run
    character(*),     intent(IN) :: kernel_type
    !> Photon packet
    type(photon),     intent(IN) :: packet

    print*,repeat("#", 20)//" Settings "//repeat("#", 20)
    print*,"# Config file: ",trim(input_file),repeat(" ", 50-16-len(trim(input_file))),"#"
    print*,"# Using: "//trim(kernel_type)//"kernel"//repeat(" ", 50-16-len(kernel_type)),"#"
    print*,"# Light source: "//trim(state%source)//repeat(" ", 50-17-len(trim(state%source))),"#"
    if(state%source == "point")then
        print*,"# Light Source Position: ["//str(packet%pos%x,4)//", "//str(packet%pos%y,4)//", "//str(packet%pos%z,4)// &
                                        "]"//repeat(" ", 6)//"#"
    else
        print*,"# Light direction: ["//str(packet%nxp,4)//", "//str(packet%nyp,4)//", "//str(packet%nzp,4)// &
                                  "]"//repeat(" ", 12)//"#"
    end if
    print*,"# Geometry: "//trim(state%experiment)//repeat(" ", 50-13-len(trim(state%experiment))),"#"
    print*,"# Seed: "//str(state%iseed,9)//repeat(" ", 32)//"#"
    if(state%tev)then
        print*,"# Tev enabled!"//repeat(" ", 35)//"#"
    end if
    if(state%render_geom)then
        print*,"# Render geometry to file enabled!"//repeat(" ", 15)//"#"
    end if
    if(state%overwrite)then
        print*,"# Overwrite Enabled!",repeat(" ", 29)//"#"
    end if
    if(state%absorb)then
        print*,"# Energy absorbed will be written to file."//repeat(" ", 7)//"#"
    end if
    print*,repeat("#", 50)
    print*,new_line("a")

end subroutine display_settings
end module kernels