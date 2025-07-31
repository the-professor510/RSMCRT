module raman_MCRTMod

    implicit none

    private
    public :: raman_MCRT

contains
    subroutine raman_MCRT(input_file)

        !Shared data
        use iarray
        use constants, only : wp

        !subroutines
        use detectors,     only : dect_array
        use historyStack,  only : history_stack_t
        use photonMod,     only : photon
        use piecewiseMod
        use random,        only : ran2, init_rng
        use sdfs,          only : sdf
        use sim_state_mod, only : state
        use kernels, only : setup, finalise, reset_detectors
        use writer_mod, only : write_raman
        use setupMod, only : setup_escapeFunction

        !external deps
        use tev_mod, only : tevipc
        use tomlf,   only : toml_table
#ifdef _OPENMP
        use omp_lib
#endif
        character(len=*), intent(in) :: input_file
        
        integer                       :: j
        type(history_stack_t)         :: history
        type(photon)                  :: packet
        type(toml_table)              :: dict
        real(kind=wp),    allocatable :: distances(:), image(:,:,:)
        type(dect_array), allocatable :: dects(:)
        type(sdf),        allocatable :: array(:)
        real(kind=wp)                 :: nscatt, start
        type(spectrum_t)              :: spectrum
        type(tevipc)                  :: tev

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

        !setup the escape function
        call setup_escapeFunction(size(dects))

        call run_RMCRT(input_file, history, packet, dict, & 
                        distances, image, dects, array, nscatt, start, & 
                        tev, spectrum)

        !store the escape funcitons for each detector
        call write_raman(dects, dict)

        call finalise(dict, dects, nscatt, start, history)
    end subroutine raman_MCRT


    subroutine run_RMCRT(input_file, history, packet, dict, & 
                        distances, image, dects, array, nscatt, start, & 
                        tev, spectrum)
        !Shared data
        use iarray
        use constants, only : wp

        !subroutines
        use detectors,     only : dect_array
        use historyStack,  only : history_stack_t
        use inttau2,       only : tauint2
        use photonMod,     only : photon
        use piecewiseMod
        use random,        only : ran2, init_rng, seq
        use sdfs,          only : sdf
        use sim_state_mod, only : state
        use utils,         only : pbar
        use writer_mod,    only : checkpoint

        !external deps
        use tev_mod, only : tevipc
        use tomlf,   only : toml_table
#ifdef _OPENMP
        use omp_lib
#endif
        character(len=*),                 intent(in) :: input_file
         
        type(history_stack_t),         intent(inout) :: history
        type(pbar)                                   :: bar
        type(photon),                  intent(inout) :: packet
        type(toml_table),              intent(inout) :: dict
        real(kind=wp),    allocatable, intent(inout) :: distances(:), image(:,:,:)
        type(dect_array), allocatable, intent(inout) :: dects(:)
        type(sdf),        allocatable, intent(inout) :: array(:)
        real(kind=wp),                 intent(inout) :: nscatt, start
        type(tevipc),                  intent(inout) :: tev
        type(seq)                                    :: seqs(2)
        type(spectrum_t),              intent(inout) :: spectrum
        real :: tic, toc
        integer :: numproc, id, j

#ifdef _OPENMP
        tic=omp_get_wtime()
        !$omp parallel default(none)& 
        !$omp& shared(dict, array, numproc, start, bar, jmean, emission, absorb, escape, input_file, phasor, tev, spectrum)& 
        !$omp& private(id, distances, image, history, seqs)& 
        !$omp& reduction(+:nscatt) firstprivate(state, packet, dects)
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

            !launch and propagate packets, either dropping full weight or partial weights (survival bias)
#ifdef survivalBias
            call survivalBiasPropagation(id, history, packet, dict, distances, image, dects, array,& 
                                        nscatt, seqs, spectrum)
#else
            call noBiasPropagation(id, history, packet, dict, distances, image, dects, array,& 
                                        nscatt, seqs, spectrum)
#endif

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
        !$OMP end  do
        
#ifdef _OPENMP
        !$OMP end parallel
        toc=omp_get_wtime()
#else
        call cpu_time(toc)
#endif
        print*,"Photons/s: ",(state%nphotons / (toc - tic))
    end subroutine run_RMCRT

    !Full weight reduction
    subroutine noBiasPropagation(id, history, packet, dict, distances, image, dects, array,& 
                                nscatt, seqs, spectrum)

        !Shared data
        use iarray
        use constants, only : wp

        !subroutines
        use detectors,     only : dect_array
        use historyStack,  only : history_stack_t
        use inttau2,       only : tauint2
        use photonMod,     only : photon
        use piecewiseMod
        use random,        only : ran2, seq
        use sdfs,          only : sdf
        use sim_state_mod, only : state
        use vec4_class,    only : vec4


        use kernels, only : reset_detectors


        !external deps
        use tomlf,   only : toml_table
        
        integer,                       intent(inout) :: id
        type(history_stack_t),         intent(inout) :: history
        type(photon),                  intent(inout) :: packet
        type(toml_table),              intent(inout) :: dict
        real(kind=wp),    allocatable, intent(inout) :: distances(:), image(:,:,:)
        type(dect_array), allocatable, intent(inout) :: dects(:)
        type(sdf),        allocatable, intent(inout) :: array(:)
        real(kind=wp),                 intent(inout) :: nscatt
        type(seq),                     intent(inout) :: seqs(2)
        type(spectrum_t),              intent(inout) :: spectrum

        real(kind=wp)   :: ran, total
        integer         :: i

        ! Release photon from source
        call packet%emit(spectrum, dict, seqs)

        do while (packet%xcell < 1 .or. packet%xcell > state%grid%nxg .or. &
                    packet%ycell < 1 .or. packet%ycell > state%grid%nyg .or. &
                    packet%zcell < 1 .or. packet%zcell > state%grid%nzg)
            call packet%emit(spectrum, dict, seqs)
        end do

        if(state%render_source)call recordEmissionLocation(packet)
        packet%step = 0
        packet%id = id
        distances = 0._wp
        do i = 1, size(distances)
            distances(i) = array(i)%evaluate(packet%pos)
        end do
        packet%layer=(maxloc(distances,dim=1, mask=(distances<0._wp)))
        
        if(state%trackHistory)call history%push(vec4(packet%pos, packet%step))
        ! Find scattering location
        call tauint2(state%grid, packet, array, dects, history)

        do while(.not. packet%tflag)
            if(state%trackHistory)call history%push(vec4(packet%pos, packet%step))
            ran = ran2()

            if(ran < array(packet%layer)%getAlbedo()) then !interacts with tissue
                ran = ran2()
                if (ran< 0.01) then !raman scatters
                    ! launch a new particle with the same weight as the current photon, 
                    ! in an isotropic direction. Record if it hits a detector in escape(:,:,:,:) 

                    !zero detectors
                    call reset_detectors(dects)

                    !propagate the Raman source
#ifdef survivalBias
                    call survivalBiasRaman(id, history, packet, dict, distances, image, dects, array,& 
                                        nscatt, seqs, spectrum)
#else
                    call noBiasRaman(id, history, packet, dict, distances, image, dects, array,& 
                                        nscatt, seqs, spectrum)
#endif


                    !we now need to loop through detectors and add it to the escape function but atomic                    
                    do i = 1, size(dects)
                        
                        total = 0._wp
                        call dects(i)%p%total_dect(total)

                        !$omp atomic
                        escape(i, packet%xcell, packet%ycell, packet%zcell)=escape(i, packet%xcell, packet%ycell, packet%zcell)&
                                                                            + total
                    end do
                    
                end if

                call packet%scatter(array(packet%layer)%gethgg(), &
                                    array(packet%layer)%getg2())
                nscatt = nscatt + 1
                packet%step = packet%step + 1
            else
                packet%tflag = .true.
                !record the fluence and absorption
                call recordWeight(packet, 1.0_wp, array(packet%layer)%getMua())
                exit
            end if
            ! Find next scattering location
            call tauint2(state%grid, packet, array, dects, history)
        end do

    end subroutine noBiasPropagation

    !Partial weight reduction with survival biasing as a variance reduction technique
    subroutine survivalBiasPropagation(id, history, packet, dict, distances, image, dects, array,& 
                                        nscatt, seqs, spectrum)

        !Shared data
        use iarray
        use constants, only : wp, CHANCE, THRESHOLD

        !subroutines
        use detectors,     only : dect_array
        use historyStack,  only : history_stack_t
        use inttau2,       only : tauint2
        use photonMod,     only : photon
        use piecewiseMod
        use random,        only : ran2, seq
        use sdfs,          only : sdf
        use sim_state_mod, only : state
        use vec4_class,    only : vec4

        use kernels, only : reset_detectors

        !external deps
        use tomlf,   only : toml_table
        
        integer,                       intent(inout) :: id
        type(history_stack_t),         intent(inout) :: history
        type(photon),                  intent(inout) :: packet
        type(toml_table),              intent(inout) :: dict
        real(kind=wp),    allocatable, intent(inout) :: distances(:), image(:,:,:)
        type(dect_array), allocatable, intent(inout) :: dects(:)
        type(sdf),        allocatable, intent(inout) :: array(:)
        real(kind=wp),                 intent(inout) :: nscatt
        type(seq),                     intent(inout) :: seqs(2)
        type(spectrum_t),              intent(inout) :: spectrum

        real(kind=wp)   :: ran, weight_absorb, total
        integer         :: i

        ! Release photon from point source
        call packet%emit(spectrum, dict, seqs)

        do while (packet%xcell < 1 .or. packet%xcell > state%grid%nxg .or. &
                    packet%ycell < 1 .or. packet%ycell > state%grid%nyg .or. &
                    packet%zcell < 1 .or. packet%zcell > state%grid%nzg)
            call packet%emit(spectrum, dict, seqs)
        end do

        if(state%render_source)call recordEmissionLocation(packet)
        packet%step = 0
        packet%id = id
        distances = 0._wp
        do i = 1, size(distances)
            distances(i) = array(i)%evaluate(packet%pos)
        end do
        packet%layer=maxloc(distances,dim=1, mask=(distances<0._wp))
        
        if(state%trackHistory)call history%push(vec4(packet%pos, packet%step))
        ! Find scattering location
        call tauint2(state%grid, packet, array, dects, history)

        do while(.not. packet%tflag)
            if(state%trackHistory)call history%push(vec4(packet%pos, packet%step))
            ran = ran2()

            !Reduce the packet weight
            weight_absorb = packet%weight * (1._wp - array(packet%layer)%getAlbedo())
            packet%weight = packet%weight - weight_absorb
        
            !record the fluence and absorption
            call recordWeight(packet, weight_absorb, array(packet%layer)%getMua())

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

            ran = ran2()
            if (ran< 0.01) then !raman scatters
                ! launch a new particle with the same weight as the current photon, 
                ! in an isotropic direction. Record if it hits a detector in escape(:,:,:,:) 

                !zero detectors
                call reset_detectors(dects)

                !propagate the Raman source
#ifdef survivalBias
                call survivalBiasRaman(id, history, packet, dict, distances, image, dects, array,& 
                                    nscatt, seqs, spectrum)
#else
                call noBiasRaman(id, history, packet, dict, distances, image, dects, array,& 
                                    nscatt, seqs, spectrum)
#endif

                
                !we now need to loop through detectors and add it to the escape function but atomic                    
                do i = 1, size(dects)
                    
                    total = 0._wp
                    call dects(i)%p%total_dect(total)

                    !$omp atomic
                    escape(i, packet%xcell, packet%ycell, packet%zcell)=escape(i, packet%xcell, packet%ycell, packet%zcell)&
                                                                        + total
                end do
                
            end if

            ! scatter the particle
            call packet%scatter(array(packet%layer)%gethgg(), array(packet%layer)%getg2())
            nscatt = nscatt + 1
            packet%step = packet%step + 1

            ! Find next scattering location
            call tauint2(state%grid, packet, array, dects, history)
        end do
    end subroutine survivalBiasPropagation






    subroutine recordEmissionLocation(packet)
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
    end subroutine recordEmissionLocation

    subroutine recordWeight(packet, weightAbsorbed, absorptionCoefficient)
        !! record energy absorbed and fluence
        use photonMod
        use iarray,     only: phasor, jmean, emission, absorb
        use constants , only : wp

        !> packet stores the photon related variables
        type(photon),    intent(IN) :: packet
        !> weight absorbed at this point in space
        real(kind=wp),   intent(IN) :: weightAbsorbed
        !> Absoption Ceofficient at this point in space
        real(kind=wp),  intent(IN) :: absorptionCoefficient
        
        integer       :: celli, cellj, cellk
        real(kind=wp) :: mua
        celli = packet%xcell
        cellj = packet%ycell
        cellk = packet%zcell

        if (absorptionCoefficient == 0.0_wp) then
            mua = 1e-15_wp
        else 
            mua = absorptionCoefficient
        end if

!$omp atomic
        absorb(celli, cellj, cellk) = absorb(celli, cellj, cellk) + weightAbsorbed
        jmean(celli, cellj, cellk) = jmean(celli, cellj, cellk) + weightAbsorbed/mua
    end subroutine recordWeight









        !Full weight reduction
    subroutine noBiasRaman(id, history, packet, dict, distances, image, dects, array,& 
                                nscatt, seqs, spectrum)

        !Shared data
        use iarray
        use constants, only : wp

        !subroutines
        use detectors,     only : dect_array
        use historyStack,  only : history_stack_t
        use inttau2,       only : tauint2
        use photonMod,     only : photon, photon_origin, set_photon
        use piecewiseMod
        use random,        only : ran2, seq
        use sdfs,          only : sdf
        use sim_state_mod, only : state
        use vec4_class,    only : vec4
        use vector_class

        !external deps
        use tomlf,   only : toml_table
        
        integer,                       intent(inout) :: id
        type(history_stack_t),         intent(inout) :: history
        type(photon),                  intent(inout) :: packet
        type(toml_table),              intent(inout) :: dict
        real(kind=wp),    allocatable, intent(inout) :: distances(:), image(:,:,:)
        type(dect_array), allocatable, intent(inout) :: dects(:)
        type(sdf),        allocatable, intent(inout) :: array(:)
        real(kind=wp),                 intent(inout) :: nscatt
        type(seq),                     intent(inout) :: seqs(2)
        type(spectrum_t),              intent(inout) :: spectrum

        real(kind=wp)   :: ran
        integer         :: i

        type(vector) :: sourcePoss, sourceDir, packetPos
        real(kind = wp) :: packetDirx, packetDiry, packetDirz, weight
        integer :: cell(3), layer


        !record the current source position and packet position and direction
        sourcePoss = photon_origin%pos
        sourceDir = vector(photon_origin%nxp, photon_origin%nyp, photon_origin%nzp)
        packetPos = packet%pos
        packetDirx = packet%nxp
        packetDiry = packet%nyp
        packetDirz = packet%nzp
        weight = packet%weight
        layer = packet%layer

        !Set photon to be a point source at current packet position
        call set_photon(packet%pos, vector(1.0_wp,0.0,0.0))
        packet = photon("point")
        packet%nxp = 1.0_wp
        packet%nyp = 0.0_wp
        packet%nzp = 0.0_wp
        packet%weight = weight

        !emit photon from point source
        call packet%emit(spectrum, dict, seqs)

        if(state%render_source)call recordEmissionLocation(packet)

        

        do while (packet%xcell < 1 .or. packet%xcell > state%grid%nxg .or. &
                    packet%ycell < 1 .or. packet%ycell > state%grid%nyg .or. &
                    packet%zcell < 1 .or. packet%zcell > state%grid%nzg)
            call packet%emit(spectrum, dict, seqs)
        end do

        
        packet%step = 0
        packet%id = id
        distances = 0._wp
        do i = 1, size(distances)
            distances(i) = array(i)%evaluate(packet%pos)
        end do
        packet%layer=(maxloc(distances,dim=1, mask=(distances<0._wp)))
        
        if(state%trackHistory)call history%push(vec4(packet%pos, packet%step))

        ! Find scattering location
        call tauint2(state%grid, packet, array, dects, history)

        do while(.not. packet%tflag)
            if(state%trackHistory)call history%push(vec4(packet%pos, packet%step))
            ran = ran2()

            if(ran < array(packet%layer)%getAlbedo()) then !interacts with tissue
                ran = ran2()
                
                call packet%scatter(array(packet%layer)%gethgg(), &
                                    array(packet%layer)%getg2())
                nscatt = nscatt + 1
                packet%step = packet%step + 1
            else
                packet%tflag = .true.
                exit
            end if

            ! Find next scattering location
            call tauint2(state%grid, packet, array, dects, history)
        end do


        
        !reset the photon source and packet position
        call set_photon(sourcePoss, sourceDir)
        packet = photon(state%source)
        packet%pos = packetPos
        packet%nxp = packetDirx
        packet%nyp = packetDiry
        packet%nzp = packetDirz
        packet%weight = weight
        packet%layer = layer

        ! update cells 
        cell = state%grid%get_voxel(packet%pos)
        packet%xcell = cell(1)
        packet%ycell = cell(2)
        packet%zcell = cell(3)

    end subroutine noBiasRaman



    !Partial weight reduction with survival biasing as a variance reduction technique
    subroutine survivalBiasRaman(id, history, packet, dict, distances, image, dects, array,& 
                                        nscatt, seqs, spectrum)

        !Shared data
        use iarray
        use constants, only : wp, CHANCE, THRESHOLD

        !subroutines
        use detectors,     only : dect_array
        use historyStack,  only : history_stack_t
        use inttau2,       only : tauint2
        use photonMod,     only : photon, photon_origin, set_photon
        use piecewiseMod
        use random,        only : ran2, seq
        use sdfs,          only : sdf
        use sim_state_mod, only : state
        use vec4_class,    only : vec4
        use vector_class

        !external deps
        use tomlf,   only : toml_table
        
        integer,                       intent(inout) :: id
        type(history_stack_t),         intent(inout) :: history
        type(photon),                  intent(inout) :: packet
        type(toml_table),              intent(inout) :: dict
        real(kind=wp),    allocatable, intent(inout) :: distances(:), image(:,:,:)
        type(dect_array), allocatable, intent(inout) :: dects(:)
        type(sdf),        allocatable, intent(inout) :: array(:)
        real(kind=wp),                 intent(inout) :: nscatt
        type(seq),                     intent(inout) :: seqs(2)
        type(spectrum_t),              intent(inout) :: spectrum

        real(kind=wp)   :: ran, weight_absorb
        integer         :: i

        type(vector) :: sourcePoss, sourceDir, packetPos
        real(kind = wp) :: packetDirx, packetDiry, packetDirz, weight
        integer :: cell(3), layer


        !record the current source position and packet position and direction
        sourcePoss = photon_origin%pos
        sourceDir = vector(photon_origin%nxp, photon_origin%nyp, photon_origin%nzp)
        packetPos = packet%pos
        packetDirx = packet%nxp
        packetDiry = packet%nyp
        packetDirz = packet%nzp
        weight = packet%weight
        layer = packet%layer

        !Set photon to be a point source at current packet position
        call set_photon(packet%pos, vector(1.0_wp,0.0,0.0))
        packet = photon("point")
        packet%nxp = 1.0_wp
        packet%nyp = 0.0_wp
        packet%nzp = 0.0_wp
        packet%weight = weight

        ! Release photon from point source
        call packet%emit(spectrum, dict, seqs)

        if(state%render_source)call recordEmissionLocation(packet)

        do while (packet%xcell < 1 .or. packet%xcell > state%grid%nxg .or. &
                    packet%ycell < 1 .or. packet%ycell > state%grid%nyg .or. &
                    packet%zcell < 1 .or. packet%zcell > state%grid%nzg)
            call packet%emit(spectrum, dict, seqs)
        end do

        if(state%render_source)call recordEmissionLocation(packet)
        packet%step = 0
        packet%id = id
        distances = 0._wp
        do i = 1, size(distances)
            distances(i) = array(i)%evaluate(packet%pos)
        end do
        packet%layer=maxloc(distances,dim=1, mask=(distances<0._wp))
        
        if(state%trackHistory)call history%push(vec4(packet%pos, packet%step))
        ! Find scattering location
        call tauint2(state%grid, packet, array, dects, history)

        do while(.not. packet%tflag)
            if(state%trackHistory)call history%push(vec4(packet%pos, packet%step))
            ran = ran2()

            !Reduce the packet weight
            weight_absorb = packet%weight * (1._wp - array(packet%layer)%getAlbedo())
            packet%weight = packet%weight - weight_absorb
        
            !record the fluence and absorption
            call recordWeight(packet, weight_absorb, array(packet%layer)%getMua())

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

        
        !reset the photon source and packet position
        call set_photon(sourcePoss, sourceDir)
        packet = photon(state%source)
        packet%pos = packetPos
        packet%nxp = packetDirx
        packet%nyp = packetDiry
        packet%nzp = packetDirz
        packet%weight = weight
        packet%layer = layer

        ! update cells 
        cell = state%grid%get_voxel(packet%pos)
        packet%xcell = cell(1)
        packet%ycell = cell(2)
        packet%zcell = cell(3)
    end subroutine survivalBiasRaman


end module raman_MCRTMod