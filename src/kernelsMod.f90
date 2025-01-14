module kernels
!! Contains the main program and scattering loop. Calls all other routine to setup, run and break down the simulation.

    implicit none
    
    private
    public :: default_MCRT, escape_Function, test_kernel

contains
!###############################################################################
!                   KERNELS

    
    subroutine default_MCRT(input_file)

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

        call run_MCRT(input_file, history, packet, dict, & 
                        distances, image, dects, array, nscatt, start, & 
                        tev, spectrum)

        call finalise(dict, dects, nscatt, start, history)
    end subroutine default_MCRT

    !calculate the escape function for each detector    
    subroutine escape_Function(input_file)

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
        use vector_class,  only : vector

        !external deps
        use tev_mod, only : tevipc
        use tomlf,   only : toml_table
#ifdef _OPENMP
        use omp_lib
#endif
        character(len=*), intent(in) :: input_file
        
        integer                       :: j, i
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
         
        integer :: m, n, o, layer, count, total
        real(kind = wp) :: x,y,z
        type(vector) :: position
        real(kind=wp),  allocatable :: escapeFunction(:,:,:,:)

        !setup the geometry and detectors
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

        !allocate(escapeFunction(size(dects), state%grid%x, state%grid%y, state%grid%z))

        !set the packet to be a isotropic source, this is an accepted assumption for either fluorescence or raman
        packet = photon("point")
        packet%nxp = 1.0_wp 
        packet%nyp = 0.0_wp 
        packet%nzp = 0.0_wp 
        
        ! loop through every position in the grid array

        !use symmetry to reduce the number of voxels we have to check

        !if there is 360rotational symmetry then we only need to consider going outwards in r and down in depth, 
        !   you will need some method to map a cartesian grid onto a cylindrical and visa versa, use bilinear interpolation
        !   this will require some work
        !   symmetry should be around the centre axis of the detector(s) to work properly
        !   first work along a know good axis that is aligned with the grid, then perform the interpolation to get along a know set   
        !
        !   
        !
        !if there is 2rotational symmetry then we have to consider half of the shape divided by a plane, remember to rotate the plane
        ! leave in cartesian to not have to perform bilinear interpolation
        !if there is flipped symmetry then we have to consider half of the shape divided by a plane
        !if there is prism symmetry then we have to consider a layer going in width and depth and then use this for each layer in the length
        !if there is none then we will have to consider all squares, which will take a very long time

        !add some warning to the tell the user that this will take a while to run

        do m = 1, state%grid%nxg
            do n = 1, state%grid%nyg
                do o = 1, state%grid%nzg

                    ! reset the arrays storing data
                    call reset(dects)

                    ! find the centre position of the voxel
                    y = (((real(n, kind = wp) - 0.5)/state%grid%nyg)*2.0_wp*state%grid%ymax) - state%grid%ymax
                    x = (((real(m, kind = wp) - 0.5)/state%grid%nxg)*2.0_wp*state%grid%xmax) - state%grid%xmax
                    z = (((real(o, kind = wp) - 0.5)/state%grid%nzg)*2.0_wp*state%grid%zmax) - state%grid%zmax
                    position = vector(x,y,z)

                    packet%pos = position   ! set the emission location

                    ! get the layer at this position
                    distances = 0._wp
                    do i = 1, size(distances)
                        distances(i) = array(i)%evaluate(position)
                    end do
                    layer=maxloc(distances,dim=1, mask=(distances<0._wp))

                    ! if the layer has a non-zero kappa then it is significant and we want to perform MCRT
                    if(array(layer)%getkappa() /= real(0, kind=wp)) then
                        call run_MCRT(input_file, history, packet, dict, & 
                                        distances, image, dects, array, nscatt, start, & 
                                        tev, spectrum)
                    end if
                    
                    ! record the efficiency for each detector and add to an array of escape functions


                end do         
            end do
        end do     

        !store the escape funcitons for each detector

    end subroutine escape_Function




    subroutine run_MCRT(input_file, history, packet, dict, & 
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
        !$omp& shared(dict, array, numproc, start, bar, jmean, emission, absorb, input_file, phasor, tev, dects, spectrum)& 
        !$omp& private(id, distances, image, history, seqs)& 
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
    end subroutine run_MCRT

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
    end subroutine survivalBiasPropagation

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

subroutine reset(dects)
    !! routine resets the values of data tracking arrays

    use constants, only: wp
    use setupMod, only : zarray
    use detectors

    type(dect_array), intent(inout) :: dects(:)

    integer :: i


    !zero jmean, absorb, emission, and phasor arrays
    call zarray()

    !zero the detectors
    do i = 1, size(dects)
        associate(x => dects(i)%p)
        call x%zero_dect()
        print*, "zeroed"
        end associate
    end do

end subroutine reset

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