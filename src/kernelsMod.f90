module kernels
!! Contains the routines to call other routines to setup and break down the simulation. Along with similar helper tools

    implicit none
    
    private
    public :: test_kernel
    public :: setup, finalise, reset_detectors

contains
!###############################################################################
!                   KERNELS

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
            packet%layer=(maxloc(distances,dim=1, mask=(distances<=0._wp)))

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


        call normalise_fluence(state%grid, jmeanGLOBAL, state%nphotons)
        call write_data(jmeanGLOBAL, trim(fileplace)//"jmean/"//state%outfile, state, dict)

        call normalise_fluence(state%grid, emissionGLOBAL, state%nphotons)
        call write_data(emissionGLOBAL, trim(fileplace)//"emission/"//state%rendersourcefile, state, dict)

        call normalise_fluence(state%grid, absorbGLOBAL, state%nphotons)
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

subroutine reset_detectors(dects)
    !! routine resets the values of data tracking arrays

    use constants, only: wp
    use setupMod, only : zarray
    use detectors

    type(dect_array), intent(inout) :: dects(:)

    integer :: i


    !zero jmean, absorb, emission, and phasor arrays
    ! this takes a long time to run, can it be sped up?
    !call zarray()

    !zero the detectors
    do i = 1, size(dects)
        call dects(i)%p%zero_dect()
    end do

end subroutine reset_detectors

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

    integer :: width

    width = max(len(trim(input_file)), len(kernel_type), len(trim(state%source)), len(trim(state%experiment)), 50)
    if (mod(width,2) == 1) then
        width = width + 21
    else
        width = width + 20
    end if


    print*,repeat("#", ((width-10)/2))//" Settings "//repeat("#", ((width-10)/2))
    print*,"# Config file: ",trim(input_file),repeat(" ", width-16-len(trim(input_file))),"#"
    print*,"# Using: "//trim(kernel_type)//"kernel"//repeat(" ", width-16-len(kernel_type)),"#"
    print*,"# Light source: "//trim(state%source)//repeat(" ", width-17-len(trim(state%source))),"#"
    if(state%source == "point")then
        print*,"# Light Source Position: ["//str(packet%pos%x,4)//", "//str(packet%pos%y,4)//", "//str(packet%pos%z,4)// &
                                        "]"//repeat(" ", width-44)//"#"
    else
        print*,"# Light direction: ["//str(packet%nxp,4)//", "//str(packet%nyp,4)//", "//str(packet%nzp,4)// &
                                  "]"//repeat(" ", width-38)//"#"
    end if
    print*,"# Geometry: "//trim(state%experiment)//repeat(" ", width-13-len(trim(state%experiment))),"#"
    print*,"# Seed: "//str(state%iseed,9)//repeat(" ", width-18)//"#"
    if(state%tev)then
        print*,"# Tev enabled!"//repeat(" ", width-15)//"#"
    end if
    if(state%render_geom)then
        print*,"# Render geometry to file enabled!"//repeat(" ", width-35)//"#"
    end if
    if(state%overwrite)then
        print*,"# Overwrite Enabled!",repeat(" ", width-21)//"#"
    end if
    if(state%absorb)then
        print*,"# Energy absorbed will be written to file."//repeat(" ", width-43)//"#"
    end if
    print*,repeat("#", width)
    print*,new_line("a")

end subroutine display_settings
end module kernels
