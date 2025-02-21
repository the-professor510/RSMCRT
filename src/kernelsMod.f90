module kernels
!! Contains the main program and scattering loop. Calls all other routine to setup, run and break down the simulation.

    implicit none
    
    private
    public :: default_MCRT, escape_Function, inverse_MCRT, test_kernel

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
        use detectors
        use historyStack,  only : history_stack_t
        use inttau2,       only : tauint2
        use photonMod,     only : photon
        use piecewiseMod
        use random,        only : ran2, init_rng
        use sdfs,          only : sdf
        use sdfHelpers,    only : rotationAlign, rotmat
        use sim_state_mod, only : state
        use vector_class
        use setupMod, only : setup_escapeFunction
        use writer_mod, only : write_escape

        !external deps
        use tev_mod, only : tevipc
        use tomlf,   only : toml_table, toml_error, get_value
#ifdef _OPENMP
        use omp_lib
#endif
        character(len=*), intent(in) :: input_file
        
        integer                       :: j, i, k
        type(history_stack_t)         :: history
        type(photon)                  :: packet
        type(toml_table)              :: dict
        real(kind=wp),    allocatable :: distances(:), image(:,:,:)
        type(dect_array), allocatable :: dects(:)
        type(sdf),        allocatable :: array(:)
        real(kind=wp)                 :: nscatt, start
        type(spectrum_t)              :: spectrum
        type(tevipc)                  :: tev
        type(toml_error), allocatable :: error

        integer :: nphotons_run,pos
        character(len=128) :: line
        character(len=:), allocatable :: checkpt_input_file
         
        integer :: m, n, o, layer
        real(kind = wp) :: x,y,z, total
        type(vector) :: position, direction, gridPos
        real(kind=wp) :: rotationOnToSym(4,4), rotationOffSym(4,4)
        real(kind=wp) :: rotationAroundZOnSym(4,4), rotationAroundZOffSym(4,4)
        character(len=:), allocatable :: symmetryType
        integer :: indices(3)
        real :: tic, toc

        call cpu_time(tic)

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

        !set the packet to be a isotropic source, this is an accepted assumption for either fluorescence or raman
        packet = photon("point")
        packet%nxp = 1.0_wp 
        packet%nyp = 0.0_wp 
        packet%nzp = 0.0_wp 

        ! Symmetries to implemented
        ! none DONE
        ! prism DONE
        ! flipped DONE
        ! uniformSlab DONE
        ! none cylindircal DONE
        ! 360rotational DONE

        !symmetries to implement at some point
        ! specified
        ! specified cylindrical?
        
        call setup_escapeFunction(size(dects))
        
        !use symmetry to reduce the number of voxels we have to calculate
        call get_value(dict, "symmetryType", symmetryType)
        select case(symmetryType)
        case("none")
            !there is no symmetry launch from every cell

            print*, "No Symmetry selected"
            print*, "User warning! This may take a long time to run"
            print*, "It is advised to try and find a geometry with symmetry or to reduce the grid size"
            print*, "Number of Monte Carlo Simmulations to run: ", (state%symmetryEscapeCartGrid%nxg* &
                                                                    state%symmetryEscapeCartGrid%nyg* & 
                                                                    state%symmetryEscapeCartGrid%nzg)
            print*, ""

            !allocate the escape symmetry grids
            allocate(escapeSymmetry(size(dects), state%symmetryEscapeCartGrid%nxg, & 
                                    state%symmetryEscapeCartGrid%nyg, & 
                                    state%symmetryEscapeCartGrid%nzg))
            escapeSymmetry = 0._wp

            !precompute the rotation vector here
            !both for going from the shifted from base
            ! and for going from base to the shifted
            direction = vector(0.0_wp, 0.0_wp, 1.0_wp)

            rotationOffSym = rotationAlign(direction, state%symGridDir)
            rotationOnToSym = rotationAlign(state%symGridDir, direction)

            rotationAroundZOffSym = rotmat(direction, -state%symGridRot)
            rotationAroundZOnSym = rotmat(direction, state%symGridRot)

            gridPos = state%symGridPos

            !loop through every cell
            do m = 1, state%symmetryEscapeCartGrid%nxg
                do n = 1, state%symmetryEscapeCartGrid%nyg
                    do o = 1, state%symmetryEscapeCartGrid%nzg

                        print*, ""
                        print*, "Running ", ((m-1)*state%symmetryEscapeCartGrid%nyg*state%symmetryEscapeCartGrid%nzg + & 
                                             (n-1)*state%symmetryEscapeCartGrid%nzg + o - 1), & 
                                " out of ", (state%symmetryEscapeCartGrid%nxg* &
                                                state%symmetryEscapeCartGrid%nyg* & 
                                                state%symmetryEscapeCartGrid%nzg)

                        !calculate the escape function
                        call cart_calc_escape_sym(m,n,o, rotationAroundZOffSym, rotationOffSym, gridPos, dects, array,& 
                                                 packet, distances, dict, history, image, input_file, nscatt, spectrum,& 
                                                 start, tev)

                    end do         
                end do
            end do

            !Go through the base grid and use some form of interpolation to figure out the best match
            call cart_map_escape_sym(dects, rotationOnToSym, rotationAroundZOnSym, gridPos)

        case("prism")
            !prism symmetry, launch from a layer of cells

            print*, "Prism symmetry selected"
            print*, "Number of Monte Carlo Simmulations to run: ", (state%symmetryEscapeCartGrid%nxg* &
                                                                    state%symmetryEscapeCartGrid%nyg)
            print*, ""

            !allocate the escape symmetry grids
            allocate(escapeSymmetry(size(dects), state%symmetryEscapeCartGrid%nxg, & 
                                    state%symmetryEscapeCartGrid%nyg, & 
                                    state%symmetryEscapeCartGrid%nzg))
            escapeSymmetry = 0._wp

            !precompute the rotation vector here
            !both for going from the shifted from base
            ! and for going from base to the shifted
            direction = vector(0.0_wp, 0.0_wp, 1.0_wp)

            !state%symGridDir is the normal to plane of the prism

            rotationOffSym = rotationAlign(direction, state%symGridDir)
            rotationOnToSym = rotationAlign(state%symGridDir, direction)

            rotationAroundZOffSym = rotmat(direction, -state%symGridRot)
            rotationAroundZOnSym = rotmat(direction, state%symGridRot)

            gridPos = state%symGridPos

            !get the position of the cell
            indices = state%symmetryEscapeCartGrid%get_voxel(vector(0.0_wp,0.0_wp,0.0_wp))
            !loop through every cell
            do m = 1, state%symmetryEscapeCartGrid%nxg
                do n = 1, state%symmetryEscapeCartGrid%nyg

                    print*, ""
                    print*, "Running ", ((m-1)*state%symmetryEscapeCartGrid%nyg + n - 1), & 
                            " out of ", (state%symmetryEscapeCartGrid%nxg* &
                                        state%symmetryEscapeCartGrid%nyg)

                    !calculate the escape function
                    call cart_calc_escape_sym(m,n,indices(3), rotationAroundZOffSym, rotationOffSym, gridPos, dects, array,& 
                                                packet, distances, dict, history, image, input_file, nscatt, spectrum,& 
                                                start, tev)
                end do
            end do

            !fill the rest of the grid
            do o = 1, state%symmetryEscapeCartGrid%nzg
                escapeSymmetry(:, :, :, o) = escapeSymmetry(:, :, :, indices(3))
            end do

            !Go through the base grid and use some form of interpolation to figure out the best match
            call cart_map_escape_sym(dects, rotationOnToSym, rotationAroundZOnSym, gridPos)

        case("flipped")
            !flipped symmetry, launch half the cells

            print*, "Flipped symmetry selected"
            print*, "Number of Monte Carlo Simmulations to run: ", (state%symmetryEscapeCartGrid%nxg* &
                                                                    state%symmetryEscapeCartGrid%nyg* &
                                                                    (state%symmetryEscapeCartGrid%nzg/2)+1)
            print*, ""
            
            !allocate the escape symmetry grids
            allocate(escapeSymmetry(size(dects), state%symmetryEscapeCartGrid%nxg, & 
                                    state%symmetryEscapeCartGrid%nyg, & 
                                    state%symmetryEscapeCartGrid%nzg))
            escapeSymmetry = 0._wp

            !precompute the rotation vector here
            !both for going from the shifted from base
            ! and for going from base to the shifted
            direction = vector(0.0_wp, 0.0_wp, 1.0_wp)

            !state%symGridDir is the normal pointing off the face to be flipped on

            rotationOffSym = rotationAlign(direction, state%symGridDir)
            rotationOnToSym = rotationAlign(state%symGridDir, direction)

            rotationAroundZOffSym = rotmat(direction, -state%symGridRot)
            rotationAroundZOnSym = rotmat(direction, state%symGridRot)

            gridPos = state%symGridPos

            !get the position of the cell
            indices = state%symmetryEscapeCartGrid%get_voxel(vector(0.0_wp,0.0_wp,0.0_wp))
            !loop through every cell
            do m = 1, state%symmetryEscapeCartGrid%nxg
                do n = 1, state%symmetryEscapeCartGrid%nyg
                    do o = 1, (state%symmetryEscapeCartGrid%nzg/2)+1

                        print*, ""
                        print*, "Running ", ((m-1)*state%symmetryEscapeCartGrid%nyg*((state%symmetryEscapeCartGrid%nzg/2)+1) + & 
                                            (n-1)*((state%symmetryEscapeCartGrid%nzg/2)+1) + o - 1), & 
                                " out of ", (state%symmetryEscapeCartGrid%nxg* &
                                            state%symmetryEscapeCartGrid%nyg* &
                                            (state%symmetryEscapeCartGrid%nzg/2)+1)

                        call cart_calc_escape_sym(m,n,o, rotationAroundZOffSym, rotationOffSym, gridPos, dects, array,& 
                                                    packet, distances, dict, history, image, input_file, nscatt, spectrum,& 
                                                    start, tev)
                    end do
                end do
            end do

            !fill the rest of the grid
            do m = 1, state%symmetryEscapeCartGrid%nxg
                do n = 1, state%symmetryEscapeCartGrid%nyg
                    do o = 1, (state%symmetryEscapeCartGrid%nzg/2)+1
                        escapeSymmetry(:, m, n, state%symmetryEscapeCartGrid%nzg - o + 1) = escapeSymmetry(:, m, n, o)
                    end do
                end do
            end do

            !Go through the base grid and use some form of interpolation to figure out the best match
            call cart_map_escape_sym(dects, rotationOnToSym, rotationAroundZOnSym, gridPos)

        case("uniformSlab")
            ! The simmulation is a slab code, light is collected uniformly

            print*, "Uniform slab symmetry selected"           
            print*, "Number of Monte Carlo Simmulations to run: ", (state%symmetryEscapeCartGrid%nzg)
            print*, ""
            
            !allocate the escape symmetry grids
            allocate(escapeSymmetry(size(dects), state%symmetryEscapeCartGrid%nxg, & 
                                    state%symmetryEscapeCartGrid%nyg, & 
                                    state%symmetryEscapeCartGrid%nzg))
            escapeSymmetry = 0._wp

            !precompute the rotation vector here
            !both for going from the shifted from base
            ! and for going from base to the shifted
            direction = vector(0.0_wp, 0.0_wp, 1.0_wp)

            !state%symGridDir is the normal pointing off the face to be flipped on

            rotationOffSym = rotationAlign(direction, state%symGridDir)
            rotationOnToSym = rotationAlign(state%symGridDir, direction)

            rotationAroundZOffSym = rotmat(direction, -state%symGridRot)
            rotationAroundZOnSym = rotmat(direction, state%symGridRot)

            gridPos = state%symGridPos

            !get the position of the cell
            indices = state%symmetryEscapeCartGrid%get_voxel(vector(0.0_wp,0.0_wp,0.0_wp))
            !loop through every cell
            do o = 1, state%symmetryEscapeCartGrid%nzg

                print*, ""
                print*, "Running ", (o - 1), & 
                        " out of ", (state%symmetryEscapeCartGrid%nzg)

                call cart_calc_escape_sym(indices(1),indices(2),o, rotationAroundZOffSym, rotationOffSym, gridPos, dects, array,& 
                                            packet, distances, dict, history, image, input_file, nscatt, spectrum,& 
                                            start, tev)
            end do

            !fill the rest of the grid
            do m = 1, state%symmetryEscapeCartGrid%nxg
                do n = 1, state%symmetryEscapeCartGrid%nyg
                    escapeSymmetry(:, m, n, :) = escapeSymmetry(:, indices(1), indices(2), :)
                end do
            end do

            !Go through the base grid and use some form of interpolation to figure out the best match
            call cart_map_escape_sym(dects, rotationOnToSym, rotationAroundZOnSym, gridPos)

        case("noneRotational")
            ! Do for all radii, theta and z values

            print*, "No Symmetry in cylindrical coordinates selected"
            print*, "User warning! This may take a long time to run"
            print*, "It is advised to try and find a geometry with symmetry or to reduce the grid size"
            print*, "Number of Monte Carlo Simmulations to run: ", (state%symmetryEscapeCylGrid%nrg* &
                                                                    state%symmetryEscapeCylGrid%ntg* & 
                                                                    state%symmetryEscapeCylGrid%nzg)
            print*, ""

            allocate(escapeSymmetry(size(dects), state%symmetryEscapeCylGrid%nrg, & 
                                    state%symmetryEscapeCylGrid%ntg, & 
                                    state%symmetryEscapeCylGrid%nzg))
            escapeSymmetry = 0._wp

            !precompute the rotation vector here
            !both for going from the shifted from base
            ! and for going from base to the shifted
            direction = vector(0.0_wp, 0.0_wp, 1.0_wp)

            rotationOffSym = rotationAlign(direction, state%symGridDir)
            rotationOnToSym = rotationAlign(state%symGridDir, direction)

            rotationAroundZOffSym = rotmat(direction, -state%symGridRot)
            rotationAroundZOnSym = rotmat(direction, state%symGridRot)

            gridPos = state%symGridPos

            !loop through every cell
            do m = 1, state%symmetryEscapeCylGrid%nrg
                do n = 1, state%symmetryEscapeCylGrid%ntg
                    do o = 1, state%symmetryEscapeCylGrid%nzg

                        print*, ""
                        print*, "Running ", ((m-1)*state%symmetryEscapeCylGrid%ntg*state%symmetryEscapeCylGrid%nzg + & 
                                             (n-1)*state%symmetryEscapeCylGrid%nzg + o - 1), & 
                                " out of ", (state%symmetryEscapeCylGrid%nrg* &
                                            state%symmetryEscapeCylGrid%ntg* & 
                                            state%symmetryEscapeCylGrid%nzg)

                        !calculate the escape function
                        call cyl_calc_escape_sym(m,n,o, rotationAroundZOffSym, rotationOffSym, gridPos, dects, array,& 
                                                 packet, distances, dict, history, image, input_file, nscatt, spectrum,& 
                                                 start, tev)

                    end do         
                end do
            end do

            !Go through the base grid and use some form of interpolation to figure out the best match
            call cyl_map_escape_sym(dects, rotationOnToSym, rotationAroundZOnSym, gridPos)

        case("360rotational")
            ! Do for all radii and z values at one theta value

            print*, "360 Rotational symmetry in cylindrical coordinates selected"
            print*, "User warning! This may take a long time to run"
            print*, "It is advised to try and find a geometry with symmetry or to reduce the grid size"
            print*, "Number of Monte Carlo Simmulations to run: ", (state%symmetryEscapeCylGrid%nrg* & 
                                                                    state%symmetryEscapeCylGrid%nzg)
            print*, ""

            allocate(escapeSymmetry(size(dects), state%symmetryEscapeCylGrid%nrg, & 
                                    state%symmetryEscapeCylGrid%ntg, & 
                                    state%symmetryEscapeCylGrid%nzg))
            escapeSymmetry = 0._wp

            !precompute the rotation vector here
            !both for going from the shifted from base
            ! and for going from base to the shifted
            direction = vector(0.0_wp, 0.0_wp, 1.0_wp)

            rotationOffSym = rotationAlign(direction, state%symGridDir)
            rotationOnToSym = rotationAlign(state%symGridDir, direction)

            rotationAroundZOffSym = rotmat(direction, -state%symGridRot)
            rotationAroundZOnSym = rotmat(direction, state%symGridRot)

            gridPos = state%symGridPos

            n=1
            do m = 1, state%symmetryEscapeCylGrid%nrg
                do o = 1, state%symmetryEscapeCylGrid%nzg

                    !calculate the escape function
                    print*, ""
                    print*, "Running ", ((m-1)*state%symmetryEscapeCylGrid%nzg + o - 1), & 
                            " out of ", (state%symmetryEscapeCylGrid%nrg*state%symmetryEscapeCylGrid%nzg)
                    call cyl_calc_escape_sym(m,n,o, rotationAroundZOffSym, rotationOffSym, gridPos, dects, array,& 
                                                packet, distances, dict, history, image, input_file, nscatt, spectrum,& 
                                                start, tev)

                end do 
            end do

            !loop through every cell
            
            do n = 1, state%symmetryEscapeCylGrid%ntg
                escapeSymmetry(:,:,n,:) = escapeSymmetry(:,:,1,:)
            end do

            !Go through the base grid and use some form of interpolation to figure out the best match
            call cyl_map_escape_sym(dects, rotationOnToSym, rotationAroundZOnSym, gridPos)
        case default                     
            print*,"Unknown symmetry type"
            stop 1
        end select
                                            
        !store the escape funcitons for each detector
        call write_escape(dects, dict)

        call finalise(dict, dects, nscatt, start, history)

        call cpu_time(toc)
        print*,"Time to Run: ",((toc - tic))

    end subroutine escape_Function


    subroutine cart_calc_escape_sym(m,n,o, rotationAroundZOffSym, rotationOffSym, gridPos, dects, array, packet, & 
                                     distances, dict, history, image, input_file, nscatt, spectrum, start, tev)

        !Calculate the cartesian symmetry escape function 
        use constants, only : wp
        use iarray
        use detectors
        use sdfs,          only : sdf
        use sdfHelpers,    only : rotationAlign, rotmat
        use sim_state_mod
        use vector_class
        use photonMod,     only : photon, set_photon
        use historyStack,  only : history_stack_t
        use piecewiseMod

        !external deps
        use tev_mod, only : tevipc
        use tomlf,   only : toml_table, toml_error, get_value

        !> indices of symmetryEscapeCartGrid
        integer, intent(in) :: m,n,o
        !> rotation matrix to unrotate around the z axis
        real(kind=wp), intent(in) :: rotationAroundZOffSym(4,4)
        !> rotation matrix to unrotate around the z axis
        real(kind=wp), intent(in) :: rotationOffSym(4,4)
        !> rotation matrix to unrotate around the z axis
        type(vector), intent(in) :: gridPos
        character(len=*), intent(in) :: input_file
        type(history_stack_t)        , intent(inout) :: history
        type(photon)                 , intent(inout) :: packet
        type(toml_table)             , intent(inout) :: dict
        real(kind=wp),    allocatable, intent(inout) :: distances(:), image(:,:,:)
        type(dect_array), allocatable, intent(inout) :: dects(:)
        type(sdf),        allocatable, intent(inout) :: array(:)
        real(kind=wp)                , intent(inout) :: nscatt, start
        type(spectrum_t)             , intent(inout) :: spectrum
        type(tevipc)                 , intent(inout) :: tev


        integer :: loopCounter, layer
        real(kind= wp) :: x,y,z, total
        type(vector) :: position

        ! reset the arrays storing data
        call reset(dects)

        ! find the centre position of the voxel
        x = (((real(m, kind = wp) - 0.5)/state%symmetryEscapeCartGrid%nxg)*& 
            2.0_wp*state%symmetryEscapeCartGrid%xmax) - state%symmetryEscapeCartGrid%xmax 
        y = (((real(n, kind = wp) - 0.5)/state%symmetryEscapeCartGrid%nyg)*& 
            2.0_wp*state%symmetryEscapeCartGrid%ymax) - state%symmetryEscapeCartGrid%ymax 
        z = (((real(o, kind = wp) - 0.5)/state%symmetryEscapeCartGrid%nzg)*& 
            2.0_wp*state%symmetryEscapeCartGrid%zmax) - state%symmetryEscapeCartGrid%zmax 
        
        position = vector(x,y,z)

        !rotate to align x and y axis after z axis alignment
        position = position .dot. rotationAroundZOffSym

        !align z axis
        position = position .dot. rotationOffSym

        !shift
        position = position + gridPos

        !set the emissin location to the centre of the voxel
        call set_photon(position, vector(0.0_wp,0.0_wp,0.0_wp))
        packet%pos = position

        ! get the layer at this position
        distances = 0._wp
        do loopCounter = 1, size(distances)
            distances(loopCounter) = array(loopCounter)%evaluate(position)
        end do
        layer=(maxloc(distances,dim=1, mask=(distances<0._wp)))

        !is this point inside the defined geometry
        if (layer == 0) then
            escapeSymmetry(loopCounter, m, n, o) = 0.0_wp
            return
        end if

        ! if the layer has a non-zero kappa then it is significant and we want to perform MCRT
        if(array(layer)%getkappa() /= real(0, kind=wp)) then
            call run_MCRT(input_file, history, packet, dict, & 
                            distances, image, dects, array, nscatt, start, & 
                            tev, spectrum)
        end if
        
        ! record the efficiency for each detector and add to an array of escape functions
        do loopCounter = 1, size(dects)
            if(array(layer)%getkappa() /= real(0, kind=wp)) then
                total = 0._wp
                call dects(loopCounter)%p%total_dect(total)
                escapeSymmetry(loopCounter, m, n, o) = total/state%nphotons


                !temporary while testing 
                !escapeSymmetry(loopCounter, m, n, o) = layer
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            else
                escapeSymmetry(loopCounter, m, n, o) = 0.0_wp

                !temporary while testing 
                !escapeSymmetry(loopCounter, m, n, o) = layer
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            end if
        end do

    end subroutine cart_calc_escape_sym

    subroutine cart_map_escape_sym(dects, rotationOnToSym, rotationAroundZOnSym, gridPos)

        use iarray
        use constants, only : wp
        use sim_state_mod, only : state
        use vector_class, only : vector
        use detectors
        use interpolate

        type(dect_array), allocatable, intent(inout) :: dects(:)
        real(kind=wp), intent(in) :: rotationOnToSym(4,4), rotationAroundZOnSym(4,4)
        !> rotation matrix to unrotate around the z axis
        type(vector), intent(in) :: gridPos

        integer :: i,j,k, m,n,o
        integer :: loopCounter, indx(3)
        integer :: xIndx(2), yIndx(2),zIndx(2)
        real(kind=wp) :: x,y,z, closestX, closestY, closestZ
        logical :: notOnXedge, notOnYedge, notOnZedge
        real(kind=wp) :: corners3D(2,2,2,4), corners2D(2,2,3), corners1D(2,2), point3D(4), point2D(3), point1D(2)
        type(vector) :: position
        

        !Go through the base grid and use some form of interpolation to figure out the best match
        do m = 1, state%grid%nxg
            do n = 1, state%grid%nyg
                do o = 1, state%grid%nzg

                    !TO DO
                    !add some function to do the part in the this for both cart and cyl

                    ! reset the arrays storing data
                    call reset(dects)

                    ! find the centre position of the voxel
                    y = (((real(n, kind = wp) - 0.5)/state%grid%nyg)*2.0_wp*state%grid%ymax) - state%grid%ymax
                    x = (((real(m, kind = wp) - 0.5)/state%grid%nxg)*2.0_wp*state%grid%xmax) - state%grid%xmax
                    z = (((real(o, kind = wp) - 0.5)/state%grid%nzg)*2.0_wp*state%grid%zmax) - state%grid%zmax

                    position = vector(x,y,z) 

                    !shift
                    position = position - gridPos

                    !rotate, there is none for this geometry
                    position = position .dot. rotationOnToSym

                    !rotate to align x and y axis after z axis alignment
                    position = position .dot. rotationAroundZOnSym
                    
                    !find the points in symmetry escape that correspond to this point?
                    !this returns the point that is closest to this
                    indx = -1
                    indx = state%symmetryEscapeCartGrid%get_voxel(position)

                    !we are out of the bounds of the escapesymmetry grid
                    if (indx(1) == -1 .or. indx(2) == -1 .or. indx(3) == -1) then
                        do loopCounter = 1, size(dects)
                            escape(loopCounter, m, n, o) = -1._wp
                        end do
                        cycle
                    end if

                    !what is the position of this square that we are in
                    closestX = (((real(indx(1), kind = wp) - 0.5)/state%symmetryEscapeCartGrid%nxg)*& 
                        2.0_wp*state%symmetryEscapeCartGrid%xmax) - state%symmetryEscapeCartGrid%xmax 
                    closestY = (((real(indx(2), kind = wp) - 0.5)/state%symmetryEscapeCartGrid%nyg)*& 
                        2.0_wp*state%symmetryEscapeCartGrid%ymax) - state%symmetryEscapeCartGrid%ymax 
                    closestZ = (((real(indx(3), kind = wp) - 0.5)/state%symmetryEscapeCartGrid%nzg)*& 
                        2.0_wp*state%symmetryEscapeCartGrid%zmax) - state%symmetryEscapeCartGrid%zmax
                    
                    !What are the indices of the closest symmetry grid cells 
                    if (closestX > position%x) then
                        xIndx(1) = indx(1) - 1
                        xIndx(2) = indx(1)
                    else 
                        xIndx(1) = indx(1)
                        xIndx(2) = indx(1) + 1
                    end if
                    if (closestY > position%y) then
                        yIndx(1) = indx(2) - 1
                        yIndx(2) = indx(2)
                    else 
                        yIndx(1) = indx(2)
                        yIndx(2) = indx(2) + 1
                    end if
                    if (closestZ > position%z) then
                        zIndx(1) = indx(3) - 1
                        zIndx(2) = indx(3)
                    else 
                        zIndx(1) = indx(3)
                        zIndx(2) = indx(3) + 1
                    end if

                    !is the current position on the edge of the symmetry grid
                    if (xIndx(1) < 1 .or. xIndx(2) > state%symmetryEscapeCartGrid%nxg) then
                        notOnXedge = .false.
                    else 
                        notOnXedge = .true.
                    end if
                    if (yIndx(1) < 1 .or. yIndx(2) > state%symmetryEscapeCartGrid%nyg) then
                        notOnYedge = .false.
                    else 
                        notOnYedge = .true.
                    end if
                    if (zIndx(1) < 1 .or. zIndx(2) > state%symmetryEscapeCartGrid%nzg) then
                        notOnZedge = .false.
                    else 
                        notOnZedge = .true.
                    end if

                    !is the point fully enclosed by 8 points of the symmetry grid or is it on a face, edge or corner
                    if ((notOnXedge) .and. (notOnYedge) .and. (notOnZedge)) then
                        !we are not on an edge or corner perform trilinear interpolation
                        !corners3D = 0._wp

                        do loopCounter = 1, size(dects)
                        do i = 1,2
                            do j = 1,2
                                do k =1,2
                                    corners3D(i,j,k,1) = (((real(xIndx(i), kind = wp) - 0.5)/state%symmetryEscapeCartGrid%nxg)*& 
                                                2.0_wp*state%symmetryEscapeCartGrid%xmax) - state%symmetryEscapeCartGrid%xmax 
                                    corners3D(i,j,k,2) = (((real(yIndx(j), kind = wp) - 0.5)/state%symmetryEscapeCartGrid%nyg)*& 
                                                2.0_wp*state%symmetryEscapeCartGrid%ymax) - state%symmetryEscapeCartGrid%ymax
                                    corners3D(i,j,k,3) = (((real(zIndx(k), kind = wp) - 0.5)/state%symmetryEscapeCartGrid%nzg)*& 
                                                2.0_wp*state%symmetryEscapeCartGrid%zmax) - state%symmetryEscapeCartGrid%zmax 
                                    corners3D(i,j,k,4) = escapeSymmetry(loopCounter, xIndx(i), yIndx(j), zIndx(k))
                                end do
                            end do
                        end do
                        point3D(1) = position%x
                        point3D(2) = position%y
                        point3D(3) = position%z
                        point3D(4) = 0._wp

                        call trilinearInterpolate(corners3D, point3D)

                        escape(loopCounter, m, n, o) = point3D(4)
                        end do
                    else if ( (notOnXedge .and. notOnYedge) .and. .not.notOnZedge) then
                        !we are on the edge of z, perform bilinear interpolation

                        !get the zindx that is inside
                        if (zIndx(1) >= 1) then
                            k = 1
                        else 
                            k = 2
                        end if

                        do loopCounter = 1, size(dects)
                            do i = 1,2
                                do j = 1,2
                                    corners2D(i,j,1) = (((real(xIndx(i), kind = wp) - 0.5)/state%symmetryEscapeCartGrid%nxg)*& 
                                                2.0_wp*state%symmetryEscapeCartGrid%xmax) - state%symmetryEscapeCartGrid%xmax 
                                    corners2D(i,j,2) = (((real(yIndx(j), kind = wp) - 0.5)/state%symmetryEscapeCartGrid%nyg)*& 
                                                2.0_wp*state%symmetryEscapeCartGrid%ymax) - state%symmetryEscapeCartGrid%ymax
                                    corners2D(i,j,3) = escapeSymmetry(loopCounter, xIndx(i), yIndx(j), zIndx(k))
                                end do
                            end do
                            point2D(1) = position%x
                            point2D(2) = position%y
                            point2D(3) = 0._wp

                            call bilinearInterpolate(corners2D, point2D)

                            escape(loopCounter, m, n, o) = point2D(3)
                        end do
                    else if ( (notOnXedge .and. notOnZedge) .and. .not.notOnYedge) then
                        !we are on the edge of y, perform bilinear interpolation
                        !get the yindx that is inside
                        if (yIndx(1) >= 1) then
                            j = 1
                        else 
                            j = 2
                        end if

                        do loopCounter = 1, size(dects)
                            do i = 1,2
                                do k = 1,2
                                    corners2D(i,k,1) = (((real(xIndx(i), kind = wp) - 0.5)/state%symmetryEscapeCartGrid%nxg)*& 
                                                2.0_wp*state%symmetryEscapeCartGrid%xmax) - state%symmetryEscapeCartGrid%xmax 
                                    corners2D(i,k,2) = (((real(zIndx(k), kind = wp) - 0.5)/state%symmetryEscapeCartGrid%nzg)*& 
                                                2.0_wp*state%symmetryEscapeCartGrid%zmax) - state%symmetryEscapeCartGrid%zmax
                                    corners2D(i,k,3) = escapeSymmetry(loopCounter, xIndx(i), yIndx(j), zIndx(k))
                                end do
                            end do
                            point2D(1) = position%x
                            point2D(2) = position%z
                            point2D(3) = 0._wp

                            call bilinearInterpolate(corners2D, point2D)

                            escape(loopCounter, m, n, o) = point2D(3)
                        end do
                    else if ( (notOnYedge .and. notOnZedge) .and. .not.notOnXedge) then
                        !we are on the edge of x, perform bilinear interpolation
                        !get the xindx that is inside
                        if (xIndx(1) >= 1) then
                            i = 1
                        else 
                            i = 2
                        end if

                        do loopCounter = 1, size(dects)
                            do j = 1,2
                                do k = 1,2
                                    corners2D(j,k,1) = (((real(yIndx(j), kind = wp) - 0.5)/state%symmetryEscapeCartGrid%nyg)*& 
                                                2.0_wp*state%symmetryEscapeCartGrid%ymax) - state%symmetryEscapeCartGrid%ymax 
                                    corners2D(j,k,2) = (((real(zIndx(k), kind = wp) - 0.5)/state%symmetryEscapeCartGrid%nzg)*& 
                                                2.0_wp*state%symmetryEscapeCartGrid%zmax) - state%symmetryEscapeCartGrid%zmax
                                    corners2D(j,k,3) = escapeSymmetry(loopCounter, xIndx(i), yIndx(j), zIndx(k))
                                end do
                            end do
                            point2D(1) = position%y
                            point2D(2) = position%z
                            point2D(3) = 0._wp

                            call bilinearInterpolate(corners2D, point2D)

                            escape(loopCounter, m, n, o) = point2D(3)
                        end do
                    else if ( (notOnXedge) .and. (.not.notOnYedge) .and. (.not.notOnZedge) ) then
                        !we are on the edge of y and z, perform linear interpolation

                        !get the yindx that is inside
                        if (yIndx(1) >= 1) then
                            j = 1
                        else 
                            j = 2
                        end if
                        !get the zindx that is inside
                        if (zIndx(1) >= 1) then
                            k = 1
                        else 
                            k = 2
                        end if

                        do loopCounter = 1, size(dects)
                            do i = 1,2
                                    corners1D(i,1) = (((real(xIndx(i), kind = wp) - 0.5)/state%symmetryEscapeCartGrid%nxg)*& 
                                                2.0_wp*state%symmetryEscapeCartGrid%xmax) - state%symmetryEscapeCartGrid%xmax 
                                    corners1D(i,2) = escapeSymmetry(loopCounter, xIndx(i), yIndx(j), zIndx(k))
                            end do
                            point1D(1) = position%x
                            point1D(2) = 0._wp

                            call linearInterpolate(corners1D, point1D)

                            escape(loopCounter, m, n, o) = point1D(2)
                        end do
                    else if ( (notOnYedge) .and. (.not.notOnXedge) .and. (.not.notOnZedge) ) then
                        !we are on the edge of x and z, perform linear interpolation
                        !get the xindx that is inside
                        if (xIndx(1) >= 1) then
                            i = 1
                        else 
                            i = 2
                        end if
                        !get the zindx that is inside
                        if (zIndx(1) >= 1) then
                            k = 1
                        else 
                            k = 2
                        end if

                        do loopCounter = 1, size(dects)
                            do j = 1,2
                                    corners1D(j,1) = (((real(yIndx(j), kind = wp) - 0.5)/state%symmetryEscapeCartGrid%nyg)*& 
                                                2.0_wp*state%symmetryEscapeCartGrid%ymax) - state%symmetryEscapeCartGrid%ymax 
                                    corners1D(j,2) = escapeSymmetry(loopCounter, xIndx(i), yIndx(j), zIndx(k))
                            end do
                            point1D(1) = position%y
                            point1D(2) = 0._wp

                            call linearInterpolate(corners1D, point1D)

                            escape(loopCounter, m, n, o) = point1D(2)
                        end do
                    else if ( (notOnZedge) .and. (.not.notOnXedge) .and. (.not.notOnYedge) ) then
                        !we are on the edge of x and y, perform linear interpolation
                        !get the xindx that is inside
                        if (xIndx(1) >= 1) then
                            i = 1
                        else 
                            i = 2
                        end if
                        !get the yindx that is inside
                        if (yIndx(1) >= 1) then
                            j = 1
                        else 
                            j = 2
                        end if

                        do loopCounter = 1, size(dects)
                            do k = 1,2
                                    corners1D(k,1) = (((real(zIndx(k), kind = wp) - 0.5)/state%symmetryEscapeCartGrid%nzg)*& 
                                                2.0_wp*state%symmetryEscapeCartGrid%zmax) - state%symmetryEscapeCartGrid%zmax 
                                    corners1D(k,2) = escapeSymmetry(loopCounter, xIndx(i), yIndx(j), zIndx(k))
                            end do
                            point1D(1) = position%z
                            point1D(2) = 0._wp

                            call linearInterpolate(corners1D, point1D)

                            escape(loopCounter, m, n, o) = point1D(2)
                        end do
                    else
                        !we are on the edge of x, y, and z, set it to be equal to the closest value
                        escape(:, m, n, o) = escapeSymmetry(:, indx(1), indx(2), indx(3))
                    end if
                end do   
            end do
        end do
    end subroutine cart_map_escape_sym

    subroutine cyl_calc_escape_sym(m,n,o, rotationAroundZOffSym, rotationOffSym, gridPos, dects, array, packet, & 
                                    distances, dict, history, image, input_file, nscatt, spectrum, start, tev)

        !Calculate the cartesian symmetry escape function 
        use constants, only : wp
        use iarray
        use detectors
        use sdfs,          only : sdf
        use sdfHelpers,    only : rotationAlign, rotmat
        use sim_state_mod
        use vector_class
        use photonMod,     only : photon, set_photon
        use historyStack,  only : history_stack_t
        use piecewiseMod

        !external deps
        use tev_mod, only : tevipc
        use tomlf,   only : toml_table, toml_error, get_value

        !> indices of symmetryEscapeCartGrid
        integer, intent(in) :: m,n,o
        !> rotation matrix to unrotate around the z axis
        real(kind=wp), intent(in) :: rotationAroundZOffSym(4,4)
        !> rotation matrix to unrotate around the z axis
        real(kind=wp), intent(in) :: rotationOffSym(4,4)
        !> rotation matrix to unrotate around the z axis
        type(vector), intent(in) :: gridPos
        character(len=*), intent(in) :: input_file
        type(history_stack_t)        , intent(inout) :: history
        type(photon)                 , intent(inout) :: packet
        type(toml_table)             , intent(inout) :: dict
        real(kind=wp),    allocatable, intent(inout) :: distances(:), image(:,:,:)
        type(dect_array), allocatable, intent(inout) :: dects(:)
        type(sdf),        allocatable, intent(inout) :: array(:)
        real(kind=wp)                , intent(inout) :: nscatt, start
        type(spectrum_t)             , intent(inout) :: spectrum
        type(tevipc)                 , intent(inout) :: tev


        integer :: loopCounter, layer
        real(kind= wp) :: rad,theta,x,y,z, total
        type(vector) :: position

        ! reset the arrays storing data
        call reset(dects)

        ! find the centre position of the voxel in radians
        rad = ((real(m, kind = wp)-0.5)/state%symmetryEscapeCylGrid%nrg)*state%symmetryEscapeCylGrid%rmax
        theta = ((real(n, kind = wp)-0.5)/state%symmetryEscapeCylGrid%ntg)*state%symmetryEscapeCylGrid%tmax
        z = (((real(o, kind = wp) - 0.5)/state%symmetryEscapeCylGrid%nzg)*& 
            2.0_wp*state%symmetryEscapeCylGrid%zmax) - state%symmetryEscapeCylGrid%zmax 

        !convert rad and theta into x and y
        x = rad * cos(theta)
        y = rad * sin(theta)

        !translate back to the main grid
        position = vector(x,y,z)

        !rotate to align x and y axis after z axis alignment
        position = position .dot. rotationAroundZOffSym

        !align z axis
        position = position .dot. rotationOffSym

        !shift
        position = position + gridPos

        !set the emissin location to the centre of the voxel
        call set_photon(position, vector(0.0_wp,0.0_wp,0.0_wp))
        packet%pos = position

        ! get the layer at this position
        distances = 0._wp
        do loopCounter = 1, size(distances)
        distances(loopCounter) = array(loopCounter)%evaluate(position)
        end do
        layer=(maxloc(distances,dim=1, mask=(distances<0._wp)))

        !is this point inside the defined geometry
        if (layer == 0) then
            escapeSymmetry(loopCounter, m, n, o) = 0.0_wp
            return
        end if
        
        ! if the layer has a non-zero kappa then it is significant and we want to perform MCRT
        if(array(layer)%getkappa() /= real(0, kind=wp)) then
            call run_MCRT(input_file, history, packet, dict, & 
                            distances, image, dects, array, nscatt, start, & 
                            tev, spectrum)
        end if

        ! record the efficiency for each detector and add to an array of escape functions
        do loopCounter = 1, size(dects)
            if(array(layer)%getkappa() /= real(0, kind=wp)) then
                total = 0._wp
                call dects(loopCounter)%p%total_dect(total)
                escapeSymmetry(loopCounter, m, n, o) = total/state%nphotons


                !temporary while testing 
                !escapeSymmetry(loopCounter, m, n, o) = layer
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            else
                escapeSymmetry(loopCounter, m, n, o) = 0.0_wp

                !temporary while testing 
                !escapeSymmetry(loopCounter, m, n, o) = layer
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            end if
        end do

    end subroutine cyl_calc_escape_sym
    
    subroutine cyl_map_escape_sym(dects, rotationOnToSym, rotationAroundZOnSym, gridPos)

        use iarray
        use constants, only : wp, PI, TWOPI
        use sim_state_mod, only : state
        use vector_class, only : vector
        use detectors
        use interpolate

        type(dect_array), allocatable, intent(inout) :: dects(:)
        real(kind=wp), intent(in) :: rotationOnToSym(4,4), rotationAroundZOnSym(4,4)
        !> rotation matrix to unrotate around the z axis
        type(vector), intent(in) :: gridPos

        integer :: i,j,k, m,n,o
        integer :: loopCounter, indx(3)
        integer :: radIndx(2), thetaIndx(2),zIndx(2)
        real(kind=wp) :: x,y,z, closestRad, closestTheta, closestZ, rad, theta
        logical :: notOnXedge, notOnYedge, notOnZedge
        real(kind=wp) :: corners3D(2,2,2,4), corners2D(2,2,3), corners1D(2,2), point3D(4), point2D(3), point1D(2)
        real(kind=wp) :: averageRad(2)
        real(kind=wp) :: a1, a2, a3, at
        real(kind=wp) :: thetaLow, thetaHigh
        type(vector) :: position
        

        !Go through the base grid and use some form of interpolation to figure out the best match
        do m = 1, state%grid%nxg
            do n = 1, state%grid%nyg
                do o = 1, state%grid%nzg

                    !TO DO
                    !add some function to do the part in the this for both cart and cyl

                    ! reset the arrays storing data
                    call reset(dects)

                    ! find the centre position of the voxel
                    y = (((real(n, kind = wp) - 0.5)/state%grid%nyg)*2.0_wp*state%grid%ymax) - state%grid%ymax
                    x = (((real(m, kind = wp) - 0.5)/state%grid%nxg)*2.0_wp*state%grid%xmax) - state%grid%xmax
                    z = (((real(o, kind = wp) - 0.5)/state%grid%nzg)*2.0_wp*state%grid%zmax) - state%grid%zmax

                    position = vector(x,y,z) 

                    !shift
                    position = position - gridPos

                    !rotate, there is none for this geometry
                    position = position .dot. rotationOnToSym

                    !rotate to align x and y axis after z axis alignment
                    position = position .dot. rotationAroundZOnSym



                    !convert to rad and theta
                    rad = sqrt(position%x**2 + position%y**2)
                    if(rad == 0)then
                        theta=0.0
                    else
                        theta=atan2(position%y,position%x)
                        if(theta < 0.0)theta=theta+2*atan2(0.0d0,-1.0d0)
                    end if
                    
                    !find the points in symmetry escape that correspond to this point?
                    !this returns the point that is closest to this
                    indx = -1
                    indx = state%symmetryEscapeCylGrid%get_voxel(position)

                    !we are out of the bounds of the escapesymmetry grid
                    if (indx(1) == -1 .or. indx(2) == -1 .or. indx(3) == -1) then
                        do loopCounter = 1, size(dects)
                            escape(loopCounter, m, n, o) = -1._wp
                        end do
                        cycle
                    end if

                    !what is the position of this square that we are in
                    closestRad = ((real(indx(1), kind = wp)-0.5)/state%symmetryEscapeCylGrid%nrg)*state%symmetryEscapeCylGrid%rmax
                    closestTheta = ((real(indx(2), kind = wp)-0.5)/state%symmetryEscapeCylGrid%ntg)*state%symmetryEscapeCylGrid%tmax
                    closestZ = (((real(indx(3), kind = wp) - 0.5)/state%symmetryEscapeCylGrid%nzg)*& 
                        2.0_wp*state%symmetryEscapeCylGrid%zmax) - state%symmetryEscapeCylGrid%zmax 
                    
                    !What are the indices of the closest symmetry grid cells 
                    if (closestRad > rad) then
                        radIndx(1) = indx(1) - 1
                        radIndx(2) = indx(1)
                    else 
                        radIndx(1) = indx(1)
                        radIndx(2) = indx(1) + 1
                    end if
                    if (closestTheta > theta) then
                        thetaIndx(1) = indx(2) - 1
                        thetaIndx(2) = indx(2)
                    else 
                        thetaIndx(1) = indx(2)
                        thetaIndx(2) = indx(2) + 1
                    end if
                    if (closestZ > position%z) then
                        zIndx(1) = indx(3) - 1
                        zIndx(2) = indx(3)
                    else 
                        zIndx(1) = indx(3)
                        zIndx(2) = indx(3) + 1
                    end if

                    !special cases
                    !rIndex(1) = - 1, weighted average between the three points that bound it between the point r = 0 and 
                    !                  r equating to rIndx of 1 and the given theta bounds
                    !rIndex(2) greater than max, on edge of r, do some form of bilinear with theta and z

                    !zIndex(1) < 1 then do bilinear cylindrical 
                    !zIndex(2) > nzg then bilinear cylindrical
                    
                    !all other cases are contained within the abilility to use trilinear interpolation

                    !store the theta values that can be used in calculations
                    thetaLow = ((real(thetaIndx(1), kind = wp)-0.5)/state%symmetryEscapeCylGrid%ntg)*& 
                                state%symmetryEscapeCylGrid%tmax
                    thetaHigh = ((real(thetaIndx(2), kind = wp)-0.5)/state%symmetryEscapeCylGrid%ntg)*& 
                                state%symmetryEscapeCylGrid%tmax

                    if (thetaIndx(1) < 1) then 
                        ! wrap around
                        thetaIndx(1) = state%symmetryEscapeCylGrid%ntg
                    end if
                    if (thetaIndx(2) > state%symmetryEscapeCylGrid%ntg) then 
                        !wrap around
                        thetaIndx(2) = 1
                    end if

                    !print*, ""
                    !print*, ""
                    !print*, rad, theta, z
                    !print*, position
                    !print*, indx
                    !print*, closestRad, closestTheta, closestZ
                    !print*, radIndx, thetaIndx, zIndx
                    !print*, thetaLow, thetaHigh
                    
                    ! radIndx(1) < 1
                    if (radIndx(1) < 1) then 
                        !we need to find the proprotion of areas

                        !find the area a1, a2, a3 and use them as weightings
                        at = PI * (((0.5_wp)/state%symmetryEscapeCylGrid%nrg)*state%symmetryEscapeCylGrid%rmax)**2 * & 
                            ((thetaHigh-thetaLow)/TWOPI)
                        a1 = (0.5_wp * (((0.5_wp)/state%symmetryEscapeCylGrid%nrg)*state%symmetryEscapeCylGrid%rmax) * & 
                            rad * sin(thetaHigh - theta))
                        a2 = (0.5_wp * (((0.5_wp)/state%symmetryEscapeCylGrid%nrg)*state%symmetryEscapeCylGrid%rmax) * & 
                            rad * sin(theta-thetaLow))
                        a3 = (at - a1 - a2)

                        a1 = a1/at
                        a2 = a2/at
                        a3 = a3/at

                        if (zIndx(1) < 1) then 
                            !we are on the bottom z edge, take average r as your point
                            do loopCounter = 1, size(dects)
                                
                                averageRad(1) = 0._wp
                                do i = 1, state%symmetryEscapeCylGrid%ntg
                                    averageRad(1) = averageRad(1) + escapeSymmetry(loopCounter, 1, i, 1)
                                end do
                                averageRad(1) = averageRad(1)/state%symmetryEscapeCylGrid%ntg

                                escape(loopCounter, m, n, o) = a1 * escapeSymmetry(loopCounter, 1, thetaIndx(1), 1) + & 
                                                               a2 * escapeSymmetry(loopCounter, 1, thetaIndx(2), 1) + & 
                                                               a3 * averageRad(1)

                            end do
                            cycle
                        else if(zIndx(2) > state%symmetryEscapeCylGrid%nzg) then 
                            !we are on the top z edge, take average r as your point
                            do loopCounter = 1, size(dects)

                                averageRad(2) = 0._wp
                                do i = 1, state%symmetryEscapeCylGrid%ntg
                                    averageRad(2) = averageRad(2) + & 
                                                    escapeSymmetry(loopCounter, 1, i, state%symmetryEscapeCylGrid%nzg)
                                end do
                                averageRad(2) = averageRad(2)/state%symmetryEscapeCylGrid%ntg

                                escape(loopCounter, m, n, o) = a1 * escapeSymmetry(loopCounter, 1, thetaIndx(1), & 
                                                                                    state%symmetryEscapeCylGrid%nzg) + & 
                                                               a2 * escapeSymmetry(loopCounter, 1, thetaIndx(2), & 
                                                                                    state%symmetryEscapeCylGrid%nzg) + & 
                                                               a3 * averageRad(2)
                            end do
                            cycle
                        else 
                            !we are in the middle so we will lineraly interpolate between z bottom and z top
                            do loopCounter = 1, size(dects)
                                averageRad(1) = 0._wp
                                averageRad(2) = 0._wp
                                do i = 1, state%symmetryEscapeCylGrid%ntg
                                    averageRad(1) = averageRad(1) + escapeSymmetry(loopCounter, 1, i, zIndx(1))
                                    averageRad(2) = averageRad(2) + escapeSymmetry(loopCounter, 1, i, zIndx(2))
                                end do
                                averageRad(1) = averageRad(1)/state%symmetryEscapeCylGrid%ntg
                                averageRad(2) = averageRad(2)/state%symmetryEscapeCylGrid%ntg


                                do k = 1,2
                                    corners1D(k,1) = (((real(zIndx(k), kind = wp) - 0.5)/state%symmetryEscapeCylGrid%nzg)*& 
                                                2.0_wp*state%symmetryEscapeCylGrid%zmax) - state%symmetryEscapeCylGrid%zmax 
                                    corners1D(k,2) = a1 * escapeSymmetry(loopCounter, 1, thetaIndx(1), zIndx(k)) + & 
                                                     a2 * escapeSymmetry(loopCounter, 1, thetaIndx(2), zIndx(k)) + & 
                                                     a3 * averageRad(k)
                                end do

                                point1D(1) = position%z
                                point1D(2) = 0._wp

                                call linearInterpolate(corners1D, point1D)

                                escape(loopCounter, m, n, o) = point1D(2)
                            end do
                            cycle
                        end if
                    end if

                    !radIndx(2) > state%symmetryEscapeCylGrid%nzg, we are on the edge of the radius
                    if (radIndx(2) > state%symmetryEscapeCylGrid%nrg) then 
                        !we need to check if we are on the edge of z
                        !zIndx(1) < 1, we are on the bottom edge, linearlyInterp using theta as your points
                        if(zIndx(1) < 1) then
                            
                            do loopCounter = 1, size(dects)
                                
                                corners1D(1,1) = thetaLow
                                corners1D(1,2) = escapeSymmetry(loopCounter, state%symmetryEscapeCylGrid%nrg, thetaIndx(1), 1)

                                corners1D(2,1) = thetaHigh
                                corners1D(2,2) = escapeSymmetry(loopCounter, state%symmetryEscapeCylGrid%nrg, thetaIndx(2), 1)
                                

                                point1D(1) = theta
                                point1D(2) = 0._wp

                                call linearInterpolate(corners1D, point1D)

                                escape(loopCounter, m, n, o) = point1D(2)
                            end do
                            cycle

                        !zIndx(2) > state%symmetryEscapeCylGrid%nzg, we are on the top edge
                        else if(zIndx(2) > state%symmetryEscapeCylGrid%nzg) then
                            do loopCounter = 1, size(dects)
                                
                                corners1D(1,1) = thetaLow
                                corners1D(1,2) = escapeSymmetry(loopCounter, state%symmetryEscapeCylGrid%nrg, thetaIndx(1), & 
                                                                state%symmetryEscapeCylGrid%nzg)

                                corners1D(2,1) = thetaHigh
                                corners1D(2,2) = escapeSymmetry(loopCounter, state%symmetryEscapeCylGrid%nrg, thetaIndx(2), & 
                                                                state%symmetryEscapeCylGrid%nzg)
                                

                                point1D(1) = theta
                                point1D(2) = 0._wp

                                call linearInterpolate(corners1D, point1D)

                                escape(loopCounter, m, n, o) = point1D(2)
                            end do
                            cycle

                        !else bilinearly interpolate between theta values and z values
                        else 
                            do loopCounter = 1, size(dects)                               
                                do k = 1,2
                                    corners2D(1,k,1) = thetaLow
                                    corners2D(1,k,2) = (((real(zIndx(k), kind = wp) - 0.5)/state%symmetryEscapeCylGrid%nzg)*& 
                                                    2.0_wp*state%symmetryEscapeCylGrid%zmax) - state%symmetryEscapeCylGrid%zmax 
                                    corners2D(1,k,3) = escapeSymmetry(loopCounter, state%symmetryEscapeCylGrid%nrg, thetaIndx(1), & 
                                                                    zIndx(k))

                                    corners2D(2,k,1) = thetaHigh
                                    corners2D(2,k,2) = (((real(zIndx(k), kind = wp) - 0.5)/state%symmetryEscapeCylGrid%nzg)*& 
                                                    2.0_wp*state%symmetryEscapeCylGrid%zmax) - state%symmetryEscapeCylGrid%zmax 
                                    corners2D(2,k,3) = escapeSymmetry(loopCounter, state%symmetryEscapeCylGrid%nrg, thetaIndx(2), & 
                                                                    zIndx(k))
                                end do

                                point2D(1) = theta
                                point2D(2) = position%z
                                point2D(3) = 0._wp

                                call bilinearInterpolate(corners2D, point2D)

                                escape(loopCounter, m, n, o) = point2D(3)
                            end do
                            cycle
                        end if
                    end if

                    !zIndx(1) < 1, we are on the bottom edge but not on the edge of r
                    if(zIndx(1) < 1) then
                        do loopCounter = 1, size(dects)

                            do i =1,2
                                corners2D(i, 1, 1) = ((real(radIndx(i), kind = wp)-0.5)/state%symmetryEscapeCylGrid%nrg) & 
                                                        *state%symmetryEscapeCylGrid%rmax
                                corners2D(i, 1, 2) = thetaLow
                                corners2D(i, 1, 3) = escapeSymmetry(loopCounter, radIndx(i), thetaIndx(1), 1)

                                corners2D(i, 2, 1) = ((real(radIndx(i), kind = wp)-0.5)/state%symmetryEscapeCylGrid%nrg) & 
                                                        *state%symmetryEscapeCylGrid%rmax
                                corners2D(i, 2, 2) = thetaHigh
                                corners2D(i, 2, 3) = escapeSymmetry(loopCounter, radIndx(i), thetaIndx(2), 1)
                            end do

                            point2D(1) = rad
                            point2D(2) = theta
                            point2D(3) = 0._wp

                            call cylBilinearInterpolate(corners2D, point2D)

                            escape(loopCounter, m, n, o) = point2D(3)
                        end do
                        cycle
                    end if

                    !zIndx(2) > state%symmetryEscapeCylGrid%nzg, we are on the top edge but not on the edge of r
                    if(zIndx(2) > state%symmetryEscapeCylGrid%nzg) then
                        do loopCounter = 1, size(dects)

                            do i =1,2
                                corners2D(i, 1, 1) = ((real(radIndx(i), kind = wp)-0.5)/state%symmetryEscapeCylGrid%nrg) & 
                                                        *state%symmetryEscapeCylGrid%rmax
                                corners2D(i, 1, 2) = thetaLow
                                corners2D(i, 1, 3) = escapeSymmetry(loopCounter, radIndx(i), thetaIndx(1), & 
                                                                    state%symmetryEscapeCylGrid%nzg)

                                corners2D(i, 2, 1) = ((real(radIndx(i), kind = wp)-0.5)/state%symmetryEscapeCylGrid%nrg) & 
                                                        *state%symmetryEscapeCylGrid%rmax
                                corners2D(i, 2, 2) = thetaHigh
                                corners2D(i, 2, 3) = escapeSymmetry(loopCounter, radIndx(i), thetaIndx(2), & 
                                                                    state%symmetryEscapeCylGrid%nzg)
                            end do

                            point2D(1) = rad
                            point2D(2) = theta
                            point2D(3) = 0._wp

                            call cylBilinearInterpolate(corners2D, point2D)

                            escape(loopCounter, m, n, o) = point2D(3)
                        end do
                        cycle
                    end if

                    do loopCounter = 1, size(dects)
                        
                        do i =1,2
                            do k = 1,2
                                corners3D(i, 1, k, 1) = ((real(radIndx(i), kind = wp)-0.5)/state%symmetryEscapeCylGrid%nrg) & 
                                                        *state%symmetryEscapeCylGrid%rmax
                                corners3D(i, 1, k, 2) = thetaLow
                                corners3D(i, 1, k, 3) = (((real(zIndx(k), kind = wp) - 0.5)/state%symmetryEscapeCylGrid%nzg)*& 
                                                    2.0_wp*state%symmetryEscapeCylGrid%zmax) - state%symmetryEscapeCylGrid%zmax 
                                corners3D(i, 1, k, 4) = escapeSymmetry(loopCounter, radIndx(i), thetaIndx(1), zIndx(k))

                                corners3D(i, 2, k, 1) = ((real(radIndx(i), kind = wp)-0.5)/state%symmetryEscapeCylGrid%nrg) & 
                                                        *state%symmetryEscapeCylGrid%rmax
                                corners3D(i, 2, k, 2) = thetaHigh
                                corners3D(i, 2, k, 3) = (((real(zIndx(k), kind = wp) - 0.5)/state%symmetryEscapeCylGrid%nzg)*& 
                                                    2.0_wp*state%symmetryEscapeCylGrid%zmax) - state%symmetryEscapeCylGrid%zmax 
                                corners3D(i, 2, k, 4) = escapeSymmetry(loopCounter, radIndx(i), thetaIndx(2), zIndx(k))
                            end do
                        end do

                        point3D(1) = rad
                        point3D(2) = theta
                        point3D(3) = position%z
                        point3D(4) = 0._wp

                        call cylTrilinearInterpolate(corners3D, point3D)

                        escape(loopCounter, m, n, o) = point3D(4)
                    end do

                end do   
            end do
        end do
    end subroutine cyl_map_escape_sym

    subroutine inverse_MCRT(input_file)

        !Shared data
        use iarray
        use constants, only : wp, fileplace

        !subroutines
        use detectors,     only : dect_array
        use historyStack,  only : history_stack_t
        use photonMod,     only : photon
        use piecewiseMod
        use random,        only : ran2, init_rng
        use sdfs,          only : sdf
        use sim_state_mod, only : state
        use opticalProperties, only : opticalProp_t, mono

        !external deps
        use tev_mod, only : tevipc
        use tomlf,   only : toml_table, get_value
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

        integer :: maxNumSteps, layer
        real(kind = wp) :: maxStepSize, gradStepSize, accuracy
        logical :: findmua, findmus, findg, findn
        real(kind=wp) :: temp, error

        real(kind=wp), allocatable :: gradDescentData(:,:)
        real(kind=wp) :: mus, mua, hgg, n
        integer :: i, SDF_array_index


        integer :: NoVariablesToOptimize
        type(opticalProp_t) :: trialOptProp


        real(kind=wp) :: probability, alpha, k
        real(kind=wp) :: musboundupper, musboundlower, muaboundupper, muaboundlower
        real(kind=wp) :: gboundupper, gboundlower, nboundupper, nboundlower
        integer :: nb_samples, indexOfMinError, indexOfMaxRatio
        real(kind=wp) :: minError, maxRatio, ratioCounter
        real(kind=wp), allocatable :: ratios(:)

        real(kind=wp) :: ran !used to temporarily store a random number
        real(kind=wp) :: leftMin, rightMax !sides of the LIPO condition
        real(kind=wp) :: runningTotal


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

        !read in the data used in gardient Descent or AdaLIPO
        call get_value(dict, "maxStepSize", maxStepSize)
        call get_value(dict, "gradStepSize", gradStepSize)
        call get_value(dict, "accuracy", accuracy)
        call get_value(dict, "maxNumSteps", maxNumSteps)
        call get_value(dict, "Findmua", findmua)
        call get_value(dict, "Findmus", findmus)
        call get_value(dict, "Findg", findg)
        call get_value(dict, "Findn", findn)
        call get_value(dict, "inverseLayer", layer)

        !check how many variables we will optimize/search for
        NoVariablesToOptimize = 0
        if (findmua) then
            NoVariablesToOptimize = NoVariablesToOptimize + 1
        end if
        if (findmus) then
            NoVariablesToOptimize = NoVariablesToOptimize + 1
        end if
        if (findg) then
            NoVariablesToOptimize = NoVariablesToOptimize + 1
        end if
        if (findn) then
            NoVariablesToOptimize = NoVariablesToOptimize + 1
        end if
        if(NoVariablesToOptimize == 0) then
            print*, "Please select at least one of mus, mua, hgg, n to find with inverse MCRT"
            return 
        end if

        !loop through the layers finding the index in the SDF array of the selected layer and its initial optical properties
        SDF_array_index = -1
        do i = 1, size(array)
            if (array(i)%getLayer() == layer) then 
                !we have found the layer, store the index of this and its optical properties
                SDF_array_index = i
                mua = array(i)%getMua()
                mus = array(i)%getKappa() - mua
                hgg = array(i)%gethgg()
                n = array(i)%getN()
                exit
            end if
        end do

        !check that the selected layer is found in the SDF array
        if (SDF_array_index == -1) then
            print*, "Selected layer not found in SDF array please choose a layer inside the SDF array"
            return
        end if


        !set the values for AdaLIPO
        probability = 1.0_wp
        alpha = 0.01
        nb_samples = 0
        k = 0.0_wp

        musboundlower = 0.0_wp
        musboundupper = 100.0_wp
        muaboundlower = 0.0_wp
        muaboundupper = 100.0_wp
        gboundlower = -1.0_wp
        gboundupper = 1.0_wp
        nboundlower = 1.0_wp
        nboundupper = 20.0_wp


        allocate(gradDescentData(maxNumSteps, 5))
        if (findmus) then
            gradDescentData(1,1) = ran2() * (musboundupper-musboundlower) + musboundlower
        else 
            gradDescentData(1,1) = mus
        end if 
        if (findmua) then
            gradDescentData(1,2) = ran2() * (muaboundupper-muaboundlower) + muaboundlower
        else
            gradDescentData(1,2) = mua
        end if 
        if (findg) then
            gradDescentData(1,3) = ran2() * (gboundupper-gboundlower) + gboundlower
        else 
            gradDescentData(1,3) = hgg
        end if 
        if (findn) then
            gradDescentData(1,4) = ran2() * (nboundupper-nboundlower) + nboundlower
        else 
            gradDescentData(1,4) = n
        end if

        !evaluate, by running MCRT
        call run_MCRT(input_file, history, packet, dict, & 
                        distances, image, dects, array, nscatt, start, & 
                        tev, spectrum)
        error = 0._wp
        call inverse_evaluate(dects, error)
        gradDescentData(1,5) = error

        ! reset the arrays storing data
        call reset(dects)  

        !store the position of minimum error
        minError = gradDescentData(1,5)
        indexOfMinError = 1

        !allocate the ratio array (it will be of the size maxNumSteps triangular number)
        allocate(ratios(int(nint(maxNumSteps * (maxNumSteps+1)/2.0_wp))))
        ratios = 0.0_wp
        ratioCounter = 1
        maxRatio = 1
        indexOfMaxRatio = 1

        do i = 2, maxNumSteps

            
            ran = ran2()
            if( ran <= 1.0_wp) then
                !we are in the explore stage
                !get the new guesses for the mua, mus, n, and g
                if (findmus) then
                    gradDescentData(i,1) = ran2() * (musboundupper-musboundlower) + musboundlower
                else 
                    gradDescentData(i,1) = mus
                end if 
                if (findmua) then
                    gradDescentData(i,2) = ran2() * (muaboundupper-muaboundlower) + muaboundlower
                else
                    gradDescentData(i,2) = mua
                end if 
                if (findg) then
                    gradDescentData(i,3) = ran2() * (gboundupper-gboundlower) + gboundlower
                else 
                    gradDescentData(i,3) = hgg
                end if 
                if (findn) then
                    gradDescentData(i,4) = ran2() * (nboundupper-nboundlower) + nboundlower
                else 
                    gradDescentData(i,4) = n
                end if

                nb_samples = nb_samples + 1

                !update the optical properties
                trialOptProp = mono(mus, mua, hgg, n)
                temp = array(SDF_array_index)%updateOptProp(trialOptProp)

                !evaluate, by running MCRT
                call run_MCRT(input_file, history, packet, dict, & 
                            distances, image, dects, array, nscatt, start, & 
                            tev, spectrum)
                error = 0._wp
                call inverse_evaluate(dects, error)
                gradDescentData(i,5) = error

                ! reset the arrays storing data
                call reset(dects) 
            else
                do while(.true.)
                    !get the new guesses for the mua, mus, n, and g
                    if (findmus) then
                        gradDescentData(i,1) = ran2() * (musboundupper-musboundlower) + musboundlower
                    else 
                        gradDescentData(i,1) = mus
                    end if 
                    if (findmua) then
                        gradDescentData(i,2) = ran2() * (muaboundupper-muaboundlower) + muaboundlower
                    else
                        gradDescentData(i,2) = mua
                    end if 
                    if (findg) then
                        gradDescentData(i,3) = ran2() * (gboundupper-gboundlower) + gboundlower
                    else 
                        gradDescentData(i,3) = hgg
                    end if 
                    if (findn) then
                        gradDescentData(i,4) = ran2() * (nboundupper-nboundlower) + nboundlower
                    else 
                        gradDescentData(i,4) = n
                    end if

                    !update the optical properties
                    trialOptProp = mono(mus, mua, hgg, n)
                    temp = array(SDF_array_index)%updateOptProp(trialOptProp)

                    !left_min = min(gradDescentData + k * (gradDescentData(i,1)**2) )
                end do
            end if
        end do
        
        
        !we have finished the gradient descent output the gradDescentData to a file and write the error
        



        !big goals, get the run_MCRT running on the GPU thus speeding up the program exponentially
        !add more SDFs, and add a tool to create a smooth SDF from a given grid mesh 

        

    end subroutine inverse_MCRT

    subroutine inverse_evaluate(dects, error)
        !Calculate the error between the detector target values and given detector actual values

        use constants, only : wp
        use detectors,     only : dect_array
        use sim_state_mod, only : state

        type(dect_array), allocatable, intent(inout) :: dects(:)
        real(kind=wp), intent(inout) :: error

        integer :: counter, loopCounter
        real(kind = wp) :: total, targetVal

        error = 0._wp
        counter = 0

        do loopCounter = 1, size(dects)
            targetVal = dects(loopCounter)%p%targetValue
            if (targetVal /= -1) then

                !get the total value of the detector
                total = 0._wp
                call dects(loopCounter)%p%total_dect(total)
                total = total / real(state%nphotons, kind=wp)

                !error is defined as average absolute difference between all detectors
                error = error + abs((total-targetVal))

                counter = counter + 1
            end if
        end do

        error = -error/counter

    end subroutine inverse_evaluate


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
        packet%layer=(maxloc(distances,dim=1, mask=(distances<0._wp)))
        
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
    ! this takes a long time to run, can it be sped up?
    !call zarray()

    !zero the detectors
    do i = 1, size(dects)
        call dects(i)%p%zero_dect()
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