module escapeFunctionMod

    implicit none

    private
    public :: escape_Function

contains
    !calculate the escape function for each detector    
    subroutine escape_Function(input_file)

        !Shared data
        use iarray
        use constants, only : wp

        !subroutines
        use detectors
        use historyStack,  only : history_stack_t
        use inttau2,       only : tauint2
        use photonmod
        use piecewiseMod
        use random,        only : ran2, init_rng
        use sdfs,          only : sdf
        use sdfHelpers,    only : rotationAlign, rotmat
        use sim_state_mod, only : state
        use vector_class
        use setupMod, only : setup_escapeFunction, zarray
        use writer_mod, only : write_escape
        use kernels, only : setup, finalise, reset_detectors

        use default_MCRTMod, only : run_MCRT

        !external deps
        use tev_mod, only : tevipc
        use tomlf,   only : toml_table, toml_error, get_value, set_value
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
        integer :: m1, n1, o1, loopCounter
        real(kind = wp) :: x,y,z, total
        type(vector) :: position, direction, gridPos
        real(kind=wp) :: rotationOnToSym(4,4), rotationOffSym(4,4)
        real(kind=wp) :: rotationAroundZOnSym(4,4), rotationAroundZOffSym(4,4)
        character(len=:), allocatable :: symmetryType
        integer :: indices(3)
        real :: tic, toc

        !temporary while testing adjoint
        real(kind=wp) :: posX, posY, posZ, dirX, dirY, dirZ, radius, acceptAngle
        character(len=:), allocatable :: dectType, dectID
        character(len=8) :: countStr
        type(vector) :: poss, dirr


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

            !setup symmetry grid
            call setup_cart_symGrid(dects)

            packet = photon("escapeCartSymmetry")

            !loop through every cell
            do m = 1, state%symmetryEscapeCartGrid%nxg
                do n = 1, state%symmetryEscapeCartGrid%nyg
                    do o = 1, state%symmetryEscapeCartGrid%nzg

                        call update_sumGridCells_inDict(dict, m, n, o)

                        print*, ""
                        print*, "Running ", ((m-1)*state%symmetryEscapeCartGrid%nyg*state%symmetryEscapeCartGrid%nzg + & 
                                             (n-1)*state%symmetryEscapeCartGrid%nzg + o - 1), & 
                                " out of ", (state%symmetryEscapeCartGrid%nxg* &
                                                state%symmetryEscapeCartGrid%nyg* & 
                                                state%symmetryEscapeCartGrid%nzg)

                        !calculate the escape function
                        call cart_calc_escape_sym(m,n,o, dects, array,& 
                                                 packet, distances, dict, history, image, input_file, nscatt, spectrum,& 
                                                 start, tev)

                    end do         
                end do
            end do

            !Go through the base grid and use some form of interpolation to figure out the best match
            call cart_map_escape_sym(dects)

        case("prism")
            !prism symmetry, launch from a layer of cells

            print*, "Prism symmetry selected"
            print*, "Number of Monte Carlo Simmulations to run: ", (state%symmetryEscapeCartGrid%nxg* &
                                                                    state%symmetryEscapeCartGrid%nyg)
            print*, ""

            !setup symmetry grid
            call setup_cart_symGrid(dects)

            packet = photon("escapeCartSymmetry")

            !get the position of the cell
            indices = state%symmetryEscapeCartGrid%get_voxel(vector(0.0_wp,0.0_wp,0.0_wp))
            !loop through every cell
            do m = 1, state%symmetryEscapeCartGrid%nxg
                do n = 1, state%symmetryEscapeCartGrid%nyg

                    call update_sumGridCells_inDict(dict, m, n, indices(3))

                    print*, ""
                    print*, "Running ", ((m-1)*state%symmetryEscapeCartGrid%nyg + n - 1), & 
                            " out of ", (state%symmetryEscapeCartGrid%nxg* &
                                        state%symmetryEscapeCartGrid%nyg)

                    !calculate the escape function
                    call cart_calc_escape_sym(m,n,indices(3), dects, array,& 
                                                packet, distances, dict, history, image, input_file, nscatt, spectrum,& 
                                                start, tev)
                end do
            end do

            !fill the rest of the grid
            do o = 1, state%symmetryEscapeCartGrid%nzg
                escapeSymmetry(:, :, :, o) = escapeSymmetry(:, :, :, indices(3))
            end do

            !Go through the base grid and use some form of interpolation to figure out the best match
            call cart_map_escape_sym(dects)

        case("flipped")
            !flipped symmetry, launch half the cells

            print*, "Flipped symmetry selected"
            print*, "Number of Monte Carlo Simmulations to run: ", (state%symmetryEscapeCartGrid%nxg* &
                                                                    state%symmetryEscapeCartGrid%nyg* &
                                                                    (state%symmetryEscapeCartGrid%nzg/2)+1)
            print*, ""
            
            !setup symmetry grid
            call setup_cart_symGrid(dects)

            packet = photon("escapeCartSymmetry")

            !get the position of the cell
            indices = state%symmetryEscapeCartGrid%get_voxel(vector(0.0_wp,0.0_wp,0.0_wp))
            !loop through every cell
            do m = 1, state%symmetryEscapeCartGrid%nxg
                do n = 1, state%symmetryEscapeCartGrid%nyg
                    do o = 1, (state%symmetryEscapeCartGrid%nzg/2)+1

                        call update_sumGridCells_inDict(dict, m, n, o)

                        print*, ""
                        print*, "Running ", ((m-1)*state%symmetryEscapeCartGrid%nyg*((state%symmetryEscapeCartGrid%nzg/2)+1) + & 
                                            (n-1)*((state%symmetryEscapeCartGrid%nzg/2)+1) + o - 1), & 
                                " out of ", (state%symmetryEscapeCartGrid%nxg* &
                                            state%symmetryEscapeCartGrid%nyg* &
                                            (state%symmetryEscapeCartGrid%nzg/2)+1)

                        call cart_calc_escape_sym(m,n,o, dects, array,& 
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
            call cart_map_escape_sym(dects)

        case("uniformSlab")
            ! The simmulation is a slab code, light is collected uniformly

            print*, "Uniform slab symmetry selected"           
            print*, "Number of Monte Carlo Simmulations to run: ", (state%symmetryEscapeCartGrid%nzg)
            print*, ""
            
            !setup symmetry grid
            call setup_cart_symGrid(dects)

            packet = photon("escapeCartSymmetry")

            !get the position of the cell
            indices = state%symmetryEscapeCartGrid%get_voxel(vector(0.0_wp,0.0_wp,0.0_wp))
            !loop through every cell
            do o = 1, state%symmetryEscapeCartGrid%nzg

                call update_sumGridCells_inDict(dict, indices(1), indices(2), o)

                print*, ""
                print*, "Running ", (o - 1), & 
                        " out of ", (state%symmetryEscapeCartGrid%nzg)

                call cart_calc_escape_sym(indices(1),indices(2),o, dects, array,& 
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
            call cart_map_escape_sym(dects)

        case("noneRotational")
            ! Do for all radii, theta and z values

            print*, "No Symmetry in cylindrical coordinates selected"
            print*, "User warning! This may take a long time to run"
            print*, "It is advised to try and find a geometry with symmetry or to reduce the grid size"
            print*, "Number of Monte Carlo Simmulations to run: ", (state%symmetryEscapeCylGrid%nrg* &
                                                                    state%symmetryEscapeCylGrid%ntg* & 
                                                                    state%symmetryEscapeCylGrid%nzg)
            print*, ""

            !setup cylindrical symmetry grid
            call setup_cyl_symGrid(dects)

            packet = photon("escapeCylSymmetry")

            !loop through every cell
            do m = 1, state%symmetryEscapeCylGrid%nrg
                do n = 1, state%symmetryEscapeCylGrid%ntg
                    do o = 1, state%symmetryEscapeCylGrid%nzg

                        call update_sumGridCells_inDict(dict, m, n, o)

                        print*, ""
                        print*, "Running ", ((m-1)*state%symmetryEscapeCylGrid%ntg*state%symmetryEscapeCylGrid%nzg + & 
                                             (n-1)*state%symmetryEscapeCylGrid%nzg + o - 1), & 
                                " out of ", (state%symmetryEscapeCylGrid%nrg* &
                                            state%symmetryEscapeCylGrid%ntg* & 
                                            state%symmetryEscapeCylGrid%nzg)

                        !calculate the escape function
                        call cyl_calc_escape_sym(m,n,o, dects, array,& 
                                                 packet, distances, dict, history, image, input_file, nscatt, spectrum,& 
                                                 start, tev)

                    end do         
                end do
            end do

            !Go through the base grid and use some form of interpolation to figure out the best match
            call cyl_map_escape_sym(dects)

        case("360rotational")
            ! Do for all radii and z values at one theta value

            print*, "360 Rotational symmetry in cylindrical coordinates selected"
            print*, "User warning! This may take a long time to run"
            print*, "It is advised to try and find a geometry with symmetry or to reduce the grid size"
            print*, "Number of Monte Carlo Simmulations to run: ", (state%symmetryEscapeCylGrid%nrg* & 
                                                                    state%symmetryEscapeCylGrid%nzg)
            print*, ""

            !setup cylindrical symmetry grid
            call setup_cyl_symGrid(dects)

            packet = photon("escapeCylSymmetry")

            n=1
            do m = 1, state%symmetryEscapeCylGrid%nrg
                do o = 1, state%symmetryEscapeCylGrid%nzg

                    call update_sumGridCells_inDict(dict, m, n, o)

                    !calculate the escape function
                    print*, ""
                    print*, "Running ", ((m-1)*state%symmetryEscapeCylGrid%nzg + o - 1), & 
                            " out of ", (state%symmetryEscapeCylGrid%nrg*state%symmetryEscapeCylGrid%nzg)
                    call cyl_calc_escape_sym(m,n,o, dects, array,& 
                                                packet, distances, dict, history, image, input_file, nscatt, spectrum,& 
                                                start, tev)

                end do 
            end do

            !loop through every cell
            
            do n = 1, state%symmetryEscapeCylGrid%ntg
                escapeSymmetry(:,:,n,:) = escapeSymmetry(:,:,1,:)
            end do

            !Go through the base grid and use some form of interpolation to figure out the best match
            call cyl_map_escape_sym(dects)
        case("adjoint")
            !Use the adjoint method to calculate the escape function

            print*, "Adjoint symmetry in cartesian symmetry"
            print*, "Using the absorption from each detector for the escape function"
            print*, symmetryType

            do n=1, size(dects)

                write(countStr, "(I8)") n

                call set_value(dict, "dectCount", countStr)
                call get_value(dict, "dect"//countStr//"type", dectType)
                call get_value(dict, "dect"//countStr//"ID", dectID)

                call get_value(dict, "dect"//countStr//"position%x", posX)
                call get_value(dict, "dect"//countStr//"position%y", posY)
                call get_value(dict, "dect"//countStr//"position%z", posZ)
                call get_value(dict, "dect"//countStr//"direction%x", dirX)
                call get_value(dict, "dect"//countStr//"direction%y", dirY)
                call get_value(dict, "dect"//countStr//"direction%z", dirZ)
            
                print*, " "
                print*, n, dectType, dectID

                !setup the new source type for the nth detector
                poss = vector(posX, posY, posZ)
                dirr = vector(-1.0_wp*dirX, -1.0_wp*dirY, -1.0_wp*dirZ)
                call set_photon(poss, dirr)
                packet%nxp = 1.0_wp 
                packet%nyp = 0.0_wp 
                packet%nzp = 0.0_wp 

                !set the detector type
                if (dectType == "circle") then
                    packet = photon("circleDect")
                else if (dectType == "annulus") then
                    packet = photon("annulusDect")
                !else if (dectType == "focus") then
                !    packet = photon("circleDect")
                else
                    print*, "source not implemented, set escape funciton to zero"
                    escape(n,:,:,:) = 0.0_wp
                    cycle
                end if               
        
                !zero all arrays
                call zarray()

                !run the adjoint MCRT
                call run_MCRT(input_file, history, packet, dict, & 
                            distances, image, dects, array, nscatt, start, & 
                            tev, spectrum)

                !store the escape function for nth detector
                escape(n,:,:,:) = jmean(:,:,:)
            end do

            !find where the layer is not part of the egg and set fluence to zero
            do m1 = 1, state%grid%nxg
                do n1 = 1, state%grid%nyg
                    do o1 = 1, state%grid%nzg

                        !get the coords at the centre of the voxel
                        x = (((real(m1, kind = wp) - 0.5)/state%grid%nxg)*& 
                            2.0_wp*state%grid%xmax) - state%grid%xmax 
                        y = (((real(n1, kind = wp) - 0.5)/state%grid%nyg)*& 
                            2.0_wp*state%grid%ymax) - state%grid%ymax 
                        z = (((real(o1, kind = wp) - 0.5)/state%grid%nzg)*& 
                            2.0_wp*state%grid%zmax) - state%grid%zmax 

                        ! get the layer at this position
                        distances = 0._wp
                        do loopCounter = 1, size(distances)
                            distances(loopCounter) = array(loopCounter)%evaluate(vector(x,y,z))
                        end do
                        layer=(maxloc(distances,dim=1, mask=(distances<0._wp)))

                        if ((layer == 0) .or. (array(layer)%getkappa() == real(0, kind=wp))) then
                            do n = 1, size(dects)
                                escape(n, m1, n1, o1) = 0.0_wp
                            end do
                        end if
                    end do
                end do 
            end do

        case default                     
            print*,"Unknown symmetry type"
            stop 1
        end select
                                            
        !store the escape funcitons for each detector
        call write_escape(dects, symmetryType, dict)

        call finalise(dict, dects, nscatt, start, history)

        call cpu_time(toc)
        print*,"Time to Run: ",((toc - tic))

    end subroutine escape_Function



    subroutine setup_cart_symGrid(dects)
        use constants, only : wp
        use iarray
        use detectors
        use vector_class
        use sdfs,       only : sdf
        use sdfHelpers, only : rotationAlign, rotmat
        use sim_state_mod

        type(dect_array), allocatable, intent(inout) :: dects(:)

        type(vector) :: direction, gridPos
        real(kind=wp) :: rotationOnToSym(4,4), rotationOffSym(4,4)
        real(kind=wp) :: rotationAroundZOnSym(4,4), rotationAroundZOffSym(4,4)

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

        !store rotation matrices in state
        state%rotationOffSym = rotationOffSym
        state%rotationOnToSym = rotationOnToSym
        state%rotationAroundZOffSym = rotationAroundZOffSym
        state%rotationAroundZOnSym = rotationAroundZOnSym
        state%gridPos = gridPos

    end subroutine setup_cart_symGrid

    subroutine setup_cyl_symGrid(dects)
        use constants, only : wp
        use iarray
        use detectors
        use vector_class
        use sdfs,       only : sdf
        use sdfHelpers, only : rotationAlign, rotmat
        use sim_state_mod

        type(dect_array), allocatable, intent(inout) :: dects(:)

        type(vector) :: direction, gridPos
        real(kind=wp) :: rotationOnToSym(4,4), rotationOffSym(4,4)
        real(kind=wp) :: rotationAroundZOnSym(4,4), rotationAroundZOffSym(4,4)

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

        state%rotationOffSym = rotationOffSym
        state%rotationOnToSym = rotationOnToSym
        state%rotationAroundZOffSym = rotationAroundZOffSym
        state%rotationAroundZOnSym = rotationAroundZOnSym
        state%gridPos = gridPos

    end subroutine setup_cyl_symGrid

    subroutine update_sumGridCells_inDict(dict, cellx, celly, cellz)
        !external deps
        use tev_mod, only : tevipc
        use tomlf,   only : toml_table, toml_error, get_value, set_value

        type(toml_table), intent(inout) :: dict
        integer, intent(in) :: cellx, celly, cellz

        call set_value(dict, "symGridCellx", cellx)
        call set_value(dict, "symGridCelly", celly)
        call set_value(dict, "symGridCellz", cellz)
    end subroutine update_sumGridCells_inDict



    subroutine cart_calc_escape_sym(m,n,o, dects, array, packet, & 
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
        use default_MCRTMod, only : run_MCRT
        use kernels, only : reset_detectors

        !external deps
        use tev_mod, only : tevipc
        use tomlf,   only : toml_table, toml_error, get_value

        !> indices of symmetryEscapeCartGrid
        integer, intent(in) :: m,n,o
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
        call reset_detectors(dects)

        ! find the centre position of the voxel
        x = (((real(m, kind = wp) - 0.5)/state%symmetryEscapeCartGrid%nxg)*& 
            2.0_wp*state%symmetryEscapeCartGrid%xmax) - state%symmetryEscapeCartGrid%xmax 
        y = (((real(n, kind = wp) - 0.5)/state%symmetryEscapeCartGrid%nyg)*& 
            2.0_wp*state%symmetryEscapeCartGrid%ymax) - state%symmetryEscapeCartGrid%ymax 
        z = (((real(o, kind = wp) - 0.5)/state%symmetryEscapeCartGrid%nzg)*& 
            2.0_wp*state%symmetryEscapeCartGrid%zmax) - state%symmetryEscapeCartGrid%zmax 
        
        position = vector(x,y,z)

        !rotate to align x and y axis after z axis alignment
        position = position .dot. state%rotationAroundZOffSym

        !align z axis
        position = position .dot. state%rotationOffSym

        !shift
        position = position + state%gridPos

        !set the emission location to the centre of the voxel
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
            do loopCounter = 1, size(dects)
                escapeSymmetry(loopCounter, m, n, o) = 0.0_wp

                !temporary while testing
                !escapeSymmetry(loopCounter, m, n, o) = layer
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            end do
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

    subroutine cart_map_escape_sym(dects)

        use iarray
        use constants, only : wp
        use sim_state_mod, only : state
        use vector_class, only : vector
        use detectors
        use interpolate

        type(dect_array), allocatable, intent(inout) :: dects(:)

        integer :: i,j,k, m,n,o
        integer :: loopCounter, indx(3)
        integer :: xIndx(2), yIndx(2),zIndx(2)
        real(kind=wp) :: x,y,z, closestX, closestY, closestZ
        logical :: notOnXedge, notOnYedge, notOnZedge
        real(kind=wp) :: corners3D(2,2,2,4), corners2D(2,2,3), corners1D(2,2), point3D(4), point2D(3), point1D(2)
        type(vector) :: position
        
        print*, " "
        print*, "Starting Interpolation"

        !Go through the base grid and use some form of interpolation to figure out the best match
        do m = 1, state%grid%nxg
            do n = 1, state%grid%nyg
                do o = 1, state%grid%nzg

                    ! find the centre position of the voxel
                    y = (((real(n, kind = wp) - 0.5)/state%grid%nyg)*2.0_wp*state%grid%ymax) - state%grid%ymax
                    x = (((real(m, kind = wp) - 0.5)/state%grid%nxg)*2.0_wp*state%grid%xmax) - state%grid%xmax
                    z = (((real(o, kind = wp) - 0.5)/state%grid%nzg)*2.0_wp*state%grid%zmax) - state%grid%zmax

                    position = vector(x,y,z) 

                    !shift
                    position = position - state%gridPos

                    !rotate, there is none for this geometry
                    position = position .dot. state%rotationOnToSym

                    !rotate to align x and y axis after z axis alignment
                    position = position .dot. state%rotationAroundZOnSym
                    
                    !find the points in symmetry escape that correspond to this point?
                    !this returns the point that is closest to this
                    indx = -1
                    indx = state%symmetryEscapeCartGrid%get_voxel(position)

                    !we are out of the bounds of the escapesymmetry grid
                    if (indx(1) == -1 .or. indx(2) == -1 .or. indx(3) == -1) then
                        do loopCounter = 1, size(dects)
                            escape(loopCounter, m, n, o) = 0._wp
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

        print*, " "
        print*, "Finished Interpolation"
    end subroutine cart_map_escape_sym

    subroutine cyl_calc_escape_sym(m,n,o, dects, array, packet, & 
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
        use default_MCRTMod, only : run_MCRT
        use kernels, only : reset_detectors

        !external deps
        use tev_mod, only : tevipc
        use tomlf,   only : toml_table, toml_error, get_value

        !> indices of symmetryEscapeCartGrid
        integer, intent(in) :: m,n,o
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
        call reset_detectors(dects)

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
        position = position .dot. state%rotationAroundZOffSym

        !align z axis
        position = position .dot. state%rotationOffSym

        !shift
        position = position + state%gridPos

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
            do loopCounter = 1, size(dects)
                escapeSymmetry(loopCounter, m, n, o) = 0.0_wp
                
                !temporary while testing 
                !escapeSymmetry(loopCounter, m, n, o) = layer
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            end do
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
                !escapeSymmetry(loopCounter, m, n, o) = x + y**2 + z**3
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            else
                escapeSymmetry(loopCounter, m, n, o) = 0.0_wp

                !temporary while testing 
                !escapeSymmetry(loopCounter, m, n, o) = layer
                !escapeSymmetry(loopCounter, m, n, o) = x + y**2 + z**3
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            end if
        end do

    end subroutine cyl_calc_escape_sym
    
    subroutine cyl_map_escape_sym(dects)

        use iarray
        use constants, only : wp, PI, TWOPI
        use sim_state_mod, only : state
        use vector_class, only : vector
        use detectors
        use interpolate

        type(dect_array), allocatable, intent(inout) :: dects(:)

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
        
        print*, " "
        print*, "Starting Interpolation"

        !Go through the base grid and use some form of interpolation to figure out the best match
        do m = 1, state%grid%nxg
            do n = 1, state%grid%nyg
                do o = 1, state%grid%nzg

                    
                    ! find the centre position of the voxel
                    y = (((real(n, kind = wp) - 0.5)/state%grid%nyg)*2.0_wp*state%grid%ymax) - state%grid%ymax
                    x = (((real(m, kind = wp) - 0.5)/state%grid%nxg)*2.0_wp*state%grid%xmax) - state%grid%xmax
                    z = (((real(o, kind = wp) - 0.5)/state%grid%nzg)*2.0_wp*state%grid%zmax) - state%grid%zmax

                    position = vector(x,y,z) 

                    !shift
                    position = position - state%gridPos

                    !rotate, there is none for this geometry
                    position = position .dot. state%rotationOnToSym

                    !rotate to align x and y axis after z axis alignment
                    position = position .dot. state%rotationAroundZOnSym



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
                            escape(loopCounter, m, n, o) = 0._wp
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

        print*, " "
        print*, "Finished Interpolation"
    end subroutine cyl_map_escape_sym


end module escapeFunctionMod