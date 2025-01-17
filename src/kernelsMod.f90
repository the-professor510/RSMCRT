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
        
        call setup_escapeFunction(size(dects))
        
        !use symmetry to reduce the number of voxels we have to calculate
        call get_value(dict, "symmetryType", symmetryType)
        select case(symmetryType)
        case("none")
            !there is no symmetry launch from every cell

            print*, "User warning! This may take a long time to run"
            print*, "It is advised to try and find a geometry with symmetry or to reduce the grid size"
            
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
            !TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !what is the normal to the plane this prism is in?


            !first allocate the escapeSymmetry, with a size that is determined by the user,
            !we will want the following settings
            !cartesian base: none, prism, flipped, specified
            !   we need a vector for num voxels
            !   we need a vector for dimensions, going out from zero
            !   we need a vector for central position for the coordinate system
            !   we need a normal vector for the prism plane and the flipped plane
            !   we need a series of points to calculate this for
            !cylindrical base: 360 rotational, n rotational
            !   we need a vector for num voxels (num in r, num in theta, num in z)
            !   we need a vector for max dimensions going out from zero (max r, max z)
            !   we need a vector for central position of the symmetry coordinate system
            !   we need a vector for the axis along which the z direction points
            !allocate(escapeSymmetry(size(dects), state%grid%nxg, state%grid%nyg, state%grid%nzg))

            !this is now done, we have a list of variables that define the system so we can allocate many two different sets of grids

            ! when working we will define the symmetry escape function looping through the x,y,z if cart or r, theta, z if cylindrical
            ! we will find what the point in the symmetry coordinate space coordinates to in the given x,y,z
            ! we will find if we have a non-zero kappa, which if we do we will then find the escape function for
            ! we will then store this escape function in the symmetry escape function array
            ! we will then loop through the base coordinate escape and convert each position into the symmetry coordinate
            ! we will find the 8 closest points to perform trilinear interpolation
            ! if we do not have 8 points we will set it to -1

            !Jobs To Do

            !Add section to input toml to read symmetry settings for the escape function, make these separate from the geometry now

            ! Add function to shift and rotate between coordinate systems

            ! begin with none, then do prism, then do flipped. To do prism, you will need to implement first a shift of coordinates and then a rotation of coordinates
            ! none
            ! prism
            ! flipped
            ! specified
            ! 360rotational
            ! nrotational

            ! Fill the escape function, if it is greater than the cartesian grid, then ignore it and don't find it

            !loop through every cell
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
                        !if(array(layer)%getkappa() /= real(0, kind=wp)) then
                        !    call run_MCRT(input_file, history, packet, dict, & 
                        !                    distances, image, dects, array, nscatt, start, & 
                        !                    tev, spectrum)
                        !end if
                        
                        ! record the efficiency for each detector and add to an array of escape functions
                        do i = 1, size(dects)
                            total = 0._wp
                            call dects(i)%p%total_dect(total)
                        !    escape(i, m, n, o) = total/state%nphotons
                        end do
                    end do         
                end do
            end do
        case("flipped")
            !TO DO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Only do half the cells, ensure the you know the plane of division
            !
            allocate(escapeSymmetry(size(dects), state%grid%nxg, state%grid%nyg, state%grid%nzg))

            !loop through every cell
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
                        !if(array(layer)%getkappa() /= real(0, kind=wp)) then
                        !    call run_MCRT(input_file, history, packet, dict, & 
                        !                    distances, image, dects, array, nscatt, start, & 
                        !                    tev, spectrum)
                        !end if
                        
                        ! record the efficiency for each detector and add to an array of escape functions
                        do i = 1, size(dects)
                            total = 0._wp
                            call dects(i)%p%total_dect(total)
                        !    escape(i, m, n, o) = total/state%nphotons
                        end do
                    end do         
                end do
            end do
        case("360rotational")
            !TO DO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! do radially and in depth only for a single angle


            ! find the axis along which we have symmetry

            ! get an axis along which we calculate the symmetry, we want to choose our radii and selection of z such that 
            ! we are garunteed to be within the grid

            ! rotate to align the z axis, assume that there is no difference due to theta, use dot product 
            !                                                       between desired axis and current axis
            ! rotate to align the theta axis, use dot product desired axis and current axis
            ! shift to align the central coordinates with one another
            ! use bilinear interprolation to convert from cylindrical to cartesian

            allocate(escapeSymmetry(size(dects), state%grid%nxg, state%grid%nyg, state%grid%nzg))

            !loop through every cell
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
                        !if(array(layer)%getkappa() /= real(0, kind=wp)) then
                        !    call run_MCRT(input_file, history, packet, dict, & 
                        !                    distances, image, dects, array, nscatt, start, & 
                        !                    tev, spectrum)
                        !end if
                        
                        ! record the efficiency for each detector and add to an array of escape functions
                        do i = 1, size(dects)
                            total = 0._wp
                            call dects(i)%p%total_dect(total)
                            escape(i, m, n, o) = total/state%nphotons
                        end do
                    end do         
                end do
            end do
        case("nrotational")
            !TO DO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! do radially and in depth only for 360/n region of angles
            allocate(escapeSymmetry(size(dects), state%grid%nxg, state%grid%nyg, state%grid%nzg))

            !loop through every cell
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
                        !if(array(layer)%getkappa() /= real(0, kind=wp)) then
                        !    call run_MCRT(input_file, history, packet, dict, & 
                        !                    distances, image, dects, array, nscatt, start, & 
                        !                    tev, spectrum)
                        !end if
                        
                        ! record the efficiency for each detector and add to an array of escape functions
                        do i = 1, size(dects)
                            total = 0._wp
                            call dects(i)%p%total_dect(total)
                            escape(i, m, n, o) = total/state%nphotons
                        end do
                    end do         
                end do
            end do
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
        packet = photon(state%source)
        packet%pos = position

        ! get the layer at this position
        distances = 0._wp
        do loopCounter = 1, size(distances)
            distances(loopCounter) = array(loopCounter)%evaluate(position)
        end do
        layer=maxloc(distances,dim=1, mask=(distances<0._wp))

        !!  if the layer has a non-zero kappa then it is significant and we want to perform MCRT
        !if(array(layer)%getkappa() /= real(0, kind=wp)) then
        !    call run_MCRT(input_file, history, packet, dict, & 
        !                    distances, image, dects, array, nscatt, start, & 
        !                    tev, spectrum)
        !end if
        
        ! record the efficiency for each detector and add to an array of escape functions
        do loopCounter = 1, size(dects)
            if(array(layer)%getkappa() /= real(0, kind=wp)) then
                !total = 0._wp
                !call dects(loopCounter)%p%total_dect(total)
                !escapeSymmetry(loopCounter, m, n, o) = total/state%nphotons


                !temporary while testing 
                escapeSymmetry(loopCounter, m, n, o) = layer
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            else
                escapeSymmetry(loopCounter, m, n, o) = 0.0_wp

                !temporary while testing 
                escapeSymmetry(loopCounter, m, n, o) = layer
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