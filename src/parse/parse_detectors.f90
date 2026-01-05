module parse_detectorsMod
    !! routine to parse the detector table from the input Toml file.
    use constants, only : wp
    use parse_HelpersMod
    use vector_class

    use tomlf
    use tomlf_error

    implicit none

    private
    public :: parse_detectors

contains
    
    subroutine parse_detectors(table, dict, dects, context, error)
        !! parse the detectors

        use detectors,     only : dect_array, circle_dect, annulus_dect, camera, fibre_dect
        use sim_state_mod, only : state

        !> Input Toml table
        type(toml_table),               intent(inout) :: table
        !> Dictonary used to store metadata
        type(toml_table),               intent(inout) :: dict
        !> Detector array to be filled.
        type(dect_array), allocatable,  intent(out)   :: dects(:)
        !> Context handle for error reporting.
        type(toml_context),             intent(in)    :: context
        !> Error message
        type(toml_error),  allocatable, intent(out)   :: error

        type(toml_array), pointer :: array
        type(toml_table), pointer :: child
        character(len=:), allocatable :: dect_type, dect_ID
        type(circle_dect), target, save, allocatable :: dect_c(:)
        type(annulus_dect), target, save, allocatable :: dect_a(:)
        type(fibre_dect), target, save, allocatable :: dect_f(:)
        type(camera), target, save, allocatable :: dect_cam(:)
        integer :: i, c_counter, a_counter, f_counter, cam_counter, j, k,origin, l
        real(kind=wp) :: targetValue

        c_counter = 0
        a_counter = 0
        f_counter = 0
        cam_counter = 0
        call get_value(table, "detectors", array)
        allocate(dects(len(array)))

        do i = 1, len(array)
            call get_value(array, i, child)
            call get_value(child, "type", dect_type, origin=origin)
            call get_value(child, "ID", dect_ID, origin=origin)
            if(.not. allocated(dect_ID)) then
                !throw an error, the user must specify a detector ID
                call make_error(error, context%report("Need to specify a detector ID", origin, &
                                              "No detector ID specified"), -1)
                return
            end if
            select case(dect_type)
            case default
                call make_error(error, &
                    context%report("Invalid detector type. Valid types are [circle, annulus, camera]", &
                    origin, "expected valid detector type"), -1)
                return
            case("circle")
                c_counter = c_counter + 1
            case("annulus")
                a_counter = a_counter + 1
            case("fibre")
                f_counter = f_counter + 1
            case("camera")
                cam_counter = cam_counter + 1
            end select
        end do

        if(c_counter > 0)then
            if(allocated(dect_c))deallocate(dect_c)
            allocate(dect_c(c_counter))
        end if
        if(a_counter > 0)then
            if(allocated(dect_a))deallocate(dect_a)
            allocate(dect_a(a_counter))
        end if
        if(f_counter > 0)then
            if(allocated(dect_f))deallocate(dect_f)
            allocate(dect_f(f_counter))
        end if
        if(cam_counter > 0)then
            if(allocated(dect_cam))deallocate(dect_cam)
            allocate(dect_cam(cam_counter))
        end if
        c_counter = 1
        a_counter = 1
        f_counter = 1
        cam_counter = 1
        state%trackHistory=.false.
        do i = 1, len(array)
            call get_value(array, i, child)
            call get_value(child, "type", dect_type)
            call get_value(child, "ID", dect_ID, "none", origin=origin)
            !call get_value(child, "historyFileName", state%historyFilename, "photPos.obj")
            call get_value(child, "inverseTarget", targetValue, -1._wp)
            select case(dect_type)
            case("circle")
                call handle_circle_dect(child, dict, dect_c, c_counter, context, error, dect_ID, targetValue)
                if(allocated(error))return
            case("annulus")
                call handle_annulus_dect(child, dict, dect_a, a_counter, context, error, dect_ID, targetValue)
                if(allocated(error))return
            case("fibre")
                call handle_fibre_collection_dect(child, dect_f, f_counter, context, error, dect_ID, targetValue)
                if(allocated(error))return
            case("camera")
                call handle_camera(child, dect_cam, cam_counter, context, error, dect_ID, targetValue)
                if(allocated(error))return
            end select
        end do

        do i = 1, c_counter-1
            allocate(dects(i)%p, source=dect_c(i))
            !dects(i)%p => dect_c(i)
        end do

        do j = 1, a_counter-1
            allocate(dects(j+i-1)%p, source=dect_a(j))
            !dects(j+i-1)%p => dect_a(j)
        end do

        do k = 1, f_counter-1
            allocate(dects(j+i+k-2)%p, source=dect_f(k))
            !dects(j+i+k-2)%p => dect_f(k)
        end do

        do l = 1, cam_counter-1
            allocate(dects(j+i+k+l-3)%p, source=dect_cam(l))
            !dects(j+i+l-3)%p => dect_cam(l)
        end do

        if(.not. allocated(state%historyFilename))state%historyFilename="photPos.obj"

    end subroutine parse_detectors

    subroutine handle_camera(child, dects, counts, context, error, dect_ID, targetValue)
        !! Read in Camera settings and initalise variable
        use detectors,     only : camera
        use sim_state_mod, only : state

        !> Detector table
        type(toml_table), pointer,     intent(in)    :: child
        !> Array of detectors
        type(camera),                  intent(inout) :: dects(:)
        !> Number dectors created
        integer,                       intent(inout) :: counts
        !> Context handle for error reporting.
        type(toml_context),            intent(in)    :: context
        !> Error message
        type(toml_error), allocatable, intent(out)   :: error
        !> Detector ID
        character(len=:), allocatable, intent(in)    :: dect_ID
        !> Target Value used for inverse MCRT
        real(kind=wp), intent(in) :: targetValue

        integer       :: layer, nbins
        real(kind=wp) :: maxval
        type(vector)  :: p1, p2, p3
        logical       :: trackHistory

        p1 = get_vector(child, "p1", default=vector(-1.0, -1.0, -1.0), context=context, error=error)
        p2 = get_vector(child, "p2", default=vector(2.0, 0.0, 0.0), context=context, error=error)
        p3 = get_vector(child, "p3", default=vector(0.0, 2.0, 0.0), context=context, error=error)

        call get_value(child, "layer", layer, 1)
        call get_value(child, "nbins", nbins, 100)
        call get_value(child, "maxval", maxval, 100._wp)
        call get_value(child, "trackHistory", trackHistory, .false.)
        if(trackHistory)state%trackHistory=.true.
#ifdef _OPENMP
        if(trackHistory)then
            call make_error(error, "Track history currently incompatable with OpenMP!", -1)
            return
        end if
#endif
        dects(counts) = camera(p1, p2, p3, layer, nbins, maxval, trackHistory, dect_ID, targetValue)
        counts = counts + 1

    end subroutine handle_camera

    subroutine handle_circle_dect(child, dict, dects, counts, context, error, dect_ID, targetValue)
        !! Read in Circle_detector settings and initalise variable
        use detectors,     only : circle_dect
        use sim_state_mod, only : state

        !> Detector table
        type(toml_table), pointer,     intent(in)    :: child
        !> Dictonary used to store metadata
        type(toml_table),               intent(inout) :: dict
        !> Array of dects
        type(circle_dect),             intent(inout) :: dects(:)
        !> Number dectors created
        integer,                       intent(inout) :: counts
        !> Context handle for error reporting.
        type(toml_context),            intent(in)    :: context
        !> Error message
        type(toml_error), allocatable, intent(out)   :: error
        !> Detector ID
        character(len=:), allocatable, intent(in)    :: dect_ID
        !> Target Value used for inverse MCRT
        real(kind=wp), intent(in) :: targetValue

        integer       :: layer, nbins
        real(kind=wp) :: maxval, radius, acceptAngle
        type(vector)  :: pos, dir
        logical       :: trackHistory



        character(len=8) :: countStr
        character(len=:), allocatable :: dectType

        pos = get_vector(child, "position", context=context, error=error)
        dir = get_vector(child, "direction", default=vector(0.0, 0.0, -1.0), context=context, error=error)
        dir = dir%magnitude()
        call get_value(child, "layer", layer, 1)
        call get_value(child, "radius", radius, 1.0_wp)

        !testing accpetanceAngle
        call get_value(child, "acceptanceAngle", acceptAngle, 90.0_wp)

        call get_value(child, "nbins", nbins, 100)
        call get_value(child, "maxval", maxval, 100._wp)
        call get_value(child, "trackHistory", trackHistory, .false.)
        if(trackHistory)state%trackHistory=.true.
#ifdef _OPENMP
        if(trackHistory)then
            call make_error(error, "Track history currently incompatable with OpenMP!", -1)
            return
        end if
#endif
#ifdef escapeFunction
        !add the detector values to the dictionary

        write(countStr, "(I8)") counts
        dectType = "circle"
        
        call set_value(dict, "dect"//countStr//"type", dectType)
        call set_value(dict, "dect"//countStr//"ID", dect_ID)
        call set_value(dict, "dect"//countStr//"position%x", pos%x)
        call set_value(dict, "dect"//countStr//"position%y", pos%y)
        call set_value(dict, "dect"//countStr//"position%z", pos%z)
        call set_value(dict, "dect"//countStr//"direction%x", dir%x)
        call set_value(dict, "dect"//countStr//"direction%y", dir%y)
        call set_value(dict, "dect"//countStr//"direction%z", dir%z)
        call set_value(dict, "dect"//countStr//"radius", radius)
        call set_value(dict, "dect"//countStr//"acceptanceAngle", acceptAngle)
#endif
        dects(counts) = circle_dect(pos, dir, layer, radius, acceptAngle, nbins, trackHistory, dect_ID, targetValue)
        counts = counts + 1

    end subroutine handle_circle_dect

    subroutine handle_fibre_collection_dect(child, dects, counts, context, error, dect_ID, targetValue)
        !! Read in handle_fibre_collection_dector settings and initalise variable
        use detectors,     only : fibre_dect
        use sim_state_mod, only : state

        !> Detector table
        type(toml_table), pointer,     intent(in)    :: child
        !> Array of dectd
        type(fibre_dect),             intent(inout) :: dects(:)
        !> Number dectors created
        integer,                       intent(inout) :: counts
        !> Context handle for error reporting.
        type(toml_context),            intent(in)    :: context
        !> Error message
        type(toml_error), allocatable, intent(out)   :: error
        !> Detector ID
        character(len=:), allocatable, intent(in)    :: dect_ID
        !> Target Value used for inverse MCRT
        real(kind=wp), intent(in) :: targetValue

        integer       :: layer, nbins
        real(kind=wp) :: maxval
        type(vector)  :: pos, dir
        logical       :: trackHistory

        real(kind=wp) :: focalLength1, focalLength2, f1Aperture, f2Aperture
        real(kind=wp) :: backOffset, frontToPinSep, pinToBackSep
        real(kind=wp) :: pinAperture, acceptAngle, coreDiameter

        pos = get_vector(child, "position", context=context, error=error)
        dir = get_vector(child, "direction", default=vector(0.0, 0.0, -1.0), context=context, error=error)
        dir = dir%magnitude()
        call get_value(child, "focalLength1", focalLength1, 1.0_wp)
        call get_value(child, "focalLength2", focalLength2, 1.0_wp)
        call get_value(child, "f1Aperture", f1Aperture, 1.0_wp)
        call get_value(child, "f2Aperture", f2Aperture, 1.0_wp)
        call get_value(child, "backOffset", backOffset, focalLength2)
        call get_value(child, "frontToPinSep", frontToPinSep, focalLength1)
        call get_value(child, "pinToBackSep", pinToBackSep, focalLength2)
        call get_value(child, "pinAperture", pinAperture, max(f1Aperture, f2Aperture))
        call get_value(child, "acceptanceAngle", acceptAngle, 90.0_wp)
        call get_value(child, "coreDiameter", coreDiameter, 0.01_wp)   
             
        call get_value(child, "layer", layer, 1)
        call get_value(child, "nbins", nbins, 1)
        call get_value(child, "maxval", maxval, 100._wp)
        call get_value(child, "trackHistory", trackHistory, .false.)
        if(trackHistory)state%trackHistory=.true.
#ifdef _OPENMP
        if(trackHistory)then
            call make_error(error, "Track history currently incompatable with OpenMP!", -1)
            return
        end if
#endif
#ifdef escapeFunction
        !add the detector values to the dictionary

        write(countStr, "(I8)") counts
        dectType = "annulus"
        
        call set_value(dict, "dect"//countStr//"type", dectType)
        call set_value(dict, "dect"//countStr//"ID", dect_ID)
        call set_value(dict, "dect"//countStr//"position%x", pos%x)
        call set_value(dict, "dect"//countStr//"position%y", pos%y)
        call set_value(dict, "dect"//countStr//"position%z", pos%z)
        call set_value(dict, "dect"//countStr//"direction%x", dir%x)
        call set_value(dict, "dect"//countStr//"direction%y", dir%y)
        call set_value(dict, "dect"//countStr//"direction%z", dir%z)
        call set_value(dict, "dect"//countStr//"focalLength1", focalLength1)
        call set_value(dict, "dect"//countStr//"focalLength2", focalLength2)
        call set_value(dict, "dect"//countStr//"f1Aperture", f1Aperture)
        call set_value(dict, "dect"//countStr//"f2Aperture", f2Aperture)
        call set_value(dict, "dect"//countStr//"backOffset", backOffset)
        call set_value(dict, "dect"//countStr//"frontToPinSep", frontToPinSep)
        call set_value(dict, "dect"//countStr//"pinToBackSep", pinToBackSep)
        call set_value(dict, "dect"//countStr//"pinAperture", pinAperture)
        call set_value(dict, "dect"//countStr//"acceptanceAngle", acceptAngle)
        call set_value(dict, "dect"//countStr//"coreDiameter", coreDiameter)
#endif
        dects(counts) = fibre_dect(pos, dir, layer, nbins, maxval, trackHistory, focalLength1, & 
                                    focalLength2, f1Aperture, f2Aperture, backOffset, & 
                                    frontToPinSep, pinToBackSep, pinAperture, acceptAngle, coreDiameter, & 
                                    dect_ID, targetValue)
        counts = counts + 1

    end subroutine handle_fibre_collection_dect

    subroutine handle_annulus_dect(child, dict, dects, counts, context, error, dect_ID, targetValue)
        !! Read in Annulus_detector settings and initalise variable
        
        use detectors,     only : annulus_dect
        use sim_state_mod, only : state
        use utils,         only : str

        !> Detector Table
        type(toml_table), pointer,     intent(in)    :: child
        !> Dictonary used to store metadata
        type(toml_table),               intent(inout) :: dict
        !> Array of dects
        type(annulus_dect),            intent(inout) :: dects(:)
        !> Number dectors created
        integer,                       intent(inout) :: counts
        !> Context handle for error reporting.
        type(toml_context),            intent(in)    :: context
        !> Error message
        type(toml_error), allocatable, intent(out)   :: error
        !> Detector ID
        character(len=:), allocatable, intent(in)    :: dect_ID
        !> Target Value used for inverse MCRT
        real(kind=wp), intent(in) :: targetValue


        integer       :: layer, nbins, origin
        real(kind=wp) :: maxval, radius1, radius2, acceptAngle
        type(vector)  :: pos, dir
        logical       :: trackHistory

        character(len=8) :: countStr
        character(len=:), allocatable :: dectType

        pos = get_vector(child, "position", context=context, error=error)
        dir = get_vector(child, "direction", default=vector(0.0, 0.0, -1.0), context=context, error=error)
        call get_value(child, "layer", layer, 1)
        call get_value(child, "radius1", radius1, 0.1_wp)
        call get_value(child, "radius2", radius2, 0.2_wp, origin=origin)
        
        !testing accpetanceAngle
        call get_value(child, "acceptanceAngle", acceptAngle, 90.0_wp)
        
        if(radius2 <= radius1)then
            call make_error(error,&
            context%report("Radii are invalid", origin, &
            "Expected radius2 ("//str(radius2,6)//") > radius 1 ("//str(radius1,6)//")"), -1)
            return
        end if
        
        call get_value(child, "nbins", nbins, 100)
        call get_value(child, "maxval", maxval, 100._wp)
        call get_value(child, "trackHistory", trackHistory, .false.)
        if(trackHistory)state%trackHistory=.true.
#ifdef _OPENMP
        if(trackHistory)then
            call make_error(error, "Track history currently incompatable with OpenMP!", -1)
            return
        end if
#endif
#ifdef escapeFunction
        !add the detector values to the dictionary

        write(countStr, "(I8)") counts
        dectType = "annulus"
        
        call set_value(dict, "dect"//countStr//"type", dectType)
        call set_value(dict, "dect"//countStr//"ID", dect_ID)
        call set_value(dict, "dect"//countStr//"position%x", pos%x)
        call set_value(dict, "dect"//countStr//"position%y", pos%y)
        call set_value(dict, "dect"//countStr//"position%z", pos%z)
        call set_value(dict, "dect"//countStr//"direction%x", dir%x)
        call set_value(dict, "dect"//countStr//"direction%y", dir%y)
        call set_value(dict, "dect"//countStr//"direction%z", dir%z)
        call set_value(dict, "dect"//countStr//"radius1", radius1)
        call set_value(dict, "dect"//countStr//"radius2", radius2)
        call set_value(dict, "dect"//countStr//"acceptanceAngle", acceptAngle)
#endif
        dects(counts) = annulus_dect(pos, dir, layer, radius1, radius2, acceptAngle, nbins, maxval, & 
                                     trackHistory, dect_ID, targetValue)
        counts = counts + 1
    end subroutine handle_annulus_dect
end module parse_detectorsMod