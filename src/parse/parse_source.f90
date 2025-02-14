module parse_sourcesMod
    !! routine to parse the source table from the input Toml file.
    use constants, only : wp
    use parse_HelpersMod
    use parse_SpectrumMod
    use vector_class

    use tomlf

    implicit none
    
    private
    public :: parse_source

contains
    
    subroutine parse_source(table, packet, dict, spectrum, context, error)
    !! Parse sources
    !! any updates here MUST be reflected in docs/config.md
        use sim_state_mod, only : state
        use vector_class,  only : length, magnitude
        use photonmod
        use piecewiseMod
        use tomlf_error
        
        !> Input Toml table
        type(toml_table),  intent(inout) :: table
        !> Dictonary used to store metadata
        type(toml_table),  intent(inout) :: dict
        !> Photon packet. Used to store information to save computation
        type(photon),      intent(out)   :: packet
        !> Spectrum type.
        type(spectrum_t),  intent(out)   :: spectrum
        !> Context handle for error reporting
        type(toml_context) :: context
        !> Error message
        type(toml_error), allocatable, intent(out) :: error

        type(toml_table), pointer :: child
        type(toml_array), pointer :: children

        type(vector) :: poss, dirr, rotation
        real(kind=wp) :: dir(3), pos(3), rot(3), corners(3, 3), radius, focalLength, rlo, rhi, sigma, beam_size, tempCoord
        integer :: i, nlen, origin
        character(len=1) :: axis(3)
        character(len=:), allocatable :: direction, annulus_type, focus_type

        axis = ["x", "y", "z"]
        pos = 0._wp
        dir = 0._wp
        tempCoord = 0._wp
        rot = 0._wp
        corners = reshape((/ -1._wp, -1._wp, 1._wp, &
                              2._wp,  0._wp, 0._wp, &
                              0._wp,  2._wp, 0._wp /), &
                           shape(corners), order=[2, 1])

        call get_value(table, "source", child, requested=.false.)
        if(associated(child))then
            call get_value(child, "name", state%source, "point")
            call get_value(child, "nphotons", state%nphotons, 1000000)

            if(state%source /= "uniform")then
                poss = get_vector(child, "position", error, context)
                if(allocated(error))return
            end if

            children => null()

            if(state%source /= "uniform" .and. state%source /= "point" .and. &
                state%source /= "circular" .and. state%source /= "pencil")then
                call get_value(child, "rotation", children, requested=.false., origin=origin)
                if(associated(children))then
                    nlen = len(children)
                    if(nlen /= 3)then
                        call make_error(error, &
                        context%report("Need a matrix row for points", origin, "expected matrix row of size 3"), -1)
                        return
                    end if
                    do i = 1, len(children)
                        call get_value(children, i, rot(i))
                        call set_value(dict, "rotation%"//axis(i), rot(i))
                    end do
                else
                    call make_error(error, &
                    context%report("Source requires rotation variable", origin, "expected rotation variable"), -1)
                    return
                end if

                
                rotation%x = rot(1)
                rotation%y = rot(2)
                rotation%z = rot(3)
                if(length(rotation) < 1e-8_wp) then
                    print'(a)',context%report(&
                    "Need to specify rotation that has length greater than 0.0", origin, level=toml_level%warning)
                    return
                else 
                    rotation = rotation%magnitude()
                end if
            end if
            
            dirr = get_vector(child, "direction", error, context)
            if(allocated(error))then
                call get_value(child, "direction", direction, origin=origin)
                if(allocated(direction))then
                    if(state%source == "point")then
                        print'(a)',context%report(&
                        "Point source needs no direction!!", origin, level=toml_level%warning)
                    end if

                    if(state%source == "annulus")then
                        print'(a)',context%report(&
                        "Annulus source needs no direction!!", origin, level=toml_level%warning)
                    end if

                    if(state%source == "focus")then
                        print'(a)',context%report(&
                        "Focus source needs no direction!!", origin, level=toml_level%warning)
                    end if
    
                    select case(direction)
                    case("x")
                        dirr = vector(1._wp, 0._wp, 0._wp)
                    case("-x")
                        dirr = vector(-1._wp, 0._wp, 0._wp)                
                    case("y")
                        dirr = vector(0._wp, 1._wp, 0._wp)
                    case("-y")
                        dirr = vector(0._wp, -1._wp, 0._wp)
                    case("z")
                        dirr = vector(0._wp, 0._wp, 1._wp)
                    case("-z")
                        dirr = vector(0._wp, 0._wp, -1._wp)
                    case default
                        call make_error(error, context%report("Direction needs a cardinal direction i.e x, y, or z", origin, &
                                                "Expected cardinal direction"), -1)
                        return 
                    end select
                elseif(state%source /= "point" .and. state%source /= "annulus" .and. state%source /= "focus")then
                    call make_error(error, context%report("Need to specify direction for source type!", origin, &
                                              "No direction specified"), -1)
                    return
                end if
            else
                if(state%source == "point")then
                    print'(a)',context%report(&
                    "Point source needs no direction!!", origin, level=toml_level%warning)
                end if
                if(state%source == "annulus")then
                    print'(a)',context%report(&
                    "Annulus source needs no direction!!", origin, level=toml_level%warning)
                end if

                if(state%source == "focus")then
                    print'(a)',context%report(&
                    "Focus source needs no direction!!", origin, level=toml_level%warning)
                end if
                return
            end if


            children => null()
            
            call get_value(child, "point1", children, requested=.false., origin=origin)
            if(associated(children))then
                nlen = len(children)
                if(nlen < 3)then
                    call make_error(error, &
                    context%report("Need a matrix row for points", origin, "expected matrix row of size 3"), -1)
                    return
                end if
                do i = 1, len(children)
                    call get_value(children, i, corners(i, 1))
                    call set_value(dict, "pos1%"//axis(i), corners(i,1))
                end do
            else
                if(state%source == "uniform")then
                    call make_error(error, &
                    context%report("Uniform source requires point1 variable", origin, "expected point1 variable"), -1)
                    return
                end if
            end if

            call get_value(child, "point2", children, requested=.false., origin=origin)
            if(associated(children))then
                nlen = len(children)
                if(nlen < 3)then
                    call make_error(error, &
                    context%report("Need a matrix row for points", origin, "expected matrix row of size 3"), -1)
                    return
                end if
                do i = 1, len(children)
                    call get_value(children, i, corners(i, 2))
                    call set_value(dict, "pos2%"//axis(i), corners(i,2))
                end do
            else
                if(state%source == "uniform")then
                    call make_error(error, &
                    context%report("Uniform source requires point2 variable", origin, "expected point2 variable"), -1)
                    return
                end if
            end if

            call get_value(child, "point3", children, requested=.false., origin=origin)
            if(associated(children))then
                nlen = len(children)
                if(nlen < 3)then
                    call make_error(error, &
                    context%report("Need a matrix row for points", origin, "expected matrix row of size 3"), -1)
                    return
                end if
                do i = 1, len(children)
                    call get_value(children, i, corners(i, 3))
                    call set_value(dict, "pos3%"//axis(i), corners(i,3))
                end do
            else
                if(state%source == "uniform")then
                    call make_error(error, &
                    context%report("Uniform source requires point3 variable", origin, "expected point3 variable"), -1)
                    return
                end if
            end if
            call get_value(child, "radius", radius, 0.5_wp)
            call set_value(dict, "radius", radius)

            ! parameters for annulus beam and focus beam types
            call get_value(child, "focalLength", focalLength, 1._wp)
            call set_value(dict, "focalLength", focalLength)

            call get_value(child, "rhi", rhi, 0.6_wp)
            call set_value(dict, "rhi", rhi)

            call get_value(child, "rlo", rlo, 0.5_wp)
            call set_value(dict, "rlo", rlo)

            call get_value(child, "sigma", sigma, 0.04_wp)
            call set_value(dict, "sigma", sigma)

            call get_value(child, "annulus_type", annulus_type, "gaussian")
            call set_value(dict, "annulus_type", annulus_type)

            call get_value(child, "focus_type", focus_type, "gaussian")
            call set_value(dict, "focus_type", focus_type)

            call get_value(child, "beam_size", beam_size, 0.5_wp)
            call set_value(dict, "beam_size", beam_size)

            ! parse spectrum
            call parse_spectrum(child, spectrum, dict, context, error)
            if(allocated(error))return
        else
            call make_error(error, context%report("Simulation needs Source table", origin, "Missing source table"), -1)
            return
        end if

        call set_photon(poss, dirr)
        packet = photon(state%source)
        packet%pos = poss
        packet%nxp = dirr%x
        packet%nyp = dirr%y
        packet%nzp = dirr%z

    end subroutine parse_source
end module parse_sourcesMod