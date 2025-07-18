module parse_symmetryMod
    !! routine to parse the source table from the input Toml file.
    use constants, only : wp
    use parse_HelpersMod
    use parse_SpectrumMod
    use vector_class

    use tomlf

    implicit none
    
    private
    public :: parse_symmetry

contains
        subroutine parse_symmetry(table, dict, error)
        !! parse symmetry information, only used in the computation of the escape function
        use constants,     only: TWOPI
        use sim_state_mod, only : state
        use gridMod,       only : init_grid_cart, init_grid_cyl 
        use vector_class,  only : vector, magnitude
        use tomlf
        use tomlf_error

        !> Input Toml table 
        type(toml_table),              intent(inout) :: table
        !> Error message
        type(toml_error), allocatable, intent(out)   :: error
        !> Dictonary used to store metadata
        type(toml_table),               intent(inout) :: dict

        type(toml_table), pointer :: child
        type(toml_array), pointer :: children

        character(len=:), allocatable :: symmetryType
        integer :: i, nlen
        integer :: nxrg, nytg, nzg
        integer :: escapenphotons
        real(kind=wp) :: xrmax, ytmax, zmax, rotation
        type(vector) :: pos, dir

        pos = vector(0._wp,0._wp,0._wp)
        dir = vector(0._wp,0._wp,1._wp)

        state%source = "point"


        call get_value(table, "symmetry", child)

        if(associated(child))then
            !symmetry, used to reduce the computation time of the escape function
            call get_value(child, "symmetryType", symmetryType, "none")
            call set_value(dict, "symmetryType", symmetryType)

            call get_value(child, "escapenphotons", escapenphotons, 100000)
            state%nphotons = escapenphotons

            call get_value(child, "GridSize", children, requested=.false.)
            if(associated(children))then
                nlen = len(children)
                if(nlen /= 3)then
                    call make_error(error, "Need a vector of size 3 for symmetry grid size.", -1)
                    return
                end if
                call get_value(children, 1, nxrg)
                call get_value(children, 2, nytg)
                call get_value(children, 3, nzg)
            else
                nxrg = 10
                nytg = 10
                nzg = 10
            end if

            call get_value(child, "maxValues", children, requested=.false.)
            if(associated(children))then
                nlen = len(children)
                if(nlen /= 3)then
                    call make_error(error, "Need a vector of size 3 for symmetry max values.", -1)
                    return
                end if
                call get_value(children, 1, xrmax)
                call get_value(children, 2, ytmax)
                call get_value(children, 3, zmax)
            else
                xrmax = 1.0
                ytmax = 1.0
                zmax = 1.0
            end if

            call get_value(child, "position", children, requested=.false.)
            if(associated(children))then
                nlen = len(children)
                if(nlen /= 3)then
                    call make_error(error, "Need a vector of size 3 for symmetry position.", -1)
                    return
                end if
                call get_value(children, 1, pos%x)
                call get_value(children, 2, pos%y)
                call get_value(children, 3, pos%z)
            else
                pos = vector(0._wp, 0._wp, 0._wp)
            end if

            call get_value(child, "direction", children, requested=.false.)
            if(associated(children))then
                nlen = len(children)
                if(nlen /= 3)then
                    call make_error(error, "Need a vector of size 3 for symmetry position.", -1)
                    return
                end if
                call get_value(children, 1, dir%x)
                call get_value(children, 2, dir%y)
                call get_value(children, 3, dir%z)
            else
                dir = vector(0._wp, 0._wp, 1._wp)
            end if

            call get_value(child, "rotation", rotation, 0._wp)
            if (rotation < 0.0_wp .or. rotation >= 360.0_wp ) then
                call make_error(error, "Must specifcy a rotation for symmetry that is between 0.0 and 360.0, inclusive of 0.0")
                return
            end if

            if (dir%x == 0._wp .and. dir%y == 0._wp .and. dir%z == 0._wp) then
                call make_error(error, "Must specify a non-zero direction for symmetry")
                return
            end if

            dir = dir%magnitude()

            if (symmetryType == "none" .or. symmetryType == "prism" .or. symmetryType == "flipped" & 
                .or. symmetryType == "uniformSlab") then
                state%symGridPos = pos
                state%symGridDir = dir
                state%symGridRot = rotation
                state%symmetryEscapeCartGrid = init_grid_cart(nxrg, nytg, nzg, xrmax, ytmax, zmax)
            else if (symmetryType == "noneRotational" .or. symmetryType == "360rotational") then
                state%symGridPos = pos
                state%symGridDir = dir
                state%symGridRot = rotation
                state%symmetryEscapeCylGrid = init_grid_cyl(nxrg, nytg, nzg, xrmax, ytmax, zmax)
            else if (symmetryType == "adjoint") then 

                !no symmetry grid is required
            else
                call make_error(error, "Unrecognised symmetry type")
                return
            end if
        else 
            !set the symmetry type to none, and set the other variables to their default values
            symmetryType = "none"
            call set_value(dict, "symmetryType", symmetryType)

            !set default number of photons to run
            state%nphotons = 100000

            !set default size of symmetry grid
            nxrg = 10
            nytg = 10
            nzg = 10

            !set max values of symmetry grid
            xrmax = 1.0
            ytmax = 1.0
            zmax = 1.0

            !set default position of symmetry grid
            state%symGridPos = pos

            !set default direction of the symmetry grid
            state%symGridDir = dir

            !define the default escape cart grid
            state%symmetryEscapeCartGrid = init_grid_cart(nxrg, nytg, nzg, xrmax, ytmax, zmax)
        end if
    end subroutine parse_symmetry

end module parse_symmetryMod