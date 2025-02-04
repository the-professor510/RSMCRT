module parse_mod
!! Module contains parses the input toml config files.
!! See [config](../|page|/config.html) for details of toml input file.
    use tomlf
    use tomlf_error, only : make_error
    use constants, only : wp
    use vector_class

    use parse_detectorsMod
    use parse_sourcesMod
    use parse_geometryMod
    
    implicit none

    private
    public :: parse_params

    contains

    subroutine parse_params(filename, packet, dects, spectrum, dict, error)
        !! entry point for parsing toml file

        use detectors,   only : dect_array
        use photonmod
        use piecewiseMod
        
        !> filename of input toml file
        character(*),      intent(IN)    :: filename
        !> dictionary that stores potential metadata to be saved with simulation output
        type(toml_table),  intent(INOUT) :: dict
        !> some input options set up data in the photon class
        type(photon),      intent(OUT)   :: packet
        !> detector array which is setup during parsing
        type(dect_array), allocatable, intent(out) :: dects(:)
        !> spectrum type which is set up during parsing
        type(spectrum_t), intent(out) :: spectrum
        !> Last error raised during parsing. Unallocated if no error raised. Need to handle this on return from parse_params.
        type(toml_error), allocatable, intent(out) :: error

        type(toml_table), allocatable :: table
        type(toml_context) :: context

        call toml_load(table, trim(filename), context=context, error=error)
        if(allocated(error))return

        call parse_source(table, packet, dict, spectrum, context, error)
        if(allocated(error))return

        call parse_grid(table, dict, error)
        if(allocated(error))return

        call parse_geometry(table, dict, context, error)
        if(allocated(error))return

        call parse_detectors(table, dects, context, error)
        if(allocated(error))return

        call parse_output(table, error)
        if(allocated(error))return

        call parse_simulation(table, error)
        if(allocated(error))return

#ifdef escapeFunction
        call parse_symmetry(table, dict, error)
        if(allocated(error))return
#elif inverseMCRT
        call parse_inverse(table, dict, error)
        if(allocated(error))return
#endif

    end subroutine parse_params
    

    subroutine parse_grid(table, dict, error)
    !! parse grid input data
        use sim_state_mod, only : state
        use gridMod,       only : init_grid_cart 
        
        !> Input Toml table
        type(toml_table),               intent(inout) :: table
        !> Dictonary used to store metadata
        type(toml_table),               intent(inout) :: dict
        !> Error message
        type(toml_error),  allocatable, intent(out)   :: error

        type(toml_table), pointer     :: child
        integer                       :: nxg, nyg, nzg
        real(kind=wp)                 :: xmax, ymax, zmax
        character(len=:), allocatable :: units

        call get_value(table, "grid", child)

        if(associated(child))then
            call get_value(child, "nxg", nxg, 200) !200 is the default 
            !value if it is no stated inside the toml
            call get_value(child, "nyg", nyg, 200)
            call get_value(child, "nzg", nzg, 200)
            call get_value(child, "xmax", xmax, 1.0_wp)
            call get_value(child, "ymax", ymax, 1.0_wp)
            call get_value(child, "zmax", zmax, 1.0_wp)
            call get_value(child, "units", units, "cm")
            call set_value(dict, "units", units)
        else
            call make_error(error, "Need grid table in input param file", -1)
            return
        end if

        !set the grid variable inside the program
        state%grid = init_grid_cart(nxg, nyg, nzg, xmax, ymax, zmax)

    end subroutine parse_grid

    subroutine parse_output(table, error)
        !! parse output file information
        use sim_state_mod, only : state

        !> Input Toml table 
        type(toml_table),              intent(inout) :: table
        !> Error message
        type(toml_error), allocatable, intent(out)   :: error

        type(toml_table), pointer :: child
        type(toml_array), pointer :: children
        integer :: i, nlen

        call get_value(table, "output", child)

        if(associated(child))then
            call get_value(child, "fluence", state%outfile, "fluence.nrrd")
            call get_value(child, "absorb", state%outfile_absorb, "absorb.nrrd")
            call get_value(child, "render_geometry_name", state%rendergeomfile, "geom_render.nrrd")
            call get_value(child, "render_geometry", state%render_geom, .false.)
            call get_value(child, "render_source_name", state%rendersourcefile, "source_render.nrrd")
            call get_value(child, "render_source", state%render_source, .false.)

            call get_value(child, "render_size", children, requested=.false.)
            if(associated(children))then
                nlen = len(children)
                if(nlen < 3)then
                    call make_error(error, "Need a vector of size 3 for render_size.", -1)
                    return
                end if
                do i = 1, len(children)
                    call get_value(children, i, state%render_size(i))
                end do
            else
                state%render_size = [200, 200, 200]
            end if

            call get_value(child, "overwrite", state%overwrite, .false.)
        else
            call make_error(error, "Need output table in input param file", -1)
            return
        end if

    end subroutine parse_output

    subroutine parse_simulation(table, error)
        !! parse simulation information
        use sim_state_mod, only : state

        !> Input Toml table 
        type(toml_table),              intent(inout) :: table
        !> Error message
        type(toml_error), allocatable, intent(out)   :: error

        type(toml_table), pointer :: child

        call get_value(table, "simulation", child)

        if(associated(child))then
            call get_value(child, "iseed", state%iseed, 123456789)
            call get_value(child, "tev", state%tev, .false.)
            call get_value(child, "absorb", state%absorb, .false.)

            call get_value(child, "load_checkpoint", state%loadckpt, .false.)
            call get_value(child, "checkpoint_file", state%ckptfile, "check.ckpt")
            call get_value(child, "checkpoint_every_n", state%ckptfreq, 1000000)

        else
            call make_error(error, "Need simulation table in input param file", -1)
            return
        end if

    end subroutine parse_simulation

    subroutine parse_symmetry(table, dict, error)
        !! parse symmetry information, only used in the computation of the escape function
        use constants,     only: TWOPI
        use sim_state_mod, only : state
        use gridMod,       only : init_grid_cart, init_grid_cyl 
        use vector_class,  only : vector, magnitude

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

    subroutine parse_inverse(table, dict, error)
        !! parse symmetry information, only used in the computation of the escape function
        use sim_state_mod, only : state
        use gridMod,       only : init_grid_cart, init_grid_cyl 
        use vector_class,  only : vector, magnitude

        !> Input Toml table 
        type(toml_table),              intent(inout) :: table
        !> Error message
        type(toml_error), allocatable, intent(out)   :: error
        !> Dictonary used to store metadata
        type(toml_table),               intent(inout) :: dict

        type(toml_table), pointer :: child
        type(toml_array), pointer :: children

        integer :: maxNumSteps, layer
        real(kind = wp) :: maxStepSize, gradStepSize, accuracy
        logical :: findmua, findmus, findg, findn

        call get_value(table, "inverse", child)

        if(associated(child))then
            
            call get_value(child, "maxStepSize", maxStepSize, 1.0_wp)
            call set_value(dict, "maxStepSize", maxStepSize)

            call get_value(child, "gradStepSize", gradStepSize, 0.0001_wp)
            call set_value(dict, "gradStepSize", gradStepSize)

            call get_value(child, "accuracy", accuracy, 0.01_wp)
            call set_value(dict, "accuracy", accuracy)

            call get_value(child, "maxNumSteps", maxNumSteps, 1000)
            call set_value(dict, "maxNumSteps", maxNumSteps)
            
            call get_value(child, "Findmua", findmua, .false.)
            call set_value(dict, "Findmua", findmua)

            call get_value(child, "Findmus", findmus, .false.)
            call set_value(dict, "Findmus", findmus)

            call get_value(child, "Findg", findg, .false.)
            call set_value(dict, "Findg", findg)

            call get_value(child, "Findn", findn, .false.)
            call set_value(dict, "Findn", findn)

            call get_value(child, "layer", layer, -985464082)
            if(layer /= -985464082) then
               call set_value(dict, "inverseLayer", layer)
            else 
                call make_error(error, "Must specifiy a layer in inverse table", -1)
                return
            end if
        else
            call make_error(error, "Need inverse table in input param file", -1)
            return
        end if

    end subroutine parse_inverse

end module parse_mod