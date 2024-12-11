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

    end subroutine parse_params
    

    subroutine parse_grid(table, dict, error)
    !! parse grid input data
        use sim_state_mod, only : state
        use gridMod,       only : init_grid 
        
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
        state%grid = init_grid(nxg, nyg, nzg, xmax, ymax, zmax)

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
            call get_value(child, "render_geomerty_name", state%rendergeomfile, "geom_render.nrrd")
            call get_value(child, "render_geomerty", state%render_geom, .false.)
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
end module parse_mod