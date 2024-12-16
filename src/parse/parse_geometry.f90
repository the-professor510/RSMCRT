module parse_geometryMod
    !! routine to parse the source table from the input Toml file.
    use constants, only : wp
    use parse_HelpersMod
    use parse_SpectrumMod
    use vector_class

    use tomlf

    implicit none
    
    private
    public :: parse_geometry

contains
    
    subroutine parse_geometry(table, dict, context, error)
    !! Parse Geometry
    !! any updates here MUST be reflected in docs/config.md
        use sim_state_mod, only : state
        use photonmod
        use piecewiseMod
        use tomlf_error
        
        !> Input Toml table
        type(toml_table),  intent(inout) :: table
        !> Dictonary used to store metadata
        type(toml_table),  intent(inout) :: dict
        !> Context handle for error reporting
        type(toml_context) :: context
        !> Error message
        type(toml_error), allocatable, intent(out) :: error
        
        type(toml_table), pointer :: child
        type(toml_array), pointer :: children
        
        real(kind=wp)             :: muaTemp, musTemp, murTemp, hggTemp, nTemp, tempCoord, tempLength, sphereRadius
        integer                   :: num_spheres, i, nlen, numOptProp, origin
        character(4) :: string 

        real(kind=wp) :: tau, musb, muab, musc, muac, hgg
        real(kind=wp) :: bottomSphereRad, topSphereRad, SphereSep, ShellThickness, YolkRadius

        child => null()
        call get_value(table, "geometry", child)
        if(associated(child))then
            call get_value(child, "geom_name", state%experiment, "sphere")

            !Left over from the old setupGeometries
            call get_value(child, "tau", tau, 10._wp)
            call set_value(dict, "tau", tau)

            call get_value(child, "num_spheres", num_spheres, 10)
            call set_value(dict, "num_spheres", num_spheres)

            call get_value(child, "musb", musb, 0.0_wp)
            call set_value(dict, "musb", musb)
            call get_value(child, "muab", muab, 0.01_wp)
            call set_value(dict, "muab", muab)
            call get_value(child, "musc", musc, 0.0_wp)
            call set_value(dict, "musc", musc)
            call get_value(child, "muac", muac, 0.01_wp)
            call set_value(dict, "muac", muac)
            call get_value(child, "hgga", hgg, 0.7_wp)
            call set_value(dict, "hgga", hgg)
            !end of old stuff

            ! "egg", "sphere", "box"

            !want to take in mus, mua, mur, and hgg for each different medium
            !for the sphere i want to take in a central position and the radius
            !for the box i want to take in a central position and the x,y,z of it
            !for the egg i want to take in a central position and the r1, r2, h, and thickness of shell
            !for anything else i will want to do similar


            call get_value(child, "numOptProp", numOptProp, 1)
            call set_value(dict, "numOptProp", numOptProp)
            if(numOptProp < 1) then
                call make_error(error, &
                context%report("Need to set an integer value of at least one or greater for numOptProp", origin, &
                                "numOptProp Incorrectly Specified"), -1)
            else if ((state%experiment == "sphere") .and. (numOptProp .NE. 1)) then
                call make_error(error, &
                context%report("For geometry of sphere must set numOptProp to one", origin, &
                                "numOptProp Incorrectly Specified"), -1)
            else if ((state%experiment == "box") .and. (numOptProp .NE. 1)) then
                call make_error(error, &
                context%report("For geometry of box must set numOptProp to one", origin, &
                                "numOptProp Incorrectly Specified"), -1)
            else if ((state%experiment == "egg") .and. (numOptProp .NE. 3)) then
                call make_error(error, &
                context%report("For geometry of egg must set numOptProp to three", origin, &
                                "numOptProp Incorrectly Specified"), -1)
            end if

            children => null()
            
            call get_value(child, "mua", children, requested=.false., origin=origin)
            if(associated(children))then
                nlen = len(children)
                if ((nlen .EQ. numOptProp)) then
                    do i = 1, numOptProp
                        write(string,'(I4)') i
                        call get_value(children, i, muaTemp)
                        call set_value(dict, "mua%"//string, muaTemp)
                    end do
                else
                    call make_error(error, &
                        context%report("length of mua must be equal to numOptProp", origin, &
                                "mua Incorrectly Specified"), -1)
                end if
            else
                do i = 1, numOptProp
                    muaTemp = 0.0_wp
                    write(string,'(I4)') i
                    call set_value(dict, "mua%"//string, muaTemp)
                end do 
            end if

            call get_value(child, "mus", children, requested=.false., origin=origin)
            if(associated(children))then
                nlen = len(children)
                if ((nlen .EQ. numOptProp)) then
                    do i = 1, numOptProp
                        write(string,'(I4)') i
                        call get_value(children, i, musTemp)
                        call set_value(dict, "mus%"//string, musTemp)
                    end do
                else
                    call make_error(error, &
                        context%report("length of mus must be equal to numOptProp", origin, &
                                "mus Incorrectly Specified"), -1)
                end if
            else
                do i = 1, numOptProp
                    musTemp = 1._wp
                    write(string,'(I4)') i
                    call set_value(dict, "mus%"//string, musTemp)
                end do 
            end if

            call get_value(child, "mur", children, requested=.false., origin=origin)
            if(associated(children))then
                nlen = len(children)
                if ((nlen .EQ. numOptProp)) then
                    do i = 1, numOptProp
                        write(string,'(I4)') i
                        call get_value(children, i, murTemp)
                        call set_value(dict, "mur%"//string, murTemp)
                    end do
                else
                    call make_error(error, &
                        context%report("length of mur must be equal to numOptProp", origin, &
                                "mur Incorrectly Specified"), -1)
                end if
            else
                do i = 1, numOptProp
                    murTemp = 0.0_wp
                    write(string,'(I4)') i
                    call set_value(dict, "mur%"//string, murTemp)
                end do 
            end if

            call get_value(child, "hgg", children, requested=.false., origin=origin)
            if(associated(children))then
                nlen = len(children)
                if ((nlen .EQ. numOptProp)) then
                    do i = 1, numOptProp
                        write(string,'(I4)') i
                        call get_value(children, i, hggTemp)
                        call set_value(dict, "hgg%"//string, hggTemp)
                    end do
                else
                    call make_error(error, &
                        context%report("length of hgg must be equal to numOptProp", origin, &
                                "hgg Incorrectly Specified"), -1)
                end if
            else
                do i = 1, numOptProp
                    hggTemp = 0._wp
                    write(string,'(I4)') i
                    call set_value(dict, "hgg%"//string, hggTemp)
                end do 
            end if

            call get_value(child, "n", children, requested=.false., origin=origin)
            if(associated(children))then
                nlen = len(children)
                if ((nlen .EQ. numOptProp)) then
                    do i = 1, numOptProp
                        write(string,'(I4)') i
                        call get_value(children, i, nTemp)
                        call set_value(dict, "n%"//string, nTemp)
                    end do
                else
                    call make_error(error, &
                        context%report("length of n must be equal to numOptProp", origin, &
                                "n Incorrectly Specified"), -1)
                end if
            else
                do i = 1, numOptProp
                    nTemp = 1._wp
                    write(string,'(I4)') i
                    call set_value(dict, "n%"//string, nTemp)
                end do 
            end if

            call get_value(child, "position", children, requested=.false., origin=origin)
            if(associated(children) )then
                nlen = len(children)
                if(nlen .NE. 3)then
                    call make_error(error, &
                    context%report("Need a matrix row for points", origin, "expected matrix row of size 3"), -1)
                    return
                end if
                do i = 1, len(children)
                    write(string,'(I4)') i
                    call get_value(children, i, tempCoord)
                    call set_value(dict, "position%"//string, tempCoord)
                end do
            else
                do i = 1, 3
                    tempCoord = 0.0_wp
                    write(string,'(I4)') i
                    call set_value(dict, "position%"//string, tempCoord)
                end do
            end if

            call get_value(child, "boundingBox", children, requested=.false., origin=origin)
            if(associated(children))then
                nlen = len(children)
                if(nlen .NE. 3)then
                    call make_error(error, &
                    context%report("Need a matrix row for points", origin, "expected matrix row of size 3"), -1)
                    return
                end if
                do i = 1, len(children)
                    write(string,'(I4)') i
                    call get_value(children, i, tempLength)
                    call set_value(dict, "boundinglength%"//string, tempLength)
                end do
            else
                do i = 1, 3
                    tempLength = 2.0_wp
                    write(string,'(I4)') i
                    call set_value(dict, "boundinglength%"//string, tempLength)
                end do
            end if

            if (state%experiment == "sphere") then 
                call get_value(child, "sphereRadius", sphereRadius, 1._wp)
                call set_value(dict, "sphereRadius", sphereRadius) 
            end if


            if (state%experiment == "box") then 
                call get_value(child, "BoxDimensions", children, requested=.false., origin=origin)
                if(associated(children))then
                    nlen = len(children)
                    if(nlen .NE. 3)then
                        call make_error(error, &
                        context%report("Need a matrix row for points", origin, "expected matrix row of size 3"), -1)
                        return
                    end if
                    do i = 1, len(children)
                        write(string,'(I4)') i
                        call get_value(children, i, tempLength)
                        call set_value(dict, "BoxDimensions%"//string, tempLength)
                    end do
                else
                    do i = 1, 3
                        tempLength = 1.0_wp
                        write(string,'(I4)') i
                        call set_value(dict, "BoxDimensions%"//string, tempLength)
                    end do
                end if
            end if

            if (state%experiment == "egg") then
                call get_value(child, "BottomSphereRadius", bottomSphereRad, 3._wp)
                call set_value(dict, "BottomSphereRadius", bottomSphereRad)
                call get_value(child, "TopSphereRadius", topSphereRad, 3._wp * sqrt(2._wp - sqrt(2._wp)))
                call set_value(dict, "TopSphereRadius", topSphereRad)
                call get_value(child, "SphereSep", SphereSep, 3._wp * sqrt(2._wp - sqrt(2._wp)))
                call set_value(dict, "SphereSep", SphereSep)
                call get_value(child, "ShellThickness", ShellThickness, 0.05_wp)
                call set_value(dict, "ShellThickness", ShellThickness)
                call get_value(child, "YolkRadius", YolkRadius, 1.5_wp)
                call set_value(dict, "YolkRadius", YolkRadius)
            end if

        else
            call make_error(error, "Need geometry table in input param file", -1)
            return
        end if


        

    end subroutine parse_geometry
end module parse_geometryMod