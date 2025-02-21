module photonMod

!! This source file contains the photon type, all the photon launch routines for different light sources, and the scattering code.

!! Below are the current types of light sources available. Check [here](config.md) for parameters needed for each light source.
!!   
!! - uniform
!! - pencil
!! - annulus
!! - focus
!! - point
!! - circular
!! - SLM (2D image source)
!! - double slit
!! - square aperture

    use constants, only : wp
    use vector_class
    use random, only : seq

    implicit none
    
    !> photon class
    type :: photon
        !> postion of photon packet in cm. (0,0,0) is the center of the grid.
        type(vector)  :: pos
        !> direction vectors
        real(kind=wp) :: nxp, nyp, nzp
        !> direction cosines
        real(kind=wp) :: sint, cost, sinp, cosp, phi
        !> Wavelength of the packet
        real(kind=wp) :: wavelength
        !> Current phase of the packet
        real(kind=wp) :: phase
        !> \[\frac{2\pi}{\lambda}\]. Used to save computational time
        real(kind=wp) :: fact
        !> Energy of the packet. TODO
        real(kind=wp) :: energy
        !> grid cell position 
        integer       :: xcell, ycell, zcell
        !> photon alive flag
        logical       :: tflag
        !> ID of the SDF the packet is in
        integer       :: layer
        !> Thread ID of the packet
        integer       :: id
        !> Debug data. Number of SDF evals
        integer       :: cnts, bounces
        !> used if photon packet weights are used
        real(kind=wp) :: weight, step!, L
        !> emission routine
        procedure(generic_emit), pointer :: emit => null()
        contains
            !> scattering routine
            procedure :: scatter => scatter
    end type photon

    interface photon
        !> assign the emission function to the photon object
        module procedure init_source
        !> intialise the photon class
        module procedure init_photon
    end interface photon

    abstract interface
        subroutine generic_emit(this, spectrum, dict, seqs)
            
            use tomlf, only : toml_table, get_value
            use random, only : seq
            use piecewiseMod
            import :: photon
            class(photon) :: this
            type(spectrum_t), intent(in) :: spectrum
            type(toml_table), optional, intent(inout) :: dict
            type(seq), optional, intent(inout) :: seqs(2)

        end subroutine generic_emit
    end interface

    !> used to save some computation time
    type(photon) :: photon_origin

    private
    public :: photon, init_source, set_photon

    contains
        
        subroutine set_photon(pos, dir)

            type(vector), intent(in) :: pos, dir

            photon_origin%pos = pos
            photon_origin%nxp = dir%x
            photon_origin%nyp = dir%y
            photon_origin%nzp = dir%z

        end subroutine set_photon

        type(photon) function init_photon(val)
        !! set up all the variables in the photon object
            !> value to assing to variables
            real(kind=wp), intent(in) :: val

            init_photon%pos = vector(val, val, val)
            init_photon%nxp = val
            init_photon%nyp = val
            init_photon%nzp = val
            init_photon%sint = val
            init_photon%cost = val
            init_photon%sinp = val
            init_photon%cosp = val
            init_photon%phi = val
            init_photon%wavelength = val
            init_photon%energy = val
            init_photon%fact = val
            init_photon%zcell = int(val)
            init_photon%ycell = int(val)
            init_photon%zcell = int(val)
            init_photon%tflag = .true.
            init_photon%layer = int(val)
            init_photon%id = int(val)
            init_photon%cnts = int(val)
            init_photon%bounces = int(val)
            init_photon%weight = val
            init_photon%step = val 

        end function init_photon

        type(photon) function init_source(choice)
        !! Bind emission function to photon object
            !> Name of light source to use
            character(*), intent(IN) :: choice

            if(choice == "uniform")then
                init_source%emit => uniform
            elseif(choice == "pencil")then
                init_source%emit => pencil
            elseif(choice == "dslit")then
                init_source%emit => dslit
            elseif(choice == "aperture")then
                init_source%emit => aperture
            elseif(choice == "annulus")then
                init_source%emit => annulus
            elseif(choice == "focus")then
                init_source%emit => focus
            elseif(choice == "point")then
                init_source%emit => point
            elseif(choice == "circular")then
                init_source%emit => circular
            elseif(choice == "slm")then
                init_source%emit => slm
            else
                error stop "No such source!"
            end if


        end function init_source

        subroutine slm(this, spectrum, dict, seqs)
            !! image source
            ! TODO fix hardcoded size of 200
            ! investigate +2 error
            use piecewiseMod
            use tomlf,         only : toml_table, get_value
            use random,        only : ran2, seq
            use sim_state_mod, only : state
            use constants,     only : TWOPI

            class(photon) :: this
            !> Input image to sample position from
            type(spectrum_t),           intent(in)    :: spectrum
            !> Metadata dictionary
            type(toml_table), optional, intent(inout) :: dict
            !> random numbers from a sequence. Quasi-Monte Carlo
            type(seq),        optional, intent(inout) :: seqs(2)
            
            integer :: cell(3)
            real(kind=wp) :: x, y

            this%pos = photon_origin%pos
            call spectrum%p%sample(x, y)
            this%pos%x = (x-100) / (state%grid%nxg / (2.*state%grid%xmax))
            this%pos%y = (y-100) / (state%grid%nyg / (2.*state%grid%ymax))

            this%nxp = photon_origin%nxp
            this%nyp = photon_origin%nyp
            this%nzp = photon_origin%nzp

            this%phi  = atan2(this%nyp, this%nxp)
            this%cosp = cos(this%phi)
            this%sinp = sin(this%phi)
            this%cost = this%nzp
            this%sint = sqrt(1._wp - this%cost**2)  

            this%tflag  = .false.
            this%cnts   = 0
            this%bounces = 0
            this%layer  = 1
            this%weight = 1.0_wp

            this%phase = 0.0_wp
            this%wavelength = 500.e-9_wp
            this%energy = 1._wp
            this%fact = TWOPI/(this%wavelength)

            ! Linear Grid 
            cell = state%grid%get_voxel(this%pos)
            this%xcell = cell(1)
            this%ycell = cell(2)
            this%zcell = cell(3)

        end subroutine slm

        subroutine circular(this, spectrum, dict, seqs)
        !! circular source

            use sim_state_mod, only : state
            use random,        only : ran2, seq
            use constants,     only : twoPI
            use tomlf,         only : toml_table, get_value
            use sdfHelpers,    only : rotationAlign, translate
            use mat_class,     only : invert
            use vector_class
            use piecewiseMod

            class(photon) :: this
            type(spectrum_t), intent(in) :: spectrum
            type(toml_table), optional, intent(inout) :: dict
            type(seq), optional, intent(inout) :: seqs(2)

            type(vector) :: a, b
            integer :: cell(3)
            real(kind=wp) :: t(4,4), radius, r, theta, tmp

            this%nxp = photon_origin%nxp
            this%nyp = photon_origin%nyp
            this%nzp = photon_origin%nzp

            call get_value(dict, "radius", radius)

            ! https://math.stackexchange.com/a/1681815
            r = radius * sqrt(ran2())
            theta = ran2() * TWOPI
            
            !set inital vector from which the source points
            a = vector(1._wp, 0._wp, 0._wp)
            a = a%magnitude()
            !set vector to rotate to. User defined.
            b = vector(this%nxp, this%nyp, this%nzp)
            b = b%magnitude()
            
            ! method fails if below condition is true. So change a vector to point down x-axis
            if(abs(a) == abs(b))then
                a = vector(0._wp, 0._wp, 1._wp)
                a = a%magnitude()
                this%pos = vector(r * cos(theta), r * sin(theta), 0._wp)
            else
                this%pos = vector(0._wp, r * cos(theta), r * sin(theta))
            end if

            ! get rotation matrix
            t = rotationAlign(a, b)
            ! get translation matrix
            t = matmul(t, invert(translate(photon_origin%pos)))
            ! transform point
            this%pos = this%pos .dot. t
            this%pos%x = - this%pos%x
            this%pos%y = - this%pos%y
            this%pos%z = - this%pos%z

            if (this%pos%x == -state%grid%xmax) then
                this%pos%x = this%pos%x + 7.9e-7_wp
            else if (this%pos%x == state%grid%xmax) then
                this%pos%x = this%pos%x - 7.9e-7_wp
            end if
            if (this%pos%y == -state%grid%ymax) then
                this%pos%y = this%pos%y + 7.9e-7_wp
            else if (this%pos%y == state%grid%ymax) then
                this%pos%y = this%pos%y - 7.9e-7_wp
            end if
            if (this%pos%z == -state%grid%zmax) then
                this%pos%z = this%pos%z + 7.9e-7_wp
            else if (this%pos%z == state%grid%zmax) then
                this%pos%z = this%pos%z - 7.9e-7_wp
            end if

            this%phi  = atan2(this%nyp, this%nxp)
            this%cosp = cos(this%phi)
            this%sinp = sin(this%phi)
            this%cost = this%nzp
            this%sint = sqrt(1._wp - this%cost**2)  

            call spectrum%p%sample(this%wavelength, tmp)
            
            this%tflag  = .false.
            this%cnts   = 0
            this%bounces = 0
            this%layer  = 1
            this%weight = 1.0_wp

            ! Linear Grid 
            cell = state%grid%get_voxel(this%pos)
            this%xcell = cell(1)
            this%ycell = cell(2)
            this%zcell = cell(3)

            !print*, cell
        end subroutine circular


        subroutine point(this, spectrum, dict, seqs)
        !! isotropic point source

            use sim_state_mod, only : state
            use random,        only : ran2, seq
            use constants,     only : twoPI
            use tomlf,         only : toml_table, get_value
            use piecewiseMod

            class(photon) :: this
            type(spectrum_t), intent(in) :: spectrum
            type(toml_table), optional, intent(inout) :: dict
            type(seq), optional, intent(inout) :: seqs(2)

            integer :: cell(3)
            real(kind=wp) :: wavelength, tmp

            this%pos = photon_origin%pos

            this%phi  = ran2()*twoPI
            this%cosp = cos(this%phi)
            this%sinp = sin(this%phi)
            this%cost = 2._wp*ran2()-1._wp
            this%sint = sqrt(1._wp - this%cost**2)

            this%nxp = this%sint * this%cosp
            this%nyp = this%sint * this%sinp
            this%nzp = this%cost

            this%tflag  = .false.
            this%cnts   = 0
            this%bounces = 0
            this%layer  = 1
            this%weight = 1.0_wp
            ! this%L = 1.0

            call spectrum%p%sample(wavelength, tmp)
            this%wavelength = wavelength

            this%energy = 1._wp
            this%fact = TWOPI/(this%wavelength)

            ! Linear Grid 
            cell = state%grid%get_voxel(this%pos)
            this%xcell = cell(1)
            this%ycell = cell(2)
            this%zcell = cell(3)
        
        end subroutine point

        subroutine focus(this, spectrum, dict, seqs)

            use random,        only : ran2, rang, seq, ranu
            use sim_state_mod, only : state
            use utils,         only : deg2rad
            use vector_class,  only : length
            use tomlf,         only : toml_table, get_value
            use sdfHelpers,    only : rotationAlign, translate
            use mat_class,     only : invert
            use constants,     only : twoPI
            use piecewiseMod

            class(photon) :: this
            type(spectrum_t), intent(in) :: spectrum
            type(toml_table), optional, intent(inout) :: dict
            type(seq), optional, intent(inout) :: seqs(2)

            type(vector)  :: pos, targ, dir, rotation, startPos, a, b
            real(kind=wp) :: dist, tmp, focalLength, radius, beam_size
            integer       :: cell(3), counter

            character(len=:), allocatable :: focus_type
            real(kind=wp) :: x, y, z, phi, sinp, cosp, t(4,4), rotationx, rotationy, rotationz
            real(kind=wp) :: xToMove, yToMove, zToMove, stepSize
            logical       :: inX, inY, inZ, triedX, triedY, triedZ

            call get_value(dict, "focalLength", focalLength)
            call get_value(dict, "focus_type", focus_type)
            call get_value(dict, "beam_size", beam_size)

            call get_value(dict, "rotation%x", rotationx)
            call get_value(dict, "rotation%y", rotationy)
            call get_value(dict, "rotation%z", rotationz)

            if (focus_type == "square") then 
                this%pos%x = ranu(-beam_size, beam_size)
                this%pos%y = ranu(-beam_size, beam_size)
                this%pos%z = 0._wp

            else if(focus_type == "circle")then
                radius = beam_size * sqrt(ran2())
                phi = TWOPI * ran2()
                cosp = cos(phi)
                sinp = sin(phi)
                x = radius * cosp
                y = radius * sinp
                z = 0._wp ! just inside surface of medium. TODO make this user configurable?
                pos = vector(x, y, z)
                this%pos = pos

            else if(focus_type == "gaussian")then
                ![ref] https://omlc.org/classroom/ece532/class4/example.html
                radius = beam_size * sqrt(-log(1-ran2()))
                ! beam_size is the 1/e radius
                phi = TWOPI * ran2()
                cosp = cos(phi)
                sinp = sin(phi)
                x = radius * cosp
                y = radius * sinp
                z = 0._wp ! just inside surface of medium. TODO make this user configurable?
                pos = vector(x, y, z)
                this%pos = pos
            else
                error stop "No such beam type!"
            end if 

            targ = vector(0._wp,0._wp,-focalLength)

            dist = length(this%pos - targ)

            dir = (-1._wp)*(this%pos-targ) / dist
            dir = dir * sign(1._wp, focalLength)
            dir = dir%magnitude()


            !set inital vector from which the source points
            a = vector(0._wp, 0._wp, -1._wp)
            a = a%magnitude()
            !set vector to rotate to. User defined.
            b = vector(rotationx, rotationy, rotationz)
            b = b%magnitude()

            startPos = photon_origin%pos
            startPos%x = -startPos%x 
            startPos%y = -startPos%y 
            startPos%z = -startPos%z 
            
            if (a == b) then 
                t(:, 1) = [1._wp, 0._wp, 0._wp, 0._wp]
                t(:, 2) = [0._wp, 1._wp, 0._wp, 0._wp]
                t(:, 3) = [0._wp, 0._wp, 1._wp, 0._wp]
                t(:, 4) = [0._wp, 0._wp, 0._wp, 1._wp]
            else if (abs(a) == abs(b)) then
                t(:, 1) = [1._wp, 0._wp, 0._wp, 0._wp]
                t(:, 2) = [0._wp, 1._wp, 0._wp, 0._wp]
                t(:, 3) = [0._wp, 0._wp, -1._wp, 0._wp]
                t(:, 4) = [0._wp, 0._wp, 0._wp, 1._wp]

                startPos%z =  startPos%z 
            else 
                t = rotationAlign(a, b)
            end if
            

            ! get rotation matrix
            dir = dir .dot. t 
            dir = dir%magnitude()

            if (abs(a) == abs(b) .and. .not.(a == b)) then
                t(3,3) = 1._wp
            end if
            
            t = matmul(t, invert(translate(startPos)))
            ! transform point
            this%pos = this%pos .dot. t

            this%nxp = dir%x
            this%nyp = dir%y
            this%nzp = dir%z
            this%phi  = atan2(this%nyp, this%nxp)
            this%cosp = cos(this%phi)
            this%sinp = sin(this%phi)
            this%cost = this%nzp
            this%sint = sqrt(1._wp - this%cost**2)

            this%tflag = .false.
            this%bounces = 0
            this%cnts = 0
            this%weight = 1.0_wp
            call spectrum%p%sample(this%wavelength, tmp)
            this%energy = 1.0_wp
            this%fact = TWOPI / this%wavelength

            inX = .false.
            inY = .false.
            inZ = .false.
            triedX = .false. 
            triedY = .false. 
            triedZ = .false.
            counter = 0

            do while (.not.(inX) .or. .not.(inY) .or. .not.(inZ))
                if(this%pos%x <= - state%grid%xmax) then
                    stepSize = ( - state%grid%xmax - this%pos%x + 9e-7_wp)/this%nxp
                    this%pos%x = this%pos%x + dir%x*stepSize
                    this%pos%y = this%pos%y + dir%y*stepSize
                    this%pos%z = this%pos%z + dir%z*stepSize
                    triedX = .true.
                else if(this%pos%x >= state%grid%xmax) then
                    stepSize = (  state%grid%xmax - this%pos%x - 9e-7_wp)/this%nxp
                    this%pos%x = this%pos%x + dir%x*stepSize
                    this%pos%y = this%pos%y + dir%y*stepSize
                    this%pos%z = this%pos%z + dir%z*stepSize
                    triedX = .true.
                else 
                    inX = .true.
                end if 

                if(this%pos%y <= - state%grid%ymax) then
                    stepSize = ( - state%grid%ymax - this%pos%y + 9e-7_wp)/this%nyp
                    this%pos%x = this%pos%x + dir%x*stepSize
                    this%pos%y = this%pos%y + dir%y*stepSize
                    this%pos%z = this%pos%z + dir%z*stepSize
                    triedY = .true.
                else if(this%pos%y >= state%grid%ymax) then
                    stepSize = (  state%grid%ymax - this%pos%y - 9e-7_wp)/this%nyp
                    this%pos%x = this%pos%x + dir%x*stepSize
                    this%pos%y = this%pos%y + dir%y*stepSize
                    this%pos%z = this%pos%z + dir%z*stepSize
                    triedY = .true.
                else 
                    inY = .true.
                end if 

                if(this%pos%z <= - state%grid%zmax) then
                    stepSize = ( - state%grid%zmax - this%pos%z + 9e-7_wp)/this%nzp
                    this%pos%x = this%pos%x + dir%x*stepSize
                    this%pos%y = this%pos%y + dir%y*stepSize
                    this%pos%z = this%pos%z + dir%z*stepSize
                    triedZ = .true.
                else if(this%pos%z >= state%grid%zmax) then
                    stepSize = (  state%grid%zmax - this%pos%z - 9e-7_wp)/this%nzp
                    this%pos%x = this%pos%x + dir%x*stepSize
                    this%pos%y = this%pos%y + dir%y*stepSize
                    this%pos%z = this%pos%z + dir%z*stepSize
                    triedZ = .true.
                else 
                    inZ = .true.
                end if 

                if ((triedX .and. triedY .and. triedZ) .or. counter >4) then
                    !This will never be within the grid
                    exit
                end if
                counter = counter +1
            end do

            ! Linear Grid 
            cell = state%grid%get_voxel(this%pos)
            this%xcell = cell(1)
            this%ycell = cell(2)
            this%zcell = cell(3)
        end subroutine focus


        subroutine uniform(this, spectrum, dict, seqs)
        !! uniformly illuminate a surface of the simulation media

            use random,        only : ranu, ran2, randint, seq
            use sim_state_mod, only : state
            use tomlf,         only : toml_table, get_value
            use constants,     only : TWOPI
            use piecewiseMod

            class(photon) :: this
            type(spectrum_t), intent(in) :: spectrum
            type(toml_table), optional, intent(inout) :: dict
            type(seq), optional, intent(inout) :: seqs(2)

            integer       :: cell(3)
            type(vector)  :: pos1, pos2, pos3
            real(kind=wp) :: rx, ry, tmp, stepSize
            logical       :: inX, inY, inZ

            this%nxp = photon_origin%nxp
            this%nyp = photon_origin%nyp
            this%nzp = photon_origin%nzp

            this%cost = this%nzp
            this%sint = sqrt(1._wp - this%cost**2)

            this%phi = atan2(this%nyp, this%nxp)
            this%cosp = cos(this%phi)
            this%sinp = sin(this%phi)

            call get_value(dict, "pos1%x", pos1%x)
            call get_value(dict, "pos1%y", pos1%y)
            call get_value(dict, "pos1%z", pos1%z)

            call get_value(dict, "pos2%x", pos2%x)
            call get_value(dict, "pos2%y", pos2%y)
            call get_value(dict, "pos2%z", pos2%z)

            call get_value(dict, "pos3%x", pos3%x)
            call get_value(dict, "pos3%y", pos3%y)
            call get_value(dict, "pos3%z", pos3%z)

            rx = ran2()!seqs(1)%next()
            ry = ran2()!seqs(2)%next()
            this%pos%x = pos1%x + rx * pos2%x + ry * pos3%x
            this%pos%y = pos1%y + rx * pos2%y + ry * pos3%y
            this%pos%z = pos1%z + rx * pos2%z + ry * pos3%z

            if (this%pos%x == -state%grid%xmax) then
                this%pos%x = this%pos%x + 7.9e-7_wp
            else if (this%pos%x == state%grid%xmax) then
                this%pos%x = this%pos%x - 7.9e-7_wp
            end if
            if (this%pos%y == -state%grid%ymax) then
                this%pos%y = this%pos%y + 7.9e-7_wp
            else if (this%pos%y == state%grid%ymax) then
                this%pos%y = this%pos%y - 7.9e-7_wp
            end if
            if (this%pos%z == -state%grid%zmax) then
                this%pos%z = this%pos%z + 7.9e-7_wp
            else if (this%pos%z == state%grid%zmax) then
                this%pos%z = this%pos%z - 7.9e-7_wp
            end if

            this%tflag = .false.
            this%cnts = 0
            this%bounces = 0
            this%weight = 1.0_wp

            !FOR PHASE
            call spectrum%p%sample(this%wavelength, tmp)
            this%energy = 1._wp
            this%fact = TWOPI/(this%wavelength)
            this%phase = 0._wp

            ! Linear Grid 
            cell = state%grid%get_voxel(this%pos)
            this%xcell = cell(1)
            this%ycell = cell(2)
            this%zcell = cell(3)
            !print*, cell
            !print*, this%pos

        end subroutine uniform


        subroutine pencil(this, spectrum, dict, seqs)
        !! pencil beam source

            use random,        only : ranu, seq
            use sim_state_mod, only : state
            use tomlf,         only : toml_table, get_value
            use piecewiseMod
            use constants,     only : TWOPI

            class(photon) :: this
            type(spectrum_t), intent(in) :: spectrum
            type(toml_table), optional, intent(inout) :: dict
            type(seq), optional, intent(inout) :: seqs(2)

            integer :: cell(3)
            real(kind=wp) :: tmp

            this%pos = photon_origin%pos

            if (this%pos%x == -state%grid%xmax) then
                this%pos%x = this%pos%x + 7.9e-7_wp
            else if (this%pos%x == state%grid%xmax) then
                this%pos%x = this%pos%x - 7.9e-7_wp
            end if
            if (this%pos%y == -state%grid%ymax) then
                this%pos%y = this%pos%y + 7.9e-7_wp
            else if (this%pos%y == state%grid%ymax) then
                this%pos%y = this%pos%y - 7.9e-7_wp
            end if
            if (this%pos%z == -state%grid%zmax) then
                this%pos%z = this%pos%z + 7.9e-7_wp
            else if (this%pos%z == state%grid%zmax) then
                this%pos%z = this%pos%z - 7.9e-7_wp
            end if

            this%nxp = photon_origin%nxp
            this%nyp = photon_origin%nyp
            this%nzp = photon_origin%nzp

            this%phi = atan2(this%nyp, this%nxp)
            this%cosp = cos(this%phi)
            this%sinp = sin(this%phi)
            this%cost = this%nzp
            this%sint = sqrt(1._wp - this%cost**2)

            this%tflag = .false.
            this%bounces = 0
            this%cnts = 0
            this%weight = 1.0_wp
            call spectrum%p%sample(this%wavelength, tmp)
            this%energy = 1.0_wp
            this%fact = TWOPI / this%wavelength

            ! Linear Grid 
            cell = state%grid%get_voxel(this%pos)
            this%xcell = cell(1)
            this%ycell = cell(2)
            this%zcell = cell(3)
        end subroutine pencil

        subroutine dslit(this, spectrum, dict, seqs)
        !!sample from double slit to produce diff pattern
        ! todo add in user defined slit widths
        ! add correct normalisation
            use random,        only : ranu, ran2, randint, seq
            use sim_state_mod, only : state
            use tomlf,         only : toml_table, get_value
            use constants,     only : TWOPI
            use piecewiseMod

            class(photon) :: this
            type(spectrum_t), intent(in) :: spectrum
            type(toml_table), optional, intent(inout) :: dict
            type(seq), optional, intent(inout) :: seqs(2)

            integer       :: cell(3)
            real(kind=wp) :: x1, y1, z1, x2, y2, z2, a, b, tmp

            call spectrum%p%sample(this%wavelength, tmp)
            this%energy = 1._wp 
            this%fact = TWOPI/(this%wavelength)

            a = 60._wp*this%wavelength !distance between slits
            b = 20._wp*this%wavelength !2 slit width

            if(ran2() > 0.5_wp)then ! pick slit and sample x, y position
                x1 = ranu(a/2._wp,  a/2._wp + b)
                y1 = ranu(-b*0.5_wp, b*0.5_wp)
            else
                x1 = ranu(-a/2._wp,  -a/2._wp - b)
                y1 = ranu(-b*0.5_wp, b*0.5_wp)
            end if

            z2 = 5.0_wp - (1.e-5_wp*(2._wp*(5.0_wp/400._wp)))

            x2 = ranu(-5.0_wp, 5.0_wp)
            y2 = ranu(-5.0_wp, 5.0_wp)
            z1 = (10000._wp * this%wavelength) - 5.0_wp !screen location

            this%pos%x = x2
            this%pos%y = y2
            this%pos%z = z2

            this%phase = sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)

            this%nxp = (x2 - x1) / this%phase
            this%nyp = (y2 - y1) / this%phase
            this%nzp = -abs((z2 - z1)) / this%phase

            this%tflag = .false.
            this%cnts = 0
            this%bounces = 0
            this%weight = 1.0_wp

            !Set direction cosine/sine
            this%cost = this%nzp
            this%sint = sqrt(1._wp - this%cost**2)

            this%phi = atan2(this%nyp, this%nxp)
            this%cosp = cos(this%phi)
            this%sinp = sin(this%phi)

            ! Linear Grid 
            cell = state%grid%get_voxel(this%pos)
            this%xcell = cell(1)
            this%ycell = cell(2)
            this%zcell = cell(3)

        end subroutine dslit

        subroutine aperture(this, spectrum, dict, seqs)
            !! sample from square aperture to produce diff pattern
            !add user defined apwid and F
            ! add correct normalisation
            use random,        only : ranu, ran2, randint, seq
            use sim_state_mod, only : state
            use tomlf,         only : toml_table, get_value
            use constants,     only : TWOPI
            use piecewiseMod

            class(photon) :: this
            type(spectrum_t), intent(in) :: spectrum
            type(toml_table), optional, intent(inout) :: dict
            type(seq), optional, intent(inout) :: seqs(2)

            integer       :: cell(3)
            real(kind=wp) :: x1, y1, z1, x2, y2, z2, b, F, apwid, tmp

            call spectrum%p%sample(this%wavelength, tmp)
            this%energy = 1._wp 
            this%fact = TWOPI/(this%wavelength)

            apwid = 200e-6_wp !aperture width
            b = apwid/2._wp !slit width
            ! Fresnel number
            F = 4.95_wp

            !sample aperture postiion
            x1 = ranu(-b,b)
            y1 = ranu(-b,b)

            z1 = (1._wp/((((F / apwid)**2) / 2._wp)*this%wavelength)) - 0.5_wp

            x2 = ranu(-0.5_wp, 0.5_wp)
            y2 = ranu(-0.5_wp, 0.5_wp)
            z2 = 0.5_wp - (1.e-5_wp*(2._wp*0.5_wp/400._wp))

            this%pos%x = x2
            this%pos%y = y2
            this%pos%z = z2

            this%phase = sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)

            this%nxp = (x2 - x1) / this%phase
            this%nyp = (y2 - y1) / this%phase
            this%nzp = -abs((z2 - z1)) / this%phase 

            this%tflag = .false.
            this%cnts = 0
            this%bounces = 0
            this%weight = 1.0_wp

            !scattering stuff - not important
            this%cost = this%nzp
            this%sint = sqrt(1._wp - this%cost**2)

            this%phi = atan2(this%nyp, this%nxp)
            this%cosp = cos(this%phi)
            this%sinp = sin(this%phi)

            ! Linear Grid 
            cell = state%grid%get_voxel(this%pos)
            this%xcell = cell(1)
            this%ycell = cell(2)
            this%zcell = cell(3)

        end subroutine aperture
        
        subroutine annulus(this, spectrum, dict, seqs)
            !! annular source
            use constants,     only : TWOPI 
            use utils,         only : deg2rad
            use tomlf,         only : toml_table, get_value
            use sdfHelpers,    only : rotationAlign, translate
            use mat_class,     only : invert
            use random,        only : ran2, rang, seq, ranu
            use sim_state_mod, only : state
            use piecewiseMod

            class(photon) :: this
            type(spectrum_t), intent(in) :: spectrum
            type(toml_table), optional, intent(inout) :: dict
            type(seq), optional, intent(inout) :: seqs(2)

            character(len=:), allocatable :: beam_type
            real(kind=wp) :: focalLength, rlo, rhi, sigma, radius, tmp, mid, x, y, z, phi, sinp, cosp, dist
            real(kind=wp) :: t(4,4), xToMove, yToMove, zToMove, stepSize, rotationx, rotationy, rotationz
            type(vector)  :: pos, dir, targ, a, b, startPos
            integer       :: cell(3), counter
            logical       :: inX, inY, inZ, triedX, triedY, triedZ

            call get_value(dict, "focalLength", focalLength)
            call get_value(dict, "rlo", rlo)
            call get_value(dict, "rhi", rhi)
            call get_value(dict, "sigma", sigma)
            call get_value(dict, "annulus_type", beam_type)

            call get_value(dict, "rotation%x", rotationx)
            call get_value(dict, "rotation%y", rotationy)
            call get_value(dict, "rotation%z", rotationz)

            if(beam_type == "tophat")then
                radius = sqrt(rlo**2 + (rhi**2 - rlo**2) * ran2())
                mid = (rhi + rlo) / 2._wp
            elseif(beam_type == "besselAnnulus")then
                radius = rlo + (rhi-rlo) * ran2()
                mid = (rhi + rlo) / 2._wp
            elseif(beam_type == "gaussian")then
                mid = (rhi + rlo) / 2._wp
                call rang(radius, tmp, mid, sigma)
            else
                error stop "No such beam type!"
            end if 

            phi = TWOPI * ran2()
            cosp = cos(phi)
            sinp = sin(phi)
            x = radius * cosp
            y = radius * sinp
            z = 0.0_wp! just inside surface of medium. TODO make this user configurable?
            pos = vector(x, y, z)
            this%pos = pos

            targ = vector(0._wp,0._wp,-focalLength)

            dist = length(vector(mid*cosp, mid*sinp, 0.0_wp) - targ)

            dir = (-1._wp)*(vector(mid*cosp, mid*sinp, 0.0_wp)-targ) / dist
            dir = dir * sign(1._wp, focalLength)
            dir = dir%magnitude()

            !print*, ""
            !print*, "Before"
            !print*, "pos", this%pos
            !print*, "dir", dir

            !set inital vector from which the source points
            a = vector(0._wp, 0._wp, -1._wp)
            a = a%magnitude()
            !set vector to rotate to. User defined.
            b = vector(rotationx, rotationy, rotationz)
            b = b%magnitude()

            startPos = photon_origin%pos
            startPos%x = -startPos%x 
            startPos%y = -startPos%y 
            startPos%z = -startPos%z 
            
            if (a == b) then 
                t(:, 1) = [1._wp, 0._wp, 0._wp, 0._wp]
                t(:, 2) = [0._wp, 1._wp, 0._wp, 0._wp]
                t(:, 3) = [0._wp, 0._wp, 1._wp, 0._wp]
                t(:, 4) = [0._wp, 0._wp, 0._wp, 1._wp]
            else if (abs(a) == abs(b)) then
                t(:, 1) = [1._wp, 0._wp, 0._wp, 0._wp]
                t(:, 2) = [0._wp, 1._wp, 0._wp, 0._wp]
                t(:, 3) = [0._wp, 0._wp, -1._wp, 0._wp]
                t(:, 4) = [0._wp, 0._wp, 0._wp, 1._wp]

                startPos%z =  startPos%z 
            else 
                t = rotationAlign(a, b)
            end if
            

            ! get rotation matrix
            dir = dir .dot. t 
            dir = dir%magnitude()

            if (abs(a) == abs(b) .and. .not.(a == b)) then
                t(3,3) = 1._wp
            end if
            
            t = matmul(t, invert(translate(startPos)))
            ! transform point
            this%pos = this%pos .dot. t

            this%nxp = dir%x
            this%nyp = dir%y
            this%nzp = dir%z
            this%phi  = atan2(this%nyp, this%nxp)
            this%cosp = cos(this%phi)
            this%sinp = sin(this%phi)
            this%cost = this%nzp
            this%sint = sqrt(1._wp - this%cost**2)

            this%tflag = .false.
            this%bounces = 0
            this%cnts = 0
            this%weight = 1.0_wp
            call spectrum%p%sample(this%wavelength, tmp)

            inX = .false.
            inY = .false.
            inZ = .false.
            triedX = .false. 
            triedY = .false. 
            triedZ = .false.
            counter = 0

            do while (.not.(inX) .or. .not.(inY) .or. .not.(inZ))
                if(this%pos%x <= - state%grid%xmax) then
                    stepSize = ( - state%grid%xmax - this%pos%x + 9e-7_wp)/this%nxp
                    this%pos%x = this%pos%x + dir%x*stepSize
                    this%pos%y = this%pos%y + dir%y*stepSize
                    this%pos%z = this%pos%z + dir%z*stepSize
                    triedX = .true.
                else if(this%pos%x >= state%grid%xmax) then
                    stepSize = (  state%grid%xmax - this%pos%x - 9e-7_wp)/this%nxp
                    this%pos%x = this%pos%x + dir%x*stepSize
                    this%pos%y = this%pos%y + dir%y*stepSize
                    this%pos%z = this%pos%z + dir%z*stepSize
                    triedX = .true.
                else 
                    inX = .true.
                end if 

                if(this%pos%y <= - state%grid%ymax) then
                    stepSize = ( - state%grid%ymax - this%pos%y + 9e-7_wp)/this%nyp
                    this%pos%x = this%pos%x + dir%x*stepSize
                    this%pos%y = this%pos%y + dir%y*stepSize
                    this%pos%z = this%pos%z + dir%z*stepSize
                    triedY = .true.
                else if(this%pos%y >= state%grid%ymax) then
                    stepSize = (  state%grid%ymax - this%pos%y - 9e-7_wp)/this%nyp
                    this%pos%x = this%pos%x + dir%x*stepSize
                    this%pos%y = this%pos%y + dir%y*stepSize
                    this%pos%z = this%pos%z + dir%z*stepSize
                    triedY = .true.
                else 
                    inY = .true.
                end if 

                if(this%pos%z <= - state%grid%zmax) then
                    stepSize = ( - state%grid%zmax - this%pos%z + 9e-7_wp)/this%nzp
                    this%pos%x = this%pos%x + dir%x*stepSize
                    this%pos%y = this%pos%y + dir%y*stepSize
                    this%pos%z = this%pos%z + dir%z*stepSize
                    triedZ = .true.
                else if(this%pos%z >= state%grid%zmax) then
                    stepSize = (  state%grid%zmax - this%pos%z - 9e-7_wp)/this%nzp
                    this%pos%x = this%pos%x + dir%x*stepSize
                    this%pos%y = this%pos%y + dir%y*stepSize
                    this%pos%z = this%pos%z + dir%z*stepSize
                    triedZ = .true.
                else 
                    inZ = .true.
                end if 

                if ((triedX .and. triedY .and. triedZ) .or. counter > 3) then
                    !This will never be within the grid
                    exit
                end if
                counter = counter +1
            end do

            ! Linear Grid 
            cell = state%grid%get_voxel(this%pos)
            this%xcell = cell(1)
            this%ycell = cell(2)
            this%zcell = cell(3)
        end subroutine annulus

        subroutine scatter(this, hgg, g2, dects)
            !! Scattering routine. Implments both isotropic and henyey-greenstein scattering
            !! taken from [mcxyz](https://omlc.org/software/mc/mcxyz/index.html)

            use constants, only : PI, TWOPI, wp
            use random,    only : ran2
            use detectors, only : dect_array

            class(photon), intent(inout) :: this
            !> g factor
            real(kind=wp), intent(in)    :: hgg
            !> g factor squared
            real(kind=wp), intent(in)    :: g2
            !> array of detectors. Only used if biased scattering is enabled.
            type(dect_array), optional, intent(in) :: dects(:)

            real(kind=wp) :: temp, uxx, uyy, uzz

            if(hgg == 0.0_wp)then
                !isotropic scattering
                this%cost = 2._wp * ran2() - 1._wp
            else
                temp = (1.0_wp - hgg*hgg) / (1.0_wp - hgg + 2._wp*hgg*ran2())
                this%cost = (1.0_wp + hgg*hgg - temp**2) / (2._wp*hgg)
            end if
            this%sint = sqrt(1._wp - this%cost**2)

            this%phi = TWOPI * ran2()
            this%cosp = cos(this%phi)
            this%sinp = sin(this%phi)

            if(this%nzp > 1._wp - 1e-12_wp)then ! near perpindicular               
                uxx = this%sint * this%cosp
                uyy = this%sint * this%sinp
                uzz = this%cost
            else if (this%nzp < -1._wp + 1e-12_wp)then ! near perpindicular               
                uxx = this%sint * this%cosp
                uyy = this%sint * this%sinp
                uzz = -this%cost
            else
                temp = sqrt(1._wp - this%nzp**2)
                uxx = this%sint * ((this%nxp * this%nzp * this%cosp - this%nyp * this%sinp) / temp) + this%nxp * this%cost
                uyy = this%sint * ((this%nyp * this%nzp * this%cosp + this%nxp * this%sinp) / temp) + this%nyp * this%cost
                uzz = -1._wp*this%sint * this%cosp * temp + this%nzp * this%cost
            end if

            temp = sqrt(uxx**2 + uyy**2 + uzz**2)
            do while(abs(temp-1) > 1e-12_wp)
                uxx = uxx/temp
                uyy = uyy/temp
                uzz = uzz/temp
                temp = sqrt(uxx**2 + uyy**2 + uzz**2)
            end do 

            this%nxp = uxx
            this%nyp = uyy
            this%nzp = uzz

        end subroutine scatter
end module photonMod