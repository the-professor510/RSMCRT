module photonMod
    
    use vector_class

    implicit none
    
    type :: photon
    
        type(vector) :: pos                   ! position
        real    :: nxp, nyp, nzp                ! direction vectors
        real    :: sint, cost, sinp, cosp, phi  ! direction cosines
        real    :: wavelength           ! Only used if tracking the phase
        integer :: xcell, ycell, zcell          ! grid cell position
        logical :: tflag                        ! Is photon dead?
        integer :: layer ! pointer to sdf packet is inside

        procedure(generic_emit), pointer :: emit => null()

    end type photon

    interface photon
        module procedure init_source
    end interface photon

    abstract interface
        subroutine generic_emit(this, grid, iseed)
            
            use gridMod

            import :: photon
            class(photon) :: this

            type(cart_grid), intent(IN)    :: grid
            integer,         intent(INOUT) :: iseed

        end subroutine generic_emit
    end interface

    private
    public :: photon, init_source

    contains
        
        type(photon) function init_source(choice)

            implicit none

            character(*), intent(IN) :: choice

            if(choice == "uniform")then
                init_source%emit => uniform
            elseif(choice == "pencil")then
                init_source%emit => pencil
            else
                init_source%emit => circular_beam
            end if


        end function init_source

        subroutine uniform(this, grid, iseed)

            use gridMod
            use random, only : ranu, ran2
            use constants, only : twoPI

            implicit none

            class(photon) :: this
            type(cart_grid), intent(IN)    :: grid
            integer,         intent(INOUT) :: iseed

            this%pos%x = ranu(-grid%xmax, grid%xmax)
            this%pos%y = .98!ranu(-grid%ymax, grid%ymax)
            this%pos%z = ranu(-grid%zmax, grid%zmax)!epsilon(1.)

            this%phi  = 0.!ran2()*twoPI
            this%cosp = 0.!cos(this%phi)
            this%sinp = -1.!sin(this%phi)
            this%cost = 0.!2.*ran2()-1.
            this%sint = 1.!sqrt(1. - this%cost**2)

            this%nxp = this%sint * this%cosp
            this%nyp = this%sint * this%sinp
            this%nzp = this%cost

            this%tflag = .false.
            this%layer = 2

            ! Linear Grid 
            this%xcell=int(grid%nxg*(this%pos%x+grid%xmax)/(2.*grid%xmax))+1
            this%ycell=int(grid%nyg*(this%pos%y+grid%ymax)/(2.*grid%ymax))+1
            this%zcell=int(grid%nzg*(this%pos%z+grid%zmax)/(2.*grid%zmax))+1

        end subroutine uniform


        subroutine pencil(this, grid, iseed)

            use gridMod
            use random, only : ranu

            implicit none

            class(photon) :: this
            type(cart_grid), intent(IN)    :: grid
            integer,         intent(INOUT) :: iseed

            this%pos%z = grid%zmax - epsilon(1.d0)
            this%pos%x = ranu(-grid%xmax/10., grid%xmax/10.)
            this%pos%y = ranu(-grid%ymax/10., grid%ymax/10.)

            this%phi = 0.
            this%cosp = 0.d0
            this%sinp = 0.d0          
            this%cost = -1.d0 
            this%sint =  0.d0

            this%nxp = this%sint * this%cosp  
            this%nyp = this%sint * this%sinp
            this%nzp = this%cost

            this%tflag = .false.

            ! Linear Grid 
            this%xcell=int(grid%nxg*(this%pos%x+grid%xmax)/(2.*grid%xmax))+1
            this%ycell=int(grid%nyg*(this%pos%y+grid%ymax)/(2.*grid%ymax))+1
            this%zcell=int(grid%nzg*(this%pos%z+grid%zmax)/(2.*grid%zmax))+1

        end subroutine pencil

        subroutine circular_beam(this, grid, iseed)

            use gridMod,   only : cart_grid
            use constants, only : PI, twoPI
            use random,    only : ranu, rang
            use surfaces, only : intersect_cone, reflect_refract
            use vector_class

            implicit none

            class(photon) :: this
            type(cart_grid), intent(IN)    :: grid
            integer,         intent(INOUT) :: iseed


            real :: seperation, beam_width, radius, height, alpha, axicon_n, base_pos
            real :: posx, posy, t, k
            type(vector) :: centre, pos, dir, normal
            logical :: flag

            seperation = 2.5d-3
            !beam paramaters
            beam_width = 100d-6
            !axicon paramaters
            radius = 12.7d-3
            height = 1.1d-3
            alpha = atan(height / radius)
            k  = (radius / height)**2
            axicon_n = 1.45
            base_pos = grid%zmax + ((seperation + beam_width) / tan(alpha * (axicon_n -1.)))
            centre = vector(0., 0., base_pos)

            call rang(posx, posy, 0., beam_width)
            pos = centre + vector(posx, posy, 2*height+base_pos)!2*height as want the upper cone
            dir = vector(0., 0., -1.)
            ! cartesian equation defines these two cones. usually want the lower cone
            ! in this case we want the upper cone as the axicon points down
            ! therefore need to change the computation of the normals and initial postion of packet 
            !
            ! \      /
            !  \    /   upper cone
            !   \  /
            !    \/
            !    /\
            !   /  \
            !  /    \   lower cone
            ! /      \

            flag = intersect_cone(pos, dir, t, centre, radius, height)
            if(flag)then
                pos = pos + t*dir
                if(pos%z >= centre%z+height)then! >= for upper cone
                    !derivative of the cartesian cone eqn
                    normal = vector(2*(pos%x-centre%x) / k, 2*(pos%y-centre%y) / k, -2*(pos%z-centre%z)+2*height)
                    normal = normal *(-1.)! upper cone so invert normals
                    normal = normal%magnitude()
                    call reflect_refract(dir, normal, axicon_n, 1., flag)
                    ! move to mediums surface and step inside
                    t = ((grid%zmax- epsilon(1.d0)) - pos%z) / dir%z
                    pos = pos + t*dir
                end if
            end if

            this%pos = pos

            this%nxp = dir%x
            this%nyp = dir%y
            this%nzp = dir%z

            this%cost = this%nzp
            this%sint = sqrt(1.d0 - this%cost**2)

            this%phi = atan2(this%nyp, this%nxp)
            this%cosp = cos(this%phi)
            this%sinp = sin(this%phi)

            this%tflag = .false.

            ! Linear Grid 
            this%xcell=int(grid%nxg*(this%pos%x+grid%xmax)/(2.*grid%xmax))+1
            this%ycell=int(grid%nyg*(this%pos%y+grid%ymax)/(2.*grid%ymax))+1
            this%zcell=int(grid%nzg*(this%pos%z+grid%zmax)/(2.*grid%zmax))+1


        end subroutine circular_beam

        ! subroutine bessel(this, grid, iseed)

        !     use gridMod,   only : cart_grid
        !     use constants, only : PI, twoPI
        !     use random,    only : ranu, rang

        !     implicit none

        !     class(photon) :: this
        !     type(cart_grid), intent(IN)    :: grid
        !     integer,         intent(INOUT) :: iseed

        !     real :: tana, r_pos, x0, y0, z0, dist, n


        !     real :: waist, d, raxi, alpha

        !     waist = .5d-1
        !     alpha = 5.
        !     raxi = 12.7
        !     d = 10.d0
        !     this%wavelength = 488d-9
        !     this%fact = twopi/this%wavelength

        !     this%phase = 0.d0
        !     tana = tan(alpha *pi/180.)

        !     call rang(this%xp, this%yp, 0.d0, sqrt(2.)*waist/4., iseed)
        !     r_pos = sqrt(this%xp**2 + this%yp**2)

        !     this%zp = (raxi - r_pos) * tana
        !     this%phase = this%zp * n

        !     x0 = ranu(-grid%xmax, grid%xmax, iseed)
        !     y0 = ranu(-grid%xmax, grid%xmax, iseed)
        !     z0 = (r_pos * tana) + d + this%zp
            

        !     dist = sqrt((x0 - this%xp)**2 + (y0 - this%yp)**2 + (z0 - this%zp)**2)

        !     this%phase = this%phase + dist

        !     this%nxp = (x0 - this%xp) / dist
        !     this%nyp = (y0 - this%yp) / dist
        !     this%nzp = -(z0 - this%zp) / dist ! -ve due to way z pos is defined

        !     this%cost = this%nzp
        !     this%sint = sqrt(1.d0 - this%cost**2)

        !     this%phi = atan2(this%nyp, this%nxp)
        !     this%cosp = cos(this%phi)
        !     this%sinp = sin(this%phi)

        !     this%xp = x0
        !     this%yp = y0
        !     this%zp = grid%zmax - grid%delta
        !     this%tflag = .false.
        !     this%xcell = int(grid%nxg * (this%xp + grid%xmax) / (2. * grid%xmax)) + 1
        !     this%ycell = int(grid%nyg * (this%yp + grid%ymax) / (2. * grid%ymax)) + 1
        !     this%zcell = int(grid%nzg * (this%zp + grid%zmax) / (2. * grid%zmax)) + 1

        ! end subroutine bessel

end module photonMod