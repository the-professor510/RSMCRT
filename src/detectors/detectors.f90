module detectors

    !! Module contains each detector type which inherits from the base detector class.
    !! detectors detect photon packets colliding with the detectors.

    use constants,    only : wp
    use detector_mod, only : detector, detector1D, detector2D, hit_t
    use vector_class, only : vector, length

    implicit none
    
    !> Circle detector
    type, extends(detector1D) :: circle_dect
        !> Radius of detector
        real(kind=wp) :: radius
    contains
        procedure :: check_hit  => check_hit_circle
    end type circle_dect

    interface circle_dect
        !> Initialise circular detector
        module procedure init_circle_dect
    end interface circle_dect

    !> Fibre detector
    type, extends(detector1D) :: fibre_dect
        !> focal length of the front lens of the 4f imaging system
        real(kind=wp) :: focalLength1
        !> focal length of the back lens of the 4f imaging system
        real(kind=wp) :: focalLength2
        !> radius/size of the front lens
        real(kind=wp) :: f1Aperture
        !> radius/size of the back lens
        real(kind=wp) :: f2Aperture
        !> distance between the detector plane and the front lens
        real(kind=wp) :: frontOffset
        !> distance between the fibre plane and the back lens
        real(kind=wp) :: backOffset
        !> distance between front lens and a pinhole in between the front and back lens
        real(kind=wp) :: frontToPinSep
        !> distance between pinhole and back lens
        real(kind=wp) :: pinToBackSep
        !> radius/size of pinhole
        real(kind=wp) :: pinAperture
        !> maximum fibre acceptance angle above the optical axis
        real(kind=wp) :: acceptAngle
        !> size of the fibre core
        real(kind=wp) :: coreDiameter
        contains
        procedure :: check_hit  => check_hit_fibre
    end type fibre_dect

    interface fibre_dect
        !> Initialise circular detector
        module procedure init_fibre_dect
    end interface fibre_dect

    !> Annuluar detector
    type, extends(detector1D) :: annulus_dect
        !> Inner radius
        real(kind=wp) :: r1
        !> Outer radius
        real(kind=wp) :: r2
        contains
        procedure :: check_hit => check_hit_annulus
    end type annulus_dect

    interface annulus_dect
        !> Initialise annuluar detector
        module procedure init_annulus_dect
    end interface annulus_dect

    !> Rectangular or "camera" detector
    type, extends(detector2D) :: camera
        !> Normal of the detector
        type(vector)  :: n
        !> Vector from pos (1st corner) to the 2nd corner of the detector
        type(vector)  :: p2
        !> Vector from pos (1st corner) to the 3rd corner of the detector
        type(vector)  :: p3
        !> Edge vector of detector
        type(vector)  :: e1
        !> Edge vector of detector
        type(vector)  :: e2
        !> Width of the detector
        real(kind=wp) :: width
        !> Height of the detector
        real(kind=wp) :: height
        contains
        procedure :: check_hit => check_hit_camera
    end type camera

    interface camera
        module procedure init_camera
    end interface camera
    
    !> Detector array
    type :: dect_array
        class(detector), pointer :: p => null()
    end type dect_array

    private
    public :: camera, annulus_dect, circle_dect, dect_array, fibre_dect

    contains
    
    function init_circle_dect(pos, dir, layer, radius, nbins, trackHistory) result(out)
        !! Initalise Circle detector
        !> Centre of detector
        type(vector),  intent(in) :: pos
        !> Normal of the detector
        type(vector),  intent(in) :: dir
        !> Layer ID
        integer,       intent(in) :: layer
        !> Number of bins in the detector
        integer,       intent(in) :: nbins
        !> Radius of the detector
        real(kind=wp), intent(in) :: radius
        !> Boolean on if to store photon's history prior to hitting the detector.
        logical,       intent(in) :: trackHistory

        type(circle_dect) :: out

        out%dir = dir
        out%pos = pos
        out%layer = layer
        !extra bin for data beyond end of array
        out%nbins = nbins + 1
        out%radius = radius
        allocate(out%data(out%nbins))
        out%data = 0.0_wp
        if(nbins == 0)then
            out%bin_wid = 1._wp
        else
            out%bin_wid = radius / real(nbins, kind=wp)
        end if
        out%trackHistory = trackHistory

    end function init_circle_dect

    logical function check_hit_circle(this, hitpoint)
        !! Check if a hitpoint is in the circle
        
        use geometry, only : intersectCircle

        class(circle_dect), intent(INOUT) :: this
        !> Hitpoint to check
        type(hit_t),        intent(inout) :: hitpoint
        
        real(kind=wp) :: t 

        check_hit_circle = .false.
        check_hit_circle = intersectCircle(this%dir, this%pos, this%radius, hitpoint%pos, hitpoint%dir, t, hitpoint%value1D)
        if(check_hit_circle)then
            if(t <= 0.0_wp .or. t > hitpoint%pointSep)check_hit_circle=.false.
            !is the interaction point outside of the packet path
        end if
    end function check_hit_circle

    function init_annulus_dect(pos, dir, layer, r1, r2, nbins, maxval, trackHistory) result(out)
        !! Initalise Annular detector

        !> Centre of detector
        type(vector),  intent(in) :: pos
        !> Normal of the detector
        type(vector),  intent(in) :: dir
        !> Layer ID
        integer,       intent(in) :: layer
        !> Inner radius
        real(kind=wp), intent(IN) :: r1
        !> Outer radius
        real(kind=wp), intent(IN) :: r2
        !> Number of bins in the detector
        integer,       intent(in) :: nbins
        !> Maximum value to store in bins
        real(kind=wp), intent(in) :: maxval
        !> Boolean on if to store photon's history prior to hitting the detector.
        logical,       intent(in) :: trackHistory

        type(annulus_dect) :: out

        out%pos = pos
        out%dir = dir
        out%layer = layer
        !extra bin for data beyond end of array
        out%nbins = nbins + 1
        out%r1 = r1
        out%r2 = r2
        allocate(out%data(out%nbins))
        out%data = 0.0_wp
        if(nbins == 0)then
            out%bin_wid = 1._wp
        else
            out%bin_wid = (r2-r1) / real(nbins, kind=wp)
        end if
        out%trackHistory = trackHistory

    end function init_annulus_dect

    logical function check_hit_annulus(this, hitpoint)

        use geometry, only : intersectCircle

        !! Check if a hitpoint is in the annulus
        class(annulus_dect), intent(INOUT) :: this
        !> Hitpoint to check
        type(hit_t),         intent(inout)    :: hitpoint

        logical :: hit_circle_r2, hit_circle_r1
        real(kind=wp) :: t

        check_hit_annulus = .false.
        !do we hit the inner void
        hit_circle_r1 = intersectCircle(this%dir, this%pos, this%r1, hitpoint%pos, hitpoint%dir, t, hitpoint%value1D)
        !do we hit the inner void and the annulus
        hit_circle_r2 = intersectCircle(this%dir, this%pos, this%r2, hitpoint%pos, hitpoint%dir, t, hitpoint%value1D)

        !if we don't hit the void but do hit the void and annulus we must be in the annulus
        if((.not. hit_circle_r1) .and. hit_circle_r2) then
            !is the interaction point outside of the packet path
            if(t <= 0.0_wp .or. t > hitpoint%pointSep) then
                !yes it is outside the packet path
                check_hit_annulus=.false.
            else 
                ! it is inside the packet path
                check_hit_annulus = .true.
            end if
        end if

        hitpoint%value1D = hitpoint%value1D - this%r1

    end function check_hit_annulus

    function init_fibre_dect(pos, dir, layer, nbins, maxval, trackHistory, & 
        focalLength1, focalLength2, f1Aperture, f2Aperture, frontOffset, backOffset, & 
        frontToPinSep, pinToBackSep, pinAperture, acceptAngle, coreDiameter) result(out)
        !! Initialise fibre detector
        
        !> Centre of detector
        type(vector),  intent(in) :: pos
        !> Normal of the detector
        type(vector),  intent(in) :: dir
        !> Layer ID
        integer,       intent(in) :: layer
        !> Number of bins in the detector
        integer,       intent(in) :: nbins
        !> Maximum value to store in bins
        real(kind=wp), intent(in) :: maxval
        !> Boolean on if to store photon's history prior to hitting the detector.
        logical,       intent(in) :: trackHistory
        !> focal length of the front lens of the 4f imaging system
        real(kind=wp), intent(in) :: focalLength1
        !> focal length of the back lens of the 4f imaging system
        real(kind=wp), intent(in) :: focalLength2
        !> radius/size of the front lens
        real(kind=wp), intent(in) :: f1Aperture
        !> radius/size of the back lens
        real(kind=wp), intent(in) :: f2Aperture
        !> distance between the detector plane and the front lens
        real(kind=wp), intent(in) :: frontOffset
        !> distance between the fibre plance and the back lens
        real(kind=wp), intent(in) :: backOffset
        !> distance between front lens and a pinhole in between the front and back lens
        real(kind=wp), intent(in) :: frontToPinSep
        !> distance between pinhole and back lens
        real(kind=wp), intent(in) :: pinToBackSep
        !> radius/size of pinhole
        real(kind=wp), intent(in) :: pinAperture
        !> maximum fibre acceptance angle above the optical axis
        real(kind=wp), intent(in) :: acceptAngle
        !> size of the fibre core
        real(kind=wp), intent(in) :: coreDiameter

        type(fibre_dect) :: out

        out%pos = pos
        out%dir = dir
        out%layer = layer
        
        out%focalLength1 = focalLength1
        out%focalLength2 = focalLength2
        out%f1Aperture = f1Aperture
        out%f2Aperture = f2Aperture
        out%frontOffset = frontOffset
        out%backOffset = backOffset
        out%frontToPinSep = frontToPinSep
        out%pinToBackSep = pinToBackSep
        out%pinAperture = pinAperture
        out%acceptAngle = acceptAngle
        out%coreDiameter = coreDiameter

        !extra bin for data beyond end of array
        out%nbins = nbins + 1
        allocate(out%data(out%nbins))
        out%data = 0.0_wp
        if(nbins == 0)then
            out%bin_wid = 1._wp
        else
            out%bin_wid = coreDiameter / 2 / real(nbins, kind=wp)
        end if
        out%trackHistory = trackHistory

    end function init_fibre_dect

    logical function check_hit_fibre(this, hitpoint)

        use constants, only : TWOPI
        use geometry, only : intersectCircle

        !! Check if a hitpoint is in the annulus
        class(fibre_dect), intent(INOUT) :: this
        !> Hitpoint to check
        type(hit_t),         intent(inout) :: hitpoint

        real(kind=wp) :: t 
        real(kind=wp) :: angle
        real(kind=wp) :: gradient
        real(kind=wp) :: radius
        real(kind=wp) :: costt, sintt

        check_hit_fibre = .false.
        check_hit_fibre = intersectCircle(this%dir, this%pos + this%dir*this%frontOffset, & 
                                            this%f1Aperture, hitpoint%pos, hitpoint%dir, t, hitpoint%value1D)
        if( check_hit_fibre) then
            if(t <= 0.0_wp .or. t > hitpoint%pointSep)check_hit_fibre=.false.
        end if

        !print*,this%pos
        !print*, this%pos + this%dir*this%frontOffset

        if(.not. check_hit_fibre) then
            ! The packet won't interact with the first lens stop following it
            return
        end if

        costt = this%dir .dot. hitpoint%dir
        if(costt>1.0_wp)costt=1.0_wp
        sintt = sqrt(1._wp - costt * costt)
        gradient = sintt/costt
        radius = hitpoint%value1D

        !print*, gradient, radius

        ! change gradient as the packet moves through the lens
        ! use thin lens approximation
        gradient = -radius/this%focalLength1 + gradient

        !print*, gradient, radius

        ! move to the pinhole
        radius = radius + gradient*this%frontToPinSep

        ! is the packet blocked by the aperture
        if(radius > this%pinAperture) then
            check_hit_fibre = .false.
            return
        end if

        ! move to the second lens
        radius = radius + gradient*this%pinToBackSep

        !print*, gradient, radius
        ! does the packet enter the second lens
        if(radius > this%f2Aperture) then
            check_hit_fibre = .false.
            return
        end if

        ! change gradient as the packet moves through the lens
        gradient = -radius/this%focalLength2 + gradient

        !move to the fibre
        radius = radius + gradient*this%backOffset

        !does the packet enter the fibre?
        angle = abs(atan(gradient))*360/TWOPI

        !print*, gradient, radius
        if(angle > this%acceptAngle .or. radius > (this%coreDiameter/2)) then
            check_hit_fibre = .false.
        end if    
        
        hitpoint%value1D = radius
    end function check_hit_fibre

    function init_camera(p1, p2, p3, layer, nbins, maxval, trackHistory) result(out)
        !! Initalise Camera detector

        !> Position of the 1st corner of the detector
        type(vector),  intent(in) :: p1
        !> Distance from p1 to the 2nd corner
        type(vector),  intent(in) :: p2
        !> Distance from p1 to the 3rd corner
        type(vector),  intent(in) :: p3
        !> Layer ID
        integer,       intent(in) :: layer
        !> Number of bins in the detector
        integer,       intent(in) :: nbins
        !> Maximum value to store in bins
        real(kind=wp), intent(in) :: maxval
        !> Boolean on if to store photon's history prior to hitting the detector.
        logical,       intent(in) :: trackHistory
        type(camera) :: out

        out%pos = p1
        out%p2 = p2
        out%p3 = p3
        out%e1 = p2 - p1
        out%e2 = p3 - p1
        out%width = length(out%e1)
        out%height = length(out%e2)
        out%n = out%e2 .cross. out%e1
        out%n = out%n%magnitude()
        out%layer = layer
        !extra bin for data beyond end of array
        out%nbinsX = nbins + 1
        out%nbinsY = nbins + 1
        allocate(out%data(out%nbinsX, out%nbinsY))
        out%data = 0.0_wp
        if(nbins == 0)then
            out%bin_wid_x = 1._wp
            out%bin_wid_y = 1._wp
        else
            out%bin_wid_x = maxval / real(out%nbinsX, kind=wp)
            out%bin_wid_y = maxval / real(out%nbinsY, kind=wp)
        end if
        out%trackHistory = trackHistory

    end function init_camera

    logical function check_hit_camera(this, hitpoint)
        !! Check if a hitpoint is in the camera detector
        !! [ref](https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-plane-and-ray-disk-intersection)
        class(camera), intent(inout) :: this
        !> Hitpoint to check
        type(hit_t),   intent(inout)    :: hitpoint

        real(kind=wp) :: t, proj1, proj2
        type(vector)  :: v

        check_hit_camera = .false.
        !if(this%layer /= hitpoint%layer)return

        t = ((this%pos - hitpoint%pos) .dot. this%n) / (hitpoint%dir .dot. this%n)
        if(t >= 0._wp)then
            v = (hitpoint%pos + t * hitpoint%dir) - this%pos
            proj1 = (v .dot. this%e1) / this%width
            proj2 = (v .dot. this%e2) / this%height
            if((proj1 < this%width .and. proj1 > 0._wp) .and. (proj2 < this%height .and. proj2 > 0._wp))then
                check_hit_camera = .true.
            end if
        end if
    end function check_hit_camera
end module detectors