module detector_mod
    !! Module contains photon detector abstract class and the derived types the inherit from it
    !! not fully implmented
    use vector_class
    use constants, only : wp
    
    implicit none
    !> Hit type, which records possible interaction information
    type :: hit_t
        !> Poition of the interaction
        type(vector)  :: pos
        !> Direction the photon came from
        type(vector)  :: dir
        !> Distance over which if -t is greater than then the particle does not interact
        real(kind=wp) :: pointSep
        !> Where to store the value, used in record_hit_1D_sub
        real(kind=wp) :: value1D
        ! !> Where to store the value, used in record_hit_2D_sub
        ! real(kind=wp) :: value2D
        !> Weight of the packet
        real(kind=wp) :: weight
    end type hit_t

    !only needed if using a stack to init with a single null value
    interface hit_t
        module procedure hit_init
    end interface hit_t

    !> abstract detector
    type, abstract :: detector
        !> position of the detector
        type(vector)  :: pos
        !> Surface normal of the detector
        type(vector)  :: dir
        !> Layer ID of the detector
        integer :: layer
        !> Boolean, if true store the history of the photon prior to detection.
        logical :: trackHistory
        !> Detector ID
        character(len=:), allocatable :: ID
        !> Target for inverse MCRT
        real(kind = wp) :: targetValue
        contains
            
            procedure(recordHitInterface), deferred, public :: record_hit
            procedure(checkHitInterface),  deferred, public :: check_hit
            procedure(zeroHitInterface), deferred, public :: zero_dect
            procedure(totalHitInterface), deferred, public :: total_dect
    end type detector

    abstract interface
        logical function checkHitInterface(this, hitpoint)
            use vector_class
            use constants, only : wp
            import detector, hit_t

            class(detector), intent(inout) :: this
            type(hit_t),     intent(inout) :: hitpoint
        end function checkHitInterface

        subroutine recordHitInterface(this, hitpoint, history)
            use constants,     only : wp
            use historyStack,  only : history_stack_t
            use vector_class
            import detector, hit_t

            class(detector),       intent(inout) :: this
            type(hit_t),           intent(inout) :: hitpoint
            type(history_stack_t), intent(inout) :: history
        end subroutine recordHitInterface

        subroutine zeroHitInterface(this)
            use vector_class
            use constants, only : wp
            import detector

            class(detector),       intent(inout) :: this
        end subroutine zeroHitInterface

        subroutine totalHitInterface(this, total)
            use vector_class
            use constants, only : wp
            import detector

            class(detector),        intent(inout) :: this
            real(kind=wp),          intent(out)   :: total
        end subroutine totalHitInterface
    end interface
    
    !> 1D detector type. Records linear information
    type, abstract, extends(detector) :: detector1D
        !> Number of bins
        integer       :: nbins
        !> Bin width
        real(kind=wp) :: bin_wid
        !> Bins
        real(kind=wp), allocatable :: data(:)
        contains
        procedure :: record_hit => record_hit_1D_sub
        procedure :: zero_dect => zero_1D_sub
        procedure :: total_dect => total_1D
    end type detector1D
    
    !> 2D detecctor type. Records spatial information
    type, abstract, extends(detector) :: detector2D
        !> Number of bins in x dimension (detector space)
        integer       :: nbinsX
        !> Number of bins in y dimension (detector space)
        integer       :: nbinsY
        !> Bin width in the x dimension
        real(kind=wp) :: bin_wid_x
        !> Bin width in the y dimension
        real(kind=wp) :: bin_wid_y
        !> Bins
        real(kind=wp), allocatable :: data(:,:)
        contains
        procedure :: record_hit => record_hit_2D_sub
        procedure :: zero_dect => zero_2D_sub
        procedure :: total_dect => total_2D
    end type detector2D

    private
    public :: detector, detector1D, detector2D, hit_t

contains
    type(hit_t) function hit_init(val)

        real(kind=wp), intent(in) :: val
        type(vector) :: tmp

        tmp = vector(val, val, val)

        hit_init = hit_t(tmp, tmp, val, val, val)

    end function hit_init
   
    subroutine record_hit_1D_sub(this, hitpoint, history)
        !! check if a hit is on the detector and record it if so

        use historyStack,  only : history_stack_t
        use sim_state_mod, only : state

        class(detector1D),     intent(inout) :: this
        !> Interaction information
        type(hit_t),           intent(inout) :: hitpoint
        !> Photon packet history
        type(history_stack_t), intent(inout) :: history

        real(kind=wp) :: value
        integer       :: idx

        if(this%check_hit(hitpoint))then
            value = hitpoint%value1D

            idx = min(nint(value / this%bin_wid) + 1, this%nbins)
            !$omp atomic
            this%data(idx) = this%data(idx) + hitpoint%weight
            if(this%trackHistory)then
                call history%write()
            end if
        end if
        if(state%trackHistory)call history%zero()
    end subroutine record_hit_1D_sub

    subroutine zero_1D_sub(this)
        class(detector1D),     intent(inout) :: this
        this%data= 0._wp
    end subroutine zero_1D_sub

    subroutine zero_2D_sub(this)
        class(detector2D),     intent(inout) :: this
        this%data = 0._wp
    end subroutine zero_2D_sub

    subroutine total_1D(this, total)

        class(detector1D),      intent(inout) :: this
        !> Sum of all values in data(:)
        real(kind=wp),          intent(out)   :: total

        integer :: i

        total = 0._wp
        do i= 1, size(this%data)
            total = total + this%data(i)
        end do
    end subroutine total_1D

    subroutine total_2D(this, total)

        class(detector2D),      intent(inout) :: this
        !> Sum of all values in data(:)
        real(kind=wp),          intent(out)   :: total

        integer :: i, j

        total = 0._wp
        do i= 1, size(this%data, dim = 1 )
            do j = 1, size(this%data, dim = 2)
                total = total + this%data(i,j)
            end do
        end do
    end subroutine total_2D


    subroutine record_hit_2D_sub(this, hitpoint, history)
        !! check if a hit is on the detector and record it if so

        use historyStack, only : history_stack_t
        use sim_state_mod, only : state

        class(detector2D),     intent(inout) :: this
        !> Interaction information
        type(hit_t),           intent(inout)    :: hitpoint
        !> Photon packet history
        type(history_stack_t), intent(inout) :: history

        real(kind=wp), volatile :: x, y
        integer       :: idx, idy

        if(this%check_hit(hitpoint))then
            x = hitpoint%pos%z + this%pos%x
            y = hitpoint%pos%y + this%pos%y
            idx = min(int(x / this%bin_wid_x) + 1, this%nbinsX)
            idy = min(int(y / this%bin_wid_y) + 1, this%nbinsY)
            if(idx < 1)idx = this%nbinsX
            if(idy < 1)idy = this%nbinsY
            !$omp atomic
            this%data(idx, idy) = this%data(idx, idy) + 1
            if(this%trackHistory)then
                call history%write()
            end if
        end if
        if(state%trackHistory)call history%zero()
        end subroutine record_hit_2D_sub
end module detector_mod
! program test
!     use detector_mod
!     use vector_class
!     use constants, only : wp
!     implicit none

!     type(hit_t) :: hit
!     type(vector) :: pos, dir
!     integer :: layer

!     type(circle_dect) :: dect_c
!     type(annulus_dect) :: dect_a

!     dect_c = circle_dect(vector(0._wp, 0._wp, 0._wp), 1, .5_wp, 100, 100._wp)
!     dect_a = annulus_dect(vector(0._wp, 0._wp, 0._wp), 1, .25_wp, .5_wp, 100, 100._wp)

!     layer = 1
!     pos = vector(0._wp, .5_wp, 0._wp)
!     dir = vector(0._wp, 0._wp, 1._wp)
    
!     hit = hit_t(pos, dir, 99._wp, layer)
!     call dect_c%record_hit(hit)
!     print*,sum(dect_c%data)

!     pos = vector(0._wp, .25_wp, 0._wp)
!     dir = vector(0._wp, 0._wp, 1._wp)
    
!     hit = hit_t(pos, dir, 99._wp, layer)
!     call dect_a%record_hit(hit)
!     print*,sum(dect_a%data)
! end program test