module testsDetectorMod

    use detectors
    use detector_mod, only : hit_t
    use testdrive,    only : new_unittest, unittest_type, error_type, check
    use constants,    only : wp

    implicit none

    private
    public :: detector_suite

    contains

    subroutine detector_suite(testsuite)

        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                new_unittest("Check_hit_circle", hit_circle), &
                new_unittest("Check_hit_annulus", hit_annulus), &
                new_unittest("Check_hit_camera", hit_camera) &
                ]

    end subroutine detector_suite

    subroutine hit_circle(error)

        use vector_class, only : vector
        use historyStack

        type(error_type), allocatable, intent(out) :: error

        type(hit_t) :: hitpoint
        type(circle_dect) :: a
        type(vector) :: pos, dir
        integer :: layer, nbins
        real(kind=wp) :: radius, maxval, pointSep, weight, targetValue
        logical :: flag
        type(history_stack_t) :: history
        character(len=:), allocatable :: dect_ID

        pos = vector(0.5_wp, 0._wp, 0._wp)
        dir = vector(1._wp, 0._wp, 0._wp)
        layer = 1
        radius = 0.5
        nbins = 100
        maxval = 100._wp
        dect_ID = "circleTest"
        targetValue = 0._wp
        a = circle_dect(pos, dir, layer, radius, nbins, .false.,dect_ID,targetValue)

        pos = vector(0._wp, 0._wp, 0._wp)
        dir = vector(1._wp, 0._wp, 0._wp)
        pointSep = 1._wp
        layer = 1
        weight = 1._wp
        hitpoint = hit_t(pos, dir, pointSep, layer, weight)

        flag = a%check_hit(hitpoint)
        call check(error, flag, .true.)
        if(allocated(error))return

        call a%record_hit(hitpoint, history)

        call check(error, sum(a%data), 1.0_wp)
        if(allocated(error))return

    end subroutine hit_circle


    subroutine hit_camera(error)

        use vector_class, only : vector
        use historyStack

        type(error_type), allocatable, intent(out) :: error

        type(hit_t) :: hitpoint
        type(camera) :: a
        type(vector) :: pos, dir, p1, p2, p3
        integer :: layer, nbins
        real(kind=wp) :: maxval, val, weight, targetValue
        logical :: flag
        type(history_stack_t) :: history
        character(len=:), allocatable :: dect_ID

        p1 = vector(-1._wp, -1._wp, -1._wp)
        p2 = vector(0._wp, 2._wp, 0._wp)
        p3 = vector(0._wp, 0._wp, 2._wp)
        layer = 1
        nbins = 100
        maxval = 100._wp
        dect_ID = "cameraTest"
        targetValue = 1.0_wp
        a = camera(p1, p2, p3, layer, nbins, maxval, .false.,dect_ID, targetValue)

        pos = vector(10._wp, 0._wp, 0._wp)
        dir = vector(-1._wp, 0._wp, 0._wp)
        val = 1._wp
        layer = 1
        weight = 1._wp
        hitpoint = hit_t(pos, dir, val, layer, weight)

        flag = a%check_hit(hitpoint)

        call check(error, flag, .true.)
        if(allocated(error))return
    
        call a%record_hit(hitpoint, history)
        call check(error, sum(a%data), 1.0_wp)
        if(allocated(error))return

        pos = vector(10._wp, 0._wp, 0._wp)
        dir = vector(1._wp, 0._wp, 0._wp)
        val = 1._wp
        layer = 1
        weight = 1._wp
        hitpoint = hit_t(pos, dir, val, layer, weight)

        flag = a%check_hit(hitpoint)

        call check(error, flag, .false.)
        if(allocated(error))return

    end subroutine hit_camera

    subroutine hit_annulus(error)

        use vector_class, only : vector
        use historyStack

        type(error_type), allocatable, intent(out) :: error

        type(hit_t) :: hitpoint
        type(annulus_dect) :: a
        type(vector) :: pos, dir
        integer :: layer, nbins
        real(kind=wp) :: maxval, pointSep, r1, r2, weight, target_value
        logical :: flag
        type(history_stack_t) :: history
        character(len=:), allocatable :: dect_ID

        layer = 1
        nbins = 100
        maxval = 100._wp
        r1 = 0.5_wp
        r2 = 1.0_wp
        pos = vector(0.5_wp, 0._wp, 0._wp)
        dir = vector(1._wp, 0._wp, 0._wp)
        dect_ID = "annulusTest"
        target_value = 1.0_wp
        a = annulus_dect(pos, dir, layer, r1, r2, nbins, maxval, .false.,dect_ID, target_value)

        pos = vector(0._wp, 0.75_wp, 0._wp)
        dir = vector(1._wp, 0._wp, 0._wp)
        pointSep = 1._wp
        layer = 1
        weight = 1._wp
        hitpoint = hit_t(pos, dir, pointSep, layer, weight)

        flag = a%check_hit(hitpoint)
        call check(error, flag, .true.)
        if(allocated(error))return

        call a%record_hit(hitpoint, history)
        call check(error, sum(a%data), 1.0_wp)
        if(allocated(error))return

        pos = vector(0._wp, 0._wp, 0._wp)
        dir = vector(1._wp, 0._wp, 0._wp)
        pointSep = 1._wp
        layer = 1
        weight = 1._wp
        hitpoint = hit_t(pos, dir, pointSep, layer, weight)

        flag = a%check_hit(hitpoint)

        call check(error, flag, .false.)
        if(allocated(error))return

    end subroutine hit_annulus
end module testsDetectorMod