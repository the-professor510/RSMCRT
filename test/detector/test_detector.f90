module testsDetectorMod

    use detector_mod
    use testdrive, only : new_unittest, unittest_type, error_type, check
    use constants, only : wp

    implicit none

    contains

    subroutine collect_suite1(testsuite)

        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                new_unittest("Check_hit_circle", hit_circle) &
                ! ! new_unittest("Vector_subtract", vector_sub), &
                ! ! new_unittest("Vector_multiply", vector_mult), &
                ! ! new_unittest("Vector_dot", vector_dot) &
                ]

    end subroutine collect_suite1

    subroutine hit_circle(error)

        use vector_class, only : vector

        type(error_type), allocatable, intent(out) :: error

        type(hit_t) :: hitpoint
        type(circle_dect) :: a
        type(vector) :: pos, dir
        integer :: layer, nbins
        real(kind=wp) :: radius, maxval, val
        logical :: flag

        pos = vector(0._wp, 0._wp, 0._wp)
        dir = vector(1._wp, 0._wp, 0._wp)
        layer = 1
        radius = 0.5
        nbins = 100
        maxval = 100._wp
        a = circle_dect(pos, dir, layer, radius, nbins, maxval, .false.)

        pos = vector(0._wp, 0._wp, 0._wp)
        dir = vector(1._wp, 0._wp, 0._wp)
        val = 1._wp
        layer = 1
        hitpoint = hit_t(pos, dir, val, layer)

        flag = a%check_hit(hitpoint)

        call check(error, flag, .true.)
        if(allocated(error))return

    end subroutine hit_circle

end module testsDetectorMod
program test_detector

    use, intrinsic :: iso_fortran_env, only: error_unit

    use constants, only : wp
    use testdrive, only : run_testsuite, new_testsuite, testsuite_type
    use testsDetectorMod

    implicit none
    
    type(testsuite_type), allocatable :: testsuites(:)
    integer :: i, stat
    character(len=*), parameter :: fmt='("#", *(1x, a))'

    stat = 0

    testsuites = [new_testsuite("Suite: Check Hits", collect_suite1) &
    !               new_testsuite("Suite: Vector .op. scalar", collect_suite2), &
    !               new_testsuite("Suite: Vector functions", collect_suite3), &
    !               new_testsuite("Suite: Vector .op. matrix", collect_suite4) &
                  ]

    do i = 1, size(testsuites)
        write(error_unit, fmt) "Testing:", testsuites(i)%name
        call run_testsuite(testsuites(i)%collect, error_unit, stat)
    end do

    if(stat > 0)then
        write(error_unit, '(i0, 1x, a)')stat, "test(s) failed"
        error stop
    end if
end program test_detector