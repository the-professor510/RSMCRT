module interpolate
    use constants, only : wp

    implicit none
    
    private
    public :: trilinearInterpolate, bilinearInterpolate, linearInterpolate

    contains

    subroutine trilinearInterpolate(corners, point)
        !take in 8 points in a regular grid and point contained within these

        !>corners (x,y,z,(pos,value))
        real(kind=wp), intent(in) :: corners(2,2,2,4)
        !> point we want to get the value you of by interpolation
        real(kind=wp), intent(inout):: point(4)

        real(kind = wp) :: c00, c01, c10, c11, c0, c1, c
        real(kind = wp) :: xd, yd, zd

        xd = (point(1) - corners(1,1,1,1))/(corners(2,1,1,1) - corners(1,1,1,1))
        yd = (point(2) - corners(1,1,1,2))/(corners(1,2,1,2) - corners(1,1,1,2))
        zd = (point(3) - corners(1,1,1,3))/(corners(1,1,2,3) - corners(1,1,1,3))

        !interpolate along x
        c00 = corners(1,1,1,4) * (1 - xd) + corners(2,1,1,4)*xd
        c01 = corners(1,1,2,4) * (1 - xd) + corners(2,1,2,4)*xd
        c10 = corners(1,2,1,4) * (1 - xd) + corners(2,2,1,4)*xd
        c11 = corners(1,2,2,4) * (1 - xd) + corners(2,2,2,4)*xd

        !interpolate along y
        c0 = c00*(1-yd) + c10*yd
        c1 = c01*(1-yd) + c11*yd

        !interpolate along z
        c = c0*(1-zd) + c1*zd

        !this is our result
        point(4) = c

    end subroutine trilinearInterpolate

    subroutine bilinearInterpolate(corners, point)
        !take in 4 points in a regular grid and point contained within these

        !>corners (x,y,(pos,value))
        real(kind=wp), intent(in) :: corners(2,2,3)
        !> point we want to get the value you of by interpolation
        real(kind=wp), intent(inout):: point(3)

        real(kind = wp) :: c0, c1, c
        real(kind = wp) :: xd, yd

        xd = (point(1) - corners(1,1,1))/(corners(2,1,1) - corners(1,1,1))
        yd = (point(2) - corners(1,1,2))/(corners(1,2,2) - corners(1,1,2))

        !interpolate along x
        c0 = corners(1,1,3)*(1-xd) + corners(2,1,3)*xd
        c1 = corners(1,2,3)*(1-xd) + corners(2,2,3)*xd

        !interpolate along y
        c = c0*(1-yd) + c1*yd

        !this is our result
        point(3) = c


    end subroutine bilinearInterpolate

    subroutine linearInterpolate(corners, point)
        !take in 2 points and a point contained wihtin these linearly interpolate between them

        !>corners (x, (pos,value))
        real(kind=wp), intent(in) :: corners(2,2)
        !> point we want to get the value you of by interpolation
        real(kind=wp), intent(inout):: point(2)

        real(kind = wp) ::  c
        real(kind = wp) :: xd

        xd = (point(1) - corners(1,1))/(corners(2,1) - corners(1,1))

        !interpolate along x
        c = corners(1,2)*(1-xd) + corners(2,2)*xd

        !this is our result
        point(2) = c

    end subroutine linearInterpolate

end module interpolate