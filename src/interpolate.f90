module interpolate
    use constants, only : wp

    implicit none
    
    private
    public :: trilinearInterpolate, bilinearInterpolate, linearInterpolate, cylTrilinearInterpolate, cylBilinearInterpolate

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

    subroutine cylTrilinearInterpolate(corners, point)

        !>corners (r,theta,z,(pos,value))
        real(kind=wp), intent(in) :: corners(2,2,2,4)
        !> point we want to get the value you of by interpolation
        real(kind=wp), intent(inout):: point(4)

        !using the fact that bilinear interpolation be seen as a weighted average, do this for cylindrical coordinates

        real(kind = wp) :: c
        real(kind = wp) :: a00, a01, a10, a11
        real(kind = wp) :: volume, v000, v001, v010, v011, v100, v101, v110, v111

        !calc total volume = 0.5 * (thetahigh - thetalow) * (rhigh^2 - rlow^2) * z
        volume = 0.5_wp * (corners(1,2,1,2) - corners(1,1,1,2)) * (corners(2,1,1,1)**2 - corners(1,1,1,1)**2) * & 
                (corners(1,1,2,3) - corners(1,1,1,3))

        !calc areas = 0.5 * (thetahigh - thetalow) * (rhigh^2 - rlow^2)
        a00 = 0.5_wp * (corners(2,2,1,2) - point(2)) * (corners(2,2,1,1)**2 - point(1)**2) !enclosed by p and h theta and h r
        a01 = 0.5_wp * (point(2) - corners(2,1,1,2)) * (corners(2,1,1,1)**2 - point(1)**2) !enclosed by p and l theta and h r
        a10 = 0.5_wp * (corners(1,2,1,2) - point(2)) * (point(1)**2 - corners(1,2,1,1)**2) !enclosed by p and h theta and l r
        a11 = 0.5_wp * (point(2) - corners(1,1,1,2)) * (point(1)**2 - corners(1,1,1,1)**2) !enclosed by p and l theta and l r

        !calc volume weightings
        v000 = a00 * (corners(1,1,2,3) - point(3)) / volume!volume between point and high z
        v001 = a00 * (point(3) - corners(1,1,1,3)) / volume!volume between point and low z
        v010 = a01 * (corners(1,1,2,3) - point(3)) / volume
        v011 = a01 * (point(3) - corners(1,1,1,3)) / volume
        v100 = a10 * (corners(1,1,2,3) - point(3)) / volume
        v101 = a10 * (point(3) - corners(1,1,1,3)) / volume
        v110 = a11 * (corners(1,1,2,3) - point(3)) / volume
        v111 = a11 * (point(3) - corners(1,1,1,3)) / volume

        !interpolate by finding the volume weighted average
        c = v000 * corners(1,1,1,4) + & 
            v001 * corners(1,1,2,4) + & 
            v010 * corners(1,2,1,4) + & 
            v011 * corners(1,2,2,4) + & 
            v100 * corners(2,1,1,4) + & 
            v101 * corners(2,1,2,4) + & 
            v110 * corners(2,2,1,4) + & 
            v111 * corners(2,2,2,4)

        !this is our result
        point(4) = c

    end subroutine cylTrilinearInterpolate

    subroutine cylBilinearInterpolate(corners, point)

        !>corners (r,theta,(pos,value))
        real(kind=wp), intent(in) :: corners(2,2,3)
        !> point we want to get the value you of by interpolation
        real(kind=wp), intent(inout):: point(3)

        !using the fact that bilinear interpolation be seen as a weighted average, do this for cylindrical coordinates

        real(kind = wp) :: c
        real(kind = wp) :: a00, a01, a10, a11
        real(kind = wp) :: area, w00, w01, w10, w11

        !calc total area = 0.5 * (thetahigh - thetalow) * (rhigh^2 - rlow^2)
        area = 0.5_wp * (corners(1,2,2) - corners(1,1,2)) * (corners(2,1,1)**2 - corners(1,1,1)**2) 

        !calc areas = 0.5 * (thetahigh - thetalow) * (rhigh^2 - rlow^2)
        a00 = 0.5_wp * (corners(2,2,2) - point(2)) * (corners(2,2,1)**2 - point(1)**2) !enclosed by p and h theta and h r
        a01 = 0.5_wp * (point(2) - corners(2,1,2)) * (corners(2,1,1)**2 - point(1)**2) !enclosed by p and l theta and h r
        a10 = 0.5_wp * (corners(1,2,2) - point(2)) * (point(1)**2 - corners(1,2,1)**2) !enclosed by p and h theta and l r
        a11 = 0.5_wp * (point(2) - corners(1,1,2)) * (point(1)**2 - corners(1,1,1)**2) !enclosed by p and l theta and l r

        !calc weightings
        w00 = a00/area
        w01 = a01/area
        w10 = a10/area
        w11 = a11/area
        
        !interpolate by finding the area weighted average
        c = w00 * corners(1,1,3) + & 
            w01 * corners(1,2,3) + & 
            w10 * corners(2,1,3) + & 
            w11 * corners(2,2,3) 
            
        !this is our result
        point(3) = c

    end subroutine cylBilinearInterpolate

end module interpolate