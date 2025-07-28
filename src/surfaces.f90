module surfaces
!! Contains the routines that handle reflection, and refraction via the Fresnel equations.
    
    use vector_class, only : vector
    use constants,    only : wp

    implicit none

    private
    public :: reflect_refract

    contains

    subroutine reflect_refract(I, N, n1, n2, rflag, ri)
    !! wrapper routine for fresnel calculation

        use random, only : ran2

        !> incident vector
        type(vector),  intent(INOUT) :: I
        !> normal vector
        type(vector),  intent(INOUT) :: N
        !> refractive indices
        real(kind=wp), intent(IN)    :: n1, n2
        real(kind=wp), intent(OUT)   :: Ri
        !> reflection flag
        logical,       intent(OUT)   :: rflag

        rflag = .FALSE.

        !draw random number, if less than fresnel coefficents, then reflect, else refract
        Ri = fresnel(I, N, n1, n2)
        if(ran2() <= Ri)then
            call reflect(I, N)
            rflag = .true.
        else
            call refract(I, N, n1/n2)
        end if

    end subroutine reflect_refract

    subroutine reflect(I, N)
    !! get vector of reflected photon

        !> incident vector
        type(vector), intent(INOUT) :: I
        !> normal vector
        type(vector), intent(IN)    :: N

        type(vector) :: R

        R = I - 2._wp * (N .dot. I) * N
        I = R

    end subroutine reflect

    subroutine refract(I, N, eta)
    !! get vector of refracted photon

        !> incident vector
        type(vector),  intent(INOUT) :: I
        !> normal vector
        type(vector),  intent(IN)    :: N
        !> \(\eta = \frac{n_1}{n_2}\)
        real(kind=wp), intent(IN)    :: eta

        type(vector)  :: T, Ntmp
        real(kind=wp) :: c1, c2

        Ntmp = N

        c1 = (Ntmp .dot. I)
        if(c1 < 0._wp)then
            c1 = -c1
        else
            Ntmp = (-1._wp) * N
        end if
        c2 = sqrt(1._wp - (eta)**2 * (1._wp-c1**2))

        T = eta*I + (eta * c1 - c2) * Ntmp 

        I = T

    end subroutine refract

    function fresnel(I, N, n1, n2) result (tir)
    !! calculates the fresnel coefficents

        use ieee_arithmetic, only : ieee_is_nan

        !> reffractive indicies
        real(kind=wp), intent(IN) :: n1, n2
        !> incident vector
        type(vector),  intent(IN) :: I
        !> Normal vector
        type(vector),  intent(IN) :: N

        real(kind=wp) :: cos_i, sin_i, sin_t, cos_t, tir, f1, f2

        cos_i = abs(I .dot. N)
        if(cos_i<0.0_wp)cos_i=-1._wp*cos_i
        if(cos_i>1.0_wp)cos_i=1.0_wp
        sin_i = sqrt(1._wp - cos_i * cos_i)
        sin_t = n1/n2 * sin_i

        if(sin_t > 1._wp)then
            !total internal reflection occurs
            tir = 1.0_wp
            return

        elseif(cos_i == 1._wp)then
            !the packet is perpendicular to the surface
            tir = ((n1-n2)**2)/((n1+n2)**2)
            return
        else
            sin_t = (n1/n2)*sin_i
            cos_t = sqrt(1._wp - sin_t * sin_t)
            f1 = abs((n1*cos_i - n2*cos_t) / (n1*cos_i + n2*cos_t))**2
            f2 = abs((n1*cos_t - n2*cos_i) / (n1*cos_t + n2*cos_i))**2
            
            tir = 0.5_wp * (f1 + f2)
        end if
        if(ieee_is_nan(tir) .or. tir > 1._wp .or. tir < 0._wp) then
            print*,'TIR: ', tir, f1, f2, cos_i,sin_i,cos_t,sin_t
            return
        end if


    end function fresnel
end module surfaces