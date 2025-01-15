module iarray

!!    The iarray module contains the variables that record the fluence. These are 3D arrays, with roughly the same dimensions as the cart_grid type.
!!    Jmean is the *local* fluence. JmeanGLOBAL is the *global* fluence grid. The global version is the one that is written to disk at the simulations end.

    use constants, only : sp

    implicit none
    !> phase data array
    complex(kind=sp), allocatable :: phasor(:,:,:), phasorGLOBAL(:,:,:)
    !> fluence data array
    real(kind=sp), allocatable :: jmean(:,:,:), jmeanGLOBAL(:,:,:)
    !> dropped packet weight absorption data array
    real(kind=sp), allocatable :: absorb(:,:,:), absorbGLOBAL(:,:,:) 
    !> emission location
    real(kind=sp), allocatable :: emission(:,:,:), emissionGLOBAL(:,:,:)
    !> escape function
    real(kind=sp), allocatable :: escape(:,:,:,:), escapeSymmetry(:,:,:,:)
end module iarray