module gridMod
    !! This module defines the cartesian grid type (cart_grid) and associated routines.

    !! The cart_grid type contains information related to the grid used to record the fluence. This includes the number of voxels in each cardinal direction (nxg, nyg, nzg), the **half** size of the grid in each direction (xmax, ymax, zmax), and the locations of the voxels walls in each direction (xface, yface, zface).
    !! The type-bound function get_voxel takes a position (vector) and returns the voxel the position falls in. 
    !! 
    !! Init_grid initialises a cart_grid instance.

    use constants, only : wp

    implicit none

    !! Grid class
    type :: cart_grid
        !> number of voxels in each cardinal direction for fluence grid
        integer       :: nxg, nyg, nzg
        !> half size of each dimension in fluence grid. 
        real(kind=wp) :: xmax, ymax, zmax
        !> Delta is the round off for near voxel cell walls
        real(kind=wp) :: delta
        !> position of each cell wall in fluence grid
        real(kind=wp), allocatable :: xface(:), yface(:), zface(:)
        contains
        procedure :: get_voxel => get_voxel_cart
    end type cart_grid

    interface cart_grid
        module procedure init_grid_cart
    end interface cart_grid

    type :: cyl_grid
        !> number of voxels in each cardinal direction for escape grid
        integer         :: nrg, ntg, nzg
        !> half size of each dimension in escape grid
        real(kind = wp) :: rmax, tmax, zmax
        
        contains 
        procedure :: get_voxel => get_voxel_cyl
    end type cyl_grid

    interface cyl_grid
        module procedure init_grid_cyl
    end interface cyl_grid


    public  :: cart_grid, cyl_grid, init_grid_cart, init_grid_cyl
    private

    contains

    function get_voxel_cart(this, pos) result(res)
        !! get current voxel the photon packet is in
        use vector_class
        
        !> grid class
        class(cart_grid)         :: this
        !> current vector position of photon packet
        type(vector), intent(IN) :: pos
    
        integer :: res(3)

        res(1) = floor(this%nxg*(pos%x + this%xmax)/(2._wp*this%xmax))+1
        res(2) = floor(this%nyg*(pos%y + this%ymax)/(2._wp*this%ymax))+1
        res(3) = floor(this%nzg*(pos%z + this%zmax)/(2._wp*this%zmax))+1

        if(res(1) < 1 .or. res(1) > this%nxg) then
            res(1) = -1
        end if
        if(res(2) < 1 .or. res(2) > this%nyg) then
            res(2) = -1
        end if
        if(res(3) < 1 .or. res(3) > this%nzg) then
            res(3) = -1
        end if

        !what if we are beyond the max or min values, need to return a different value

    end function get_voxel_cart

    function get_voxel_cyl(this, pos) result(res)
        !! get current voxel the photon packet is in
        use vector_class
        use constants, only: PI
        
        !> grid class
        class(cyl_grid)         :: this
        !> current vector position of photon packet
        type(vector), intent(IN) :: pos
    
        integer :: res(3)

        res(1) = floor(this%nrg*(sqrt(pos%x**2 + pos%y**2)/this%rmax))+1
        if (pos%y ==0 .and. pos%x == 0) then
            res(2) = 1
        else
            res(2) = floor(this%ntg*((atan2(pos%y, pos%x)+ PI)/this%tmax)) + 1
        end if
        res(3) = floor(this%nzg*(pos%z+this%zmax)/(2._wp*this%zmax)) + 1

        if(res(1) < 1 .or. res(1) > this%nrg) then
            res(1) = -1
        end if
        if(res(2) < 1 .or. res(2) > this%ntg) then
            res(2) = -1
        end if
        if(res(3) < 1 .or. res(3) > this%nzg) then
            res(3) = -1
        end if

    end function get_voxel_cyl

    type(cart_grid) function init_grid_cart(nxg, nyg, nzg, xmax, ymax, zmax)
    !! setup grid
        !> number of voxels in each cardinal direction for fluence grid
        integer,       intent(IN) :: nxg, nyg, nzg
        !> half size of each dimension in fluence grid. 
        real(kind=wp), intent(IN) :: xmax, ymax, zmax
        
        integer :: i

        init_grid_cart%nxg = nxg
        init_grid_cart%nyg = nyg
        init_grid_cart%nzg = nzg

        init_grid_cart%xmax = xmax
        init_grid_cart%ymax = ymax
        init_grid_cart%zmax = zmax

        allocate(init_grid_cart%xface(nxg + 1), init_grid_cart%yface(nyg + 1), init_grid_cart%zface(nzg + 2))

        init_grid_cart%xface = 0._wp
        init_grid_cart%yface = 0._wp
        init_grid_cart%zface = 0._wp

        ! Set small distance for use in optical depth integration routines 
        ! for roundoff effects when crossing cell walls
        init_grid_cart%delta = 1.e-8_wp * min(((2._wp*xmax)/nxg), ((2._wp*ymax)/nyg), ((2._wp*zmax)/nzg))


        do i = 1, nxg + 1
            init_grid_cart%xface(i) = (i - 1) * 2._wp * xmax/nxg
        end do

        do i = 1, nyg + 1
            init_grid_cart%yface(i) = (i - 1) * 2._wp * ymax/nyg
        end do

        do i = 1, nzg + 2
            init_grid_cart%zface(i) = (i - 1) * 2._wp * zmax/nzg
        end do

    end function init_grid_cart

    type(cyl_grid) function init_grid_cyl(nrg, ntg, nzg, rmax, tmax, zmax)
        !! setup grid
        !> number of voxels in each cardinal direction for fluence grid
        integer,       intent(IN) :: nrg, ntg, nzg
        !> half size of each dimension in fluence grid. 
        real(kind=wp), intent(IN) :: rmax, tmax, zmax

        init_grid_cyl%nrg = nrg
        init_grid_cyl%ntg = ntg
        init_grid_cyl%nzg = nzg

        init_grid_cyl%rmax = rmax
        init_grid_cyl%tmax = tmax
        init_grid_cyl%zmax = zmax

    end function init_grid_cyl
end module gridMod