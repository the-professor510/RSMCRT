program mcpolar
!! Entry point for program
    use kernels, only : weight_scatter, pathlength_scatter, pathlength_scatter2

    integer :: num_args, i
    character(len=64), allocatable :: args(:)

    num_args = command_argument_count()
    if(num_args == 0)then
        allocate(args(1))
        !args(1) = "scat_test.toml"
        !args(1) = "validation1.toml"   !Validate hgg scattering, absorption and scattering
        !args(1) = "validation2.toml"   !Validate refractive index mismatch is working
        args(1) = "validation3.toml"   !Validate refractive index mismatch
    else
        allocate(args(num_args))
        do i = 1, num_args
            call get_command_argument(i, args(i))
        end do
    end if
    
    ! call weight_scatter(trim(args(1)))
#ifdef survivalBias
    call pathlength_scatter2(trim(args(1)))
#else
    call pathlength_scatter(trim(args(1)))
#endif

end program