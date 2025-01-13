program mcpolar
!! Entry point for program
    use kernels, only : run_MCRT_Default, run_MCRT_Survival_Bias

    integer :: num_args, i
    character(len=64), allocatable :: args(:)

    num_args = command_argument_count()
    if(num_args == 0)then
        allocate(args(1))
        args(1) = "scat_test.toml"
    else
        allocate(args(num_args))
        do i = 1, num_args
            call get_command_argument(i, args(i))
        end do
    end if
    
#ifdef survivalBias
    call run_MCRT_Survival_Bias(trim(args(1)))
#else
    call run_MCRT_Default(trim(args(1)))
#endif

end program