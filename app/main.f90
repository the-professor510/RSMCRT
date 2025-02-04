program mcpolar
!! Entry point for program
    use kernels, only : default_MCRT, escape_Function, inverse_MCRT

    integer :: num_args, i
    character(len=64), allocatable :: args(:)

    num_args = command_argument_count()
    if(num_args == 0)then
        allocate(args(1))
        args(1) = "default.toml"
    else
        allocate(args(num_args))
        do i = 1, num_args
            call get_command_argument(i, args(i))
        end do
    end if
    
#ifdef escapeFunction
    call escape_Function(trim(args(1)))
#elif inverseMCRT
    call inverse_MCRT(trim(args(1)))
#else
    call default_MCRT(trim(args(1)))
#endif


end program