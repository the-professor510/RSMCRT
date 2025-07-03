module Interfaces

    use constants, only : sp
    
    implicit none
    
    interface stdlib_sposv
        subroutine sposv(uplo, n, nrhs, A, Lda, B, Ldb, info)
            character, intent(in) :: uplo
            integer, intent(in) :: n
            integer, intent(in) :: nrhs
            real(kind=4), intent(inout) :: A(1:Lda,*)
            integer, intent(in) :: Lda
            real(kind=4), intent(inout) :: B(1:Lda,*)
            integer, intent(in) :: Ldb
            integer, intent(out) :: info

        end subroutine sposv
    end interface stdlib_sposv
end module Interfaces