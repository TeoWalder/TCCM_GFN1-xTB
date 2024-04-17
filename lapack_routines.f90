subroutine diagonalize_matrix(N,A,e)

  ! Diagonalize a square matrix

  implicit none

  ! Input variables
  integer, intent(in)    :: N
  real(8), intent(inout) :: A(N,N)
  real(8), intent(out)   :: e(N)
  ! Local variables
  integer                :: lwork,info
  real(8), allocatable   :: work(:)

  ! Memory allocation

  allocate(work(3*N))
  lwork = size(work)

  call dsyev('V','U',N,A,N,e,work,lwork,info)
 
  if(info.ne.0) then 
    write(*,*) 'Problem in diagonalize_matrix (dsyev)!!'
    stop
  endif

end subroutine diagonalize_matrix
