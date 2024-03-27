subroutine inv_sqrt_S(nBas,S,X)

  implicit none

  integer, intent(in)   :: nBas
  real(8), intent(in)   :: S(nBas,nBas)
  real(8), intent(out)  :: X(nBas,nBas)

  real(8)               :: U(nBas,nBas)
  real(8)               :: e(nBas)
  real(8)               :: thresh
  integer               :: i, j, k

  ! diagonalize S

  U(:,:) = S(:,:)
  call diagonalize_matrix(nBas,U,e)

  ! take the inverse of the square root of e

  thresh = 1.d-6
  do i = 1,nBas
    if (e(i).gt.thresh) then
      e(i) = 1.d0/sqrt(e(i))
    else
      write(*,"(20a,f10.6,20a)") "S matrix eigenvalue", e(i), "below the trashold"
      stop
    end if
  end do

  ! transform back

  do i = 1,nBas
    do j = 1,nBas
      X(i,j) = 0.d0
      do k = 1,nBas
        X(i,j) = X(i,j) + U(i,k)*e(k)*U(j,k)
      end do
    end do
  end do

end subroutine inv_sqrt_S
