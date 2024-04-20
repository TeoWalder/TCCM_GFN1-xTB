subroutine repulsion_energy(n,atype,dist,Z,alpha,kf,E)

  implicit none

  integer, intent(in)  :: n, atype(n)
  real(8), intent(in)  :: dist(n,n)
  real(8), intent(in)  :: Z(4), alpha(4), kf
  real(8), intent(out) :: E

  integer :: i, j, A, B

  ! Interatomic atomic repulsion energy (eq.9)

  E = 0.d0
  
  do i = 1,n-1
    do j = i+1,n
      A = atype(i)
      B = atype(j)
      E = E + Z(A)*Z(B)/dist(i,j)*exp(-sqrt(alpha(A)*alpha(B))*dist(i,j)**kf)
    end do
  end do

  return

end subroutine repulsion_energy
