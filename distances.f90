subroutine distances(n,pos,dist)

  implicit none

  integer, intent(in)  :: n
  real(8), intent(in)  :: pos(n,3)
  real(8), intent(out) :: dist(n,n)

  integer :: i, j

  ! Compute euclidean distance sqrt(x^2 + y^2 + z^2)
  dist = 0.d0

  do i = 1,n-1
    do j = i+1,n
      dist(i,j) = dsqrt((pos(i,1) - pos(j,1))*(pos(i,1) - pos(j,1)) + &
                        (pos(i,2) - pos(j,2))*(pos(i,2) - pos(j,2)) + &
                        (pos(i,3) - pos(j,3))*(pos(i,3) - pos(j,3)))
    end do
  end do

  return

end subroutine distances
