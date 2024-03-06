module module_chem

  implicit none

contains

!---------- COMPUTE DISTANCE MATRIX -------------------------------------------!
  subroutine distance(n,pos,dist)

    implicit none

    integer, intent(in)  :: n
    real(8), intent(in)  :: pos(n,3)
    real(8), intent(out) :: dist(n,n)

    integer     :: i, j

    dist(:,:) = 0.d0

    do i = 1,n-1
      do j = i+1,n
        dist(i,j) = dsqrt((pos(i,1) - pos(j,1))*(pos(i,1) - pos(j,1)) + &
                          (pos(i,2) - pos(j,2))*(pos(i,2) - pos(j,2)) + &
                          (pos(i,3) - pos(j,3))*(pos(i,3) - pos(j,3)))
      end do
    end do

    return

  end subroutine distance

!---------- REPULSION ENERGY --------------------------------------------------!
  subroutine repulsion_energy(n,atype,dist,Z,alpha,kf,E)

    implicit none

    integer, intent(in)  :: n, atype(n)
    real(8), intent(in)  :: dist(n,n)
    real(8), intent(in)  :: Z(4), alpha(4), kf
    real(8), intent(out) :: E

    integer :: i, j, A, B

    E = 0.d0
    
    do i = 1,n-1
      do j = i+1,n
        A = atype(i)
        B = atype(j)
        E = E + Z(A)*Z(B)/dist(i,j)*dexp(-dsqrt(alpha(A)*alpha(B))*dist(i,j)**kf)
      end do
    end do

    return

  end subroutine repulsion_energy

end module module_chem
