module module_chem

  implicit none

contains

!---------- COMPUTE DISTANCE MATRIX -------------------------------------------!
  subroutine distance(n,nc,pos,dist)

    implicit none

    integer    , intent(in)  :: n
    integer(16), intent(in)  :: nc
    real(8)    , intent(in)  :: pos(n,3)
    real(8)    , intent(out) :: dist(nc,3)

    integer     :: i, j
    integer(16) :: k

    k = 1_16
    do i = 1,n
      do j = i+1,n
        dist(k,1) = dble(i)
        dist(k,2) = dble(j)
        dist(k,3) = dsqrt((pos(i,1) - pos(j,1))*(pos(i,1) - pos(j,1)) + &
                          (pos(i,2) - pos(j,2))*(pos(i,2) - pos(j,2)) + &
                          (pos(i,3) - pos(j,3))*(pos(i,3) - pos(j,3)))
        k = k + 1_16
      end do
    end do

    return

  end subroutine distance

!---------- REPULSION ENERGY --------------------------------------------------!
  subroutine repulsion_energy(n,nc,atype,dist,Z,alpha,kf,E)

    implicit none

    integer(16), intent(in)  :: nc
    integer    , intent(in)  :: n, atype(n)
    real(8)    , intent(in)  :: dist(nc,3)
    real(8)    , intent(in)  :: Z(4), alpha(4), kf
    real(8)    , intent(out) :: E

    integer :: i, A, B

    E = 0.d0
    
    do i = 1,nc
      A = atype(int(dist(i,1)))
      B = atype(int(dist(i,2)))
      E = E + Z(A)*Z(B)/dist(i,3)*dexp(-dsqrt(alpha(A)*alpha(B))*dist(i,3)**kf)
    end do

    return

  end subroutine repulsion_energy

end module module_chem
