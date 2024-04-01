subroutine CK_matrix(nat,dist,eta,ckm)

  implicit none

  integer, intent(in) :: nat
  real(8), intent(in) :: eta(nat,2)
  real(8), intent(in) :: dist(nat,nat)

  real(8), intent(out) :: ckm(nat,nat,2,2)

  real(8) :: inv_etaAB(nat,nat,2,2)
  integer :: i, j
  
  ! Average hardness between atom i shell l and atom j shell l'
  ! (eq.32)
  inv_etaAB = 0.d0

  do i = 1,nat-1
    do j = i+1,nat
      inv_etaAB(i,j,1,1) = 1/eta(i,1) + 1/eta(j,1)
      inv_etaAB(i,j,1,2) = 1/eta(i,1) + 1/eta(j,2)
      inv_etaAB(i,j,2,1) = 1/eta(i,2) + 1/eta(j,1)
      inv_etaAB(i,j,2,2) = 1/eta(i,2) + 1/eta(j,2)
    end do
  end do

  inv_etaAB(:,:,:,:) = 0.5d0*inv_etaAB(:,:,:,:)

  ! Mataga-Nishimoto-Ohno-Klopman function
  ! with kg = 2 (eq.31)
  ckm = 0.d0

  do i = 1,nat-1 
    do j = i+1,nat
      ckm(i,j,1,1) = dist(i,j)*dist(i,j) + inv_etaAB(i,j,1,1)*inv_etaAB(i,j,1,1)
      ckm(i,j,1,2) = dist(i,j)*dist(i,j) + inv_etaAB(i,j,1,2)*inv_etaAB(i,j,1,2)
      ckm(i,j,2,1) = dist(i,j)*dist(i,j) + inv_etaAB(i,j,2,1)*inv_etaAB(i,j,2,1)
      ckm(i,j,2,2) = dist(i,j)*dist(i,j) + inv_etaAB(i,j,2,2)*inv_etaAB(i,j,2,2)

      ckm(i,j,:,:) = 1/dsqrt(ckm(i,j,:,:))
    end do
  end do

  print *, 'CKM:'
  do i = 1,nat-1
    do j = i+1,nat
      print *, ckm(i,j,:,:)
    end do
  end do

  stop
  return

end subroutine CK_matrix
