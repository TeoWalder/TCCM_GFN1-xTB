subroutine CK_matrix()

  implicit none

  

  inv_etaAB(:,:,:,:) = 0.5d0*inv_etaAB(:,:,:,:)

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

  return

subroutine CK_matrix
