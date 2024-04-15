subroutine CK_matrix(nAt,dist,eta,ckm)

  implicit none

  integer, intent(in)  :: nAt
  real(8), intent(in)  :: eta(nAt,2)
  real(8), intent(in)  :: dist(nAt,nAt)

  real(8), intent(out) :: ckm(nAt,nAt,2,2)

  real(8)              :: inv_etaAB(nAt,nAt,2,2)
  integer              :: A, B
  
  ! Average hardness between atom i shell l and atom j shell l'
  ! (eq.32)

  inv_etaAB = 0.d0

  do A = 1,nAt
    do B = A,nAt
      inv_etaAB(A,B,1,1) = 1.d0/eta(A,1) + 1.d0/eta(B,1)
      inv_etaAB(A,B,1,2) = 1.d0/eta(A,1) + 1.d0/eta(B,2)
      inv_etaAB(A,B,2,1) = 1.d0/eta(A,2) + 1.d0/eta(B,1)
      inv_etaAB(A,B,2,2) = 1.d0/eta(A,2) + 1.d0/eta(B,2)
    end do
  end do

  inv_etaAB(:,:,:,:) = 0.5d0*inv_etaAB(:,:,:,:)

  ! Mataga-Nishimoto-Ohno-Klopman function
  ! with kg = 2 (eq.31)

  ckm = 0.d0

  do A = 1,nAt 
    do B = A,nAt
      ckm(A,B,1,1) = dist(A,B)*dist(A,B) + inv_etaAB(A,B,1,1)*inv_etaAB(A,B,1,1)
      ckm(A,B,1,2) = dist(A,B)*dist(A,B) + inv_etaAB(A,B,1,2)*inv_etaAB(A,B,1,2)
      ckm(A,B,2,1) = dist(A,B)*dist(A,B) + inv_etaAB(A,B,2,1)*inv_etaAB(A,B,2,1)
      ckm(A,B,2,2) = dist(A,B)*dist(A,B) + inv_etaAB(A,B,2,2)*inv_etaAB(A,B,2,2)

      ckm(A,B,:,:) = 1.d0/dsqrt(ckm(A,B,:,:))
    end do
  end do

  return

end subroutine CK_matrix
