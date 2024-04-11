subroutine SCF(nBas,nOcc,nAt,S,H0,X,shell,CKM,Gamm,E1,E2,E3,Etot,qA)

  implicit none

  ! Input variables
  integer, intent(in)  :: nBas
  integer, intent(in)  :: nOcc
  integer, intent(in)  :: nAt
  real(8), intent(in)  :: S(nBas,nBas)
  real(8), intent(in)  :: H0(nBas,nBas)
  real(8), intent(in)  :: X(nBas,nBas)
  real(8), intent(in)  :: shell(nBas,5)
  real(8), intent(in)  :: CKM(nAt,nAt,2,2)
  real(8), intent(in)  :: Gamm(nAt)
  ! Output variables
  real(8), intent(out) :: Etot
  real(8), intent(out) :: E1, E2, E3
  real(8), intent(out) :: qA(nAt)
  ! Local variables
  integer, parameter   :: maxSCF = 60
  real(8), parameter   :: thresh = 1.d-5
  real(8), parameter   :: qthresh = 1.d-3
  real(8)              :: Conv
  real(8)              :: DqA, DqS
  real(8), allocatable :: C(:,:)
  real(8), allocatable :: Cp(:,:)
  real(8), allocatable :: P(:,:)
  real(8), allocatable :: F(:,:), Fp(:,:)
  real(8), allocatable :: qS_old(:,:), qS_new(:,:)
  real(8), allocatable :: qA_old(:,:), qA_new(:,:)
  ! Counters
  integer              :: nSCF
  integer              :: mu, nu, si, la
  integer              :: A, B, l, lp

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|       Self-Consistent Field Calculation      |'
  write(*,*)'************************************************'
  write(*,*)

  ! Memory allocation

  allocate(C(nBas,nBas), Cp(nBas,nBas), P(nBas,nBas), F(nBas,nBas), Fp(nBas,nBas))
  allocate()

  ! Guess coefficients and eigenvalues

  F(:,:) = H0(:,:)

  qA_old = 0.d0
  qA = 0.d0

  ! Initialization

  nSCF = 0
  Conv = 1d0

!-----------------------------------------------------------
! SCF loop
!-----------------------------------------------------------

  write(*,*)
  write(*,*)'----------------------------------------------------'
  write(*,*)'| SCF calculation                                  |'
  write(*,*)'----------------------------------------------------'
  write(*,'(1x,a1,1x,a3,1x,a1,1x,a16,1x,a1,1x,a10,1X,a1,1x,a10,1x,a1,1x)') &
            '|','#','|','DqS change','|','DqS change','|'
  write(*,*)'----------------------------------------------------'

  do while(conv.gt.thresh.and.nSCF.lt.maxSCF)

    ! Increment 

    nSCF = nSCF + 1

    ! F' and C'

    Fp = matmul(transpose(X),matmul(F,X))

    Cp = Fp

    call diagonalize_matrix(nBas,Cp,e)

    ! Coefficients Matrix

    C = matmul(X,Cp)

    ! Update Density Matrix

    do mu = 1,nBas
      do nu = 1,nBas
        P(mu,nu) = 0.d0
        do si = 1,nOcc
          P(mu,nu) = P(mu,nu) + 2.d0*C(mu,si)*C(nu,si)
        end do
      end do
    end do

    ! Compute Shell Charges (eq.29)

    do A = 1,nAt
      do l = 1,2
        qS_new(A,l) = shell(2*(A-1)+l,3)
        do nu = 1,nBas
          qS_new(A,l) = qS_new(A,l) - S(2*(A-1)+l,nu)*P(2*(A-1)+l,nu)
        end do
      end do
    end do
    
    ! Compute Atomic Charges (eq.20)

    do A = 1,nAt
      qA_new(A) = sum(qS_new(A,:))
    end do

    ! Update F

    F(:,:) = H0(:,:)

    do A = 1,nAt-1
      do l = 1,2

        shift_sh_A = sum(CKM(A,:,l,:)*qS_new(:,:))
        shift_at_A = Gamm(A)*qA_new(A)*qA_new(A)
        do B = A+1,nAt
          do lp = 1,2
            shift_sh_B = sum(CKM(B,:,lp,:)*qS_new(:,:))
            shift_at_B = Gamm(B)*qA_new(B)
    
            mu = 2*(A-1) + l
            nu = 2*(B-1) + lp
  
            F(mu,nu) = F(mu,nu) - 0.5d0*S(mu,nu)  &
                       *(shift_sh_A + shift_sh_B + shift_at_A + shift_at_B)
          end do
        end do

      end do
    end do

    ! Convergency Criterium

    DqS = abs(max(qS_new - qS_old))
    DqA = abs(max(qA_new - qA_old))
    
    conv = max(Dqs,DqA)
       
    ! New Charges
    
    if (DqS.ge.qthresh) then
      qS_old = qS_new + 0.4d0*DqS   ! dumped
    else 
      qS_old = qS_new               ! not dumped
    end if
    
    if (DqA.ge.qthresh) then
      qA_old = qA_new + 0.4d0*DqA   ! dumped
    else
      qA_old = qA_new               ! not dumped
    end if

    ! Compute energies



    ! Dump results

    write(*,'(1x,a1,1x,i3,1x,a1,1x,f16.10,1x,a1,1x,f10.6,1x,a1,1x)') &
            '|',nSCF,'|',DqS,'|',DqA,'|'
 
  enddo
  write(*,*)'----------------------------------------------------'
!-----------------------------------------------------------
! End of SCF loop
!-----------------------------------------------------------

! Did it actually converge?

  if(nSCF.eq.maxSCF) then

    write(*,*)
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)'                 Convergence failed                 '
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)

    stop

  endif

end subroutine SCF
