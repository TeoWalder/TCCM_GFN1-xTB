subroutine SCF(nBas,nO,S,T,V,H0,ERI,X,ENuc,EHF,e,c)

  implicit none

  ! Input variables
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nO
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: H0(nBas,nBas) 
  double precision,intent(in)   :: X(nBas,nBas) 
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ENuc
  ! Local variables
  integer,parameter             :: maxSCF = 64
  double precision,parameter    :: thresh = 1d-5
  integer                       :: nSCF
  double precision              :: Conv
  double precision              :: Gap
  double precision              :: ET,EV,EJ
  double precision              :: EK
  double precision,allocatable  :: cp(:,:)
  double precision,allocatable  :: P(:,:)
  double precision,allocatable  :: J(:,:)
  double precision,allocatable  :: K(:,:)
  double precision,allocatable  :: F(:,:),Fp(:,:)
  double precision,allocatable  :: error(:,:)
  double precision              :: trace_matrix(nBas,nBas)
  ! Counters
  integer                       :: mu, nu, si, la
  ! Output variables
  double precision,intent(out)  :: EHF
  double precision,intent(out)  :: e(nBas)
  double precision,intent(out)  :: c(nBas,nBas)

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|       Self-Consistent Field Calculation      |'
  write(*,*)'************************************************'
  write(*,*)

  ! Memory allocation

  allocate(cp(nBas,nBas),P(nBas,nBas),      &
           J(nBas,nBas),K(nBas,nBas),F(nBas,nBas),Fp(nBas,nBas), &
           error(nBas,nBas))

  ! Guess coefficients and eigenvalues

  F(:,:) = H0(:,:)
  c(:,:) = 0.d0

  ! Initialization

  nSCF = 0
  Conv = 1d0

!------------------------------------------------------------------------
! Main SCF loop
!------------------------------------------------------------------------

  write(*,*)
  write(*,*)'----------------------------------------------------'
  write(*,*)'| SCF calculation                                  |'
  write(*,*)'----------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A16,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X)') &
            '|','#','|','Energies','|','Conv','|','HL Gap','|'
  write(*,*)'----------------------------------------------------'

  do while(Conv.gt.thresh.and.nSCF.lt.maxSCF)

    ! Increment 

    nSCF = nSCF + 1

    ! Update Density Matrix

    do mu = 1,nBas
      do nu = 1,nBas
        P(mu,nu) = 0.d0
        do si = 1,nO
          P(mu,nu) = P(mu,nu) + 2.d0*c(mu,si)*c(nu,si)
        end do
      end do
    end do

    
    ! J, K. and F

    do mu = 1,nBas
      do nu = 1,nBas
        J(mu,nu) = 0.d0
        do la = 1,nBas
          do si = 1,nBas
            J(mu,nu) = J(mu,nu) + P(la,si)*ERI(mu,la,nu,si)
          end do
        end do
      end do
    end do

    do mu = 1,nBas
      do nu = 1,nBas
        K(mu,nu) = 0.d0
        do la = 1,nBas
          do si = 1,nBas
            K(mu,nu) = K(mu,nu) - 0.5d0*P(la,si)*ERI(mu,nu,la,si)
          end do
        end do
      end do
    end do

    F(:,:) = Hc(:,:) + J(:,:) + K(:,:)

    ! F' and C'

    Fp = matmul(transpose(X),matmul(F,X))

    cp = Fp

    call diagonalize_matrix(nBas,cp,e)

    ! Coefficients Matrix

    c = matmul(X,cp)

    ! Convergency Criterium

    if (nSCF.eq.1) cycle
    error = matmul(F,matmul(P,S)) - matmul(S,matmul(P,F))
    Conv = maxval(dabs(error))

    ! Step Energy and HL gap

    trace_matrix = matmul(P,F+Hc)
    EHF = 0.d0
    do mu = 1,nBas
      EHF = EHF + trace_matrix(mu,mu)
    end do
    EHF = EHF*0.5d0

    Gap = e(nO+1) - e(nO)

    ! Dump results

    write(*,'(1X,A1,1X,I3,1X,A1,1X,F16.10,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X)') &
      '|',nSCF,'|',EHF+ENuc,'|',Conv,'|',Gap,'|'
 
  enddo
  write(*,*)'----------------------------------------------------'
!------------------------------------------------------------------------
! End of SCF loop
!------------------------------------------------------------------------

! Did it actually converge?

  if(nSCF.eq.maxSCF) then

    write(*,*)
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)'                 Convergence failed                 '
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)

    stop

  endif

  ! Compute final HF energy

  call print_RHF(nBas,nO,e,C,ENuc,ET,EV,EJ,EK,EHF)

end subroutine SCF
