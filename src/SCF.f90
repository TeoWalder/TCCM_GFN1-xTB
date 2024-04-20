subroutine SCF(nBas,nOcc,nAt,nSh,atype,S,H0,X,shell,CKM,Gamm,Eel,qA,ev_hartree)

  implicit none

  ! Input variables
  integer, intent(in)  :: nBas
  integer, intent(in)  :: nOcc
  integer, intent(in)  :: nAt
  integer, intent(in)  :: nSh
  integer, intent(in)  :: atype(nAt)
  real(8), intent(in)  :: ev_hartree
  real(8), intent(in)  :: S(nBas,nBas)
  real(8), intent(in)  :: H0(nBas,nBas)
  real(8), intent(in)  :: X(nBas,nBas)
  real(8), intent(in)  :: shell(nSh,5)
  real(8), intent(in)  :: CKM(nAt,nAt,2,2)
  real(8), intent(in)  :: Gamm(nAt)
  ! Output variables
  real(8), intent(out) :: Eel
  real(8), intent(out) :: qA(nAt)
  ! Local variables
  integer, parameter   :: maxSCF = 60
  real(8), parameter   :: thresh = 1.d-7
  real(8), parameter   :: dump_thresh = 1.d-3
  real(8)              :: Conv
  real(8)              :: DqA, DqS
  real(8)              :: Eel_old, E1, E2, E3
  real(8), allocatable :: C(:,:)
  real(8), allocatable :: Cp(:,:)
  real(8), allocatable :: P(:,:)
  real(8), allocatable :: F(:,:), Fp(:,:)
  real(8), allocatable :: e(:)
  real(8), allocatable :: qS_old(:,:), qS(:,:)
  real(8), allocatable :: qA_old(:)
  real(8)              :: shift_at_A, shift_at_B
  real(8)              :: shift_sh_A, shift_sh_B
  ! Counters
  integer              :: nSCF
  integer              :: mu, nu, si
  integer              :: ish_A, ish_B
  integer              :: A, B, Ca, l, lp, lpp

  write(*,'(x,a74)') '___________________________________________________________________________'
  write(*,*)
  write(*,'(x,a74)') '                         SELF-CONSISTENT FIELD                             '
  write(*,'(x,a74)') '___________________________________________________________________________'
  write(*,*)

  ! Memory allocation

  allocate(C(nBas,nBas), Cp(nBas,nBas), P(nBas,nBas), F(nBas,nBas), Fp(nBas,nBas), e(nBas))
  allocate(qS_old(nAt,2), qS(nAt,2), qA_old(nAt))

  ! Initial Guess

  F(:,:) = H0(:,:)
  qS_old = 0.d0
  qA_old = 0.d0
  Eel_old = 0.d0

  ! Initialization

  nSCF = 0
  Conv = 1.d0

!------ SCF loop ------------------------------------------!

  write(*,'(x,a1,x,a3,x,a1,x,a13,x,a1,x,a10,x,a1,x,a10,x,a1,x,a10,x,a1,x,a10,x,a1,x,a10,x,a1,x)') &
             '|','#','|','Elec. Energy','|','E(1)','|','E(2)','|','E(3)','|','DE','|'
  write(*,*) '---------------------------------------------------------------------------'   

  do while(conv.gt.thresh.and.nSCF.lt.maxSCF)

    ! Increment 

    nSCF = nSCF + 1

    ! F' and C'

    Fp = matmul(transpose(X),matmul(F,X))

    Cp = Fp

    call diagonalize_matrix(nBas,Cp,e)

    ! Coefficients Matrix

    C = matmul(X,Cp)

    ! Cmpute Density Matrix

    P(:,:) = 0.d0
    do mu = 1,nBas
      do nu = 1,nBas 
        do si = 1,nOcc
          P(mu,nu) = P(mu,nu) + 2.d0*C(mu,si)*C(nu,si)
        end do
      end do
    end do

    ! Compute Charges

    ! (eq.29)
    do A = 1,nAt
      do l = 1,2
        qS(A,l) = shell(2*(A-1)+l,3)

        ish_A = 0
        if (l.eq.2.and.atype(A).ne.1) ish_A = 2

        do mu = int(shell(2*(A-1)+l,4)), int(shell(2*(A-1)+l,4)) + ish_A
          do nu = 1,nBas

            qS(A,l) = qS(A,l) - S(mu,nu)*P(mu,nu)

          end do
        end do

      end do
    end do

    ! (eq.20)
    do A = 1,nAt
      qA(A) = sum(qS(A,:))
    end do

    ! Compute energies

    ! (eq.37) 
    E1 = 0.d0
    do mu = 1,nBas
      do nu = 1,nBas
        E1 = E1 + P(mu,nu)*H0(mu,nu)
      end do
    end do

    ! (eq.19)
    E2 = 0.d0
    do A = 1,nAt
      do B = 1,nAt
        do l = 1,2
          do lp = 1,2
            E2 = E2 + CKM(A,B,l,lp)*qS(A,l)*qS(B,lp)
          end do
        end do
      end do
    end do
    E2 = E2*0.5d0

    ! (eq.20)
    E3 = 0.d0
    do A = 1,nAt
      E3 = E3 + Gamm(A)*qA(A)*qA(A)*qA(A)
    end do
    E3 = E3/3.d0

    ! Total Electronic Energy
    Eel = E1 + E2 + E3

    ! Damp Charges

    DqS = maxval(abs(qS - qS_old))
    DqA = maxval(abs(qA - qA_old))

    if (DqS.ge.dump_thresh.or.DqA.ge.dump_thresh) then
      qS = qS_old + 0.4d0*(qS - qS_old)
      qA = qA_old + 0.4d0*(qA - qA_old)
    end if

    print*, qA(:)

    ! Compute Fock Matrix

    F(:,:) = H0(:,:)

    do A = 1,nAt
      shift_at_A = Gamm(A)*qA(A)*qA(A)
      do l = 1,2

        shift_sh_A = sum(ckm(A,:,l,:)*qS(:,:))

        do B = 1,nAt
          shift_at_B = Gamm(B)*qA(B)*qA(B)
          do lp = 1,2

            shift_sh_B = sum(ckm(B,:,lp,:)*qS(:,:))

            ish_A = 0
            ish_B = 0
            if (l.eq.2.and.atype(A).ne.1) ish_A = 2
            if (lp.eq.2.and.atype(B).ne.1) ish_B = 2

            do mu = int(shell(2*(A-1)+l,4)), int(shell(2*(A-1)+l,4)) + ish_A
              do nu = int(shell(2*(B-1)+lp,4)), int(shell(2*(B-1)+lp,4)) + ish_B

                F(mu,nu) = F(mu,nu) - 0.5d0*S(mu,nu)  &
                           *(shift_sh_A + shift_sh_B + shift_at_A + shift_at_B)
              end do
            end do

          end do
        end do

      end do
    end do

    ! Convergency Criterium
   
    conv = abs(Eel - Eel_old)
       
    ! New Charges
    
    qS_old = qS
    qA_old = qA
    Eel_old = Eel

    ! Write cycle results

    write(*,'(x,a1,x,i3,x,a1,x,f13.8,x,a1,x,f10.6,x,a1,x,f10.6,x,a1,x,f10.6,x,a1,x,f10.8,x,a1,x)') &
            '|',nSCF,'|',Eel,'|',E1,'|',E2,'|',E3,'|',conv,'|'

  enddo
  write(*,*) '---------------------------------------------------------------------------'

!------ End of SCF loop -------------------------------------------------------!

  ! Did it converge?

  if(nSCF.eq.maxSCF) then ! No :(

    write(*,*)
    write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*) '                 Convergence failed                 '
    write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)

    stop

  else ! Yes, write SCF results

    write(*,*)
    write(*,'(x,a20,x,i2,x,a6)') '### SCF CONVERGED in', nSCF, 'steps!'
    write(*,*)

    write(*,'(x,a17)') 'Orbital Energies:'
    write(*,'(x,a3,x,a10,2x,a10)') '#', 'Hartree','eV'
    do mu = 1,nBas
      write(*,'(x,i3,x,f10.6,2x,f10.6)') mu, e(mu), e(mu)*ev_hartree
    end do
    write(*,*)

    write(*,'(x,a20,f10.6,2x,a7)') '1st order Energy:   ', E1, 'Hartree'
    write(*,'(x,a20,f10.6,2x,a7)') '2nd order Energy:   ', E2, 'Hartree'
    write(*,'(x,a20,f10.6,2x,a7)') '3rd order Energy:   ', E3, 'Hartree'
    write(*,*)
    write(*,'(x,a20,f10.6,2x,a7)') 'Electronic Energy:  ', Eel, 'Hartree'
    write(*,*)

  end if

  ! Deallocate variables

  deallocate(C, Cp, P, F, Fp, e)
  deallocate(qS_old, qS, qA_old)

end subroutine SCF
