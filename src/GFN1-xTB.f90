program GFN1_xTB

  implicit none

! PARAMETERS  
  real(8) :: ang_bohr, ev_hartree, kf
  real(8) :: alpha(4), Zeff(4)
  real(8) :: Gamm(4)
! PHYSICAL VARIABLES
  integer                :: nAt, nel
  integer                :: nBas, nSh, nOcc
  integer  , allocatable :: atype(:)
  character, allocatable :: symbol(:)
  real(8)  , allocatable :: shell(:,:)
  real(8)  , allocatable :: pos(:,:), dist(:,:)
  real(8)                :: Erep, E0, Eel, Etot!, Edisp
  real(8)  , allocatable :: H0(:,:), S(:,:), X(:,:)
  real(8)  , allocatable :: eta(:,:)
  real(8)  , allocatable :: CKM(:,:,:,:)
  real(8)  , allocatable :: G(:)
  real(8)  , allocatable :: q(:)
! TECHNICAL VARIABLES
  integer :: i

!---------- HELLO WORLD -------------------------------------------------------!

  write(*,*) '                  ______________________________________ '
  write(*,*) '                 |                                      |' 
  write(*,*) '                 |      Online Programming Project      |'
  write(*,*) '                 |      GFN1-xTB Semiempirical DFT      |'
  write(*,*) '                 |______________________________________|'
  write(*,*)

!---------- INITIALIZATION ----------------------------------------------------!

  ! Open parameters.dat
  open(10, file='parameters.dat', status='old')
  
  call read_parameters(ang_bohr, ev_hartree, kf, alpha, Zeff, Gamm)

  do i = 1,5
    read(5,*)
  end do
 
  ! Read number of atoms
  read(5,*) nAt

  allocate(pos(nAt,3), symbol(nat), dist(nat,nat), atype(nat), q(nat), G(nat))

  read(5,*)
  read(5,*)

  write(*,'(x,a74)') '___________________________________________________________________________'
  write(*,*)
  write(*,'(x,a74)') '                          MOLECULAR STRUCTURE                              '
  write(*,'(x,a74)') '___________________________________________________________________________'
  write(*,*)
  ! Read CoordinAtes
  call read_coordinAtes(nat, ang_bohr, atype, pos, symbol, Gamm, G)

  ! Compute Distances
  call distances(nAt, pos, dist)

  ! Print Distances
  write(*,*) 'Distances (Bohr):'
  write(*,*)
  call print_matrix(nAt, dist)

!---------- BASIS SET ---------------------------------------------------------!

  ! Read number of shells and basis functions
  read(5,*) nSh, nBas
  read(5,*)

  allocate(shell(nSh,5), eta(nAt,2))
  allocate(H0(nBas,nBas), S(nBas,nBas), X(nBas,nBas))
  allocate(CKM(nAt,nat,2,2))

  ! Read number of electrons and occupied orbitals
  read(5,*) nel, nOcc

  read(5,*)
  read(5,*)
  read(5,*) 

  ! Basis Set, Electronegativity, H0 & S
  call read_basis(nAt, nSh, nBas, H0, S, eta, shell)

  ! Print Basis Set
  write(*,'(x,a74)') '___________________________________________________________________________'
  write(*,*)
  write(*,'(x,a74)') '                              BASIS SET                                    '
  write(*,'(x,a74)') '___________________________________________________________________________'
  write(*,*)
  write(*,'(x,a27,i3)') 'Total number of electrons: ', nel
  write(*,'(x,a27,i3)') 'Occupied orbitals:         ', nOcc
  write(*,'(x,a27,i3)') 'Number of basis functions: ', nBas
  write(*,*)

!---------- ZEROth ORDER ENERGY -----------------------------------------------!

  call repulsion_energy(nAt, atype, dist, Zeff, alpha, kf, Erep)
  ! call dispersion_energy()

  E0 = Erep !+ Edis

!---------- HIGHER ORDER ENERGIES ---------------------------------------------!

  ! take the inverse square root of the Overlap
  call inv_sqrt_S(nBas, S, X)

  ! compute the Coulomb-Kernel Matrix
  call CK_matrix(nAt, dist, eta, ckm)

  ! Self-Consistent Field
  call SCF(nBas, nOcc, nAt, nSh, atype, S, H0, X, shell, CKM, G, Eel, q, ev_hartree)

!--------- TOTAL ENERGY & CHARGES ---------------------------------------------!

  Etot = E0 + Eel

  write(*,'(x,a74)') '___________________________________________________________________________'
  write(*,*)
  write(*,'(x,a74)') '                               RESULTS                                     '
  write(*,'(x,a74)') '___________________________________________________________________________'
  write(*,*)
  write(*,'(x,a20,f10.6,2x,a7)') 'Repulsion  Energy:  ', Erep, 'Hartree'
  write(*,'(x,a20,a10,2x,a7)')   'Dispersion Energy:  ', '--', 'Hartree'
  write(*,'(x,a20,f10.6,2x,a7)') 'Electronic Energy:  ', Eel , 'Hartree'
  write(*,'(x,a20,f10.6,2x,a7)') 'Total      Energy:  ', Etot, 'Hartree'
  write(*,*)
  write(*,'(x,a25)') 'Mulliken Charges:        '
  write(*,*)
  do i = 1,nAt
    write(*,'(x,a1,f10.6)') symbol(i), q(i)
  end do

!---------- END PROGRAM -------------------------------------------------------!

  ! Close parameters.dat
  close(10)
  ! Deallocate variables
  deallocate(pos, symbol, dist, atype)
  deallocate(shell, eta, q, G)
  deallocate(H0, S, X, CKM)

  stop

end program GFN1_xTB
