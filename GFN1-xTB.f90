program GFN1_xTB

  implicit none

! PARAMETERS  
  real(8) :: ang_bohr, ev_hartree, kf
  real(8) :: alpha(4), Zeff(4)
  real(8) :: Gamm(4)
! PHYSICAL VARIABLES
  integer                :: nat, nel
  integer                :: nbasis, nshell, nocc
  integer  , allocatable :: atype(:)
  integer  , allocatable :: shell(:,:)
  character, allocatable :: symbol(:)
  real(8)  , allocatable :: pos(:,:), dist(:,:)
  real(8)                :: Erep, E0, E1, E2, E3, Etot
  real(8)  , allocatable :: H0(:,:), S(:,:), X(:,:)
  real(8)  , allocatable :: eta(:), ev(:)
  real(8)  , allocatable :: CKM(:,:,:,:)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real(8)  , allocatable :: qorb(:,:), qat(:)
! TECHNICAL VARIABLES
  integer :: i

!---------- HELLO WORLD -------------------------------------------------------!

  write(*,*) ' ______________________________________ '
  write(*,*) '~                                      ~' 
  write(*,*) '~      Online Programming Project      ~'
  write(*,*) '~      GFN1-xTB Semiempirical DFT      ~'
  write(*,*) '~______________________________________~'
  write(*,*)

!---------- INITIALIZATION ----------------------------------------------------!

  ! Open parameters.dat
  open(10, file='parameters.dat', status='old')
  
  call read_parameters(ang_bohr, ev_hartree, kf, alpha, Zeff, Gamm)

  do i = 1,5
    read(5,*)
  end do
 
  ! Read number of atoms
  read(5,*) nat

  allocate(pos(nat,3), symbol(nat), dist(nat,nat), atype(nat))

  read(5,*)
  read(5,*)

  write(*,'(a74)') '___________________________________________________________________________'
  write(*,*)
  write(*,'(a74)') '                          MOLECULAR STRUCTURE                              '
  write(*,'(a74)') '___________________________________________________________________________'
  write(*,*)
  ! Read Coordinates
  call read_coordinates(nat, ang_bohr, atype, pos, symbol)

  ! Compute Distances
  call distances(nat, pos, dist)

  ! Print Distances
  write(*,*) 'Distances (Bohr):'
  write(*,*)
  call print_matrix(nat, dist)

!---------- BASIS SET ---------------------------------------------------------!

  ! Read number of shells and basis functions
  read(5,*) nshell, nbasis
  read(5,*)

  allocate(shell(nshell,4), eta(nshell), ev(nbasis))
  allocate(H0(nbasis,nbasis), S(nbasis,nbasis), X(nbasis,nbasis))
  allocate(CKM(nat,nat,2,2), qat(), qorb(nat,n))

  ! Read number of electrons and occupied orbitals
  read(5,*) nel, nocc

  read(5,*)
  read(5,*)
  read(5,*) 

  ! Basis Set, Electronegativity, H0 & S
  call read_basis(nshell, nbasis, H0, S, eta, shell)

  ! Print Basis Set
  write(*,'(a74)') '___________________________________________________________________________'
  write(*,*)
  write(*,'(a74)') '                              BASIS SET                                    '
  write(*,'(a74)') '___________________________________________________________________________'
  write(*,*)
  write(*,'(x,a27,i3)') 'Total number of electrons: ', nel
  write(*,'(x,a27,i3)') 'Occupied orbitals:         ', nocc
  write(*,'(x,a27,i3)') 'Number of basis functions: ', nbasis
  write(*,*)

!---------- ZEROth ORDER ENERGY -----------------------------------------------!

  call repulsion_energy(nat, atype, dist, Zeff, alpha, kf, Erep)
  ! call dispersion_energy()

  E0 = Erep !+ Edis

  write(*,'(a74)') '___________________________________________________________________________'
  write(*,*)
  write(*,'(a74)') '                          0th ORDER ENERGY                                 '
  write(*,'(a74)') '___________________________________________________________________________'
  write(*,*)
  write(*,'(x,a25,f15.10,2x,a7)') 'Repulsion Energy:        ', Erep, 'Hartree'
  write(*,'(x,a25,a15,2x,a7)')    'Dispersion Energy:       ', '--', 'Hartree'
  write(*,*)
  write(*,'(x,a25,f15.10,2x,a7)') 'Total 0th Order Energy:  ', E0  , 'Hartree'

!---------- HIGHER ORDER ENERGIES ---------------------------------------------!

  write(*,'(a74)') '___________________________________________________________________________'
  write(*,*)

  ! take the inverse square root of the Overlap
  call inv_sqrt_S(nbasis, S, X)

  ! compute the Coulomb-Kernel Matrix
!  call CK_matrix()

  ! Self-Consistent Field
!  call SCF()

!---------- END PROGRAM -------------------------------------------------------!

  ! Close parameters.dat
  close(10)
  ! Deallocate variables
  deallocate(pos, symbol, dist, atype)
  deallocate(shell, eta)
  deallocate(H0, S, X, ev, CKM)

  stop

end program GFN1_xTB
