program GFN1_xTB

  implicit none

! PARAMETERS  
  real(8) :: ang_bohr, ev_hartree, kf
  real(8) :: alpha(4), Zeff(4)
! PHYSICAL VARIABLES
  integer                :: nat, nel
  integer                :: nbasis, nshell, nocc
  integer  , allocatable :: atype(:)
  integer  , allocatable :: shell(:,:)
  character, allocatable :: symbol(:)
  real(8)                :: Erep, E0, E1, E2, E3, Etot
  real(8)  , allocatable :: H0(:,:), S(:,:)
  real(8)  , allocatable :: pos(:,:), dist(:,:)
  real(8)  , allocatable :: eta(:)
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
  
  ! Read parameters.dat
  call read_parameters(ang_bohr, ev_hartree, kf, alpha, Zeff)

  do i = 1,5
    read(5,*)
  end do
 
  ! Read number of atoms
  read(5,*) nat

  allocate(pos(nat,3),symbol(nat),dist(nat,nat),atype(nat))

  read(5,*)
  read(5,*)

  ! Read Coordinates
  call read_coordinates(nat,ang_bohr,atype,pos,symbol)

  ! Read number of shells and basis functions
  read(5,*) nshell, nbasis
  read(5,*)

  allocate(H0(nbasis,nbasis),S(nbasis,nbasis),shell(nshell,4),eta(nshell))

  ! Read number of electrons and occupied orbitals
  read(5,*) nel, nocc

  read(5,*)
  read(5,*)
  read(5,*) 

  ! Basis Set, Electronegativity, H0 & S
  call read_basis(nshell,nbasis,H0,S,eta,shell)

  ! Compute Distances
  call distances(nat,pos,dist)

  ! Print Distances
  write(*,*) 'Distances (Bohr):'
  write(*,*)
  call print_matrix(nat,dist)


!---------- ZEROTH ORDER ENERGY -----------------------------------------------!

  call repulsion_energy(nat,atype,dist,Zeff,alpha,kf,Erep)

  E0 = Erep

  write(*,'(a74)') '___________________________________________________________________________'
  write(*,*)
  write(*,'(x,a25)')              '0th ORDER ENERGY         '
  write(*,*)
  write(*,'(x,a25,f15.10,2x,a7)') 'Repulsion Energy:        ', Erep, 'Hartree'
  write(*,'(x,a25,f15.10,2x,a7)') 'Total 0th Order Energy:  ', E0  , 'Hartree'

!---------- FIRST ORDER ENERGY ------------------------------------------------!

  write(*,'(a74)') '___________________________________________________________________________'
  write(*,*)
 

!---------- END PROGRAM -------------------------------------------------------!

  deallocate(pos,symbol,dist,atype)
  deallocate(shell,eta)
  deallocate(H0,S)
  stop

end program GFN1_xTB
