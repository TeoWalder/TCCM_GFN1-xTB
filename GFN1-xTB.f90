!-------------------------------------------------------!
!---------- TIGHT BINDING SEMI-EMPIRICAL DFT -----------!
!-------------------------------------------------------!

program xTB

  use module_chem
  use module_utils
  implicit none

! PARAMETERS  
  integer                :: natoms
  real(8)                :: ang_bohr, ev_hartree, kf
  real(8)                :: alpha(4), Zeff(4)
  character, allocatable :: symbol(:)
! PHYSICAL VARIABLES
  integer, allocatable  :: atype(:)
  real(8)               :: Erep, Edisp, E0, E1, E2, E3, Etot
  real(8), allocatable  :: pos(:,:), dist(:,:)
! TECHNICAL VARIABLES
  integer       :: i, j, ios
  character(30) :: filename

!------------------------------------------------------------------------------!

!---------- INITIALIZATION ----------------------------------------------------!
  
  ! Reading parameters.dat
  call parameters(ang_bohr, ev_hartree, kf, alpha, Zeff)

  ! Reading input file
  open(10, file='input', status='old', iostat=ios)

    if (ios.ne.0) then
      write(*,*) '*** ERROR opening the input file'
      stop
    end if
    read(10,*)
    read(10,*) filename

  close(10)

  ! Reading xyz file
  open(11, file=filename, status='old', iostat=ios)

    if (ios.ne.0) then
      write(*,*) '*** ERROR opening', filename
      stop
    end if
    read(11,*) natoms
    read(11,*)

  allocate(pos(natoms,3), symbol(natoms), dist(natoms,natoms), atype(natoms))

    do i = 1,natoms
      read(11,*) symbol(i), pos(i,:)
    end do

  close(11)

  do i = 1,natoms
    if (symbol(i).eq.'H') then 
      atype(i) = 1
    else if (symbol(i).eq.'C') then
      atype(i) = 2
    else if (symbol(i).eq.'N') then
      atype(i) = 3
    else if (symbol(i).eq.'O') then
      atype(i) = 4
    else
      stop '*** ERROR: wrong atom type in xyz file'
    end if
  end do

  write(*,*) 'Molecule:      ', filename
  write(*,*)
  write(*,*) 'Positions (Angstrom):'
  do i = 1,natoms
    write(*,'(3f10.3)') pos(i,:)
  end do
  write(*,*)

  ! Convert Angstrom to a.u.
  pos(:,:) = pos(:,:)/ang_bohr

  write(*,*) 'Positions (bohr):'
  do i = 1,natoms
    write(*,'(3f10.3)') pos(i,:)
  end do
  write(*,*)

!---------- COMPUTE DISTANCES -------------------------------------------------!

  call distance(natoms, pos, dist)

  write(*,*) 'Distances (bohr):'
  write(*,'(2x,25i10)') (i, i=1,natoms)
  do i = 1,natoms
    write(*,'(i2,25f10.6)') i, dist(1:i,i)
  end do
  write(*,*)

!---------- ZEROTH ORDER ENERGY -----------------------------------------------!

  call repulsion_energy(natoms, atype, dist, Zeff, alpha, kf, Erep)
!  call dispersion_energy()

  E0 = Erep !+ Edisp

  write(*,*) '0th ORDER ENERGY'
  write(*,'(a25,f15.10)') 'Repulsion Energy:        ', Erep
!  write(*,'(a25,f15.10)') 'Dispersion Energy:       ', Edisp
  write(*,'(a25,f15.10)') 'Total 0th Order Energy:  ', E0
  write(*,*)

!---------- FIRST ORDER ENERGY ------------------------------------------------!

  

!---------- END PROGRAM -------------------------------------------------------!

  deallocate(pos,symbol,dist,atype)
  stop

end program xTB
