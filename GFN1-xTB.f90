!-------------------------------------------------------!
!---------- TIGHT BINDING SEMI-EMPIRICAL DFT -----------!
!-------------------------------------------------------!

program xTB

  use chem_module
  use utilities_module
  implicit none

! PARAMETERS  
  integer                :: natoms
  real(8)                :: ang_bohr, ev_hartree, kf
  real(8)                :: alpha(4), Zeff(4)
  character, allocatable :: symbol(:)
! PHYSICAL VARIABLES
  integer    , allocatable  :: atype(:)
  integer(16)               :: ncouples
  real(8)                   :: Erep, Edisp, E0, E1, E2, E3, Etot
  real(8)    , allocatable  :: pos(:,:), dist(:,:)
! TECHNICAL STUFF
  integer       :: i, j, A, B, ios
  character(20) :: filename

!------------------------------------------------------------------------------!

!---------- INITIALIZATION ----------------------------------------------------!
  
  ! Reading parameters.dat
  call parameters(ang_bohr, ev_hartree, kf, alpha, Zeff)

  ! Reading input file
  open(10, file="input", status='old', iostat=ios)

    if (ios.ne.0) then
      write(*,*) "*** ERROR opening the input file ***"
      stop
    end if
    read(10,*)
    read(10,*) filename

  close(10)

  ! Reading xyz file
  open(11, file=filename, status='old', iostat=ios)

    if (ios.ne.0) then
      write(*,*) "*** ERROR opening ", filename, "***"
      stop
    end if
    read(11,*) natoms
    read(11,*)

  ncouples = fact(natoms)/(2_16*fact(natoms-2))
  allocate(pos(natoms,3), symbol(natoms), dist(ncouples,3), atype(natoms))

    do i = 1,natoms
      read(11,*) symbol(i), pos(i,:)
    end do

  close(11)

  do i = 1,natoms
    if (symbol(i).eq.'H') atype(i) = 1
    if (symbol(i).eq.'C') atype(i) = 2
    if (symbol(i).eq.'N') atype(i) = 3
    if (symbol(i).eq.'O') atype(i) = 4
  end do

  print*, "Molecule:      ", filename
  print*, "********** POSITIONS **************"
  do i = 1,natoms
    write(*,"(3f10.3)") pos(i,:)
  end do

!---------- COMPUTE DISTANCES -------------------------------------------------!

  call distance(natoms, ncouples, pos, dist)

  print*, "********* DISTANCES (A) ***********"
  do i = 1,ncouples
    write(*,"(2f10.0,f10.6)") dist(i,:)
  end do

  call dist_conversion(ncouples, dist, ang_bohr)

  print*, "********* DISTANCES (b) ***********"
  do i = 1,ncouples
    write(*,"(2f10.0,f10.6)") dist(i,:)
  end do  

!---------- ZEROTH ORDER ENERGY -----------------------------------------------!

  call repulsion_energy(natoms, ncouples, atype, dist, Zeff, alpha, kf, Erep)
!  call dispersion_energy()

  E0 = Erep !+ Edisp

  print*, "********* 0th ORDER ENERGY (a.u.) ***********"
  write(*,"(a25,f15.10)") "Repulsion Energy:        ", Erep
!  write(*,"(a25,f15.10)") "Dispersion Energy:       ", Edisp
  write(*,"(a25,f15.10)") "Total 0th Order Energy:  ", E0

!---------- FIRST ORDER ENERGY ------------------------------------------------!

  

!---------- END PROGRAM -------------------------------------------------------!

  deallocate(pos,symbol,dist,atype)
  stop

end program xTB
