subroutine read_coordinates(nat,ang_bohr,atype,pos,symbol)

  implicit none

  integer, intent(in) :: nat
  real(8), intent(in) :: ang_bohr

  integer  , intent(out) :: atype(nat)
  real(8)  , intent(out) :: pos(nat,3)
  character, intent(out) :: symbol(nat)

  integer :: i

  ! Atom coordinates (angstrom)
  write(*,*) 'Positions (Angstrom):'
  write(*,*)
  do i = 1,nat
    read(5,*) symbol(i), pos(i,:)
    write(*,'(x,a,3f10.3)') symbol(i), pos(i,:)
  end do
  write(*,*)

  ! Convert Angstrom to a.u.
  pos(:,:) = pos(:,:)/ang_bohr

  ! Atom coordinates (bohr)
  write(*,*) 'Positions (Bohr):'
  write(*,*)
  do i = 1,nat
    write(*,'(x,a,3f10.3)') symbol(i), pos(i,:)
  end do
  write(*,*)

  ! Associate symbols to numbers
  do i = 1,nat
    if (symbol(i).eq.'H') then 
      atype(i) = 1
    else if (symbol(i).eq.'C') then
      atype(i) = 2
    else if (symbol(i).eq.'N') then
      atype(i) = 3
    else if (symbol(i).eq.'O') then
      atype(i) = 4
    end if
  end do

  read(5,*)
  read(5,*)

  return

end subroutine read_coordinates
