subroutine read_coordinates(nat,ang_bohr,atype,pos,symbol,Gamm,G)

  implicit none

  integer, intent(in) :: nat
  real(8), intent(in) :: ang_bohr
  real(8), intent(in) :: Gamm(4)

  integer  , intent(out) :: atype(nat)
  real(8)  , intent(out) :: pos(nat,3)
  real(8)  , intent(out) :: G(nat)
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
    if (symbol(i).eq.'H'.or.symbol(i).eq.'h') then 
      atype(i) = 1
      G(i) = Gamm(1)
    else if (symbol(i).eq.'C'.or.symbol(i).eq.'c') then
      atype(i) = 2
      G(i) = Gamm(2)
    else if (symbol(i).eq.'N'.or.symbol(i).eq.'n') then
      atype(i) = 3
      G(i) = Gamm(3)
    else if (symbol(i).eq.'O'.or.symbol(i).eq.'o') then
      atype(i) = 4
      G(i) = Gamm(4)
    end if
  end do

  read(5,*)
  read(5,*)

  return

end subroutine read_coordinates
