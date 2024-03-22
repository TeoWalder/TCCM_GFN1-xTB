subroutine read_basis(nshell,nbasis,H0,S,eta,shell)

  implicit none

  integer, intent(in) :: nshell, nbasis

  real(8), intent(out) :: H0(nbasis,nbasis)
  real(8), intent(out) :: S(nbasis,nbasis)
  real(8), intent(out) :: eta(nshell)
  integer, intent(out) :: shell(nshell,4)

  integer      :: i, j, cnt, dummy
  character(1) :: cdummy

  ! Shells and electronegativity
  do i = 1,nshell
    read(5,*) dummy, shell(i,1), cdummy, shell(i,2:4), eta(i)
  end do

  read(5,*)
  read(5,*)

  ! S overlap matrix elements
  do i = 1,nbasis
    read(5,*)
    cnt = 0
    do j = 1,nbasis/7
      read(5,*) S(i,cnt+1:cnt+7)
      cnt = cnt + 7
    end do
    if (mod(nbasis,7).ne.0) read(5,*) S(i,nbasis-mod(nbasis,7)+1:nbasis) 
  end do

  read(5,*)
  read(5,*)

  ! H0 hamiltonian elements
  do i = 1,nbasis
    read(5,*)
    cnt = 0
    do j = 1,nbasis/7
      read(5,*) H0(i,cnt+1:cnt+7)
      cnt = cnt + 7
    end do
    if (mod(nbasis,7).ne.0) read(5,*) H0(i,nbasis-mod(nbasis,7)+1:nbasis) 
  end do



end subroutine
