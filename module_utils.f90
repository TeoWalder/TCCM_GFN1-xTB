module module_utils

  implicit none

  real(8), parameter, public :: pi = 4.d0*atan(1.d0)

contains

 !---------- READ FROM PARAMETERS.DAT -------------------------------!
  subroutine parameters(ang_bohr, ev_hartree, kf, alpha, Zeff)

    real(8), intent(inout) :: ang_bohr, ev_hartree, kf
    real(8), intent(inout) :: alpha(4), Zeff(4)

    integer :: i 

    ! Read parameters for a.u. conversion and
    !to compute Repulsion Energy (eq. 9)

                                  ! Lines:
    do i = 1,6
      read(10,*)                  ! 1-6
    end do
    read(10,*) ang_bohr           ! 7
    read(10,*)                    ! 8
    read(10,*) ev_hartree         ! 9
    do i = 1,5
      read(10,*)                  ! 10-14
    end do
    read(10,*) kf                 ! 15
    read(10,*)                    ! 16
    read(10,*) alpha(:)           ! 17
    read(10,*)                    ! 18
    read(10,*) Zeff(:)            ! 19

    return

  end subroutine parameters

 !---------- READ FROM PARAMETERS.DAT -------------------------------!

  subroutine basis_set()

    implicit none

    ! Read wavefunction parameters (eq. 39-40)

                                         ! Lines:
    do i = 1,70
      read(10,*)                         ! 20-89
    end do

    do i = 1,8                           ! 70-113
    read(10,*) atom, shell, np           ! Basis: (a)
    if (atom) then
      allocate(zeta(np),d(np))
      read(10,*) zeta(:)                 ! Basis: (b)
      read(10,*) contraction(:)          ! Basis: (c)
    else
      read(10,*)
      read(10,*)
    end if

end module module_utils
