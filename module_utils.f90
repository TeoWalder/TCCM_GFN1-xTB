module module_utils

  implicit none

  public
  real(8), parameter :: pi = 4.d0*atan(1.d0)

contains

!---------- READ FROM PARAMETERS.DAT ------------------------------------------!
  subroutine parameters(ang_bohr, ev_hartree, kf, alpha, Zeff)

    real(8), intent(inout) :: ang_bohr, ev_hartree, kf
    real(8), intent(inout) :: alpha(4), Zeff(4)

    integer :: i

    open(10, file="parameters.dat")

                              ! Records
    do i = 1,6
      read(10,*)              ! 1-6
    end do
    read(10,*) ang_bohr       ! 7
    read(10,*)                ! 8
    read(10,*) ev_hartree     ! 9
    do i = 1,5
      read(10,*)              ! 10-14
    end do
    read(10,*) kf             ! 15
    read(10,*)                ! 16
    read(10,*) alpha(:)       ! 17
    read(10,*)                ! 18
    read(10,*) Zeff(:)        ! 19

    return

  end subroutine parameters

end module module_utils
