module utilities_module

  implicit none

  public
  real(8), parameter :: pi = 3.141592653589793d0

contains

!---------- READ FROM PARAMETERS.DAT ------------------------------------------!
  subroutine parameters(ang_bohr, ev_hartree, kf, alpha, Zeff)

    real(8), intent(inout) :: ang_bohr, ev_hartree, kf
    real(8), intent(inout) :: alpha(:), Zeff(:)

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

!---------- COMPUTE FACTORIAL -------------------------------------------------!
  recursive integer(16) function fact(n) result(f)

    integer, intent(in) :: n
    
    if (n.eq.0) then
      f = 1_16
    else
      f = int(n,16)*fact(n-1)
    end if

  end function fact

!------------ CONVERSIONS -----------------------------------------------------!
  subroutine dist_conversion(nc,d,a_b)

    implicit none

    integer(16), intent(in)    :: nc
    real(8)    , intent(inout) :: d(nc,3)
    real(8)    , intent(in)    :: a_b

    integer :: i

    do i = 1,nc
      d(i,3) = d(i,3)/a_b     ! Distance in a.u.
    end do

    return

  end subroutine dist_conversion

  subroutine E_conversion(E, ev_h)

    implicit none

    real(8), intent(inout) :: E
    real(8), intent(in)    :: ev_h

    E = E*ev_h                ! Energy in eV

    return

  end subroutine E_conversion

end module utilities_module
