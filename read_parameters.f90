subroutine read_parameters(ang_bohr, ev_hartree, kf, alpha, Zeff)

    real(8), intent(inout) :: ang_bohr, ev_hartree, kf
    real(8), intent(inout) :: alpha(4), Zeff(4)

    integer :: i 

    ! Read parameters for a.u. conversion and to compute 
    !repulsion energy (eq. 9) from file 'parameters.dat'

    open(10, file='parameters.dat', status='old')

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

    close(10)

    return

end subroutine read_parameters
