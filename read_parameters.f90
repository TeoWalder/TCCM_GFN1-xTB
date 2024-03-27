subroutine read_parameters(ang_bohr,ev_hartree,kf,alpha,Zeff,Gamm)

    real(8), intent(inout) :: ang_bohr, ev_hartree, kf
    real(8), intent(inout) :: alpha(4), Zeff(4), Gamm(4)

    integer :: i 

    ! Read parameters for a.u. conversion and to compute 
    !repulsion energy (eq. 9) from file 'parameters.dat'

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

    do i = 1,115
      read(10,*)                  ! 20-134
    end do

    read(10,*) Gamm(:)            ! 135

    return

end subroutine read_parameters
