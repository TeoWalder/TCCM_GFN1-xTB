subroutine print_matrix(n,mat)

  implicit none

  integer, intent(in) :: n
  real(8), intent(in) :: mat(n,n)

  integer :: i, j, cnt

  ! Print a matrix grouped by 7 columns
  
  cnt = 0

  do j = 1,n/7

    write(*,'(2x,7i10)') (i, i=cnt+1,cnt+7)
    do i = cnt+1,n
      if (i.le.cnt+7) then
        write(*,'(i2,7f10.6)') i, mat(cnt+1:i,i)
      elseif (i.gt.cnt+7) then
        write(*,'(i2,7f10.6)') i, mat(cnt+1:cnt+7,i)
      end if
    end do

    cnt = cnt + 7

    write(*,*)

  end do

  if (mod(n,7).ne.0) then
    write(*,'(2x,7i10)') (i, i=n-mod(n,7)+1,n)
    do i = n-mod(n,7)+1,n
      write(*,'(i2,7f10.6)') i, mat(n-mod(n,7)+1:i,i)
    end do
  end if

  write(*,*)

  return

end subroutine print_matrix
