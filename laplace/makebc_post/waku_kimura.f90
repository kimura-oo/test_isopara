!
!
program waku
	implicit none
!
	integer :: i, j, ibc, n, m, nx, mx, i1, i2, i3
	integer, allocatable :: nc(:,:), kk(:), nbc(:)
	double precision, allocatable :: xx(:), yy(:), zz(:)
	double precision :: xxmin, yymin, xxmax, yymax
	integer :: min_x_index, min_y_index, max_x_index, max_y_index
      integer, allocatable :: max_x_indices(:), min_y_indices(:), other_indices(:)
      integer :: max_x_count, min_y_count, other_count
      logical :: is_max_x, is_min_y
!
	character(128) :: mesfile, resfile
!
	open(9, file = 'file_bc_waku.dat', status = 'unknown')
	read(9, '(a)') mesfile
	read(9, '(a)') resfile
	close(9)
!
	open(20, file = mesfile, status = 'old',     action = 'read' )
	open(50, file = resfile, status = 'replace', action = 'write')
!
	read(20,*) nx, mx
      write(*,*) "nx =", nx, "mx =", mx
	allocate( xx(nx), yy(nx), nc(3,mx) )
	read(20,*) ( n, xx(n), yy(n), i = 1, nx )
	read(20,*) ( m, nc(1,m), nc(2,m), nc(3,m), j = 1, mx )
	close(20)
!
	allocate( kk(nx), nbc(nx) )
!
      kk(i) = 0
	do i = 1, mx
		i1 = nc(1, i)
		i2 = nc(2, i)
		i3 = nc(3, i)
		kk(i1) = kk(i1) + i2 - i3
		kk(i2) = kk(i2) + i3 - i1
		kk(i3) = kk(i3) + i1 - i2
            
            if (i1 == 13 ) then 
                  write(*,*) i2, i3
            else if (i2 == 13) then
                  write(*,*) i3, i1
            else if (i3 == 13) then
                  write(*,*) i1, i2
            end if

	end do
!     
      write(*,*) kk(13)


	ibc = 0
!
	do i = 1, nx
		if( kk(i) .ne. 0 ) then
			ibc = ibc + 1
			nbc(ibc) = i
		end if
	end do
!
      write(50,501) ( i, nbc(i), xx( nbc(i) ), yy( nbc(i) ), i = 1, ibc )

      xxmax = -1.0d30
      max_x_count = 0
      allocate( max_x_indices(ibc), min_y_indices(ibc) )
!     
      do i = 1, ibc
            if (xx(nbc(i)) > xxmax) then
                  xxmax = xx(nbc(i))
                  max_x_count = 1
                  max_x_indices(max_x_count) = nbc(i)
            else if (xx(nbc(i)) == xxmax) then
                  max_x_count = max_x_count + 1
                  max_x_indices(max_x_count) = nbc(i)
            end if
      end do
      
!
      yymin = 1.0d30
      min_y_count = 0
!
      do i = 1, ibc
            if (yy(nbc(i)) < yymin) then
                  yymin = yy(nbc(i))
                  min_y_count = 1
                  min_y_indices(min_y_count) = nbc(i)
            else if (yy(nbc(i)) == yymin) then
                  min_y_count = min_y_count + 1
                  min_y_indices(min_y_count) = nbc(i)
            end if
      end do
!      
!
      other_count = 0
      allocate(other_indices(ibc))
!
      do i = 1, ibc
            is_max_x = .false.
            is_min_y = .false.
! max_x_indices に含まれているか確認
      do j = 1, max_x_count
            if (nbc(i) == max_x_indices(j)) then
                  is_max_x = .true.
                  exit
            end if
      end do
      ! min_y_indices に含まれているか確認
      do j = 1, min_y_count
            if (nbc(i) == min_y_indices(j)) then
                  is_min_y = .true.
                  exit
            end if
      end do
    ! どちらにも該当しない場合
            if (.not. is_max_x .and. .not. is_min_y) then
                  other_count = other_count + 1
                  other_indices(other_count) = nbc(i)
            end if
      end do
!
            write(50,*) max_x_count
            do i = 1, max_x_count
                  write(50, *) i, max_x_indices(i), xx(max_x_indices(i)), yy(max_x_indices(i))
            end do
                  write(50,*) min_y_count
            do i = 1, min_y_count
                  write(50, *) i, min_y_indices(i), xx(min_y_indices(i)), yy(min_y_indices(i))
            end do
            write(50,*) other_count
            do i = 1, other_count
                  write(50, *) i, other_indices(i), xx(other_indices(i)), yy(other_indices(i))
            end do

	501 format( i5, 1x, i5, 1x, f10.6, 1x, f10.6 )
end program waku
