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
	kk(:) = 0
	allocate( kk(nx), nbc(nx) )
!
	do i = 1, mx
		i1 = nc(1, i)
		i2 = nc(2, i)
		i3 = nc(3, i)
		kk(i1) = kk(i1) + i2 - i3
		kk(i2) = kk(i2) + i3 - i1
		kk(i3) = kk(i3) + i1 - i2
	end do
!
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
	501 format( i0, x, i0, x, f0.6, x, f0.6 )
end program waku
