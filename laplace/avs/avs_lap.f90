!************************************************
!    2次元Laplace方程式のavs可視化用プログラム  *
!    　　　　　　　　　　　　　　　　　　　     *
!         　　　　　2022. 3. 2. coded by Yasui　*
!         　　　2025. 2. 28. revised by Miyake　*
!************************************************
!
program main
	implicit none
	integer :: i, j, k, node, nelm
	integer , allocatable :: nc(:,:)
	double precision , allocatable :: xx(:), yy(:) , uu(:), vv(:), phi(:)
	character(50) :: mesfile, resfile, inpfile
!
	open(9, file = 'file_avs.dat', status = 'unknown')
!
	read(9,'(a)') mesfile
	read(9,'(a)') resfile
	read(9,'(a)') inpfile
	close(9)
!
	open(11, file = mesfile, status = 'unknown')
	open(50, file = resfile, status = 'unknown')
	open(60, file = inpfile, status = 'unknown')
!
!====== Mesh Data ======
!
	read(11,*) node, nelm
	allocate( xx(node), yy(node), nc(3,nelm) )
	read(11,*) ( i, xx(i), yy(i), j = 1, node )
	read(11,*) ( i, (nc(j,i), j = 1, 3 ), k = 1, nelm )
	close(11)
!
	allocate( uu(nelm), vv(nelm), phi(node) )
!
	read(50,*)( i, phi(i), j = 1, node )
	read(50,*)( i, uu(i), vv(i), j = 1, nelm )
	close(50)
!
!-------- ステップ数
	write(60,'(i0)') 1
!
!-------- データの繰り返しタイプ
	write(60,'(a)') 'data_geom'
!
!-----ステップ番号
	write(60,'(a)') 'step1'
!
!---- mesh data
	write(60,601) node, nelm
	write(60,602) (i, xx(i), yy(i), 0.0d0, i = 1, node)
	write(60,603) (i, (nc(j,i), j = 1, 3), i = 1, nelm)
!
!---- 節点のデータ数，要素のデータ数
	write(60,'(a)') '1  2'
!
!---- 節点データ成分数，第1成分ベクトル長
	write(60,'(a)') '1  1'
!
!---- 第1成分名
	write(60,'(a)') 'phi,'
!
!---- 節点データ
	write(60,604) (i, phi(i),  i = 1, node)
!
!---- 要素データ成分数，第1成分ベクトル長，第2成分ベクトル長
	write(60,'(a)') '2  1  1'
!
!---- 第1成分名
	write(60,'(a)') 'vel x,'
!
!---- 第2成分名
	write(60,'(a)') 'vel y,'
!
!---- 要素データ
	write(60,606) (i, uu(i), vv(i), i = 1, nelm)
!
	601 format(2i7)
	602 format(i7,3e15.6)
	603 format(i7, 5x, '1', 5x, 'tri', 3i7)
	604 format(i7,e15.6)
	605 format(3i7)
	606 format(i7, 2e15.6)
!
	end program main
