program makebc
      implicit none 
      double precision :: eps
      integer :: node , nelm
      double precision, allocatable :: xx(:),yy(:)
      integer :: ifbcu,sfbcu
      integer , allocatable :: nfbcu(:)
      double precision , allocatable :: ffbcu(:)
      double precision :: xmax,xmin,ymax,ymin

      integer :: i,j

      character (50) :: mesfile, bdcfile

      open(9,file='makebc.file',status='unknown')
      read(9,'(a)') mesfile
      read(9,'(a)') bdcfile

      open(10,file=mesfile, status='unknown')
      open(50,file=bdcfile, status='unknown')

      read(10,*) node, nelm
      allocate(xx(node),yy(node))
      read(10,*) (i,xx(i),yy(i),j=1,node)

      eps = 1.0d-10
      sfbcu = 0

      xmax = maxval(xx)
      xmin = minval(xx)
      ymax = maxval(yy)
      ymin = minval(yy)

!write(*,*) xmax,xmin,ymax,ymin
!stop

      do i = 1, node
            if((dabs(xx(i)-xmin).le.eps).or.(dabs(xx(i)-xmax).le.eps).or.(dabs(yy(i)-ymin).le.eps).or.(dabs(yy(i)-ymax).le.eps))then
                  sfbcu = sfbcu + 1
            end if
      end do
!write(6,*) sfbcu
! stop
      allocate( nfbcu(sfbcu), ffbcu(sfbcu) )

      ifbcu = 0

      do i = 1, node
            if((dabs(xx(i) - xmin) .le. eps) .or. (dabs(xx(i) - xmax) .le. eps).or.(dabs(yy(i) - ymin) .le. eps) .or. (dabs(yy(i) - ymax) .le. eps)) then
                  ifbcu = ifbcu + 1
                  nfbcu(ifbcu) = i
                  ffbcu(ifbcu) = 0.0d0
	      end if
      end do

      write(50,120) sfbcu
      write(50,121) (i, nfbcu(i), ffbcu(i), i = 1, sfbcu)

      120 format (3i9)
      121 format (2i9,e15.6)
      stop
end program makebc
