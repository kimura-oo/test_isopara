program makebc_func
      implicit none
      double precision :: eps
      integer :: node, nelm
      double precision, allocatable :: xx(:), yy(:)
      integer :: ifbcu, sfbcu
      integer, allocatable :: nfbcu(:)
      double precision, allocatable :: ffbcu(:)
      double precision :: xmax, xmin, ymax, ymin
      integer :: i, j
      character(50) :: mesfile, bdcfile
    
      open(9, file='makebc.file', status='unknown')
      read(9, '(a)') mesfile
      read(9, '(a)') bdcfile
      close(9)

      open(10, file=mesfile, status='old')
      open(50, file=bdcfile, status='replace')

      read(10, *) node, nelm
      allocate(xx(node), yy(node))

      do i = 1, node
            read(10, *) j, xx(i), yy(i)
      end do
      close(10)

      eps = 1.0d-10
      sfbcu = 0

      xmax = maxval(xx)
      xmin = minval(xx)
      ymax = maxval(yy)
      ymin = minval(yy)

      do i = 1, node
            if ((dabs(xx(i) - xmin) .le. eps) .or. (dabs(xx(i) - xmax) .le. eps) .or. (dabs(yy(i) - ymin) .le. eps) .or. (dabs(yy(i) - ymax) .le. eps)) then
                  sfbcu = sfbcu + 1
            end if
      end do
    
      allocate(nfbcu(sfbcu), ffbcu(sfbcu))

      ifbcu = 0

      do i = 1, node
            if ((dabs(xx(i) - xmin) .le. eps) .or. (dabs(xx(i) - xmax) .le. eps) .or. (dabs(yy(i) - ymin) .le. eps) .or. (dabs(yy(i) - ymax) .le. eps)) then
                  ifbcu = ifbcu + 1
                  nfbcu(ifbcu) = i
                  ffbcu(ifbcu) = f(xx(i), yy(i))
            end if
      end do

      write(50, 120) sfbcu
      write(50, 121) (i, nfbcu(i), ffbcu(i), i = 1, sfbcu)
      close(50)

120 format (i9)
121 format (2i9, e15.6)

contains
    function f(x, y)
        implicit none
        double precision, intent(in) :: x, y
        double precision :: f, pi

        pi = acos(-1.0d0)

        !f = x**3 - 3*x*y**2 + 2*x**2 + y
        f = 2*pi**2*sin(pi*x)*sin(pi*y)
    end function f

end program makebc_func