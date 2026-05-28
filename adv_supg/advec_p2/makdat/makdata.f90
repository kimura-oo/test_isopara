program makdata
      implicit none

      integer :: node, nelm
      double precision, allocatable :: xx(:), yy(:)
      character(80) :: mesfile, uvfile, inifile, bcfile
      integer :: i, j

      open(9, file='file.dat', status='old', action='read')
      read(9, '(a)') mesfile
      read(9, '(a)') uvfile
      read(9, '(a)') inifile
      read(9, '(a)') bcfile
      close(9)

      write(*,*) 'Input mesh file     : ', trim(mesfile)
      write(*,*) 'Output velocity file: ', trim(uvfile)
      write(*,*) 'Output initial file : ', trim(inifile)
      write(*,*) 'Output boundary file: ', trim(bcfile)

      open(10, file=mesfile, status='old', action='read')
      read(10, *) node, nelm
      allocate(xx(node), yy(node))
      do i = 1, node
            read(10, *) j, xx(i), yy(i) 
      end do
      close(10)

      write(*,*) 'Number of nodes       : ', node
      write(*,*) 'Number of elements    : ', nelm
      
      call make_uv(node, xx, yy, uvfile)
      call make_initial(node, xx, yy, inifile)
      call make_boundary(node, xx, yy, bcfile)
      
      deallocate(xx, yy)

            write(*,*) 'All data files have been created successfully.'
  
      stop

end program makdata

subroutine make_uv(node, xx, yy, filename)
      implicit none
      integer, intent(in) :: node
      double precision, intent(in) :: xx(node), yy(node)
      character(*), intent(in) :: filename
  
      integer :: i
      double precision, allocatable :: uu(:), vv(:)
  
      allocate(uu(node), vv(node))

      do i = 1, node
            uu(i) = 0.0d0
            vv(i) = -0.5d0
      end do

      ! do i = 1, node
      !       uu(i) = -yy(i)
      !       vv(i) = xx(i)
      ! end do

      open(50, file=filename, status='unknown', action='write')
      do i = 1, node
            write(50, '(i9, 2e15.6)') i, uu(i), vv(i)
      end do
      close(50)
  
      deallocate(uu, vv)
  
end subroutine make_uv

subroutine make_initial(node, xx, yy, filename)
      implicit none
      integer, intent(in) :: node
      double precision, intent(in) :: xx(node), yy(node)
      character(*), intent(in) :: filename

      integer :: i
      double precision :: pi
      double precision, allocatable :: r(:), u(:)
  
      pi = 4.0d0 * atan(1.0d0) 
  
      allocate(r(node), u(node))
      
      do i = 1, node
            if ((yy(i) >= 0.25d0 .and. yy(i) <= 0.75d0) .and. &
                (xx(i) >= -0.25d0 .and. xx(i) <= 0.25d0)) then
                  u(i) = 1.0d0
            else
                  u(i) = 0.0d0
            end if
      end do

      ! do i = 1, node
      !       r(i) = sqrt(xx(i)**2 + (yy(i) - 0.5d0)**2)
      !             if (r(i) <= 0.5d0) then
      !                   u(i) = 0.5d0 * (cos(2.0d0 * pi * r(i)) + 1.0d0)
      !             else
      !                   u(i) = 0.0d0
      !             end if
      ! end do
  
      open(50, file=filename, status='unknown', action='write')
      do i = 1, node
            write(50, '(i9, e15.6)') i, u(i)
      end do
      close(50)
  
      deallocate(r, u)
  
end subroutine make_initial

subroutine make_boundary(node, xx, yy, filename)
      implicit none
      integer, intent(in) :: node
      double precision, intent(in) :: xx(node), yy(node)
      character(*), intent(in) :: filename
  
      double precision :: eps
      integer :: sfbcu
      integer, allocatable :: nfbcu(:)
      double precision, allocatable :: ffbcu(:)
      double precision :: xmax, xmin, ymax, ymin
      integer :: i, ifbcu
  
      eps = 1.0d-10
  
      xmax = maxval(xx)
      xmin = minval(xx)
      ymax = maxval(yy)
      ymin = minval(yy)
  
      sfbcu = 0
      do i = 1, node
            if ((dabs(xx(i) - xmin) <= eps) .or. &
            (dabs(xx(i) - xmax) <= eps) .or. &
            (dabs(yy(i) - ymin) <= eps) .or. &
            (dabs(yy(i) - ymax) <= eps)) then
                  sfbcu = sfbcu + 1
            end if
      end do
      
      allocate(nfbcu(sfbcu), ffbcu(sfbcu))
  
      ifbcu = 0
      do i = 1, node
            if ((dabs(xx(i) - xmin) <= eps) .or. &
                  (dabs(xx(i) - xmax) <= eps) .or. &
                  (dabs(yy(i) - ymin) <= eps) .or. &
                  (dabs(yy(i) - ymax) <= eps)) then
                  ifbcu = ifbcu + 1
                  nfbcu(ifbcu) = i
                  ffbcu(ifbcu) = 0.0d0
            end if
      end do

      open(50, file=filename, status='unknown', action='write')
      write(50, '(i9)') sfbcu
      do i = 1, sfbcu
            write(50, '(2i9, e15.6)') i, nfbcu(i), ffbcu(i)
      end do
      close(50)

      deallocate(nfbcu, ffbcu)
  
end subroutine make_boundary