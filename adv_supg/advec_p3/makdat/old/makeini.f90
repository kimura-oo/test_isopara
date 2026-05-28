implicit none

integer :: node , nelm
double precision, allocatable :: xx(:),yy(:)
integer :: i,j
double precision ::pi
double precision,allocatable :: r(:),u(:)
character (50) :: mesfile,bdcfile
pi = 3.141592

open(9,file='makeini.file',status='unknown')

read(9,'(a)') mesfile
read(9,'(a)') bdcfile

open(10,file=mesfile, status='unknown')
open(50,file=bdcfile, status='unknown')

read(10,*) node, nelm
allocate(xx(node),yy(node),r(node),u(node))
read(10,*) (i,xx(i),yy(i),j=1,node)

do i=1,node
	r(i) = sqrt((xx(i)-0.5)**2+(yy(i)-0.75)**2)
	if (r(i)<=0.5) then
            u(i) = 0.5 * (cos(2 * pi * r(i)) + 1)
	else
	u(i) = 0
	end if
end do
!                                     
!
!do i=1,node
!   if ((xx(i) <= 0.25 ) .and. (xx(i) >=-0.25) .and. (yy(i)<=0.75) .and. (yy(i)>=0.25)) then
!     u(i) = 1.0d0
!  else 
!    u(i) =0
!  end if
!end do
!

write(50,121) (i,u(i),i=1,node)


121 format (i9,e15.6)
stop
end