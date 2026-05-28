implicit none

integer :: node , nelm
double precision, allocatable :: xx(:),yy(:)
integer :: i,j
double precision, allocatable :: uu(:),vv(:)

character (50) :: mesfile,bdcfile

open(9,file='makeuv.file',status='unknown')

read(9,'(a)') mesfile
read(9,'(a)') bdcfile

open(10,file=mesfile, status='unknown')
open(50,file=bdcfile, status='unknown')

read(10,*) node, nelm
allocate(xx(node),yy(node))
read(10,*) (i,xx(i),yy(i),j=1,node)



allocate(uu(node),vv(node))

do i=1,node
	uu(i)=-yy(i)
	vv(i)=xx(i)
end do


write(50,121) (i,uu(i),vv(i),i=1,node)


121 format (i9,2e15.6)
stop
end

