!  ==================================================================================
!
!              2-Dimensional Finite Element Analysis(Numerical Integration)
!                 of Steady Ideal Flows Governed by Laplace Eq.
!                 (Modified for Edge-based Boundary Conditions)
!
!  ==================================================================================
module indata
	implicit none
	! --- mesh
	integer :: node, nelm
	integer, allocatable :: nc(:,:)
	double precision, allocatable :: xx(:,:)
	! Potential (Dirichlet)
	integer :: n_potential_edges
	integer, allocatable :: potential_nodes(:,:) ! [2, n_potential_edges]
	double precision, allocatable :: potential_values(:)
	! Velocity (Neumann)
	integer :: n_velocity_edges
	integer, allocatable :: velocity_nodes(:,:) ! [2, n_velocity_edges]
	double precision, allocatable :: velocity_values(:)

	! --- Laplace
	double precision, allocatable :: phi(:)      ! Potential values on node (Solution)
	double precision, allocatable :: rhs(:)      ! R.H.S. Vector
	double precision, allocatable :: ue(:), ve(:) ! Velocity ( on Element )
	double precision, allocatable :: ea(:,:)     ! L.H.S. Matrix (Stiffness Matrix)
	integer :: mgp                               ! Number of integration points
	double precision :: area_total
	double precision, allocatable :: shpx(:,:,:), shpy(:,:,:) ! Shape Function derivatives
	double precision, allocatable :: jac(:,:)                 ! Jacobian
	double precision, allocatable :: xi_hat(:), eta_hat(:), sg(:) ! Gauss points and weights
end module indata
!
!------------------------				 MAIN PROGRAM 			------------------------------
!
program main
	use indata
	implicit none
	integer :: i
	double precision :: rr
	double precision :: begin_time, end_time

	write(*,*) 'main program start!'
	call CPU_TIME(begin_time)
	
      call datain
	call makeSSF      ! Make Shape Function for Numerical Integration
	call makeLHS      ! Make L.H.S. Matrix
	call makeRHS      ! Make R.H.S. Vector from Neumann conditions
	call boundc       ! Impose Dirichlet Boundary Condition
	call sweep        ! Solve the linear system
	call calvel( node, nelm, nc, phi, ue, ve, shpx, shpy, mgp )
	call output( node, nelm, phi, ue, ve )


	call CPU_TIME(end_time)
	write(*,'(a,f15.6,a)') 'Time of operation was ' , end_time - begin_time, ' seconds'
	write(*,*) 'main program end!'

end program main
!
!----------------------       MAIN PROGRAM END       ----------------------------
!
!-- SUBROUTINES --
! -------------------------------------------------------------------------------
subroutine datain
! -------------------------------------------------------------------------------
	use indata
	implicit none
	integer :: i, j, k, n, m
	character(100) :: mesfile, bdcfile, outfile, out2file, integfile
	character(20) :: bc_type_str ! àÍéûìIÇ»ï∂éöóÒïœêî
	integer :: n_vel_in, n_vel_wall, temp_n1, temp_n2

	mgp = 3 ! Number of integration points

	! FILE I/O
	open(9, file = 'file_FEM.dat', status = 'old' )
	read(9,'(a)') mesfile
	read(9,'(a)') bdcfile
	read(9,'(a)') outfile
	read(9,'(a)') out2file
	read(9,'(a)') integfile
	close(9)

	open(10,file = mesfile,   status = 'old',action='read')
	open(11,file = bdcfile,   status = 'old',action='read')
	open(60,file = outfile,   status = 'replace')

	! READ MESH
	read(10,*) node, nelm
	allocate( xx(2,node), nc(3,nelm) )
	read(10,*) ( n, xx(1,n), xx(2,n), i = 1, node )
	read(10,*) ( m, ( nc(i,m), i = 1, 3 ), j = 1, nelm)
	close(10)

	! READ NEW BOUNDARY CONDITION FILE (mesh.bc)
	! Potential
	read(11,*) n_potential_edges, bc_type_str
	allocate(potential_nodes(2, n_potential_edges), potential_values(n_potential_edges))
	do i = 1, n_potential_edges
		read(11,*) j, potential_nodes(1,i), potential_nodes(2,i), potential_values(i)
	end do

	! Velocity (Inlet and Wall)
	read(11,*) ! Skip empty line
	read(11,*) n_vel_in, bc_type_str
	read(11,*) ! Skip empty line
	read(11,*) n_vel_wall, bc_type_str
	
	n_velocity_edges = n_vel_in + n_vel_wall
	allocate(velocity_nodes(2, n_velocity_edges), velocity_values(n_velocity_edges))

	! Read inlet velocity
	rewind(11)
	read(11,*) j, bc_type_str; do i = 1, j; read(11,*); end do; read(11,*)
	read(11,*) n_vel_in, bc_type_str
	do i = 1, n_vel_in
		read(11,*) j, velocity_nodes(1,i), velocity_nodes(2,i), velocity_values(i)
	end do

	! Read wall velocity
	read(11,*) ! Skip empty line
	read(11,*) n_vel_wall, bc_type_str
	do i = 1, n_vel_wall
		read(11,*) j, velocity_nodes(1, i + n_vel_in), velocity_nodes(2, i + n_vel_in), velocity_values(i + n_vel_in)
	end do

	close(11)

	! Dynamic Memory Allocation
	allocate( phi(node), rhs(node), ea(node,node) )
	allocate( xi_hat(mgp), eta_hat(mgp), sg(mgp) )
	allocate( jac(mgp,nelm) )
	allocate( shpx(3,mgp,nelm), shpy(3,mgp,nelm) )
	allocate( ue(nelm), ve(nelm) )

end subroutine datain
!
! -------------------------------------------------------------------------------
subroutine makeRHS   ! Make R.H.S. Vector from Neumann conditions
! -------------------------------------------------------------------------------
	use indata
	implicit none
	integer :: i, n1, n2
	double precision :: x1, y1, x2, y2, edge_length, value

	! Initialize RHS vector
	rhs(:) = 0.0d0

	! Loop over all velocity (Neumann) boundary edges
	do i = 1, n_velocity_edges
		n1 = velocity_nodes(1, i)
		n2 = velocity_nodes(2, i)
		value = velocity_values(i)

		! Get node coordinates
		x1 = xx(1, n1)
		y1 = xx(2, n1)
		x2 = xx(1, n2)
		y2 = xx(2, n2)

		! Calculate the length of the edge
		edge_length = sqrt((x2 - x1)**2 + (y2 - y1)**2)

		! Add contribution to the global RHS vector
		! The contribution is (value * Length / 2) for each node of the edge
		rhs(n1) = rhs(n1) + value * edge_length / 2.0d0
		rhs(n2) = rhs(n2) + value * edge_length / 2.0d0
	end do

end subroutine makeRHS
!
! -------------------------------------------------------------------------------
subroutine boundc      ! Impose Dirichlet Boundary Condition
! -------------------------------------------------------------------------------
	use indata
	implicit none
	integer :: i, j, k, n
	double precision :: f
	logical, allocatable :: is_dirichlet_node(:)

	allocate(is_dirichlet_node(node))
	is_dirichlet_node = .false.

	! First, mark all nodes that have a Dirichlet condition
	do i = 1, n_potential_edges
		is_dirichlet_node(potential_nodes(1, i)) = .true.
		is_dirichlet_node(potential_nodes(2, i)) = .true.
	end do

	! Apply Dirichlet conditions using a large number (penalty method)
	do i = 1, n_potential_edges
		! Apply to the first node of the edge
		n = potential_nodes(1, i)
		f = potential_values(i)
		if (is_dirichlet_node(n)) then
			do j = 1, node
				rhs(j) = rhs(j) - ea(j, n) * f
			end do
			ea(n, :) = 0.0d0
			ea(:, n) = 0.0d0
			ea(n, n) = 1.0d0
			rhs(n) = f
			is_dirichlet_node(n) = .false. ! Mark as processed
		end if

		! Apply to the second node of the edge
		n = potential_nodes(2, i)
		f = potential_values(i)
		if (is_dirichlet_node(n)) then
			do j = 1, node
				rhs(j) = rhs(j) - ea(j, n) * f
			end do
			ea(n, :) = 0.0d0
			ea(:, n) = 0.0d0
			ea(n, n) = 1.0d0
			rhs(n) = f
			is_dirichlet_node(n) = .false. ! Mark as processed
		end if
	end do
	
	deallocate(is_dirichlet_node)

end subroutine boundc
!
! ------------------------------------------------------------------------------
subroutine sweep  ! Solve Linear-Systems (Gauss Elimination)
! ------------------------------------------------------------------------------
	use indata
	implicit none
	integer :: i, j, k, i1, l
	double precision :: ai, cc

	! Initialize solution vector
	phi = rhs

	! Forward elimination
	do i = 1, node
		if( dabs(ea(i,i)) <= 1.0d-12 ) then
			write(6,*) 'Diagonal Value',i,'is ZERO or too small!'
			stop
		end if

		ai = 1.0d0 / ea(i,i)
		phi(i) = phi(i) * ai
		ea(i,i+1:node) = ea(i,i+1:node) * ai ! Vectorization
		!do j = i + 1, node
		!	ea(i,j) = ea(i,j) * ai
		!end do

		if( i < node ) then
			i1 = i + 1
			do k = i1, node
				cc = ea(k,i)
				ea(k,i1:node) = ea(k,i1:node) - cc * ea(i,i1:node) ! Vectorization
				!do j = i1, node
				!	ea(k,j) = ea(k,j) - cc * ea(i,j)
				!end do
				phi(k) = phi(k) - cc * phi(i)
			end do
		end if
	end do

	! Backward substitution
	do i = node - 1, 1, -1
		do k = i + 1, node
			phi(i) = phi(i) - ea(i,k) * phi(k)
		end do
	end do

end subroutine sweep
!
! (The rest of the subroutines: makeSSF, set_integration_points, makeLHS, calvel, output
!  do not need significant changes and are omitted here for brevity.
!  Please use the versions from your original file.)
!
! ===============================================================================
!  NOTE: The following subroutines are assumed to be present and correct from
!        your original file. They are included here for completeness.
! ===============================================================================
!
! -------------------------------------------------------------------------------
	subroutine makeSSF
! -------------------------------------------------------------------------------
	use indata
	implicit none
	integer :: m, ig, i, j
	integer :: n1, n2, n3
	double precision :: DetJ, tmp
	double precision :: xx_int, yy_int
	double precision :: pi
	double precision, allocatable :: aa(:)
	double precision , allocatable :: dxdxi(:,:), dxidx(:,:), sq(:,:), sqd(:,:,:)
	
	allocate( dxidx(2,2), dxdxi(2,2), sq(3,mgp), sqd(2,3,mgp) )
	allocate( aa(nelm) )
	call set_integration_points(mgp, xi_hat, eta_hat, sg)
	aa(:) = 0.0d0
	do m = 1, nelm
		n1 = nc(1,m); n2 = nc(2,m); n3 = nc(3,m)
		do ig = 1, mgp

                  sq(1,ig) = xi_hat(ig); sq(2,ig) = eta_hat(ig); sq(3,ig) = 1.0d0 - xi_hat(ig) - eta_hat(ig) !shape function for 3-nodes

                  !sq(1,ig) = (1-xi_hat(ig)-eta_hat(ig))*(1-2*xi_hat(ig)-2*eta_hat(ig))
                  !sq(2,ig) = (xi_hat(ig))*(2*xi_hat(ig)-1)
                  !sq(3,ig) = (eta_hat(ig))*(2*eta_hat(ig)-1)
                  !sq(4,ig) = 4*xi_hat(ig)*(1-xi_hat(ig)-eta_hat(ig))
                  !sq(5,ig) = 4*xi_hat(ig)*eta_hat(ig)
                  !sq(6,ig) = 4*eta_hat(ig)*(1-xi_hat(ig)-eta_hat(ig))  !shape function for 6-nodes

                  sqd(1,1,ig) = 1.0d0; sqd(1,2,ig) = 0.0d0; sqd(1,3,ig) = -1.0d0 !shape function derivative by xi for 3-nodes
			sqd(2,1,ig) = 0.0d0; sqd(2,2,ig) = 1.0d0; sqd(2,3,ig) = -1.0d0 !shape function derivative by eta for 3-nodes

                  !sqd(1,1,ig) = 4*(xi_hat(ig)+eta_hat(ig))-3; sqd(1,2,ig) = 4*(xi_hat(ig))-1; sqd(1,3,ig) = 4*(eta_hat(ig))-1
                  !sqd(1,4,ig) = 4*(1-2*xi_hat(ig)-eta_hat(ig)); sqd(1,5,ig) = 4*eta_hat(ig); sqd(1,6,ig) = -4*eta_hat(ig) 
                  !shape function derivative by xi for 6-nodes
			
                  !sqd(2,1,ig) = 4*(xi_hat(ig)+eta_hat(ig))-3; sqd(2,2,ig) = 0; sqd(2,3,ig) = 4*(eta_hat(ig))-1
                  !sqd(2,4,ig) = -4*(xi_hat(ig)); sqd(2,5,ig) = 4*xi_hat(ig); sqd(2,6,ig) = 4*(1-xi_hat(ig)-2*eta_hat(ig)) 
                  !shape function derivative by eta for 6-nodes
			
                  dxdxi(1,1) = sqd(1,1,ig)*xx(1,n1) + sqd(1,2,ig)*xx(1,n2) + sqd(1,3,ig)*xx(1,n3) 
                               !sqd(1,4,ig)*xx(1,n4) + sqd(1,5,ig)*xx(1,n5) + sqd(1,6,ig)*xx(1,n6)
			dxdxi(1,2) = sqd(2,1,ig)*xx(1,n1) + sqd(2,2,ig)*xx(1,n2) + sqd(2,3,ig)*xx(1,n3) 
                               !sqd(2,4,ig)*xx(1,n4) + sqd(2,5,ig)*xx(1,n5) + sqd(2,6,ig)*xx(1,n6)
			dxdxi(2,1) = sqd(1,1,ig)*xx(2,n1) + sqd(1,2,ig)*xx(2,n2) + sqd(1,3,ig)*xx(2,n3) 
                               !sqd(1,4,ig)*xx(2,n4) + sqd(1,5,ig)*xx(2,n5) + sqd(1,6,ig)*xx(2,n6)
			dxdxi(2,2) = sqd(2,1,ig)*xx(2,n1) + sqd(2,2,ig)*xx(2,n2) + sqd(2,3,ig)*xx(2,n3) 
                               !sqd(2,4,ig)*xx(2,n4) + sqd(2,5,ig)*xx(2,n5) + sqd(2,6,ig)*xx(2,n6)

                  DetJ = dxdxi(1,1) * dxdxi(2,2) - dxdxi(1,2) * dxdxi(2,1)

			if ( DetJ <= 0.0d0 ) then
				write(*,*) 'error in J' , m ,ig, DetJ; stop
			else
				tmp = 1.0d0 / DetJ
                        jac(ig,m) = DetJ
			end if
			
                  aa(m) = aa(m) + DetJ * sg(ig)
			dxidx(1,1) = dxdxi(2,2) * tmp; dxidx(1,2) = -dxdxi(1,2) * tmp
			dxidx(2,1) = -dxdxi(2,1) * tmp; dxidx(2,2) = dxdxi(1,1) * tmp
			
                  do i = 1, 3
				shpx(i,ig,m) = sqd(1,i,ig) * dxidx(1,1) + sqd(2,i,ig) * dxidx(2,1)
				shpy(i,ig,m) = sqd(1,i,ig) * dxidx(1,2) + sqd(2,i,ig) * dxidx(2,2)
			end do
		end do
	end do
	area_total = sum(aa)
	deallocate( dxdxi, dxidx, sq, sqd, aa )
	end subroutine makeSSF
!
! -------------------------------------------------------------------------------
	subroutine set_integration_points(mgp, xi_hat, eta_hat, sg)
! -------------------------------------------------------------------------------
	implicit none
	integer, intent(in) :: mgp
	double precision, intent(out) :: xi_hat(mgp), eta_hat(mgp), sg(mgp)
	if (mgp==1) then
		xi_hat(1)=1.d0/3.d0; eta_hat(1)=1.d0/3.d0; sg(1)=0.5d0
	else if (mgp==3) then
		xi_hat(1)=1.d0/6.d0; eta_hat(1)=1.d0/6.d0; sg(1)=1.d0/6.d0
		xi_hat(2)=2.d0/3.d0; eta_hat(2)=1.d0/6.d0; sg(2)=1.d0/6.d0
		xi_hat(3)=1.d0/6.d0; eta_hat(3)=2.d0/3.d0; sg(3)=1.d0/6.d0
	end if
	end subroutine set_integration_points
!
! -------------------------------------------------------------------------------
	subroutine makeLHS
! -------------------------------------------------------------------------------
	use indata
	implicit none
	integer :: i, j, m, ig, ni, nj
	ea(:,:) = 0.0d0
	do m = 1, nelm
		do j = 1, 3
			do i = 1, 3
				ni = nc(i,m)
                        nj = nc(j,m)
				do ig = 1, mgp
					ea(ni,nj) = ea(ni,nj) + (shpx(i,ig,m) * shpx(j,ig,m) &
					                      +  shpy(i,ig,m) * shpy(j,ig,m)) * jac(ig,m) * sg(ig)
				end do
			end do
		end do
	end do
	end subroutine makeLHS
!
! ------------------------------------------------------------------------------
	subroutine calvel( node, nelm, nc, phi, ue, ve, dx, dy, mgp )
! ------------------------------------------------------------------------------
	implicit none
	integer, intent(in) :: node, nelm, nc(3,nelm), mgp
	double precision, intent(in) :: phi(node), dx(3,mgp,nelm), dy(3,mgp,nelm)
	double precision, intent(out) :: ue(nelm), ve(nelm)
	integer :: i, m, ig
	ue(:) = 0.0d0; ve(:) = 0.0d0
	do m = 1, nelm
		do i = 1, 3
			ue(m) = ue(m) - (dx(i,1,m)+dx(i,2,m)+dx(i,3,m))/3.0d0 * phi(nc(i,m))
			ve(m) = ve(m) - (dy(i,1,m)+dy(i,2,m)+dy(i,3,m))/3.0d0 * phi(nc(i,m))
		end do
	end do
	end subroutine calvel
!
! -------------------------------------------------------------------------------
	subroutine output( node, nelm, phi, ue, ve )
! -------------------------------------------------------------------------------
	implicit none
	integer, intent(in) :: node, nelm
	double precision, intent(in) :: phi(node), ue(nelm), ve(nelm)
	integer :: n
	write(60,600) ( n, phi(n), n =1, node )
	write(60,601) ( n, ue(n), ve(n), n =1, nelm )
	close(60)
	600 format(i7, d15.6)
	601 format(i7,2d15.6)
	end subroutine output
