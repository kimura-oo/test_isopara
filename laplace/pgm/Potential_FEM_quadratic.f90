! ==================================================================================
!
!              2-Dimensional Finite Element Analysis(Numerical Integration)
!                 of Steady Ideal Flows Governed by Laplace Eq.
!                 for Quadratic Triangular Elements (6-node Tria)
!
! ==================================================================================
module indata
	implicit none
	! --- mesh
	integer :: node, nelm
	integer, allocatable :: nc(:,:) ! MODIFIED: (6, nelm)
	double precision, allocatable :: xx(:,:)

	! --- boundary condition (Edge-based for Quadratic Elements)
	! Potential (Dirichlet)
	integer :: n_potential_edges
	integer, allocatable :: potential_nodes(:,:) ! MODIFIED: [3, n_potential_edges]
	double precision, allocatable :: potential_values(:)
	! Velocity (Neumann)
	integer :: n_velocity_edges
	integer, allocatable :: velocity_nodes(:,:) ! MODIFIED: [3, n_velocity_edges]
	double precision, allocatable :: velocity_values(:)

	! --- Laplace
	double precision, allocatable :: phi(:)
	double precision, allocatable :: rhs(:)
	double precision, allocatable :: ue(:), ve(:)
	double precision, allocatable :: ea(:,:)
	integer :: mgp
	double precision :: area_total
	double precision, allocatable :: shpx(:,:,:), shpy(:,:,:) ! MODIFIED: (6, mgp, nelm)
	double precision, allocatable :: jac(:,:)
	double precision, allocatable :: xi_hat(:), eta_hat(:), sg(:)
end module indata
!
!------------------------				 MAIN PROGRAM 			------------------------------
!
program main
	use indata
	implicit none
	double precision :: begin_time, end_time

	write(*,*) 'main program start! (Quadratic Elements)'
	call CPU_TIME(begin_time)

    call datain
	call makeSSF
	call makeLHS
	call makeRHS
	call boundc
	call sweep
	call calvel
	call output

	call CPU_TIME(end_time)
	write(*,'(a,f15.6,a)') 'Time of operation was ' , end_time - begin_time, ' seconds'
	write(*,*) 'main program end!'
end program main
!
!----------------------       MAIN PROGRAM END       ----------------------------
!
!-------------------------------------------------------------------------------
subroutine datain
! -------------------------------------------------------------------------------
	use indata
	implicit none
	integer :: i, j, n, m
	character(100) :: mesfile, bdcfile, outfile, out2file, integfile
	character(20) :: bc_type_str
	integer :: n_vel_in, n_vel_wall

	mgp = 3 ! Number of integration points for quadratic elements

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
	open(61,file = out2file,  status = 'replace')
	open(62,file = integfile, status = 'replace')

	! READ MESH
	read(10,*) node, nelm
	allocate( xx(2,node), nc(6,nelm) ) ! MODIFIED for 6 nodes
	read(10,*) ( n, xx(1,n), xx(2,n), i = 1, node )
	read(10,*) ( m, ( nc(i,m), i = 1, 6 ), j = 1, nelm) ! MODIFIED for 6 nodes
	close(10)

	! READ BOUNDARY CONDITION FILE (3 nodes per edge for quadratic)
	read(11,*) n_potential_edges, bc_type_str
	allocate(potential_nodes(3, n_potential_edges), potential_values(n_potential_edges)) ! MODIFIED
	do i = 1, n_potential_edges
		read(11,*) j, potential_nodes(1,i), potential_nodes(2,i), potential_nodes(3,i), potential_values(i) ! MODIFIED
	end do

      read(11,*) ! Skip empty line
	read(11,*) n_vel_in, bc_type_str
      read(11,*) ! Skip empty line
	read(11,*) n_vel_wall, bc_type_str
	
	n_velocity_edges = n_vel_in + n_vel_wall
	allocate(velocity_nodes(3, n_velocity_edges), velocity_values(n_velocity_edges)) ! MODIFIED

	rewind(11)
	read(11,*) j, bc_type_str; do i = 1, j; read(11,*); end do; read(11,*)
	read(11,*) n_vel_in, bc_type_str
	do i = 1, n_vel_in
		read(11,*) j, velocity_nodes(1,i), velocity_nodes(2,i), velocity_nodes(3,i), velocity_values(i) ! MODIFIED
	end do

      read(11,*) ! Skip empty line
	read(11,*) n_vel_wall, bc_type_str
	do i = 1, n_vel_wall
		read(11,*) j, velocity_nodes(1, i + n_vel_in), velocity_nodes(2, i + n_vel_in), &
                     velocity_nodes(3, i + n_vel_in), velocity_values(i + n_vel_in) ! MODIFIED
	end do

	close(11)

	! Dynamic Memory Allocation
	allocate( phi(node), rhs(node), ea(node,node) )
	allocate( xi_hat(mgp), eta_hat(mgp), sg(mgp) )
	allocate( jac(mgp,nelm) )
	allocate( shpx(6,mgp,nelm), shpy(6,mgp,nelm) ) ! MODIFIED
	allocate( ue(nelm), ve(nelm) )

end subroutine datain
!
! -------------------------------------------------------------------------------
subroutine makeRHS   ! Make R.H.S. Vector from Neumann conditions
! -------------------------------------------------------------------------------
	use indata
	implicit none
	integer :: i, n1, n2, n3
	double precision :: x1, y1, x2, y2, edge_length, value

	! Initialize RHS vector
	rhs(:) = 0.0d0

	! Loop over all velocity (Neumann) boundary edges
	do i = 1, n_velocity_edges
		n1 = velocity_nodes(1, i) ! Corner node 1
		n2 = velocity_nodes(2, i) ! Corner node 2
        n3 = velocity_nodes(3, i) ! Midside node
		value = velocity_values(i)

		x1 = xx(1, n1); y1 = xx(2, n1)
		x2 = xx(1, n2); y2 = xx(2, n2)

		edge_length = sqrt((x2 - x1)**2 + (y2 - y1)**2)

		! MODIFIED: Distribute flux for quadratic element (1/6, 4/6, 1/6 rule)
		rhs(n1) = rhs(n1) + value * edge_length / 6.0d0
		rhs(n2) = rhs(n2) + value * edge_length / 6.0d0
        rhs(n3) = rhs(n3) + value * edge_length * (4.0d0 / 6.0d0)
	end do

end subroutine makeRHS
!
! -------------------------------------------------------------------------------
subroutine boundc      ! Impose Dirichlet Boundary Condition
! -------------------------------------------------------------------------------
	use indata
	implicit none
	integer :: i, j, n
	double precision :: f
	logical, allocatable :: is_dirichlet_node(:)

	allocate(is_dirichlet_node(node))
	is_dirichlet_node = .false.

	! ALGORITHM IMPROVEMENT: First, mark all unique nodes that have a Dirichlet condition
	do i = 1, n_potential_edges
		is_dirichlet_node(potential_nodes(1, i)) = .true.
		is_dirichlet_node(potential_nodes(2, i)) = .true.
            is_dirichlet_node(potential_nodes(3, i)) = .true.
	end do

	! Apply Dirichlet conditions
	do i = 1, n_potential_edges
        f = potential_values(i)
        do j = 1, 3 ! Loop over all 3 nodes of the edge
		    n = potential_nodes(j, i)
		    if (is_dirichlet_node(n)) then
                ! Modify RHS vector
			    rhs(:) = rhs(:) - ea(:, n) * f
                ! Modify LHS matrix
			    ea(n, :) = 0.0d0
			    ea(:, n) = 0.0d0
			    ea(n, n) = 1.0d0
			    rhs(n) = f
			    is_dirichlet_node(n) = .false. ! Mark as processed to avoid re-application
		    end if
        end do
	end do
	
	deallocate(is_dirichlet_node)

end subroutine boundc
!
! ------------------------------------------------------------------------------
subroutine sweep  ! Solve Linear-Systems (Gauss Elimination)
! (This subroutine does not need to be changed)
! ------------------------------------------------------------------------------
	use indata
	implicit none
	integer :: i, j, k, i1
	double precision :: ai, cc

	phi = rhs

	! Forward elimination
	do i = 1, node
		if( dabs(ea(i,i)) <= 1.0d-12 ) then
			write(6,*) 'Diagonal Value',i,'is ZERO or too small!'
			stop
		end if

		ai = 1.0d0 / ea(i,i)
		phi(i) = phi(i) * ai
		ea(i,i+1:node) = ea(i,i+1:node) * ai

		if( i < node ) then
			i1 = i + 1
			do k = i1, node
				cc = ea(k,i)
				ea(k,i1:node) = ea(k,i1:node) - cc * ea(i,i1:node)
				phi(k) = phi(k) - cc * phi(i)
			end do
		end if
	end do

	! Backward substitution
	do i = node - 1, 1, -1
		phi(i) = phi(i) - dot_product(ea(i,i+1:node), phi(i+1:node))
	end do

end subroutine sweep
!
! -------------------------------------------------------------------------------
subroutine makeSSF ! MODIFIED for Quadratic Elements
! -------------------------------------------------------------------------------
	use indata
	implicit none
	integer :: m, ig, i
	integer, dimension(6) :: n_local
	double precision :: DetJ, tmp
	double precision, allocatable :: dxdxi(:,:), dxidx(:,:)
      double precision, allocatable :: sqd(:,:,:), sq(:,:)
      double precision :: L1, L2, L3

	allocate( dxdxi(2,2), dxidx(2,2), sqd(2,6,mgp) )
	call set_integration_points(mgp, xi_hat, eta_hat, sg)

	do m = 1, nelm
        n_local = nc(:,m)

		do ig = 1, mgp
            L2 = xi_hat(ig)
            L3 = eta_hat(ig)
            L1 = 1.0d0 - L2 - L3

            sq(1,ig) = L1 * (2.0d0 * L1 - 1.0d0)
            sq(2,ig) = L2 * (2.0d0 * L2 - 1.0d0)
            sq(3,ig) = L3 * (2.0d0 * L3 - 1.0d0)

            sq(4,ig) = 4.0d0 * L1 * L2
            sq(5,ig) = 4.0d0 * L2 * L3
            sq(6,ig) = 4.0d0 * L3 * L1

            ! Derivatives of shape functions w.r.t. xi (L2) and eta (L3)
            ! dN/d_xi = (dN/dL1)*(-1) + (dN/dL2)*(1)
            sqd(1,1,ig) = (4.0d0*L1 - 1.0d0) * (-1.0d0)
            sqd(1,2,ig) = (4.0d0*L2 - 1.0d0) * ( 1.0d0)
            sqd(1,3,ig) = 0.0d0
            sqd(1,4,ig) = 4.0d0 * (L2*(-1.0d0) + L1*(1.0d0))
            sqd(1,5,ig) = 4.0d0 * L3
            sqd(1,6,ig) = 4.0d0 * L3 * (-1.0d0)

            ! dN/d_eta = (dN/dL1)*(-1) + (dN/dL3)*(1)
            sqd(2,1,ig) = (4.0d0*L1 - 1.0d0) * (-1.0d0)
            sqd(2,2,ig) = 0.0d0
            sqd(2,3,ig) = (4.0d0*L3 - 1.0d0) * ( 1.0d0)
            sqd(2,4,ig) = 4.0d0 * L2 * (-1.0d0)
            sqd(2,5,ig) = 4.0d0 * L2
            sqd(2,6,ig) = 4.0d0 * (L3*(-1.0d0) + L1*(1.0d0))
			
            ! Calculate Jacobian matrix [d(x,y)/d(xi,eta)]
            dxdxi(:,:) = 0.0d0
            do i = 1, 6
                dxdxi(1,1) = dxdxi(1,1) + sqd(1,i,ig)*xx(1,n_local(i)) ! dx/d_xi
                dxdxi(1,2) = dxdxi(1,2) + sqd(2,i,ig)*xx(1,n_local(i)) ! dx/d_eta
                dxdxi(2,1) = dxdxi(2,1) + sqd(1,i,ig)*xx(2,n_local(i)) ! dy/d_xi
                dxdxi(2,2) = dxdxi(2,2) + sqd(2,i,ig)*xx(2,n_local(i)) ! dy/d_eta
            end do
			
            DetJ = dxdxi(1,1) * dxdxi(2,2) - dxdxi(1,2) * dxdxi(2,1)
			if ( DetJ <= 0.0d0 ) then
				write(*,*) 'error in J' , m ,ig, DetJ; stop
			else
				tmp = 1.0d0 / DetJ
                jac(ig,m) = DetJ
			end if

			dxidx(1,1) =  dxdxi(2,2) * tmp; dxidx(1,2) = -dxdxi(1,2) * tmp
			dxidx(2,1) = -dxdxi(2,1) * tmp; dxidx(2,2) =  dxdxi(1,1) * tmp
			
            ! Map derivatives to physical coordinates (x,y)
            do i = 1, 6
				shpx(i,ig,m) = sqd(1,i,ig) * dxidx(1,1) + sqd(2,i,ig) * dxidx(2,1)
				shpy(i,ig,m) = sqd(1,i,ig) * dxidx(1,2) + sqd(2,i,ig) * dxidx(2,2)
			end do
		end do
	end do
	deallocate( dxdxi, dxidx, sqd )
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
		do j = 1, 6 ! MODIFIED
			do i = 1, 6 ! MODIFIED
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
subroutine calvel
! ------------------------------------------------------------------------------
	use indata
	implicit none
	integer :: i, m, ig
    double precision :: u_int, v_int, w_sum

	do m = 1, nelm
        u_int = 0.0d0
        v_int = 0.0d0
        w_sum = 0.0d0
        do ig = 1, mgp
            do i = 1, 6 ! MODIFIED
                u_int = u_int - shpx(i,ig,m) * phi(nc(i,m)) * jac(ig,m) * sg(ig)
                v_int = v_int - shpy(i,ig,m) * phi(nc(i,m)) * jac(ig,m) * sg(ig)
            end do
            w_sum = w_sum + jac(ig,m) * sg(ig)
        end do
        ! ALGORITHM IMPROVEMENT: Calculate area-averaged velocity
        if (w_sum > 1.0d-12) then
            ue(m) = u_int / w_sum
            ve(m) = v_int / w_sum
        else
            ue(m) = 0.0d0
            ve(m) = 0.0d0
        end if
	end do
end subroutine calvel
!
! -------------------------------------------------------------------------------
subroutine output
! -------------------------------------------------------------------------------
    use indata
	implicit none
	integer :: n
	write(60,600) ( n, phi(n), n =1, node )
	write(60,601) ( n, ue(n), ve(n), n =1, nelm )
	close(60)
	600 format(i7, d15.6)
	601 format(i7,2d15.6)
end subroutine output