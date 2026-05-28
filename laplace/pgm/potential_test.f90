! ==================================================================================
!
!              2-Dimensional Finite Element Analysis(Numerical Integration)
!                 of Steady Ideal Flows Governed by Laplace Eq.
!                 for Quadratic Triangular Elements (6-node Tria) - Corrected
!
! ==================================================================================
module indata
	implicit none
	integer :: node, nelm
	integer, allocatable :: nc(:,:)
	double precision, allocatable :: xx(:,:)

	! --- boundary condition (Edge-based for Quadratic Elements)
	integer :: n_potential_edges
	integer, allocatable :: potential_nodes(:,:)
	double precision, allocatable :: potential_values(:)
	integer :: n_velocity_edges
	integer, allocatable :: velocity_nodes(:,:)
	double precision, allocatable :: velocity_values(:)

	! --- Laplace
	double precision, allocatable :: phi(:)
	double precision, allocatable :: rhs(:)
	double precision, allocatable :: ue(:), ve(:)
	double precision, allocatable :: ea(:,:)
	integer :: mgp
	double precision, allocatable :: shpx(:,:,:), shpy(:,:,:)
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
      call set_integration_points
	call makeSSF !OK
	call makeLHS !OK
	call makeRHS
      call add_source_term 
	call boundc
	call sweep !OK
	call calvel !OK
	call output !OK

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

	mgp = 3 ! 3-point quadrature is suitable for quadratic elements

	! FILE I/O
	open(9, file = 'file_FEM.dat', status = 'old' )
	read(9,'(a)') mesfile
	read(9,'(a)') bdcfile
	read(9,'(a)') outfile
	close(9)

	open(10,file = mesfile,   status = 'old',action='read')
	open(11,file = bdcfile,   status = 'old',action='read')
	open(60,file = outfile,   status = 'replace')

	! READ MESH
	read(10,*) node, nelm
	allocate( xx(2,node), nc(6,nelm) )
	do i = 1, node
		read(10,*) n, xx(1,i), xx(2,i)
	end do
	do m = 1, nelm
		read(10,*) j, (nc(i,m), i = 1, 6)
	end do
	close(10)

	! READ BOUNDARY CONDITION FILE (3 nodes per edge for quadratic)
	read(11,*) n_potential_edges, bc_type_str
      if (n_potential_edges > 0) then
	      allocate(potential_nodes(3, n_potential_edges), potential_values(n_potential_edges))
	      do i = 1, n_potential_edges
		      read(11,*) j, potential_nodes(1,i), potential_nodes(2,i), potential_nodes(3,i), potential_values(i)
	      end do
      else ! Allocate dummy array if no conditions are present
            allocate(potential_nodes(3, 1), potential_values(1))
      end if

      read(11,*, end=111) ! Read until EOF or next line
111 continue
	read(11,*, end=112) n_vel_in, bc_type_str
112 continue
      read(11,*, end=113)
113 continue
	read(11,*, end=114) n_vel_wall, bc_type_str
114 continue
	
	n_velocity_edges = n_vel_in + n_vel_wall
      if (n_velocity_edges > 0) then
	      allocate(velocity_nodes(3, n_velocity_edges), velocity_values(n_velocity_edges))
      else
            allocate(velocity_nodes(3, 1), velocity_values(1))
      end if

	rewind(11)
	read(11,*, end=211) j, bc_type_str
      do i = 1, n_potential_edges
            read(11,*)
      end do
211 continue
      read(11,*, end=212)
	read(11,*, end=212) j, bc_type_str
	      do i = 1, n_vel_in
		      read(11,*) j, velocity_nodes(1,i), velocity_nodes(2,i), velocity_nodes(3,i), velocity_values(i)
	      end do
212 continue
      read(11,*, end=213)
	read(11,*, end=213) j, bc_type_str
	      do i = 1, n_vel_wall
		      read(11,*) j, velocity_nodes(1, i + n_vel_in), velocity_nodes(2, i + n_vel_in), &
                              velocity_nodes(3, i + n_vel_in), velocity_values(i + n_vel_in)
	      end do
213 continue
	close(11)

	! Dynamic Memory Allocation
	allocate( phi(node), rhs(node), ea(node,node) )
	allocate( xi_hat(mgp), eta_hat(mgp), sg(mgp) )
	allocate( jac(mgp,nelm) )
	allocate( shpx(6,mgp,nelm), shpy(6,mgp,nelm) )
	allocate( ue(nelm), ve(nelm) )

end subroutine datain
! -------------------------------------------------------------------------------
! 2次元三角形要素の積分点（ガウス点）と重みを設定する
! -------------------------------------------------------------------------------
subroutine set_integration_points
    use indata
    implicit none
    if (mgp == 3) then
      xi_hat(1)  = 0.5d0
      eta_hat(1) = 0.0d0
      sg(1)      = 1.0d0 / 6.0d0

      xi_hat(2)  = 0.5d0
      eta_hat(2) = 0.5d0
      sg(2)      = 1.0d0 / 6.0d0

      xi_hat(3)  = 0.0d0
      eta_hat(3) = 0.5d0
      sg(3)      = 1.0d0 / 6.0d0
    elseif (mgp == 2) then
      xi_hat(1) = 0.57735
      sg(1)     = 1.0
      xi_hat(2) = 0.57735
      sg(2)     = 1.0
    end if
end subroutine set_integration_points

!
! -------------------------------------------------------------------------------
subroutine makeRHS
	use indata
	implicit none
	integer :: i, k, ig
	integer :: n(3) ! n(1), n(2) are end nodes, n(3) is the midpoint node
	double precision :: val
      double precision :: x_coords(3), y_coords(3)
      double precision :: Js, dx_ds, dy_ds
    
    ! --- Variables for 1D Gauss Quadrature ---
    integer, parameter :: mgp_line = 2 ! Number of integration points for the line
    double precision :: s_hat(mgp_line), w_hat(mgp_line)
    double precision :: s, w
    double precision :: N_1D(3) ! Shape functions for the 3 nodes on the edge

	rhs(:) = 0.0d0
    call set_line_integration_points(mgp_line, s_hat, w_hat)

	do k = 1, n_velocity_edges
        ! Get nodes and value for the current boundary edge
		n(1) = velocity_nodes(1, k)
		n(2) = velocity_nodes(2, k)
        n(3) = velocity_nodes(3, k)
		val  = velocity_values(k)

        ! Get coordinates of the 3 nodes
        x_coords(1) = xx(1, n(1)); y_coords(1) = xx(2, n(1))
        x_coords(2) = xx(1, n(2)); y_coords(2) = xx(2, n(2))
        x_coords(3) = xx(1, n(3)); y_coords(3) = xx(2, n(3))

        ! Loop over integration points (implements the summation Sigma)
		do ig = 1, mgp_line
            s = s_hat(ig)
            w = w_hat(ig)

            ! Calculate 1D quadratic shape functions at the integration point 's'
            ! Node 1 corresponds to s=-1, Node 2 to s=+1, Node 3 to s=0
            N_1D(1) = 0.5d0 * s * (s - 1.0d0)
            N_1D(2) = 0.5d0 * s * (s + 1.0d0)
            N_1D(3) = 1.0d0 - s**2

            ! Calculate Jacobian for the 1D mapping (for a straight line)
            ! Assumes n(1) and n(2) are the corner nodes of the edge
            dx_ds = (x_coords(2) - x_coords(1)) / 2.0d0
            dy_ds = (y_coords(2) - y_coords(1)) / 2.0d0
            Js = sqrt(dx_ds**2 + dy_ds**2) ! Corresponds to edge_length / 2

            ! Add the contribution from this integration point to the RHS vector
            ! This is the core of the numerical integration
            do i = 1, 3
                rhs(n(i)) = rhs(n(i)) + val * N_1D(i) * Js * w
            end do
		end do
	end do
end subroutine makeRHS
!
! -------------------------------------------------------------------------------
! Subroutine to add the source term integral to the RHS vector
! Calculates integral of (f * N_i) over the domain
! -------------------------------------------------------------------------------
subroutine add_source_term
    use indata
    implicit none
    integer :: m, i, ig, ni
    double precision :: f_val, DetJ, w
    double precision :: L(3), N_vals(6)

    ! Loop over all elements in the mesh
    do m = 1, nelm
        ! Loop over the integration points within the element
        do ig = 1, mgp
            ! --- At each integration point, do the following ---

            ! (1) Get Jacobian and integration weight (pre-calculated in makeSSF)
            DetJ = jac(ig, m)
            w    = sg(ig)

            ! (2) Define the source term value at this point
            f_val = 1.0d0

            ! (3) Calculate shape function values N_i at this point
            L(2) = xi_hat(ig)  ! L2 = xi
            L(3) = eta_hat(ig)  ! L3 = eta
            L(1) = 1.0d0 - L(2) - L(3)

            N_vals(1) = L(1) * (2.0d0 * L(1) - 1.0d0)
            N_vals(2) = L(2) * (2.0d0 * L(2) - 1.0d0)
            N_vals(3) = L(3) * (2.0d0 * L(3) - 1.0d0)
            N_vals(4) = 4.0d0 * L(1) * L(2)
            N_vals(5) = 4.0d0 * L(2) * L(3)
            N_vals(6) = 4.0d0 * L(3) * L(1)

            ! (4) Calculate the contribution and add to the global RHS vector
            do i = 1, 6
                ni = nc(i, m) ! Get the global node index
                rhs(ni) = rhs(ni) + f_val * N_vals(i) * DetJ * w
            end do
        end do
    end do

end subroutine add_source_term
! -------------------------------------------------------------------------------
subroutine boundc
	use indata
	implicit none
	integer :: i, j, n
	double precision :: f
	logical, allocatable :: is_dirichlet_node(:)

	allocate(is_dirichlet_node(node))
	is_dirichlet_node = .false.

	do i = 1, n_potential_edges
		is_dirichlet_node(potential_nodes(1, i)) = .true.
		is_dirichlet_node(potential_nodes(2, i)) = .true.
            is_dirichlet_node(potential_nodes(3, i)) = .true.
	end do

	do i = 1, n_potential_edges
            f = potential_values(i)
            do j = 1, 3
		    n = potential_nodes(j, i)
		    if (is_dirichlet_node(n)) then
			    rhs(:) = rhs(:) - ea(:, n) * f
			    ea(n, :) = 0.0d0
			    ea(:, n) = 0.0d0
			    ea(n, n) = 1.0d0
			    rhs(n) = f
			    is_dirichlet_node(n) = .false.
		    end if
            end do
	end do
	
	deallocate(is_dirichlet_node)

end subroutine boundc
!
! ------------------------------------------------------------------------------
subroutine sweep
	use indata
	implicit none
	integer :: i, k
	double precision :: ai, cc

	phi = rhs

	do i = 1, node
		if( dabs(ea(i,i)) <= 1.0d-12 ) then
			write(6,*) 'Diagonal Value',i,'is ZERO or too small!'; stop
		end if
		ai = 1.0d0 / ea(i,i)
		phi(i) = phi(i) * ai
		if (i < node) then
			ea(i,i+1:node) = ea(i,i+1:node) * ai
			do k = i + 1, node
				cc = ea(k,i)
				ea(k,i+1:node) = ea(k,i+1:node) - cc * ea(i,i+1:node)
				phi(k) = phi(k) - cc * phi(i)
			end do
		end if
	end do

	do i = node - 1, 1, -1
		phi(i) = phi(i) - dot_product(ea(i,i+1:node), phi(i+1:node))
	end do
end subroutine sweep
!
! -------------------------------------------------------------------------------
subroutine makeSSF
	use indata
	implicit none
	integer :: m, ig, i
	integer, dimension(6) :: n_local
	double precision :: DetJ, tmp
	double precision, allocatable :: dxdxi(:,:), dxidx(:,:)
    double precision, allocatable :: sqd(:,:,:)
    
    ! --- Variables for Chain Rule Calculation ---
    double precision :: L(3)         ! (L1, L2, L3)
    double precision :: N_vals(6)    ! Shape function values N(L1,L2,L3)
    double precision :: dNdL(6,3)    ! dN/dL: Derivatives of N w.r.t. Area Coords
    double precision, parameter :: dLdxi(3,2) = reshape([ &
        -1.0d0, -1.0d0, & ! dL1/d_xi, dL1/d_eta
         1.0d0,  0.0d0, & ! dL2/d_xi, dL2/d_eta
         0.0d0,  1.0d0], & ! dL3/d_xi, dL3/d_eta
         shape=[3,2])
    double precision :: dNdxi(6,2)   ! dN/d_xi, dN/d_eta: Final derivatives in reference coords

	allocate( dxdxi(2,2), dxidx(2,2), sqd(2,6,mgp) )
	call set_integration_points(mgp, xi_hat, eta_hat, sg)

	do m = 1, nelm
        n_local = nc(:,m)
		do ig = 1, mgp
            ! Step 1: Define Area Coordinates (L1, L2, L3) from reference coordinates
            L(2) = xi_hat(ig)  ! L2 = xi
            L(3) = eta_hat(ig)  ! L3 = eta
            L(1) = 1.0d0 - L(2) - L(3)

            ! ==========================================================
            !  *** ADDED PART: Calculate and store shape function values ***
            ! ==========================================================
            ! For Corner Nodes (Ni = Li * (2*Li - 1))
            N_vals(1) = L(1) * (2.0d0 * L(1) - 1.0d0)
            N_vals(2) = L(2) * (2.0d0 * L(2) - 1.0d0)
            N_vals(3) = L(3) * (2.0d0 * L(3) - 1.0d0)
            ! For Midpoint Nodes
            N_vals(4) = 4.0d0 * L(1) * L(2) ! N4 = 4*L1*L2
            N_vals(5) = 4.0d0 * L(2) * L(3) ! N5 = 4*L2*L3
            N_vals(6) = 4.0d0 * L(3) * L(1) ! N6 = 4*L3*L1
            ! (Note: N_vals is not used later in this specific subroutine,
            !  but storing it is a clear and logical step.)

            ! Step 2: Calculate derivatives of shape functions w.r.t. Area Coordinates (dN/dL)
            dNdL(1,1) = 4.0d0 * L(1) - 1.0d0; dNdL(1,2) = 0.0d0;                dNdL(1,3) = 0.0d0
            dNdL(2,1) = 0.0d0;                dNdL(2,2) = 4.0d0 * L(2) - 1.0d0; dNdL(2,3) = 0.0d0
            dNdL(3,1) = 0.0d0;                dNdL(3,2) = 0.0d0;                dNdL(3,3) = 4.0d0 * L(3) - 1.0d0
            dNdL(4,1) = 4.0d0 * L(2);         dNdL(4,2) = 4.0d0 * L(1);         dNdL(4,3) = 0.0d0
            dNdL(5,1) = 0.0d0;                dNdL(5,2) = 4.0d0 * L(3);         dNdL(5,3) = 4.0d0 * L(2)
            dNdL(6,1) = 4.0d0 * L(3);         dNdL(6,2) = 0.0d0;                dNdL(6,3) = 4.0d0 * L(1)

            ! Step 3: Apply Chain Rule using matrix multiplication to get dN/d(xi,eta)
            dNdxi = matmul(dNdL, dLdxi)

            ! Store results into the sqd array for subsequent calculations
            sqd(1, :, ig) = dNdxi(:, 1) ! dN/d_xi
            sqd(2, :, ig) = dNdxi(:, 2) ! dN/d_eta
			
            ! --- The rest of the routine remains the same ---
            dxdxi(:,:) = 0.0d0
            do i = 1, 6
                dxdxi(1,1) = dxdxi(1,1) + sqd(1,i,ig)*xx(1,n_local(i))
                dxdxi(1,2) = dxdxi(1,2) + sqd(2,i,ig)*xx(1,n_local(i))
                dxdxi(2,1) = dxdxi(2,1) + sqd(1,i,ig)*xx(2,n_local(i))
                dxdxi(2,2) = dxdxi(2,2) + sqd(2,i,ig)*xx(2,n_local(i))
            end do
			
            DetJ = dxdxi(1,1) * dxdxi(2,2) - dxdxi(1,2) * dxdxi(2,1)
			if ( DetJ <= 1.0d-12 ) then
				write(*,*) 'error: Jacobian is zero or negative at element, igp:', m, ig, DetJ; stop
			else
				tmp = 1.0d0 / DetJ
                jac(ig,m) = DetJ
			end if

			dxidx(1,1) =  dxdxi(2,2) * tmp; dxidx(1,2) = -dxdxi(1,2) * tmp
			dxidx(2,1) = -dxdxi(2,1) * tmp; dxidx(2,2) =  dxdxi(1,1) * tmp
			
            do i = 1, 6
				shpx(i,ig,m) = sqd(1,i,ig) * dxidx(1,1) + sqd(2,i,ig) * dxidx(2,1)
				shpy(i,ig,m) = sqd(1,i,ig) * dxidx(1,2) + sqd(2,i,ig) * dxidx(2,2)
			end do
		end do
	end do
	deallocate( dxdxi, dxidx, sqd )
end subroutine makeSSF

! -------------------------------------------------------------------------------
! Helper subroutine for 1D Gauss Quadrature points and weights
! -------------------------------------------------------------------------------
subroutine set_line_integration_points(mgp_line, s_hat, w_hat)
    implicit none
    integer, intent(in) :: mgp_line
    double precision, intent(out) :: s_hat(mgp_line), w_hat(mgp_line)
    
    if (mgp_line == 2) then ! 2-point quadrature is exact for polynomials up to degree 3
        s_hat(1) = -1.0d0 / sqrt(3.0d0)
        s_hat(2) =  1.0d0 / sqrt(3.0d0)
        w_hat(1) =  1.0d0
        w_hat(2) =  1.0d0
    else if (mgp_line == 1) then
        s_hat(1) = 0.0d0
        w_hat(1) = 2.0d0
    else
        write(*,*) 'Error: Unsupported number of line integration points:', mgp_line
        stop
    end if
end subroutine set_line_integration_points
!
! -------------------------------------------------------------------------------
subroutine makeLHS
	use indata
	implicit none
	integer :: i, j, m, ig, ni, nj
	ea(:,:) = 0.0d0
	do m = 1, nelm
		do j = 1, 6
			do i = 1, 6
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
	use indata
	implicit none
	integer :: i, m, ig
    double precision :: u_int, v_int, area_e

	do m = 1, nelm
        u_int = 0.0d0
        v_int = 0.0d0
        area_e = 0.0d0
        do ig = 1, mgp
            do i = 1, 6
                u_int = u_int - shpx(i,ig,m) * phi(nc(i,m)) * jac(ig,m) * sg(ig)
                v_int = v_int - shpy(i,ig,m) * phi(nc(i,m)) * jac(ig,m) * sg(ig)
            end do
            area_e = area_e + jac(ig,m) * sg(ig)
        end do
        if (area_e > 1.0d-12) then
            ue(m) = u_int / area_e
            ve(m) = v_int / area_e
        else
            ue(m) = 0.0d0
            ve(m) = 0.0d0
        end if
	end do
end subroutine calvel
!
! -------------------------------------------------------------------------------
subroutine output
    use indata
	implicit none
	integer :: n
	write(60,600) ( n, phi(n), n =1, node )
	write(60,601) ( n, ue(n), ve(n), n =1, nelm )
	close(60)
	600 format(i7, d15.6)
	601 format(i7,2d15.6)
end subroutine output