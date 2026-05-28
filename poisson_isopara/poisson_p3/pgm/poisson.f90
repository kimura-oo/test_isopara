module indata
	implicit none
	integer :: node, nelm
	integer, allocatable :: nc(:,:)
	double precision, allocatable :: xx(:,:)

	integer :: n_potential_nodes
	integer, allocatable :: potential_nodes(:)
	double precision, allocatable :: potential_values(:)
      integer :: n_velocity_edges
	integer, allocatable :: velocity_nodes(:,:)
	double precision, allocatable :: velocity_values(:)

	double precision, allocatable :: phi(:)
	double precision, allocatable :: rhs(:)
	double precision, allocatable :: ue(:), ve(:)
	double precision, allocatable :: ea(:,:)
	integer :: mgp
	double precision, allocatable :: shpx(:,:,:), shpy(:,:,:)
	double precision, allocatable :: jac(:,:)
	double precision, allocatable :: xi_hat(:), eta_hat(:), sg(:)
      integer :: mgp_line
      double precision, allocatable :: s_hat(:), w_hat(:)
end module indata
!
!------------------------				 MAIN PROGRAM 			------------------------------
!
program main
	use indata
	implicit none
	double precision :: begin_time, end_time

	write(*,*) 'main program start! (Cubic Elements)'
	call CPU_TIME(begin_time)
      
      call datain
      call set_integration_points
	call makeSSF
	call makeLHS
      call set_line_integration_points
	call makeRHS
      call add_source_term 
	call boundc
	call sweep
	call calvel
	call output
      call output2 
      call output_x05_nodes
      call output_y05_nodes

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
	character(100) :: mesfile, bdcfile, outfile, out2file, integfile, outfile_vtk, outfile_x05, outfile_y05

	mgp = 7 

	! FILE I/O
	open(9, file = 'file_FEM.dat', status = 'old' )
	read(9,'(a)') mesfile
	read(9,'(a)') bdcfile
	read(9,'(a)') outfile
      read(9,'(a)') outfile_vtk
      read(9,'(a)') outfile_x05
      read(9,'(a)') outfile_y05
	close(9)

	open(10,file = mesfile,   status = 'old',action='read')
	open(11,file = bdcfile,   status = 'old',action='read')
	open(60,file = outfile,   status = 'replace')
      open(70,file = outfile_vtk, status = 'replace')
      open(80,file = outfile_x05, status = 'replace')
      open(90,file = outfile_y05, status = 'replace')

	! READ MESH
	read(10,*) node, nelm
	allocate( xx(2,node), nc(10,nelm) )
	do i = 1, node
		read(10,*) n, xx(1,i), xx(2,i)
	end do
	do m = 1, nelm
		read(10,*) j, (nc(i,m), i = 1, 10)
	end do
	close(10)

      ! READ BOUNDARY CONDITION FILE
      read(11,*) n_potential_nodes
	      allocate(potential_nodes(n_potential_nodes), potential_values(n_potential_nodes))
	      do i = 1, n_potential_nodes
		      read(11,*) j, potential_nodes(i), potential_values(i)
	      end do
      close(11)

	allocate( phi(node), rhs(node), ea(node,node) )
	allocate( xi_hat(mgp), eta_hat(mgp), sg(mgp) )
	allocate( jac(mgp,nelm) )
	allocate( shpx(10,mgp,nelm), shpy(10,mgp,nelm) )
	allocate( ue(nelm), ve(nelm) )

      mgp_line = 2
      allocate(s_hat(mgp_line), w_hat(mgp_line))
end subroutine datain
! -------------------------------------------------------------------------------
! Set integration points and weights for triangular elements
! -------------------------------------------------------------------------------
!
subroutine set_integration_points
      use indata
      implicit none
      if (mgp == 1) then
            xi_hat(1)  = 0.3333333333d0
            eta_hat(1) = 0.3333333333d0
            sg(1)      = 1.0d0
      elseif (mgp == 3) then
            xi_hat(1)  = 0.1666666667d0
            eta_hat(1) = 0.1666666667d0
            sg(1)      = 0.1666666667d0

            xi_hat(2)  = 0.6666666667d0
            eta_hat(2) = 0.1666666667d0
            sg(2)      = 0.1666666667d0

            xi_hat(3)  = 0.1666666667d0
            eta_hat(3) = 0.6666666667d0
            sg(3)      = 0.1666666667d0
      elseif (mgp == 4) then
            xi_hat(1) = 0.3333333333d0
            eta_hat(1) = 0.3333333333d0
            sg(1)      = -0.28125d0

            xi_hat(2) = 0.6d0
            eta_hat(2) = 0.2d0
            sg(2) = 3.84d0

            xi_hat(3) = 0.2d0
            eta_hat(3) = 0.6d0
            sg(3) = 3.84d0

            xi_hat(4) = 0.2d0
            eta_hat(4) = 0.2d0
            sg(4) = 3.84d0
      elseif (mgp == 7) then
            xi_hat(1)  = 0.0d0
            eta_hat(1) = 0.0d0
            sg(1)      = 1.0d0/40.0d0  ! 0.025

            xi_hat(2)  = 0.5d0
            eta_hat(2) = 0.0d0
            sg(2)      = 1.0d0/15.0d0  ! ~0.0667

            xi_hat(3)  = 1.0d0
            eta_hat(3) = 0.0d0
            sg(3)      = 1.0d0/40.0d0  ! 0.025

            xi_hat(4)  = 0.5d0
            eta_hat(4) = 0.5d0
            sg(4)      = 1.0d0/15.0d0  ! ~0.0667

            xi_hat(5)  = 0.0d0
            eta_hat(5) = 1.0d0
            sg(5)      = 1.0d0/40.0d0  ! 0.025

            xi_hat(6)  = 0.0d0
            eta_hat(6) = 0.5d0
            sg(6)      = 1.0d0/15.0d0  ! ~0.0667

            xi_hat(7)  = 1.0d0/3.0d0
            eta_hat(7) = 1.0d0/3.0d0
            sg(7)      = 9.0d0/40.0d0  ! 0.225
      end if

end subroutine set_integration_points
! -------------------------------------------------------------------------------
! MODIFIED for CUBIC (10-node) TRIANGULAR ELEMENTS
! -------------------------------------------------------------------------------
subroutine makeSSF
	use indata
	implicit none
	integer :: m, ig, i
	integer, dimension(10) :: n_local
	double precision :: DetJ, tmp
	double precision, allocatable :: dxdxi(:,:), dxidx(:,:)
      double precision, allocatable :: sqd(:,:,:) 
      
      ! --- Variables for Chain Rule Calculation ---
      double precision :: L(3)         ! (L1, L2, L3)
      double precision :: N_vals(10)    ! Shape function values N(L1,L2,L3)
      double precision :: dNdL(10,3)    ! dN/dL: Derivatives of N w.r.t. Area Coords
      double precision :: dLdxi(3,2)
      double precision :: dNdxi(10,2)   ! dN/d_xi, dN/d_eta: Final derivatives in reference coords
	
      dLdxi(1,1) = -1.0d0; dLdxi(1,2) = -1.0d0
      dLdxi(2,1) =  1.0d0; dLdxi(2,2) =  0.0d0
      dLdxi(3,1) =  0.0d0; dLdxi(3,2) =  1.0d0
      allocate( dxdxi(2,2), dxidx(2,2), sqd(2,10,mgp) )
      
	do m = 1, nelm
            n_local = nc(1:10,m)
		do ig = 1, mgp
                  L(2) = xi_hat(ig)  ! L2 = xi
                  L(3) = eta_hat(ig)  ! L3 = eta
                  L(1) = 1.0d0 - L(2) - L(3)
                  ! Corner Nodes (1-3)
                  N_vals(1) = 0.5d0 * L(1) * (3.0d0*L(1) - 1.0d0) * (3.0d0*L(1) - 2.0d0)
                  N_vals(2) = 0.5d0 * L(2) * (3.0d0*L(2) - 1.0d0) * (3.0d0*L(2) - 2.0d0)
                  N_vals(3) = 0.5d0 * L(3) * (3.0d0*L(3) - 1.0d0) * (3.0d0*L(3) - 2.0d0)
                  ! Side Nodes (4-9)
                  N_vals(4) = 4.5d0 * L(1) * L(2) * (3.0d0*L(1) - 1.0d0) ! Edge 1-2, near V1
                  N_vals(5) = 4.5d0 * L(1) * L(2) * (3.0d0*L(2) - 1.0d0) ! Edge 1-2, near V2
                  N_vals(6) = 4.5d0 * L(2) * L(3) * (3.0d0*L(2) - 1.0d0) ! Edge 2-3, near V2
                  N_vals(7) = 4.5d0 * L(2) * L(3) * (3.0d0*L(3) - 1.0d0) ! Edge 2-3, near V3
                  N_vals(8) = 4.5d0 * L(3) * L(1) * (3.0d0*L(3) - 1.0d0) ! Edge 3-1, near V3
                  N_vals(9) = 4.5d0 * L(3) * L(1) * (3.0d0*L(1) - 1.0d0) ! Edge 3-1, near V1
                  ! Internal Node (10)
                  N_vals(10) = 27.0d0 * L(1) * L(2) * L(3)

                  ! Step 3: Calculate derivatives of CUBIC shape functions w.r.t. Area Coordinates (dN/dL)
                  dNdL = 0.0d0
                  ! Corner Nodes (1-3)
                  dNdL(1,1) = 0.5d0 * (27.0d0*L(1)**2 - 18.0d0*L(1) + 2.0d0)
                  dNdL(2,2) = 0.5d0 * (27.0d0*L(2)**2 - 18.0d0*L(2) + 2.0d0)
                  dNdL(3,3) = 0.5d0 * (27.0d0*L(3)**2 - 18.0d0*L(3) + 2.0d0)
                  ! Side Nodes (4-9)
                  dNdL(4,1) = 4.5d0 * L(2) * (6.0d0*L(1) - 1.0d0); dNdL(4,2) = 4.5d0 * L(1) * (3.0d0*L(1) - 1.0d0)
                  dNdL(5,1) = 4.5d0 * L(2) * (3.0d0*L(2) - 1.0d0); dNdL(5,2) = 4.5d0 * L(1) * (6.0d0*L(2) - 1.0d0)
                  dNdL(6,2) = 4.5d0 * L(3) * (6.0d0*L(2) - 1.0d0); dNdL(6,3) = 4.5d0 * L(2) * (3.0d0*L(2) - 1.0d0)
                  dNdL(7,2) = 4.5d0 * L(3) * (3.0d0*L(3) - 1.0d0); dNdL(7,3) = 4.5d0 * L(2) * (6.0d0*L(3) - 1.0d0)
                  dNdL(8,3) = 4.5d0 * L(1) * (6.0d0*L(3) - 1.0d0); dNdL(8,1) = 4.5d0 * L(3) * (3.0d0*L(3) - 1.0d0)
                  dNdL(9,3) = 4.5d0 * L(1) * (3.0d0*L(1) - 1.0d0); dNdL(9,1) = 4.5d0 * L(3) * (6.0d0*L(1) - 1.0d0)
                  ! Internal Node (10)
                  dNdL(10,1) = 27.0d0 * L(2) * L(3)
                  dNdL(10,2) = 27.0d0 * L(1) * L(3)
                  dNdL(10,3) = 27.0d0 * L(1) * L(2)

            ! Step 4: Apply Chain Rule
                  dNdxi = matmul(dNdL, dLdxi)

            ! Step 5: Store results into the sqd array
                  sqd(1, :, ig) = dNdxi(:, 1) ! dN/d_xi
                  sqd(2, :, ig) = dNdxi(:, 2) ! dN/d_eta
		
            ! Step 6: Assemble Jacobian Matrix and calculate its determinant
                  dxdxi = 0.0d0
                  do i = 1, 10 ! LOOP BOUND UPDATED TO 10
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
			
                  do i = 1, 10 ! LOOP BOUND UPDATED TO 10
				shpx(i,ig,m) = sqd(1,i,ig) * dxidx(1,1) + sqd(2,i,ig) * dxidx(2,1)
				shpy(i,ig,m) = sqd(1,i,ig) * dxidx(1,2) + sqd(2,i,ig) * dxidx(2,2)
			end do
		end do
	end do
deallocate( dxdxi, dxidx, sqd )
end subroutine makeSSF

! -------------------------------------------------------------------------------

subroutine makeRHS
      use indata
	implicit none
	integer :: i, k, ig
	integer :: n(4) ! 4 nodes per edge (2 ends + 2 internal)
	double precision :: val
      double precision :: x_coords(4), y_coords(4) 
      double precision :: Js, dx_ds, dy_ds
      double precision :: s, w
      double precision :: N_1D(4) 
      integer :: n_velocity_nodes
	rhs(:) = 0.0d0

	do k = 1, n_velocity_edges
      ! Get nodes and value for the current boundary edge (now 4 nodes)
	      n(1) = velocity_nodes(1, k) ! Start node (s=-1)
		n(2) = velocity_nodes(2, k) ! End node   (s=+1)
            n(3) = velocity_nodes(3, k) ! Internal node (s=-1/3)
            n(4) = velocity_nodes(4, k) ! Internal node (s=+1/3)
		val  = velocity_values(k)

      ! Get coordinates of the 4 nodes
            x_coords(1) = xx(1, n(1)); y_coords(1) = xx(2, n(1))
            x_coords(2) = xx(1, n(2)); y_coords(2) = xx(2, n(2))
            x_coords(3) = xx(1, n(3)); y_coords(3) = xx(2, n(3))
            x_coords(4) = xx(1, n(4)); y_coords(4) = xx(2, n(4))

		do ig = 1, mgp_line
                  s = s_hat(ig)
                  w = w_hat(ig)

            ! UPDATED: Calculate 1D CUBIC shape functions at integration point 's'
            ! Assumes nodes are at s = -1, +1, -1/3, +1/3
                  N_1D(1) = (1.0d0/16.0d0) * (-9.0d0*s**3 + 9.0d0*s**2 + s - 1.0d0)
                  N_1D(2) = (1.0d0/16.0d0) * ( 9.0d0*s**3 + 9.0d0*s**2 - s - 1.0d0)
                  N_1D(3) = (1.0d0/16.0d0) * ( 27.0d0*s**3 - 9.0d0*s**2 - 27.0d0*s + 9.0d0)
                  N_1D(4) = (1.0d0/16.0d0) * (-27.0d0*s**3 - 9.0d0*s**2 + 27.0d0*s + 9.0d0)

            ! Jacobian for 1D mapping (assumes a straight edge between end nodes)
                  dx_ds = (x_coords(2) - x_coords(1)) / 2.0d0
                  dy_ds = (y_coords(2) - y_coords(1)) / 2.0d0
                  Js = sqrt(dx_ds**2 + dy_ds**2)

            ! Add the contribution to the RHS vector
                  do i = 1, 4 ! UPDATED: Loop over all 4 nodes
                        rhs(n(i)) = rhs(n(i)) + val * N_1D(i) * Js * w
                  end do
		end do
	end do
end subroutine makeRHS
!
! ===================================================================
!   SUBROUTINE: add_source_term (MODIFIED)
!   PURPOSE: Add source term contribution to the RHS vector
!            for the analytical solution u = sin(pi*x)*sin(pi*y)
! ===================================================================
!
subroutine add_source_term
      use indata
      implicit none
      integer :: m, i, ig, ni
      double precision :: f_val_int, x_int, y_int, pi
      double precision :: N_vals(10)
      double precision :: L(3)

      pi = acos(-1.0d0)

    ! Loop over all elements to compute the domain integral
      do m = 1, nelm
            do ig = 1, mgp
      ! --- Step 1: Calculate shape function VALUES at the integration point ---
                  L(2) = xi_hat(ig)  ! L2 = xi
                  L(3) = eta_hat(ig) ! L3 = eta
                  L(1) = 1.0d0 - L(2) - L(3)
      ! Corner Nodes (1-3)
                  N_vals(1) = 0.5d0 * L(1) * (3.0d0*L(1) - 1.0d0) * (3.0d0*L(1) - 2.0d0)
                  N_vals(2) = 0.5d0 * L(2) * (3.0d0*L(2) - 1.0d0) * (3.0d0*L(2) - 2.0d0)
                  N_vals(3) = 0.5d0 * L(3) * (3.0d0*L(3) - 1.0d0) * (3.0d0*L(3) - 2.0d0)
      ! Side Nodes (4-9)
                  N_vals(4) = 4.5d0 * L(1) * L(2) * (3.0d0*L(1) - 1.0d0)
                  N_vals(5) = 4.5d0 * L(1) * L(2) * (3.0d0*L(2) - 1.0d0)
                  N_vals(6) = 4.5d0 * L(2) * L(3) * (3.0d0*L(2) - 1.0d0)
                  N_vals(7) = 4.5d0 * L(2) * L(3) * (3.0d0*L(3) - 1.0d0)
                  N_vals(8) = 4.5d0 * L(3) * L(1) * (3.0d0*L(3) - 1.0d0)
                  N_vals(9) = 4.5d0 * L(3) * L(1) * (3.0d0*L(1) - 1.0d0)
      ! Internal Node (10)
                  N_vals(10) = 27.0d0 * L(1) * L(2) * L(3)

            ! --- Step 2: Calculate the physical coordinates (x,y) of the integration point ---
                  x_int = 0.0d0
                  y_int = 0.0d0
                  do i = 1, 10
                        x_int = x_int + N_vals(i) * xx(1, nc(i,m))
                        y_int = y_int + N_vals(i) * xx(2, nc(i,m))
                  end do

            ! --- Step 3: Calculate the source term f value at (x_int, y_int) ---
            f_val_int = 2.0d0 * pi**2 * sin(pi * x_int) * sin(pi * y_int)
            !f_val_int = -4.0d0

            ! --- Step 4: Add contribution to each node's RHS value ---
            do i = 1, 10
                  ni = nc(i,m)
                ! This corresponds to the integral: integral( N_i * f * d_Omega )
                  rhs(ni) = rhs(ni) + f_val_int * N_vals(i) * jac(ig,m) * sg(ig)
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

	do i = 1, n_potential_nodes
            n = potential_nodes(i)
            f = potential_values(i)
		rhs(:) = rhs(:) - ea(:, n) * f
		ea(n, :) = 0.0d0
		ea(:, n) = 0.0d0
		ea(n, n) = 1.0d0
		rhs(n) = f
      end do

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
! -------------------------------------------------------------------------------
! Helper subroutine for 1D Gauss Quadrature points and weights
! -------------------------------------------------------------------------------
subroutine set_line_integration_points
      use indata
      implicit none
    
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
		do j = 1, 10
			do i = 1, 10
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
                  do i = 1, 10
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
	write(60,600) ( n, xx(1,n), xx(2,n), phi(n), n = 1, node )
	write(60,601) ( n, ue(n), ve(n), n = 1, nelm )
	close(60)
	600 format(i7, 3d15.6)
	601 format(i7,2d15.6)
end subroutine output

! -------------------------------------------------------------------------------
! Subroutine to output results in VTK legacy format for ParaView (for CUBIC Elements)
! CORRECTED FORMAT SPECIFIERS
! -------------------------------------------------------------------------------
subroutine output2
      use indata
      implicit none
      integer :: i, m, cell_list_size

      write(70, '(a)') '# vtk DataFile Version 3.0'
      write(70, '(a)') '2D FEM Cubic Analysis Results'
      write(70, '(a)') 'ASCII'
      write(70, '(a)') 'DATASET UNSTRUCTURED_GRID'
      write(70, *)

      write(70, '(a, i8, a)') 'POINTS', node, ' double'
      do i = 1, node
            write(70, '(3E16.8)') xx(1, i), xx(2, i), 0.0d0
      end do
      write(70, *)

      cell_list_size = nelm * 11
      write(70, '(a, 2i8)') 'CELLS', nelm, cell_list_size
      do m = 1, nelm
            write(70, '(11I8)') 10, (nc(i, m) - 1, i = 1, 10)
      end do
      write(70, *)

      write(70, '(a, i8)') 'CELL_TYPES', nelm
      do m = 1, nelm
            write(70, '(i3)') 69
      end do
      write(70, *)

      write(70, '(a, i8)') 'POINT_DATA', node
      write(70, '(a)') 'SCALARS potential double 1'
      write(70, '(a)') 'LOOKUP_TABLE default'
      do i = 1, node
            write(70, '(E16.8)') phi(i)
      end do
      write(70, *)

      write(70, '(a, i8)') 'CELL_DATA', nelm
      write(70, '(a)') 'VECTORS velocity double'
      do m = 1, nelm
            write(70, '(3E16.8)') ue(m), ve(m), 0.0d0
      end do

      close(70)

end subroutine output2

! -------------------------------------------------------------------------------
! -------------------------------------------------------------------------------
subroutine output_x05_nodes
    use indata
    implicit none
    integer :: n, count
    double precision :: tol

    tol = 1.0d-1

    do n = 1, node
        if (abs(xx(1,n) - 0.5d0) < tol) then
            write(80,'(i6,2f12.6)') n, xx(2,n), phi(n)
        end if
    end do
end subroutine output_x05_nodes

subroutine output_y05_nodes
    use indata
    implicit none
    integer :: n, count
    double precision :: tol

    tol = 1.0d-1

    do n = 1, node
        if (abs(xx(2,n) - 0.5d0) < tol) then
            write(90,'(i6,2f12.6)') n, xx(1,n), phi(n)
        end if
    end do
end subroutine output_y05_nodes