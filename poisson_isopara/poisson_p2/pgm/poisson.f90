
module indata
	implicit none
	integer :: node, nelm
	integer, allocatable :: nc(:,:)
	double precision, allocatable :: xx(:,:)
      
      ! --- boundary condition (Edge-based for Quadratic Elements)
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

	write(*,*) 'main program start! (Quadratic Elements)'
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
	character(20) :: bc_type_str

	mgp = 3

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
	allocate( xx(2,node), nc(6,nelm) )
	do i = 1, node
		read(10,*) n, xx(1,i), xx(2,i)
	end do
	do m = 1, nelm
		read(10,*) j, (nc(i,m), i = 1, 6)
	end do
	close(10)

      ! READ BOUNDARY CONDITION FILE
	read(11,*) n_potential_nodes
	      allocate(potential_nodes(n_potential_nodes), potential_values(n_potential_nodes))
	      do i = 1, n_potential_nodes
		      read(11,*) j, potential_nodes(i), potential_values(i)
	      end do
	close(11)

	! Dynamic Memory Allocation
	allocate( phi(node), rhs(node), ea(node,node) )
	allocate( xi_hat(mgp), eta_hat(mgp), sg(mgp) )
	allocate( jac(mgp,nelm) )
	allocate( shpx(6,mgp,nelm), shpy(6,mgp,nelm) )
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
     if (mgp == 3) then
            xi_hat(1)  = 0.1666666667d0
            eta_hat(1) = 0.1666666667d0
            sg(1)      = 0.1666666667d0

            xi_hat(2)  = 0.6666666667d0
            eta_hat(2) = 0.1666666667d0
            sg(2)      = 0.1666666667d0

            xi_hat(3)  = 0.1666666667d0
            eta_hat(3) = 0.6666666667d0
            sg(3)      = 0.1666666667d0
      elseif (mgp == 2) then
            xi_hat(1) = 0.5773502692d0
            sg(1)     = 1.0d0
            xi_hat(2) = 0.5773502692d0
            sg(2)     = 1.0d0      
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
!
! -------------------------------------------------------------------------------
!
subroutine makeSSF
	use indata
	implicit none
	integer :: m, ig, i
	integer, dimension(6) :: n_local
	double precision :: DetJ, tmp
	double precision, allocatable :: dxdxi(:,:), dxidx(:,:)
      double precision, allocatable :: sqd(:,:,:)
    
      double precision :: L(3)         ! (L1, L2, L3)
      double precision :: N_vals(6)    ! Shape function values N(L1,L2,L3)
      double precision :: dNdL(6,3)    ! dN/dL: Derivatives of N w.r.t. Area Coords
      double precision :: dLdxi(3,2)
      double precision :: dNdxi(6,2)   ! dN/d_xi, dN/d_eta: Final derivatives in reference coords
	
      dLdxi(1,1) = -1.0d0; dLdxi(1,2) = -1.0d0
      dLdxi(2,1) =  1.0d0; dLdxi(2,2) =  0.0d0
      dLdxi(3,1) =  0.0d0; dLdxi(3,2) =  1.0d0
      
      allocate( dxdxi(2,2), dxidx(2,2), sqd(2,6,mgp) )
      
	do m = 1, nelm
            n_local = nc(:,m)
		      do ig = 1, mgp
                        L(2) = xi_hat(ig)  ! L2 = xi
                        L(3) = eta_hat(ig)  ! L3 = eta
                        L(1) = 1.0d0 - L(2) - L(3)
                  ! Corner Nodes (Ni = Li * (2*Li - 1))
                        N_vals(1) = L(1) * (2.0d0 * L(1) - 1.0d0)
                        N_vals(2) = L(2) * (2.0d0 * L(2) - 1.0d0)
                        N_vals(3) = L(3) * (2.0d0 * L(3) - 1.0d0)
                  ! Midpoint Nodes
                        N_vals(4) = 4.0d0 * L(1) * L(2) ! N4 = 4*L1*L2
                        N_vals(5) = 4.0d0 * L(2) * L(3) ! N5 = 4*L2*L3
                        N_vals(6) = 4.0d0 * L(3) * L(1) ! N6 = 4*L3*L1
                  ! (Note: N_vals is not used later in here, but storing it is a clear and logical step.)

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
                  dxdxi(1,1) = dxdxi(1,1) + sqd(1,i,ig)*xx(1,n_local(i)) !dx/d_xi = dN/d_xi xx
                  dxdxi(1,2) = dxdxi(1,2) + sqd(2,i,ig)*xx(1,n_local(i))
                  dxdxi(2,1) = dxdxi(2,1) + sqd(1,i,ig)*xx(2,n_local(i))
                  dxdxi(2,2) = dxdxi(2,2) + sqd(2,i,ig)*xx(2,n_local(i))
            end do
            DetJ = dxdxi(1,1) * dxdxi(2,2) - dxdxi(1,2) * dxdxi(2,1) !calculate Jacobian for each element
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
      double precision :: s, w
      double precision :: N_1D(3)

      ! --- Variables for 1D Gauss Quadrature ---
      !integer, parameter :: mgp_line = 2 ! Number of integration points for the line
      !double precision :: s_hat(mgp_line), w_hat(mgp_line)
      !double precision :: s, w
      !double precision :: N_1D(3) ! Shape functions for the 3 nodes on the edge

	rhs(:) = 0.0d0
      !call set_line_integration_points(mgp_line, s_hat, w_hat)

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
!
subroutine add_source_term
    use indata
    implicit none
    integer :: m, i, ig, ni
    double precision :: f_val_int, x_int, y_int, pi
    double precision :: N_vals(6)
    double precision :: L(3)

    pi = acos(-1.0d0)

    ! Loop over all elements to compute the domain integral
    do m = 1, nelm
        do ig = 1, mgp
            ! --- Step 1: Calculate shape function VALUES at the integration point ---
            L(2) = xi_hat(ig)  ! L2 = xi
            L(3) = eta_hat(ig) ! L3 = eta
            L(1) = 1.0d0 - L(2) - L(3)

            ! Corner Nodes
            N_vals(1) = L(1) * (2.0d0 * L(1) - 1.0d0)
            N_vals(2) = L(2) * (2.0d0 * L(2) - 1.0d0)
            N_vals(3) = L(3) * (2.0d0 * L(3) - 1.0d0)
            ! Midside Nodes
            N_vals(4) = 4.0d0 * L(1) * L(2)
            N_vals(5) = 4.0d0 * L(2) * L(3)
            N_vals(6) = 4.0d0 * L(3) * L(1)

            ! --- Step 2: Calculate the physical coordinates (x,y) of the integration point ---
            x_int = 0.0d0
            y_int = 0.0d0
            do i = 1, 6
                x_int = x_int + N_vals(i) * xx(1, nc(i,m))
                y_int = y_int + N_vals(i) * xx(2, nc(i,m))
            end do

            ! --- Step 3: Calculate the source term f value at (x_int, y_int) ---
            f_val_int = -4.0d0
            !f_val_int = 2.0d0 * pi**2 * sin(pi * x_int) * sin(pi * y_int)

            ! --- Step 4: Add contribution to each node's RHS value ---
            do i = 1, 6
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

! -------------------------------------------------------------------------------
! Subroutine to output results in VTK legacy format for ParaView
! -------------------------------------------------------------------------------
subroutine output2
      use indata
      implicit none
      integer :: i, m, cell_list_size

      ! 1. VTKヘッダー
      write(70, '(a)') '# vtk DataFile Version 3.0'
      write(70, '(a)') '2D FEM Laplace Analysis Results'
      write(70, '(a)') 'ASCII'
      write(70, '(a)') 'DATASET UNSTRUCTURED_GRID'
      write(70, *)

      ! 2. 節点座標 (POINTS)
      !    2次元データを3次元として出力するため、Z座標には0.0を設定
      write(70, '(a, i8, a)') 'POINTS', node, ' double'
      do i = 1, node
            write(70, '(3e20.10)') xx(1, i), xx(2, i), 0.0d0
      end do
      write(70, *)

      ! 3. 要素の接続情報 (CELLS)
      !    各要素は6節点なので、リストのサイズは nelm * (1+6) となる
      cell_list_size = nelm * 7
      write(70, '(a, 2i8)') 'CELLS', nelm, cell_list_size
      do m = 1, nelm
            ! VTKは0ベースのインデックスなので、-1 する
            write(70, '(i2, 6i8)') 6, (nc(i, m) - 1, i = 1, 6)
      end do
      write(70, *)

      ! 4. 要素のタイプ (CELL_TYPES)
      !    6節点二次三角形要素のVTKタイプは '22'
      write(70, '(a, i8)') 'CELL_TYPES', nelm
      do m = 1, nelm
            write(70, '(i3)') 22
      end do
      write(70, *)

      ! 5. 節点データ (POINT_DATA) - スカラーポテンシャル
      write(70, '(a, i8)') 'POINT_DATA', node
      write(70, '(a)') 'SCALARS potential double 1'
      write(70, '(a)') 'LOOKUP_TABLE default'
      do i = 1, node
            write(70, '(e20.10)') phi(i)
      end do
      write(70, *)

      ! 6. 要素データ (CELL_DATA) - 速度ベクトル
      write(70, '(a, i8)') 'CELL_DATA', nelm
      write(70, '(a)') 'VECTORS velocity double'
      do m = 1, nelm
            ! 2次元ベクトルを3成分で書き出す（Z成分は0）
            write(70, '(3e20.10)') ue(m), ve(m), 0.0d0
      end do

      close(70)

end subroutine output2

subroutine output_x05_nodes
    use indata
    implicit none
    integer :: n, count
    double precision :: tol

    tol = 1.0d-3   

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

    tol = 1.0d-3   

    do n = 1, node
        if (abs(xx(2,n) - 0.5d0) < tol) then
            write(90,'(i6,2f12.6)') n, xx(1,n), phi(n)
        end if
    end do
end subroutine output_y05_nodes

