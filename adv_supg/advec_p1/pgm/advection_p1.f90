module dimdata
	implicit none
	double precision, allocatable :: phi(:)
	double precision, allocatable :: rhs(:)
	double precision, allocatable :: ue(:), ve(:)
      integer :: ista, iend, istep
      integer :: kx, kuvp
      double precision :: dt, dti, th, th2
      double precision :: eps, time
      character(100) :: vtkfile_base 
! --- Mesh 
      integer :: node, nelm
      integer :: nx, mx
      integer, allocatable :: nc(:,:)
      double precision, allocatable :: xx(:,:) 
! --- Boundary Condition
      integer :: ifbc
      integer, allocatable :: nfbc(:)
      double precision, allocatable :: ffbc(:)
! --- Isoparametric Element Data ---
      integer :: mgp 
      double precision, allocatable :: shp(:,:,:), shpx(:,:,:), shpy(:,:,:)
      double precision, allocatable :: jac(:,:)
      double precision, allocatable :: xi_hat(:), eta_hat(:), sg(:)
! --- tau
      double precision, allocatable :: tau(:,:)
! --- velocity, phi
      double precision, allocatable :: us(:), vs(:), ff(:)
! --- For Matrix
      double precision, allocatable :: ea(:,:,:), ad(:), buvp(:)
! --- Step phi
      double precision, allocatable :: duvp(:)
      double precision, allocatable :: aa(:)
! --- Initial Mass
      double precision :: initial_mass, total_mass, conservation_rate, calculate_total_mass
      external calculate_total_mass
!
end module dimdata

! --- Main program ---
!
program advection

      use dimdata
      implicit none
      integer :: n, i, k
      double precision :: epsw, sum1, sum2
      double precision :: begin_time, end_time

      call CPU_TIME(begin_time)
!
      call datain 
      write(6,*)'Finish datain'

      if (ista == 1) then
            write(6,*) 'Writing initial condition (step 0)...'
            istep = 0
            call output
            call write_vtk_seq
            write(6,*) 'Finished writing initial condition.'
      end if

      call set_integration_points 
      write(6,*)'Finish set_integration_points'

      call makeSSF
      write(6,*)'Finish makeSSF'

      initial_mass = calculate_total_mass()
      write(6,'(A, E20.12)') 'Calculated Initial Mass (Integral):', initial_mass

      allocate(tau(mgp,nelm))
      allocate(ad(node), ea(3,3,nelm), buvp(node))

      call maktau 
      write(6,*)'Finish maktau'
      ! stop

      do 1000 istep = ista, iend

            time = dble(istep) * dt

            call makmat
            !write(6,*)'Finish makmat'

            call mkbuvp

            call dscale ! diagonal scaling

            call bicgst(eps, kx, node, ad, kuvp, duvp, buvp, nelm, nc, ea, ifbc, nfbc) ! bi-conjugate gradient solver

            do n = 1, node
                  ff(n) = duvp(n)
            end do

            call mass_conservation

            if (istep == 1000) then
                  call output_100
            end if

            call output
            call write_vtk_seq

            write(6,*), 'Finish timestep: ', istep
            
      1000 continue

      call CPU_TIME(end_time)

      write(*,'(a,f15.6,a)') 'Time of operation was ' , end_time - begin_time, ' seconds'
	write(*,*) 'main program end!'

end program Advection

subroutine datain
      use dimdata
      implicit none
      integer :: i, j, k, n, m
      character(100) :: inpfile, mesfile, uvsfile, bdcfile, finfile, resfile, outfile100, mass_conservation
      
      open(9, file='file.txt', status='unknown')
      read(9,'(a)') inpfile
      read(9,'(a)') mesfile
      read(9,'(a)') uvsfile
      read(9,'(a)') bdcfile
      read(9,'(a)') finfile
      read(9,'(a)') resfile
      read(9,'(a)') vtkfile_base
      read(9,'(a)') outfile100
      read(9,'(a)') mass_conservation
!
      open(10, file = inpfile, status = 'unknown')
      open(11, file = mesfile, status = 'unknown')
      open(12, file = uvsfile, status = 'unknown')
      open(13, file = bdcfile, status = 'unknown')
      open(14, file = finfile, status = 'unknown')
      open(50, file = resfile, status = 'replace')
      open(100, file = outfile100, status = 'replace')
      open(60, file = mass_conservation, status = 'replace')
!----------------------------------------------------------------
      write(6,608)
608  format(/,'*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*',/)

      write(6,*)' [inpfile] ; ', inpfile
      write(6,*)' [mesfile] ; ', mesfile
      write(6,*)' [uvsfile] ; ', uvsfile
      write(6,*)' [bdcfile] ; ', bdcfile
      write(6,607)
607  format(/,'*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*',/)
!----------------------------------------------------------------
!
      read(10,*) ista, iend
      read(10,*) dt
      read(10,*) eps, kx
!
      dti = 1.0d0 / dt
      th  = 0.5d0
      mgp = 3
!
! --- Mesh Data
      read(11,*) node, nelm
      allocate( xx(2,node), nc(3,nelm) )
      do i = 1, node
            read(11,*) n,  xx(1,i),  xx(2,i)
      end do
      
      do m = 1, nelm
            read(11,*) j, (nc(i,m), i = 1, 3)
      end do
      close(11)
!
! --- Advection Velocity Data
      allocate( us(node), vs(node) )
      do i = 1, node         
            read(12,*) j, us(i), vs(i)
      enddo
      close(12)

! --- Boundary Condition Data
      read(13,*) ifbc
      allocate( nfbc(ifbc), ffbc(ifbc) )
      do i = 1, ifbc
            read(13,*) j, nfbc(i), ffbc(i)
      enddo
      close(13)
      allocate( duvp(node), ff(node) )
!
! --- Initial Condition Data
      if( ista == -1 ) then
            ista = 1
!
            do j = 1, node
                  read(14,*) i, ff(i)
            enddo
            close(14)
!
	      do n = 1, node
	            duvp(n) = ff(n)   ! Initial condition for solver
            end do
      end if
!
      do j = 1, ifbc
            ff(nfbc(j)) = ffbc(j)
      end do

      ! initial_mass = 0.0d0
      ! do n = 1, node
      !       initial_mass = initial_mass + ff(n)
      ! end do
      
      write(6,600) ista, iend, dt, eps, kx
600 format(' istart  :',i12,'    iend :',i12,/,                 &                          
             ' delta t :',f12.5,/,                              &
             ' eps     :',d12.5,'    kx   :',i12)
!
      write(6,601) node, nelm
 601  format(' node      :',i12,'    nelm   :',i12,/                &
           /,'*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*')
!----------------------------------------------------------------
!
      allocate( xi_hat(mgp), eta_hat(mgp), sg(mgp) )
      allocate( jac(mgp,nelm) )
      allocate( shp(3,mgp,nelm), shpx(3,mgp,nelm), shpy(3,mgp,nelm) ) !set allocation in end as possible as u can
!
end subroutine datain

subroutine set_integration_points
      use dimdata
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
!
subroutine makeSSF
	use dimdata
	implicit none
	integer :: m, ig, i, j, k 
	integer, dimension(3) :: n_local
	double precision :: DetJ, tmp
	double precision, allocatable :: dxdxi(:,:), dxidx(:,:)
      double precision, allocatable :: sqd(:,:,:)
      double precision, allocatable :: L(:), N_vals(:), dNdL(:,:), dLdxi(:,:), dNdxi(:,:)
      
      allocate(dxdxi(2,2), dxidx(2,2), sqd(2,3,mgp), L(3), N_vals(3), dNdL(3,3), dLdxi(3,2), dNdxi(3,2))

      dLdxi(1,1) = -1.0d0; dLdxi(1,2) = -1.0d0
      dLdxi(2,1) =  1.0d0; dLdxi(2,2) =  0.0d0
      dLdxi(3,1) =  0.0d0; dLdxi(3,2) =  1.0d0

	do m = 1, nelm
            n_local = nc(:,m)
		do ig = 1, mgp
                  ! Step 1: Define Area Coordinates (L1, L2, L3) from reference coordinates
                  L(2) = xi_hat(ig)  ! L2 = xi
                  L(3) = eta_hat(ig)  ! L3 = eta
                  L(1) = 1.0d0 - L(2) - L(3)

                  N_vals(1) = L(1) 
                  N_vals(2) = L(2) 
                  N_vals(3) = L(3)

                  shp(:,ig,m) = N_vals

                  ! do i = 1, 3
                  !     write(6,'(A,I1,A,F20.10)') 'N_vals(', i, ') =', N_vals(i)
                  ! end do

                  ! Step 2: Calculate derivatives of shape functions w.r.t. Area Coordinates (dN/dL)
                  dNdL(1,1) = 1.0d0; dNdL(1,2) = 0.0d0;  dNdL(1,3) = 0.0d0
                  dNdL(2,1) = 0.0d0; dNdL(2,2) = 1.0d0;  dNdL(2,3) = 0.0d0
                  dNdL(3,1) = 0.0d0; dNdL(3,2) = 0.0d0;  dNdL(3,3) = 1.0d0

                  ! Step 3: Apply Chain Rule using matrix multiplication to get dN/d(xi,eta)
                  dNdxi = matmul(dNdL, dLdxi)

                  ! do i = 1, 3
                  !     do j = 1, 2
                  !         write(6,'(A,I1,A,I1,A,F20.10)') 'dNdxi(', i, ',', j, ') =', dNdxi(i,j)
                  !     end do
                  ! end do
                  
                  ! Store results into the sqd array for subsequent calculations
                  sqd(1, :, ig) = dNdxi(:, 1) ! dN/d_xi
                  sqd(2, :, ig) = dNdxi(:, 2) ! dN/d_eta

                  dxdxi(:,:) = 0.0d0
                  
                  do i = 1, 3
                        dxdxi(1,1) = dxdxi(1,1) + sqd(1,i,ig)*xx(1,n_local(i))
                        dxdxi(1,2) = dxdxi(1,2) + sqd(2,i,ig)*xx(1,n_local(i))
                        dxdxi(2,1) = dxdxi(2,1) + sqd(1,i,ig)*xx(2,n_local(i))
                        dxdxi(2,2) = dxdxi(2,2) + sqd(2,i,ig)*xx(2,n_local(i))
                  end do

                  ! do i = 1, 2
                  !     do j = 1, 2
                  !         write(6,'(A,I1,A,I1,A,F20.10)') 'dxdxi(', i, ',', j, ') =', dxdxi(i,j)
                  !     end do
                  ! end do

                  DetJ = dxdxi(1,1) * dxdxi(2,2) - dxdxi(1,2) * dxdxi(2,1)
                  ! if (m == 1) then
                  !     write(6,*) 'makeSSF: m=1, ig=', ig, ' >> DetJ calculated:', DetJ
                  !     call flush(6)
                  ! end if
                  
                  ! if (mod(m,400) == 0) then
                  !       write(6, '(A, I6, A, I3, A, E22.14)') 'Element, IGP, DetJ: ', m, ',', ig, ' , ', DetJ
                  ! end if

                  if ( DetJ <= 1.0d-12 ) then
			      write(*,*) 'error: Jacobian is zero or negative at element, igp:', m, ig, DetJ; stop
                        call flush(6)
		      else
			      tmp = 1.0d0 / DetJ
                        jac(ig,m) = DetJ
		      end if

                  dxidx(1,1) =  dxdxi(2,2) * tmp; dxidx(1,2) = -dxdxi(1,2) * tmp
			dxidx(2,1) = -dxdxi(2,1) * tmp; dxidx(2,2) =  dxdxi(1,1) * tmp
			
                  ! do i = 1, 2
                  !       do j = 1, 2
                  !       write(6,'(A,I1,A,I1,A,F20.10)') 'dxidx(', i, ',', j, ') =', dxidx(i,j)
                  !     end do
                  ! end do
                  
                  do i = 1, 3
			      shpx(i,ig,m) = sqd(1,i,ig) * dxidx(1,1) + sqd(2,i,ig) * dxidx(2,1)
			      shpy(i,ig,m) = sqd(1,i,ig) * dxidx(1,2) + sqd(2,i,ig) * dxidx(2,2)
                        !write(6,*) 'shpx(i,ig,m) =', shpx(i,ig,m)
                        !write(6,*) 'shpy(i,ig,m) =', shpy(i,ig,m)
		      end do
		end do
	end do

      deallocate(dxdxi, dxidx, sqd)

end subroutine makeSSF
!
subroutine maktau       !Stabilization Parameter
      use dimdata
      implicit none
      integer, dimension(3) :: n_local
      double precision :: u, v, u_norm
      double precision :: he, sum_abs_G, dot_product_val
      double precision :: dt05, t1, t2
      integer :: m, ig, i
!
      do m = 1, nelm
            n_local = nc(:,m)

            u = 0.0d0
            v = 0.0d0
            do i = 1, 3
                  u = u + us(n_local(i))
                  v = v + vs(n_local(i))
            end do
            u = u / 3.0d0
            v = v / 3.0d0
            
            u_norm = dsqrt(u*u + v*v)
            
            ! Initialize tau for the element
            tau(:,m) = 0.0d0

            do ig = 1, mgp
                  ! --- Calculate h_e at this integration point (ig) ---
                  sum_abs_G = 0.0d0
                  do i = 1, 3
                        ! dot_product_val = u_e * (dN_i/dx) + v_e * (dN_i/dy)
                        dot_product_val = u * shpx(i,ig,m) + v * shpy(i,ig,m)
                        sum_abs_G = sum_abs_G + dabs(dot_product_val)
                  end do
                  he = (2.0d0 * u_norm) / sum_abs_G
                                    
                  ! --- Calculate tau at this integration point (ig) ---
                  dt05 = 2.0d0 / dt
                  t1 = dt05**2

                  t2 = (2.0d0 * u_norm / he)**2

                  ! tau(ig,m) = 1.0d0 / dsqrt(t1 + t2)
                  tau(ig,m) = 1.0d0 / dsqrt(t1)

                  ! tau(ig,m) = 0.1d0

                  if (m <= 10) then 
                        write(6, '(A, I6, A, I2)') 'DEBUG maktau: Elem=', m, ' IGP=', ig
                        write(6, '(A, 2E15.6)')   '  Avg Vel (u, v) =', u, v
                        write(6, '(A, E15.6)')    '  |u| (u_norm)   =', u_norm
                        write(6, '(A, E15.6)')    '  sum_abs_G      =', sum_abs_G
                        write(6, '(A, E15.6)')    '  he             =', he
                        write(6, '(A, E15.6)')    '  t1 ( (2/dt)^2 )  =', t1
                        write(6, '(A, E15.6)')    '  t2 ( (2|u|/he)^2 )=', t2
                        write(6, '(A, E15.6)')    '  tau(ig,m)      =', tau(ig,m)
                        call flush(6)
                  end if
            end do
      end do
end subroutine maktau
!
subroutine makmat
      use dimdata
      implicit none
      integer :: m, i, j, ig
      integer, dimension(3) :: n_local
      double precision :: u_int, v_int, wj, t, adv_i, adv_j
      double precision, dimension(3,3) :: Me, Ke, Ke_supg, Me_supg, E_elem
      ad(:) = 0.0d0
      ea(:, :, :) = 0.0d0

      do m = 1, nelm
            n_local = nc(:,m)
            Me = 0.0d0; Ke = 0.0d0; Me_supg = 0.0d0; Ke_supg = 0.0d0; u_int = 0.0d0; v_int = 0.0d0
            
            do ig = 1, mgp      
                  wj = jac(ig,m) * sg(ig)
                  
                  u_int = dot_product(shp(:,ig,m), us(n_local))
                  v_int = dot_product(shp(:,ig,m), vs(n_local))
                  
                  do j = 1, 3
                        do i = 1, 3
                              ! Advection term components (u . grad(N))
                              adv_i = u_int * shpx(i,ig,m) + v_int * shpy(i,ig,m)
                              adv_j = u_int * shpx(j,ig,m) + v_int * shpy(j,ig,m)
                              
                              ! 1. Standard Mass Matrix
                              Me(i,j) = Me(i,j) + shp(i,ig,m) * shp(j,ig,m) * wj 
                              
                              ! 2. Standard Advection Matrix
                              Ke(i,j) = Ke(i,j) + shp(i,ig,m) * &
                                                (u_int * shpx(j,ig,m) + v_int * shpy(j,ig,m)) * wj 
                              ! 3. SUPG Mass Matrix 
                              Me_supg(i,j) = Me_supg(i,j) + adv_i * tau(ig,m) * shp(j,ig,m) * wj
                              
                              ! 4. SUPG Advection (Stabilization) Matrix
                              Ke_supg(i,j) = Ke_supg(i,j) + adv_i * tau(ig,m) * adv_j * wj
                        end do
                  end do
            end do

            E_elem = dti * (Me + Me_supg) + th * (Ke + Ke_supg)

            do j = 1, 3
                  do i = 1, 3
                        ea(i, j, m) = E_elem(i,j)
                  end do
                  ad(n_local(j)) = ad(n_local(j)) + E_elem(j,j) 
            end do
      end do

end subroutine makmat
!
subroutine mkbuvp ! Right-hand side vector for Crank-Nicolson
      use dimdata
      implicit none
      integer :: m, i, j, ig
      integer, dimension(3) :: n_local
      double precision :: u_int, v_int, wj, t, th0, adv_i, adv_j
      
      ! Element-level matrices and vectors
      double precision, dimension(3,3) :: Me, Ke, Me_supg, Ke_supg
      double precision, dimension(3) :: phi_e_n, be_M, be_K, be

      buvp(:) = 0.0d0
      th0 = 1.0d0 - th

      do m = 1, nelm
            
            n_local = nc(:,m)
            phi_e_n = ff(n_local) ! Get previous time step solution for the element

            ! Initialize element matrices (must be recalculated here for the RHS)
            Me = 0.0d0; Ke = 0.0d0; Me_supg = 0.0d0; Ke_supg = 0.0d0
            
            do ig = 1, mgp
                  wj = jac(ig,m) * sg(ig)
                  u_int = dot_product(shp(:,ig,m), us(n_local))
                  v_int = dot_product(shp(:,ig,m), vs(n_local))
                  
                  do j = 1, 3
                        do i = 1, 3
                              adv_i = u_int * shpx(i,ig,m) + v_int * shpy(i,ig,m)
                              adv_j = u_int * shpx(j,ig,m) + v_int * shpy(j,ig,m)
                              
                              Me(i,j) = Me(i,j) + shp(i,ig,m) * shp(j,ig,m) * wj
                              Ke(i,j) = Ke(i,j) + shp(i,ig,m) * adv_j * wj
                              Me_supg(i,j) = Me_supg(i,j) + adv_i * tau(ig,m) * shp(j,ig,m) * wj
                              Ke_supg(i,j) = Ke_supg(i,j) + adv_i * tau(ig,m) * adv_j * wj
                        end do
                  end do
            end do
      
            ! Calculate element RHS vector contribution
            ! be = (dti * (Me + Me_supg) - th0 * (Ke + Ke_supg)) * phi_e_n
            be_M = matmul(Me + Me_supg, phi_e_n)
            be_K = matmul(Ke + Ke_supg, phi_e_n)
            be = dti * be_M - th0 * be_K
            
            ! Assemble into global RHS vector
            do i = 1, 3
                  buvp(n_local(i)) = buvp(n_local(i)) + be(i)
            end do
      end do
end subroutine mkbuvp

subroutine dscale
      use dimdata
      implicit none
      integer :: n, m, i, j
      double precision, allocatable :: d(:)
      nx = node
      mx = nelm
      allocate( d(3) )

      do n = 1, nx
            ad(n) = 1.0d0 / dsqrt(ad(n))
      end do

      do m = 1, mx
            do i = 1, 3
                  d(i) = ad(nc(i,m))
            end do
            do i = 1, 3
                  do j = 1, 3
                        ea(i, j, m) = ea(i, j, m) * d(i) * d(j)
                  end do
            end do
      enddo

end subroutine dscale
!
subroutine bicgst(eps, kx, nx, ad, k, xx, bb, mx, nc, ea, ifbc, nfbc)
      
      implicit none
      integer, intent(in) :: kx, nx, mx, ifbc
      integer, intent(in) :: nfbc(ifbc), nc(3,mx)
      double precision, intent(inout) :: ad(nx), ea(3,3,mx)
      double precision, intent(inout) ::xx(nx), bb(nx)
      double precision, intent(in) :: eps
      integer, intent(out) :: k
!
      double precision, allocatable :: rr(:), pp(:), tt(:)
      double precision, allocatable :: rs(:), Ap(:), At(:) 
      double precision, allocatable :: p(:)
      double precision :: brb, rtr
      double precision :: b, r
      double precision :: rsr, btb
      double precision :: epsbtb
      double precision :: rAp
      double precision :: rsr0, alph
      double precision :: Att, At2
      double precision :: zeta, beta
      integer :: i, j, n, m, ib

      allocate( rr(nx), pp(nx), tt(nx), rs(nx), Ap(nx), At(nx), p(3) )
!
      do n = 1, nx
            bb(n) = bb(n) * ad(n)
            xx(n) = xx(n) / ad(n)
      end do
      k = 0
      do n = 1, nx
            rr(n) = bb(n)
      end do
      do ib = 1, ifbc
            rr(nfbc(ib)) = 0.0d0
      end do
      do n = 1, nx
            Ap(n) = 0.0d0
      end do
      do m = 1, mx
            do i = 1, 3
                  p(i)   = xx(nc(i,m))
            end do
            do i = 1, 3
                  do j = 1, 3
                        Ap(nc(i,m)) = Ap(nc(i,m)) + ea(i, j, m) * p(j)
                  end do
            end do
      end do
      do ib = 1, ifbc
            Ap(nfbc(ib)) = 0.0d0
      end do
      btb = 0.0d0
      rtr = 0.0d0
      do n = 1, nx
            b = rr(n)
            r = b - Ap(n)
            rr(n) = r
            rs(n) = r
            pp(n) = r
            btb = btb + b * b
            rtr = rtr + r * r
      end do
      rsr = rtr
! check    btb = rtr
      epsbtb = eps * btb
!--------------------------------------------
      if(btb .lt. 1.d-30)  goto 1000
      if(rtr .lt. epsbtb)  goto 1000
!--------------------------------------------
      do  k = 1, kx
            do  n = 1, nx
                  Ap(n) = 0.0d0
            enddo
            do  m = 1, mx
                  do i = 1, 3
                        p(i) = pp(nc(i,m))
                  end do
                  do i = 1, 3
                        do j = 1, 3
                              Ap(nc(i,m)) = Ap(nc(i,m)) + ea(i, j, m) * p(j)
                        end do
                  end do
            enddo
            do  ib = 1, ifbc
            Ap(nfbc(ib)) = 0.0d0
            enddo
            rAp = 0.0d0
            do n = 1, nx
                  rAp = rAp + rs(n) * Ap(n)
            enddo
            rsr0 = rsr
            alph = rsr0 / rAp

            ! if (k <= 5) then
            !       write(6,*) '  -> k=', k, ' rsr0=', rsr0, ' rAp=', rAp, ' alph=', alph 
            ! endif

            do n = 1, nx
                  tt(n) = rr(n) - alph * Ap(n)
            enddo
            do n = 1, nx
                  At(n) = 0.0d0
            enddo
            do m = 1, mx
                  do i = 1, 3
                        p(i) = tt(nc(i,m))
                  end do
                  do i = 1, 3
                        do j = 1, 3
                              At(nc(i,m)) = At(nc(i,m)) + ea(i, j, m) * p(j)
                        end do
                  end do
            enddo
            do ib = 1, ifbc
                  At(nfbc(ib)) = 0.0d0
            enddo
            Att = 0.0d0
            At2 = 0.0d0
            do n = 1, nx
                  Att = Att + At(n) * tt(n)
                  At2 = At2 + At(n) * At(n)
            enddo
            zeta = Att / At2
            ! if (k <= 5) then
            !       write(6,*) ' ', 'Att=', Att, ' At2=', At2, ' zeta=', zeta 
            ! endif
            rsr = 0.0d0
            rtr = 0.0d0
            do n = 1, nx
                  xx(n) = xx(n) + alph * pp(n) + zeta * tt(n)
                  rr(n) = tt(n) - zeta * At(n)
                  rsr = rsr + rs(n) * rr(n)
                  rtr = rtr + rr(n) * rr(n)
            enddo
!--------------------------------------------
      if(rtr .lt. epsbtb)  goto 1000
!--------------------------------------------
            beta = alph / zeta * rsr / rsr0
!
            do n = 1, nx
                  pp(n) = rr(n) + beta * (pp(n) - zeta * Ap(n))
            enddo
      enddo
      write(6,*)   'rtr =',  rtr
      write(6,*)   'k > kx'
      stop
1000 continue
      do n = 1, nx
            xx(n) = xx(n) * ad(n)
            bb(n) = bb(n) / ad(n)
      end do
      deallocate( rr, pp, tt, rs, Ap, At, p )
end subroutine bicgst

subroutine output
      use dimdata
      implicit none
      integer :: n, i, j, k
!
      write(50,*) istep, dt
      do n = 1, node
            write(50,*) n, ff(n)
      end do
end subroutine output

subroutine output_100
      use dimdata
      implicit none
      integer :: n, i, j
      integer :: count1, count2
      real(8), allocatable :: data1(:,:), data2(:,:)
      real(8) :: temp_coord, temp_val

      count1 = 0
      do n = 1, node
            if (abs(xx(2, n)-0.5d0) < 1.0d-9) then
                  count1 = count1 + 1
            end if
      end do

      if (count1 > 0) then
            allocate(data1(count1, 2))
            i = 0
            do n = 1, node
                  if (abs(xx(2, n)-0.5d0) < 1.0d-9) then
                        i = i + 1
                        data1(i, 1) = xx(1, n)  
                        data1(i, 2) = ff(n)    
                  end if
            end do

            do i = 1, count1 - 1
                  do j = i + 1, count1
                        if (data1(i, 1) > data1(j, 1)) then
                              temp_coord  = data1(i, 1)
                              temp_val    = data1(i, 2)
                              data1(i, 1) = data1(j, 1)
                              data1(i, 2) = data1(j, 2)
                              data1(j, 1) = temp_coord
                              data1(j, 2) = temp_val
                        end if
                  end do
            end do

            do i = 1, count1
                  write(100,'(2(1x, G20.9))') data1(i, 1), data1(i, 2)
            end do
            deallocate(data1)
      end if

end subroutine output_100


subroutine write_vtk(filename, xx, ff, nc, node, nelm, us, vs)
      implicit none
      character(len=*), intent(in) :: filename
      double precision, intent(in) :: xx(2, node)
      double precision, intent(in) :: ff(node)
      integer, intent(in) :: nc(3, nelm)
      integer, intent(in) :: node, nelm
      double precision, intent(in) :: us(node)
      double precision, intent(in) :: vs(node)
  
      integer :: i, m
      integer :: vtk_unit

      open(newunit=vtk_unit, file=filename, status='replace')

      write(vtk_unit,'(A)') '# vtk DataFile Version 3.0'
      write(vtk_unit,'(A)') 'Advection result'
      write(vtk_unit,'(A)') 'ASCII'
      write(vtk_unit,'(A)') 'DATASET UNSTRUCTURED_GRID'

      !
      write(vtk_unit,'(A, 1X, I0, 1X, A)') 'POINTS', node, 'double'
      do i = 1, node
            write(vtk_unit,'(3F16.8)') xx(1,i), xx(2,i), 0.0d0
      end do

      !
      write(vtk_unit,'(A, 1X, I0, 1X, I0)') 'CELLS', nelm, nelm*4
      do m = 1, nelm
            write(vtk_unit,'(4I8)') 3, nc(1,m)-1, nc(2,m)-1, nc(3,m)-1
      end do

      !
      write(vtk_unit,'(A, 1X, I0)') 'CELL_TYPES', nelm
      do m = 1, nelm
            write(vtk_unit,'(I8)') 5
      end do

      !
      write(vtk_unit,'(A, 1X, I0)') 'POINT_DATA', node
      write(vtk_unit,'(A)') 'SCALARS phi double 1'
      write(vtk_unit,'(A)') 'LOOKUP_TABLE default'
      do i = 1, node
            write(vtk_unit,'(F18.9)') ff(i)
      end do
      
      write(vtk_unit,'(A)') 'VECTORS velocity double'   
      do i = 1, node
            write(vtk_unit,'(3F18.9)') us(i), vs(i), 0.0d0  
      end do
    close(vtk_unit)
end subroutine write_vtk

subroutine write_vtk_seq
      use dimdata
      implicit none
      character(len=200) :: filename

      write(filename, '(A, I4.4, ".vtk")') trim(vtkfile_base), istep
      call write_vtk(filename, xx, ff, nc, node, nelm, us, vs)
end subroutine write_vtk_seq

subroutine mass_conservation
      use dimdata
      implicit none
      integer :: m, ig, i

      ! total_mass = 0.0d0
      ! do i = 1, node
      !       total_mass = total_mass + ff(i)
      ! end do
      total_mass = calculate_total_mass()

      conservation_rate = total_mass / initial_mass * 100.0d0

      write(60,*) istep, conservation_rate

end subroutine mass_conservation

function calculate_total_mass() result(mass)
      use dimdata
      implicit none
      double precision :: mass
      
      integer :: m, ig, i, j
      integer, dimension(3) :: n_local
      double precision, dimension(3,3) :: Me
      double precision, dimension(3) :: phi_e
      double precision, dimension(3) :: be_M
      double precision :: wj
      
      mass = 0.0d0

      do m = 1, nelm

            n_local = nc(:,m)
            phi_e = ff(n_local) 
            
            Me = 0.0d0
            
            do ig = 1, mgp
                  wj = jac(ig,m) * sg(ig)
                  do j = 1, 3
                        do i = 1, 3
                              ! Me(i,j) = Me(i,j) + shp(i,ig,m) * shp(j,ig,m) * wj  [cite: 122, 164]
                              Me(i,j) = Me(i,j) + shp(i,ig,m) * shp(j,ig,m) * wj
                        end do
                  end do
            end do

            be_M = matmul(Me, phi_e)
            
            mass = mass + sum(be_M)
            
      end do
      
end function calculate_total_mass
