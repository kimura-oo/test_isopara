module dimdata
	implicit none
	double precision, allocatable :: phi(:)
      integer :: ista, iend, iout, istep
      integer :: kx, kuvp, alpha
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
!
end module dimdata

program main
      use dimdata

      call CPU_TIME(begin_time)
      write(6,*)'Finish datain'

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

      allocate(tau(mgp,nelm))
      allocate(ad(node), ea(3,3,nelm), buvp(node))

      call maketau
      write(6,*)'Finish maketau'
      
      do 1000 istep = ista, iend

            time = dble(istep) * dt

            call makemat
            write(6,*)'Finish makmat'

            call makebuvp
            write(6,*)'Finish makebuvp'

            call dscale
            write(6,*)'Finish dscale'

            do n = 1, node
                  ff(n) = duvp(n)
            end do

            call bicgst(eps, kx, node, ad, kuvp, duvp, buvp, nelm, nc, ea, ifbc, nfbc)
            write(6,*)'Finish bicgst'

            call output
            
            call write_vtk_seq

            write(6,*), 'Finish timestep: ', istep

      1000 continue

      call CPU_TIME(end_time)
      
      write(*,'(a,f15.6,a)') 'Time of operation was ' , end_time - begin_time, ' seconds'
	write(*,*) 'main program end!'

end program main

subroutine datain
      use dimdata
      implicit none
      integer :: i, j, k, n, m
      character(100) :: inpfile, mesfile, uvsfile, bdcfile, finfile, resfile
      !character(100) :: vtkfile_base
      
      open(9, file='file.txt', status='unknown')
      read(9,'(a)') inpfile
      read(9,'(a)') mesfile
      read(9,'(a)') uvsfile
      read(9,'(a)') bdcfile
      read(9,'(a)') finfile
      read(9,'(a)') resfile
      read(9,'(a)') vtkfile_base

      open(10, file = inpfile, status = 'unknown')
      open(11, file = mesfile, status = 'unknown')
      open(12, file = uvsfile, status = 'unknown')
      open(13, file = bdcfile, status = 'unknown')
      open(14, file = finfile, status = 'unknown')
      open(50, file = resfile, status = 'replace')
!
!----------------------------------------------------------------
      write(6,608)
608  format(/,'*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*',/)

      write(6,*)' [inpfile] ; ', inpfile
      write(6,*)' [mesfile] ; ', mesfile
      write(6,*)' [uvsfile] ; ', uvsfile
      write(6,*)' [bdcfile] ; ', bdcfile
      write(6,*)' [finfile] ; ', finfile

      write(6,607)
607  format(/,'*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*',/)
!----------------------------------------------------------------
!
      read(10,*) ista, iend
      read(10,*) dt
      read(10,*) eps, kx
      read(10,*) alpha
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
      
      write(6,600) ista, iend, iout, dt, eps, kx
600 format(' istart  :',i12,'    iend :',i12,/,                 &
             ' iout    :',i12,/,                                &
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
      allocate( shp(3,mgp,nelm), shpx(3,mgp,nelm), shpy(3,mgp,nelm) ) 
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
            sg(1)      = 1.0d0/40.0d0 

            xi_hat(2)  = 0.5d0
            eta_hat(2) = 0.0d0
            sg(2)      = 1.0d0/15.0d0  

            xi_hat(3)  = 1.0d0
            eta_hat(3) = 0.0d0
            sg(3)      = 1.0d0/40.0d0  

            xi_hat(4)  = 0.5d0
            eta_hat(4) = 0.5d0
            sg(4)      = 1.0d0/15.0d0  

            xi_hat(5)  = 0.0d0
            eta_hat(5) = 1.0d0
            sg(5)      = 1.0d0/40.0d0  

            xi_hat(6)  = 0.0d0
            eta_hat(6) = 0.5d0
            sg(6)      = 1.0d0/15.0d0  

            xi_hat(7)  = 1.0d0/3.0d0
            eta_hat(7) = 1.0d0/3.0d0
            sg(7)      = 9.0d0/40.0d0  
      end if
end subroutine set_integration_points

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

                  ! Step 2: Calculate derivatives of shape functions w.r.t. Area Coordinates (dN/dL)
                  dNdL(1,1) = 1.0d0; dNdL(1,2) = 0.0d0;  dNdL(1,3) = 0.0d0
                  dNdL(2,1) = 0.0d0; dNdL(2,2) = 1.0d0;  dNdL(2,3) = 0.0d0
                  dNdL(3,1) = 0.0d0; dNdL(3,2) = 0.0d0;  dNdL(3,3) = 1.0d0

                  ! Step 3: Apply Chain Rule using matrix multiplication to get dN/d(xi,eta)
                  dNdxi = matmul(dNdL, dLdxi)

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

                  DetJ = dxdxi(1,1) * dxdxi(2,2) - dxdxi(1,2) * dxdxi(2,1)

                  if ( DetJ <= 1.0d-12 ) then
			      write(*,*) 'error: Jacobian is zero or negative at element, igp:', m, ig, DetJ; stop
                        call flush(6)
		      else
			      tmp = 1.0d0 / DetJ
                        jac(ig,m) = DetJ
		      end if

                  dxidx(1,1) =  dxdxi(2,2) * tmp; dxidx(1,2) = -dxdxi(1,2) * tmp
			dxidx(2,1) = -dxdxi(2,1) * tmp; dxidx(2,2) =  dxdxi(1,1) * tmp
			                  
                  do i = 1, 3
			      shpx(i,ig,m) = sqd(1,i,ig) * dxidx(1,1) + sqd(2,i,ig) * dxidx(2,1)
			      shpy(i,ig,m) = sqd(1,i,ig) * dxidx(1,2) + sqd(2,i,ig) * dxidx(2,2)
		      end do
		end do
	end do
      deallocate(dxdxi, dxidx, sqd)

end subroutine makeSSF

subroutine maketau      
      use dimdata
      implicit none
      integer, dimension(3) :: n_local
      double precision :: u, v, u_norm
      double precision :: he, sum_abs_G, dot_product_val
      double precision :: dt05, t1, t2, t3
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
                  t1 = dt05 * dt05
                  t2 = (2.0d0 * u_norm / he)**2
                  t3 = (4.0d0 * alpha / (he * he))**2

                  ! tau(ig,m) = 1.0d0 / dsqrt(t1 + t2 + t3)
                  tau(ig,m) = 1.0d0 / dsqrt(t2 + t3)


            end do
      end do
end subroutine maketau

subroutine makemat
      use dimdata
      implicit none
      integer :: m, i, j, ig
      integer, dimension(3) :: n_local
      double precision :: u_int, v_int, wj, adv_i, adv_j
      double precision, dimension(3,3) :: Me, Se, Ke
      double precision, dimension(3,3) :: Me_supg, Se_supg, E_elem
      
      ad(:) = 0.0d0
      ea(:, :, :) = 0.0d0

      do m = 1, nelm
            n_local = nc(:,m)
            
            ! Initialize element matrices
            Me = 0.0d0; Se = 0.0d0; Ke = 0.0d0
            Me_supg = 0.0d0; Se_supg = 0.0d0
            E_elem = 0.0d0
            
            do ig = 1, mgp      
                  wj = jac(ig,m) * sg(ig)
                  
                  ! Velocity at integration point
                  u_int = dot_product(shp(:,ig,m), us(n_local))
                  v_int = dot_product(shp(:,ig,m), vs(n_local))
                  
                  do j = 1, 3
                        do i = 1, 3
                              
                              ! Advection term components (v . grad(N))
                              adv_i = u_int * shpx(i,ig,m) + v_int * shpy(i,ig,m)
                              adv_j = u_int * shpx(j,ig,m) + v_int * shpy(j,ig,m)
                   
                              ! 1. Standard Mass Matrix (from advection_p1.f90 [cite: 267])
                              Me(i,j) = Me(i,j) + shp(i,ig,m) * shp(j,ig,m) * wj 
                              
                              ! 2. Standard Advection Matrix (from advection_p1.f90 [cite: 268])
                              Se(i,j) = Se(i,j) + shp(i,ig,m) * adv_j * wj 
 
                              ! 3. Standard Diffusion Matrix (from poisson.f90, added 'alpha')
                              Ke(i,j) = Ke(i,j) + alpha * (shpx(i,ig,m) * shpx(j,ig,m) &
                                                   + shpy(i,ig,m) * shpy(j,ig,m)) * wj
                              
                              ! 4. SUPG Mass Matrix (from advection_p1.f90 [cite: 270])
                              Me_supg(i,j) = Me_supg(i,j) + adv_i * tau(ig,m) * shp(j,ig,m) * wj
                              
                              ! 5. SUPG Advection (Stabilization) Matrix (from advection_p1.f90 [cite: 271])
                              Se_supg(i,j) = Se_supg(i,j) + adv_i * tau(ig,m) * adv_j * wj

                        end do
                  end do
            end do

            ! Combine all matrices for the LHS
            ! E_elem = (1/dt * M_total) + (theta * K_total)
            E_elem = dti * (Me + Me_supg) + th * (Se + Se_supg + Ke)
       
            ! Assemble into global element-wise matrix array 'ea'
            do j = 1, 3
                  do i = 1, 3
                        ea(i, j, m) = E_elem(i,j)
                  end do
                  ! Assemble diagonal for scaling 
                  ad(n_local(j)) = ad(n_local(j)) + E_elem(j,j) 
            end do
      end do

end subroutine makemat

subroutine makebuvp
      use dimdata
      implicit none
      integer :: m, i, j, ig
      integer, dimension(3) :: n_local
      double precision :: u_int, v_int, wj, adv_i, adv_j, th0
      double precision :: x_int, y_int, pi, t_n
      double precision :: f_val_n, f_prt1, f_prt2, f_prt3
      double precision, dimension(3) :: Fe, Fe_supg
      
      ! Element-level matrices
      double precision, dimension(3,3) :: Me, Se, Ke, Me_supg, Se_supg
      ! Element-level vectors
      double precision, dimension(3) :: phi_e_n, be_M, be_SCK, be_u_n
      double precision, dimension(3) :: be_F, be

      
      pi = acos(-1.0d0)
      buvp(:) = 0.0d0
      th0 = 1.0d0 - th
      t_n = time - dt

      do m = 1, nelm
            n_local = nc(:,m)
            phi_e_n = ff(n_local) ! Get solution {u}^n for the element [cite: 283]

            ! Initialize element matrices and vectors
            Me = 0.0d0; Se = 0.0d0; Ke = 0.0d0
            Me_supg = 0.0d0; Se_supg = 0.0d0
            Fe = 0.0d0; Fe_supg = 0.0d0
            
            do ig = 1, mgp
                  wj = jac(ig,m) * sg(ig)
                  u_int = dot_product(shp(:,ig,m), us(n_local))
                  v_int = dot_product(shp(:,ig,m), vs(n_local))
                  
                  ! === 1. Calculate matrices for the {u}^n part ===
                  do j = 1, 3
                        do i = 1, 3
                              adv_i = u_int * shpx(i,ig,m) + v_int * shpy(i,ig,m)
                              adv_j = u_int * shpx(j,ig,m) + v_int * shpy(j,ig,m)
                              
                              Me(i,j) = Me(i,j) + shp(i,ig,m) * shp(j,ig,m) * wj
                              Se(i,j) = Se(i,j) + shp(i,ig,m) * adv_j * wj
                              Ke(i,j) = Ke(i,j) + alpha * (shpx(i,ig,m) * shpx(j,ig,m) &
                                                   + shpy(i,ig,m) * shpy(j,ig,m)) * wj

                              Me_supg(i,j) = Me_supg(i,j) + adv_i * tau(ig,m) * shp(j,ig,m) * wj
                              Se_supg(i,j) = Se_supg(i,j) + adv_i * tau(ig,m) * adv_j * wj
                        end do
                  end do
                  
                  x_int = dot_product(shp(:,ig,m), xx(1, n_local))
                  y_int = dot_product(shp(:,ig,m), xx(2, n_local))
                  
                  ! f at t=t_n
                  f_prt1 = cos(t_n) * sin(pi*x_int) * sin(pi*y_int) 
                  f_prt2 = 2.0d0 * alpha * pi**2 * sin(t_n) * sin(pi*x_int) * sin(pi*y_int)
                  f_prt3 = pi * sin(t_n) * (u_int * cos(pi*x_int)*sin(pi*y_int) &
                                               + v_int * sin(pi*x_int)*cos(pi*y_int))
                  
                  f_val_n = f_prt1 + f_prt2 + f_prt3
                  
                  do i = 1, 3
                        Fe(i) = Fe(i) + shp(i,ig,m) * f_val_n * wj
                        adv_i = u_int * shpx(i,ig,m) + v_int * shpy(i,ig,m)
                        Fe_supg(i) = Fe_supg(i) + (adv_i * tau(ig,m)) * f_val_n * wj
                  end do
            end do
            
            ! === 3. Combine parts and assemble ===
      
            ! Part 1: Contribution from {u}^n (from advection_p1.f90 [cite: 291])
            be_M = matmul(Me + Me_supg, phi_e_n)
            be_SCK = matmul(Se + Se_supg + Ke, phi_e_n)
            be_u_n = dti * be_M - th0 * be_SCK

            ! Part 2: Contribution from source {F_theta}
            be_F = Fe + Fe_supg
            
            ! Total element RHS vector
            be = be_u_n + be_F
            ! be = be_u_n

            ! Assemble into global RHS vector
            do i = 1, 3
                  buvp(n_local(i)) = buvp(n_local(i)) + be(i)
            end do
      end do

end subroutine makebuvp

subroutine dscale

      use dimdata
      implicit none
      integer :: n, m, i, j
      double precision, allocatable :: d(:)
      nx = node
      mx = nelm
      allocate( d(3) )

      do n = 1, nx
            !ad(n) = 1.0d0 / dsqrt(ad(n))
            ad(n) = 1.0d0 / (ad(n))
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

subroutine write_vtk(filename, xx, ff, nc, node, nelm)
      implicit none
      character(len=*), intent(in) :: filename
      double precision, intent(in) :: xx(2, node)
      double precision, intent(in) :: ff(node)
      integer, intent(in) :: nc(3, nelm)
      integer, intent(in) :: node, nelm
  
      integer :: i, m
      integer :: vtk_unit

      open(newunit=vtk_unit, file=filename, status='replace')

      write(vtk_unit,'(A)') '# vtk DataFile Version 3.0'
      write(vtk_unit,'(A)') 'Advection result'
      write(vtk_unit,'(A)') 'ASCII'
      write(vtk_unit,'(A)') 'DATASET UNSTRUCTURED_GRID'

      write(vtk_unit,'(A, 1X, I0, 1X, A)') 'POINTS', node, 'double'
      do i = 1, node
            write(vtk_unit,'(3F16.8)') xx(1,i), xx(2,i), 0.0d0
      end do

      write(vtk_unit,'(A, 1X, I0, 1X, I0)') 'CELLS', nelm, nelm*4
      do m = 1, nelm
            write(vtk_unit,'(4I8)') 3, nc(1,m)-1, nc(2,m)-1, nc(3,m)-1
      end do

      write(vtk_unit,'(A, 1X, I0)') 'CELL_TYPES', nelm
      do m = 1, nelm
            write(vtk_unit,'(I8)') 5
      end do

      write(vtk_unit,'(A, 1X, I0)') 'POINT_DATA', node
      write(vtk_unit,'(A)') 'SCALARS phi double 1'
      write(vtk_unit,'(A)') 'LOOKUP_TABLE default'
      do i = 1, node
            write(vtk_unit,'(F18.9)') ff(i)
      end do

    close(vtk_unit)
end subroutine write_vtk

subroutine write_vtk_seq

      use dimdata
      implicit none
      character(len=200) :: filename

      write(filename, '(A, I4.4, ".vtk")') trim(vtkfile_base), istep
      call write_vtk(filename, xx, ff, nc, node, nelm)

end subroutine write_vtk_seq





