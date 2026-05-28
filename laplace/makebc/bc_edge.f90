! =======================================================================
!
!       メッシュデータから境界辺を探索し、物理条件を自動で割り当てる境界条件生成プログラム
!
! =======================================================================
program bc_edge
    implicit none

    ! --- ファイル名 ---
    character(128) :: mesfile, bcfile

    ! --- メッシュデータ ---
    integer :: n_nodes, n_elems
    double precision, allocatable :: xx(:), yy(:)
    integer, allocatable :: nc(:,:)

    ! --- 境界辺探索用 ---
    integer, allocatable :: boundary_edges(:,:)
    integer :: n_boundary_edges

    ! --- 境界分類用 ---
    double precision :: xxmax, yymin
    double precision, parameter :: TOL = 1.0d-6 ! 浮動小数点数比較用の許容誤差
    integer :: n1, n2
    integer, allocatable :: potential_edges(:,:), velocity_in_edges(:,:), velocity_wall_edges(:,:)
    integer :: n_potential, n_velocity_in, n_velocity_wall

    ! --- ループカウンタ等 ---
    integer :: i, j

    ! ===================================================================
    !   STEP 1: ファイル名とメッシュデータの読み込み
    ! ===================================================================
    open(9, file = 'file_bc_gene.dat', status = 'old', action = 'read')
    read(9, '(a)') mesfile
    read(9, '(a)') bcfile
    close(9)
    write(*,*) 'Input mesh file : ', trim(mesfile)
    write(*,*) 'Output BC file  : ', trim(bcfile)

    open(10, file = mesfile, status = 'old', action = 'read')
    read(10, *) n_nodes, n_elems
    write(*,*) 'Nodes:', n_nodes, ' / Elements:', n_elems

    allocate(xx(n_nodes), yy(n_nodes), nc(3, n_elems))
    do i = 1, n_nodes
        read(10, *) j, xx(i), yy(i) ! 節点番号はインデックスと一致する
    end do
    do i = 1, n_elems
        read(10, *) j, nc(1, i), nc(2, i), nc(3, i)
    end do
    close(10)

    ! ===================================================================
    !   STEP 2: 境界辺の探索（出現回数=1で境界辺と判定）
    ! ===================================================================
    call find_boundary_edges(n_elems, nc, n_boundary_edges, boundary_edges)

    ! ===================================================================
    !   STEP 3: 座標の最大値・最小値を計算
    ! ===================================================================
    xxmax = maxval(xx)
    yymin = minval(yy)
    write(*,*) 'X_max = ', xxmax
    write(*,*) 'Y_min = ', yymin

    ! ===================================================================
    !   STEP 4: 境界辺を物理条件に基づいて分類
    ! ===================================================================
    ! 分類結果を格納する配列を、最大サイズ(全境界辺数)で確保
    allocate(potential_edges(2, n_boundary_edges))
    allocate(velocity_in_edges(2, n_boundary_edges))
    allocate(velocity_wall_edges(2, n_boundary_edges))
    n_potential = 0
    n_velocity_in = 0
    n_velocity_wall = 0

    do i = 1, n_boundary_edges
        n1 = boundary_edges(1, i)
        n2 = boundary_edges(2, i)

        ! --- ポテンシャル境界 (出口) の判定 ---
        if (abs(xx(n1) - xxmax) < TOL .and. abs(xx(n2) - xxmax) < TOL) then
            n_potential = n_potential + 1
            potential_edges(:, n_potential) = boundary_edges(:, i)

        ! --- 流入境界 (入口) の判定 ---
        else if (abs(yy(n1) - yymin) < TOL .and. abs(yy(n2) - yymin) < TOL) then
            n_velocity_in = n_velocity_in + 1
            velocity_in_edges(:, n_velocity_in) = boundary_edges(:, i)

        ! --- その他の境界 (壁) ---
        else
            n_velocity_wall = n_velocity_wall + 1
            velocity_wall_edges(:, n_velocity_wall) = boundary_edges(:, i)
        end if
    end do
    
    write(*,*) 'Potential edges (outlet) : ', n_potential
    write(*,*) 'Velocity=5.0 edges (inlet): ', n_velocity_in
    write(*,*) 'Velocity=0.0 edges (wall) : ', n_velocity_wall

    ! ===================================================================
    !   STEP 5: 分類結果をファイルに書き出し
    ! ===================================================================
    open(20, file = bcfile, status = 'replace', action = 'write')

    ! --- ポテンシャル条件の書き出し ---
    write(20, '(i5, 3x, a)') n_potential, 'potential'
    do i = 1, n_potential
        write(20, '(i5, 2i8, f12.6)') i, potential_edges(1, i), potential_edges(2, i), 1.0d0
    end do
    write(20, *) ! 空行

    ! --- 流入条件の書き出し ---
    write(20, '(i5, 3x, a)') n_velocity_in, 'velocity_in'
    do i = 1, n_velocity_in
        write(20, '(i5, 2i8, f12.6)') i, velocity_in_edges(1, i), velocity_in_edges(2, i), 5.0d0
    end do
    write(20, *) ! 空行

    ! --- 壁条件の書き出し ---
    write(20, '(i5, 3x, a)') n_velocity_wall, 'velocity_wall'
    do i = 1, n_velocity_wall
        write(20, '(i5, 2i8, f12.6)') i, velocity_wall_edges(1, i), velocity_wall_edges(2, i), 0.0d0
    end do

    close(20)
    write(*,*) 'Processing complete. BC file generated: ', trim(bcfile)

    ! ===================================================================
    !   後処理
    ! ===================================================================
    deallocate(xx, yy, nc, boundary_edges)
    deallocate(potential_edges, velocity_in_edges, velocity_wall_edges)

contains

! =======================================================================
!   境界辺を探索するサブルーチン
! =======================================================================
subroutine find_boundary_edges(n_elems, nc, n_boundary_edges, boundary_edges)
    implicit none
    integer, intent(in)  :: n_elems
    integer, intent(in)  :: nc(3, n_elems)
    integer, intent(out) :: n_boundary_edges
    integer, allocatable, intent(out) :: boundary_edges(:,:)

    integer, allocatable :: all_edges(:,:), unique_edges(:,:), counts(:)
    integer :: n_unique_edges
    integer :: n1, n2, n3, edge(2)
    integer :: i, j, m
    logical :: found    

    allocate(all_edges(2, 3 * n_elems))
    do m = 1, n_elems
        n1 = nc(1, m); n2 = nc(2, m); n3 = nc(3, m)
        all_edges(:, 3*m-2) = [min(n1, n2), max(n1, n2)]
        all_edges(:, 3*m-1) = [min(n2, n3), max(n2, n3)]
        all_edges(:, 3*m-0) = [min(n3, n1), max(n3, n1)]
    end do

    allocate(unique_edges(2, 3 * n_elems), counts(3 * n_elems))
    n_unique_edges = 0
    do i = 1, 3 * n_elems
        edge(:) = all_edges(:, i)
        found = .false.
        do j = 1, n_unique_edges
            if (unique_edges(1, j) == edge(1) .and. unique_edges(2, j) == edge(2)) then
                counts(j) = counts(j) + 1; found = .true.; exit
            end if
        end do
        if (.not. found) then
            n_unique_edges = n_unique_edges + 1
            unique_edges(:, n_unique_edges) = edge(:)
            counts(n_unique_edges) = 1
        end if
    end do
    deallocate(all_edges)

    n_boundary_edges = count(counts(1:n_unique_edges) == 1)
    allocate(boundary_edges(2, n_boundary_edges))
    j = 0
    do i = 1, n_unique_edges
        if (counts(i) == 1) then
            j = j + 1
            boundary_edges(:, j) = unique_edges(:, i)
        end if
    end do
    deallocate(unique_edges, counts)
    write(*,*) 'Found ', n_boundary_edges, ' boundary edges.'
end subroutine find_boundary_edges

end program bc_edge
