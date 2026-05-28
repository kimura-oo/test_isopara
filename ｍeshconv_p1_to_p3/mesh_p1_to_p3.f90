program mesh_p1_to_p3
      implicit none
      integer, parameter :: IUNIT = 10, OUNIT = 20
      character(len=256), parameter :: IN_FILE = "input.txt"
      character(len=256), parameter :: OUT_FILE = "output.p3.txt"
      character(100) :: prefile, postfile

      integer :: num_nodes_p1, num_elements
      real(8), allocatable :: coords_p1(:,:)
      integer, allocatable :: elems_p1(:,:)

      integer :: num_nodes_p3
      real(8), allocatable :: coords_p3(:,:)
      integer, allocatable :: elems_p3(:,:)

      integer, allocatable :: edges(:,:) ! (n_small, n_large, new_node_id_1, new_node_id_2)
      integer, allocatable :: internal_nodes(:) ! (element_id) -> new_node_id
      integer :: num_unique_edges, current_new_node_id
      integer :: i, j, k, stat, temp_id
      integer :: n1, n2, n3, n_pair(3, 2)
      integer, allocatable :: edge_node_ids(:,:)
      integer :: u, v, temp
      logical :: found

      open(9, file='file.txt', status='unknown', iostat=stat)
      if (stat == 0) then
            read(9,'(a)') prefile
            read(9,'(a)') postfile
            close(9)
      else
            prefile = IN_FILE
            postfile = OUT_FILE
      end if

      open(IUNIT, file=prefile, status='old', iostat=stat)
      if (stat /= 0) then
            write(6,*) "Error: Input file not found: ", trim(prefile)
            stop
      end if

      read(IUNIT,*) num_nodes_p1, num_elements

      allocate(coords_p1(num_nodes_p1, 2))
      allocate(elems_p1(num_elements, 3))

      do i = 1, num_nodes_p1
            read(IUNIT, *) temp_id, coords_p1(i, 1), coords_p1(i, 2)
      end do
      do i = 1, num_elements
            read(IUNIT, *) temp_id, elems_p1(i, 1), elems_p1(i, 2), elems_p1(i, 3)
      end do
      close(IUNIT)
      write(6,*) "Input file read successfully: ", trim(prefile)
      write(6,*) "  Number of nodes (P1)   : ", num_nodes_p1
      write(6,*) "  Number of elements     : ", num_elements

      allocate(coords_p3(num_nodes_p1 + num_elements * 3 * 2 + num_elements, 2))
      allocate(edges(num_elements * 3, 4))
      allocate(internal_nodes(num_elements))
      allocate(elems_p3(num_elements, 10))
      allocate(edge_node_ids(3, 2))

      coords_p3(1:num_nodes_p1, :) = coords_p1(:,:)
      num_unique_edges = 0
      current_new_node_id = num_nodes_p1 + 1

      do i = 1, num_elements
            n_pair(1, :) = [elems_p1(i, 1), elems_p1(i, 2)]
            n_pair(2, :) = [elems_p1(i, 2), elems_p1(i, 3)]
            n_pair(3, :) = [elems_p1(i, 3), elems_p1(i, 1)]

            do j = 1, 3
                  u = n_pair(j, 1)
                  v = n_pair(j, 2)

                  if (u > v) then
                        temp = u; u = v; v = temp
                  end if

                  found = .false.
                  do k = 1, num_unique_edges
                        if (edges(k, 1) == u .and. edges(k, 2) == v) then
                              found = .true.
                              exit
                        end if
                  end do

                  if (.not. found) then
                        num_unique_edges = num_unique_edges + 1
                        edges(num_unique_edges, 1) = u
                        edges(num_unique_edges, 2) = v
                        edges(num_unique_edges, 3) = current_new_node_id
                        edges(num_unique_edges, 4) = current_new_node_id + 1

                        coords_p3(current_new_node_id, :) = &
                             (2.0d0 * coords_p1(u, :) + 1.0d0 * coords_p1(v, :)) / 3.0d0
                        coords_p3(current_new_node_id + 1, :) = &
                             (1.0d0 * coords_p1(u, :) + 2.0d0 * coords_p1(v, :)) / 3.0d0

                        current_new_node_id = current_new_node_id + 2
                  end if
            end do
      end do
      write(6,*) "Edge nodes generated."

      do i = 1, num_elements
            n1 = elems_p1(i, 1)
            n2 = elems_p1(i, 2)
            n3 = elems_p1(i, 3)
            coords_p3(current_new_node_id, :) = &
                  (coords_p1(n1, :) + coords_p1(n2, :) + coords_p1(n3, :)) / 3.0d0
            internal_nodes(i) = current_new_node_id
            current_new_node_id = current_new_node_id + 1
      end do
      write(6,*) "Internal nodes (centroids) generated."

      num_nodes_p3 = current_new_node_id - 1
      write(6,*) "  Number of nodes (P3)   : ", num_nodes_p3

      do i = 1, num_elements
            n1 = elems_p1(i, 1)
            n2 = elems_p1(i, 2)
            n3 = elems_p1(i, 3)
            elems_p3(i, 1:3) = [n1, n2, n3]

            call find_edge_node_ids(n1, n2, edges(1:num_unique_edges,:), edge_node_ids(1,:))
            call find_edge_node_ids(n2, n3, edges(1:num_unique_edges,:), edge_node_ids(2,:))
            call find_edge_node_ids(n3, n1, edges(1:num_unique_edges,:), edge_node_ids(3,:))
            elems_p3(i, 4:5) = edge_node_ids(1, :)
            elems_p3(i, 6:7) = edge_node_ids(2, :)
            elems_p3(i, 8:9) = edge_node_ids(3, :)

            elems_p3(i, 10) = internal_nodes(i)
      end do
      write(6,*) "P3 element connectivity created."

      open(OUNIT, file=postfile, status='replace')

      write(OUNIT, '(I10, I10)') num_nodes_p3, num_elements
      do i = 1, num_nodes_p3
            write(OUNIT, '(I8, 2F14.6)') i, coords_p3(i, 1), coords_p3(i, 2)
      end do
      do i = 1, num_elements
            write(OUNIT, '(I8, 10I8)') i, elems_p3(i, :)
      end do
      close(OUNIT)
      write(6,*) "Output file written successfully: ", trim(postfile)

      deallocate(coords_p1, elems_p1)
      deallocate(coords_p3, elems_p3, edges, internal_nodes, edge_node_ids)

contains

subroutine find_edge_node_ids(node_a, node_b, edge_db, node_ids)
      integer, intent(in) :: node_a, node_b
      integer, intent(in) :: edge_db(:, :)
      integer, intent(out) :: node_ids(2)
      integer :: u, v, i
      integer :: db_size

      db_size = size(edge_db, dim=1)

      ! êﬂì_î‘çÜÇÃè¨Ç≥Ç¢ï˚ÇuÇ∆Ç∑ÇÈ
      if (node_a < node_b) then
            u = node_a; v = node_b
      else
            u = node_b; v = node_a
      end if

      node_ids = -1 ! Error code

      do i = 1, db_size
            if (edge_db(i, 1) == u .and. edge_db(i, 2) == v) then
                  ! óvëfÇÃï”ÇÃå¸Ç´Ç…çáÇÌÇπÇƒêﬂì_IDÇÃèáèòÇåàíËÇ∑ÇÈ
                  if (node_a == u) then
                      node_ids(1) = edge_db(i, 3)
                      node_ids(2) = edge_db(i, 4)
                  else
                      node_ids(1) = edge_db(i, 4)
                      node_ids(2) = edge_db(i, 3)
                  end if
                  return
            end if
      end do
end subroutine find_edge_node_ids

end program mesh_p1_to_p3