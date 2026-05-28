program mesh_p1_to_p2
      implicit none

      integer, parameter :: IUNIT = 10, OUNIT = 20
      ! character(len=256), parameter :: IN_FILE = "input.txt"
      ! character(len=256), parameter :: OUT_FILE = "output.txt"
      character(100) :: prefile, postfile



      integer :: num_nodes_p1, num_elements
      real(8), allocatable :: coords_p1(:,:)
      integer, allocatable :: elems_p1(:,:)

      integer :: num_nodes_p2
      real(8), allocatable :: coords_p2(:,:)
      integer, allocatable :: elems_p2(:,:)

      integer, allocatable :: edges(:,:) ! (n_small, n_large, mid_node_id)
      integer :: num_unique_edges, current_new_node_id
      integer :: i, j, k, stat, temp_id
      integer :: n1, n2, n3, n_pair(3, 2)
      integer :: u, v, temp
      logical :: found

      open(9, file='file.txt', status='unknown')
      read(9,'(a)') prefile
      read(9,'(a)') postfile
      close(9)

      open(10, file = prefile, status = 'unknown')
      open(11, file = postfile, status = 'unknown')

      ! open(IUNIT, file=IN_FILE, status='old', iostat=stat)

      read(10,*) num_nodes_p1, num_elements

      allocate(coords_p1(num_nodes_p1, 2))
      allocate(elems_p1(num_elements, 3))

      do i = 1, num_nodes_p1
            read(10, *) temp_id, coords_p1(i, 1), coords_p1(i, 2)
      end do
      do i = 1, num_elements
            read(10, *) temp_id, elems_p1(i, 1), elems_p1(i, 2), elems_p1(i, 3)
      end do
      close(10)
      write(6,*) "Input file read successfully."
      write(6,*) "  Number of nodes (P1)   : ", num_nodes_p1
      write(6,*) "  Number of elements     : ", num_elements

      allocate(coords_p2(num_nodes_p1 + num_elements * 3, 2))
      allocate(edges(num_elements * 3, 3))

      coords_p2(1:num_nodes_p1, :) = coords_p1(:,:)
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
                        coords_p2(current_new_node_id, 1) = (coords_p1(u, 1) + coords_p1(v, 1)) / 2.0_8
                        coords_p2(current_new_node_id, 2) = (coords_p1(u, 2) + coords_p1(v, 2)) / 2.0_8
                        current_new_node_id = current_new_node_id + 1
                  end if
            end do
      end do

      num_nodes_p2 = current_new_node_id - 1
      write(6,*) "Mid-point nodes generated."
      write(6,*) "  Number of nodes (P2)   : ", num_nodes_p2

      allocate(elems_p2(num_elements, 6))

      do i = 1, num_elements
            n1 = elems_p1(i, 1)
            n2 = elems_p1(i, 2)
            n3 = elems_p1(i, 3)
            elems_p2(i, 1:3) = [n1, n2, n3]

            if (num_unique_edges > 0) then
                  elems_p2(i, 4) = find_mid_node_id(n1, n2, edges(1:num_unique_edges, :))
                  elems_p2(i, 5) = find_mid_node_id(n2, n3, edges(1:num_unique_edges, :))
                  elems_p2(i, 6) = find_mid_node_id(n3, n1, edges(1:num_unique_edges, :))
            else
                  elems_p2(i, 4:6) = -1
            end if
      end do

      write(6,*) "P2 element connectivity created."

      ! open(OUNIT, file=OUT_FILE, status='replace')
      open(11, file = postfile, status = 'unknown')

      write(11, '(I10, I10)') num_nodes_p2, num_elements
      do i = 1, num_nodes_p2
            write(11, '(I8, 2F14.6)') i, coords_p2(i, 1), coords_p2(i, 2)
      end do
      do i = 1, num_elements
            write(11, '(I8, 6I8)') i, elems_p2(i, :)
      end do
      close(11)
      write(6,*) "Output file written successfully: ", trim(postfile)

      deallocate(coords_p1, elems_p1)
      deallocate(coords_p2, elems_p2, edges)

      contains

function find_mid_node_id(node_a, node_b, edge_db) result(mid_node_id)
      integer, intent(in) :: node_a, node_b
      integer, intent(in) :: edge_db(:, :)
      integer :: mid_node_id
      integer :: u, v, i
      integer :: db_size

      db_size = size(edge_db, dim=1)

      if (node_a < node_b) then
            u = node_a; v = node_b
      else
            u = node_b; v = node_a
      end if

      mid_node_id = -1 ! Error code
      do i = 1, db_size
            if (edge_db(i, 1) == u .and. edge_db(i, 2) == v) then
                  mid_node_id = edge_db(i, 3)
                  return
            end if
      end do
end function find_mid_node_id

end program mesh_p1_to_p2