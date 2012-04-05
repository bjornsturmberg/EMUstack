cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     P2 basis function over the unit Tetrahedron
c
c         Quadradic basis function = P2 * Grad P1
c
c          basis_list(1,j,i) = k : number on data to be stored: if k=3 only one gradient will be used; k=4 => 2 gradients
c          basis_list(2,j,i) = m : corresponds to the P2 Lagrange polynomial phi_m
c          basis_list(3,j,i) = n : corresponds to the gradient of the P1 Lagrange polynomial phi_n
c          basis_list(4,j,i)     : it will be used only if k=4  
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine basis_ls (nod_el, basis_list)
c
      implicit none
      integer*8 nnodes, nddl_t
      parameter (nnodes = 6, nddl_t=4)
      integer*8 nod_el(nnodes), basis_list(4,3,nddl_t)
c     Local variables
      integer*8 i, j, j1, j2, j3, list_end(2,3)
      integer*8 ls_e(3), ls_e_sorted(3)
      integer*8 ls_n(3), ls_n_sorted(3)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Endpoints of the 6 edges (mid-point) of the reference tetrahedron
c
      i = 1
      list_end(1,i) = 1
      list_end(2,i) = 2
      i = 2
      list_end(1,i) = 2
      list_end(2,i) = 3
      i = 3
      list_end(1,i) = 1
      list_end(2,i) = 3
c
      do i=1,1 ! scan the element face
c       The mid-edge nodes of the face
        do j=1,3 ! scan the mid-edge nodes of the face
          basis_list(1,j,i) = 3  ! number on data to be stored
          basis_list(2,j,i) = j+3 ! the mid-edge number
          j2 = mod(j+2,3) ! give the node opposite to the mid-edge node (j+3)
          if( j2 .eq. 0 ) j2 = 3
          basis_list(3,j,i) = j2
          basis_list(4,j,i) = 0  !  actually, it will not be used 
        enddo
      enddo
      do i=2,4 ! scan the 3 element edges
c       2 end-point basis vectors are attached to the edge i
        do j=1,2 ! scan the end nodes of the edge
          j1 = list_end(j,i-1)
          ls_n(j) = nod_el(j1)
        enddo
        ls_n(3) = 0
        j3 = 2
        call sort_n(j3, ls_n, ls_n_sorted)
        j1 = ls_n_sorted(1)
        j1 = list_end(j1,i-1)
        j2 = ls_n_sorted(2)
        j2 = list_end(j2,i-1)
        j = 1
        basis_list(1,j,i) = 3  ! number on data to be stored
        basis_list(2,j,i) = j1
        basis_list(3,j,i) = j2
        basis_list(4,j,i) = 0
        j = 2
        basis_list(1,j,i) = 3  ! number on data to be stored
        basis_list(2,j,i) = j2
        basis_list(3,j,i) = j1
        basis_list(4,j,i) = 0
        j = 3
        basis_list(1,j,i) = 4  ! number on data to be stored
        basis_list(2,j,i) = i+2 ! add 2 to get the correct edge number
        basis_list(3,j,i) = j1
        basis_list(4,j,i) = j2
          if (j1 .eq. j2) then
            write(*,*) "basis_ls: j1 = j2:"
            write(*,*) "basis_ls: ", i, j, j1, j2
            write(*,*) "basis_ls: Aborting..."
            stop
          endif
      enddo
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      return
      end
