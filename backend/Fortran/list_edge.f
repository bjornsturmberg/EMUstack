c
      subroutine list_edge (nel, npt, nnodes, 
     *    n_edge, type_nod, table_nod, 
     *    table_edge, table_edge_face, visite)
c
      implicit none
      integer*8 nel, npt, nnodes
      integer*8 n_edge
      integer*8 type_nod(npt), visite(npt)
      integer*8 table_nod(nnodes,nel), table_edge_face(14,nel)
      integer*8 table_edge(4,npt)
      integer*8 i, j, k, j1, j2, n_face, debug
      integer*8 list_end(2,3)

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

      n_face = nel
      debug = 0

      do i=1,npt
        visite(i) = 0
      enddo
      n_edge = 0
      do i=1,nel
        do j=4,nnodes
          if (type_nod(table_nod(j,i)) .ne. 0) then
            j1 = list_end(1,j-3)
            j2 = list_end(2,j-3)
            if (type_nod(table_nod(j1,i)) .eq. 0 .or. 
     *        type_nod(table_nod(j2,i)) .eq. 0) then
              write(*,*) "list_edge: table_nod = ", 
     *        type_nod(table_nod(j1,i)), type_nod(table_nod(j2,i)), 
     *        type_nod(table_nod(j,i)) 
              write(*,*) "type_nod(j1) = ", table_nod(j1,i)
              write(*,*) "type_nod(j2) = ", table_nod(j2,i)
              write(*,*) "type_nod(j) = ", table_nod(j,i)
              write(*,*) "list_edge: Aborting..."
              stop
            endif
          endif
        enddo
        do j=4,nnodes ! scan the element edge
          j1 = table_nod(j,i)
          k = visite(j1)
          if (k .eq. 0) then
            n_edge = n_edge + 1
            visite(j1) = n_edge
            j2 = list_end(1,j-3)
            table_edge(1,n_edge) = table_nod(j2,i)
            j2 = list_end(2,j-3)
            table_edge(2,n_edge) = table_nod(j2,i)
            table_edge(3,n_edge) = j1
c           Table of connectivity for the face (with respect to the triangle element)
            table_edge_face(j-2,i) = n_edge + n_face
            table_edge(4,n_edge) = n_edge + n_face
          else
            table_edge_face(j-2,i) = k + n_face
            table_edge(4,k) = k + n_face
          endif
        enddo
      enddo
      if (debug .eq. 1) then
        write(*,*) "list_edge: npt, n_edge, nel = ", npt,n_edge,nel
      endif
c
      return
      end
