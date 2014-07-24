c
c***********************************************************************
c
      subroutine type_node_edge_face (nel, npt, nnodes, n_ddl, 
     *      type_nod, table_nod, table_N_E_F, 
     *      visite, type_N_E_F, x, x_E_F)
c
c***********************************************************************
c
      implicit none
      integer*8 nel, npt, nnodes, n_ddl
      integer*8 type_nod(npt)
      integer*8 table_nod(nnodes,nel), table_N_E_F(14,nel)
      integer*8 visite(n_ddl), type_N_E_F(2,n_ddl)
      double precision x(2,npt), x_E_F(2,n_ddl)
c     Local variables

      integer*8 i, j, j1
      integer*8 type_n(10)
c      integer*8 list_point_F(6,4)
      integer*8 nddl_0
      parameter (nddl_0 = 14)
      double precision xel(2,6)
c
c***********************************************************************
c
c
      if ( nnodes .ne. 6 ) then
        write(*,*) "type_node_edge_face: problem nnodes = ", nnodes
        write(*,*) "type_node_edge_face: nnodes should be equal to 6 !"
        write(*,*) "type_node_edge_face: Aborting..."
        stop
      endif
c
      if ( nnodes .ne. 6 ) then
        write(*,*) "type_node_edge_face: problem nnodes = ", nnodes
        write(*,*) "type_node_edge_face: nnodes should be equal to 6 !"
        write(*,*) "type_node_edge_face: Aborting..."
        stop
      endif
c
cccccccccccccccccccccccc
c
c     Initialisation
      do j=1,n_ddl
        type_N_E_F(1,j) = 0
        type_N_E_F(2,j) = 0
      enddo
      do i=1,n_ddl
        visite(i) = 0
      enddo
c
      do i=1,nel
        do j=1,nnodes
          j1 = table_nod(j,i)
          type_n(j) = type_nod(j1)
          xel(1,j) = x(1,j1)
          xel(2,j) = x(2,j1)
        enddo
        j=1 ! an element is a face
        j1 = table_N_E_F(j,i)
        x_E_F(1,j1) = (xel(1,1) + xel(1,2) + xel(1,3))/3.0d0 ! centre of the elements
        x_E_F(2,j1) = (xel(2,1) + xel(2,2) + xel(2,3))/3.0d0
        type_N_E_F(1,j1) = 0  ! Topologically, a face is an interior domain
        type_N_E_F(2,j1) = 2 ! Face => dimension two
        do j=1,3 ! scan the 3 element edges
          j1 = table_N_E_F(j+1,i)
          x_E_F(1,j1) = xel(1,j+3)
          x_E_F(2,j1) = xel(2,j+3)
          if (visite(j1) .eq. 0) then
            visite(j1) = 1
            type_N_E_F(1,j1) = type_n(j+3)
            type_N_E_F(2,j1) = 1 ! Edge => dimension one
          endif
        enddo
      enddo
c
      return
      end
