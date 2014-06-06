c
c

c***********************************************************************
c***********************************************************************
c
c     Input : initial P2 FEM mesh
c     Output: Set the coordinates and note type of the P3 node
c
c***********************************************************************
c
c
      subroutine get_coord_p3(nel, npt, nnodes, n_ddl, 
     *      table_nod, type_nod, table_N_E_F, 
     *      type_N_E_F, x, x_N_E_F, visite)

c
      implicit none
      integer*8 nel, npt, nnodes, n_ddl
      integer*8 table_nod(nnodes,nel), table_N_E_F(14,nel)
      integer*8 type_nod(npt), type_N_E_F(2,n_ddl)
      integer*8 visite(n_ddl)

      double precision x(2,npt), x_N_E_F(2,n_ddl)

c     Variable local
      integer*8 nddl_0
      parameter (nddl_0 = 14)
      integer*8 j, k, k1, n, ind, ip(2,3)
      integer*8 iel, inod, inod1, inod2, mm
      integer*8 nut0(6), nut_N_E_F(nddl_0)

      double precision xx1, xx2, xx3, yy1, yy2, yy3
      double precision dx1, dy1
      double precision tmp1, tmp2, tmp3, tmp_x, tmp_y
c
c     debut de la programmation de la sous-routine maillage
c     -----------------------------------------------------
c
      if ( nnodes .ne. 6 ) then
        write(*,*) "get_coord_p3: problem nnodes = ", nnodes
        write(*,*) "get_coord_p3: nnodes should be equal to 6 !"
        write(*,*) "get_coord_p3: Aborting..."
        stop
      endif
c
c     ip(1,i) = i+1 MOD 3 ; Number of the next vertex to vertex i
c     ip(2,i) = i+2 MOD 3 ; Number of the second next vertex to vertex i
      ip(1,1) = 2
      ip(1,2) = 3
      ip(1,3) = 1
c
      ip(2,1) = 3
      ip(2,2) = 1
      ip(2,3) = 2
c
cccccccccccccccccccccccc
c
c     Initialisation
      do j=1,n_ddl
        visite(j) = 0
      enddo
c
cccccccccccccccccccccccc
c
      mm = 4 ! The first 4 entries of table_N_E_F(*,i) correspond to face and edges
      do iel=1,nel
        do inod=1,nnodes
          nut0(inod) = table_nod(inod,iel)
        enddo
        do inod=1,10  ! the 10 node of a P3 element
          nut_N_E_F(inod) = table_N_E_F(inod+mm,iel)
        enddo
        do inod=1,3  ! scan the vertices ############
          k = nut0(inod)
          if(visite(k) .eq. 0) then
            visite(k) = iel
            inod1 = nut0(inod)
            inod2 = nut_N_E_F(inod)
            x_N_E_F(1,inod2) = x(1,inod1)
            x_N_E_F(2,inod2) = x(2,inod1)
            type_N_E_F(1,inod2) = type_nod(inod1)
            type_N_E_F(2,inod2) = 0 ! Vertex => dimension zero
          endif
        enddo
        do inod=4,nnodes ! scan the nodes located on the edges ############
          k=nut0(inod)
          if(k .lt. 1) then
            print*, 'k = ', k
            print*, 'iel, inod = ', iel, inod
            print*, 'nut0 = ', (nut0(inod2), inod2=1,nnodes)
            stop
          endif
          if(visite(k) .eq. 0) then
            visite(k) = iel
c           Endpoints of the edge
            k1 = nut0(inod-3)
            xx1 = x(1,k1)
            yy1 = x(2,k1)
            k1 = nut0(ip(1,inod-3))
            xx2 = x(1,k1)
            yy2 = x(2,k1)
            dx1 = (xx2-xx1)/3.0d0
            dy1 = (yy2-yy1)/3.0d0
            ind = type_nod(nut0(inod))  ! type of the mid-edge node of the initial P2 mesh
            do inod2=1,2  ! 2 nodes per edge (for P3 element)
              k1 = nut_N_E_F(inod2+2*(inod-4)+3)
              x_N_E_F(1,k1) = xx1 + inod2*dx1
              x_N_E_F(2,k1) = yy1 + inod2*dy1
              type_N_E_F(1,k1) = ind
              type_N_E_F(2,k1) = 0 ! Node => dimension zero
            enddo
          endif
        enddo
c       Coordonate of the vertices
        k1 = nut0(1)
        xx1 = x(1,k1)
        yy1 = x(2,k1)
        k1 = nut0(2)
        xx2 = x(1,k1)
        yy2 = x(2,k1)
        k1 = nut0(3)
        xx3 = x(1,k1)
        yy3 = x(2,k1)
c       The tenth node is a at the center of the triangle
        n = 10 ! dimension(P3) = 10
        k1 = nut_N_E_F(n) ! this node is an interior node of the triangle ############
        tmp1 = 1.0d0/3.0d0
        tmp2 = 1.0d0/3.0d0
        tmp3 = 1.0d0/3.0d0
        tmp_x = xx1*tmp1+xx2*tmp2+xx3*tmp3
        tmp_y = yy1*tmp1+yy2*tmp2+yy3*tmp3
        x_N_E_F(1,k1) = tmp_x
        x_N_E_F(2,k1) = tmp_y
        type_N_E_F(1,k1) = 0  ! interior node
        type_N_E_F(2,k1) = 0 ! Node => dimension zero
      enddo
c
      return
      end
