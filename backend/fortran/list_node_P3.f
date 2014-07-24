
      subroutine list_node_P3 (nel, npt, nnodes, n_edge, 
     *    npt_p3, table_nod, table_N_E_F, visite)
c
      implicit none
      integer*8 nel, npt, nnodes
      integer*8 n_edge, npt_p3
      integer*8 visite(npt)

      integer*8 table_nod(nnodes,nel), table_N_E_F(14,nel)

c     Local variables
      integer*8 nnodes_0
      parameter (nnodes_0 = 6)
      integer*8 nod_el_p(nnodes_0)
      integer*8 j, k, j1, m1, m2
      integer*8 iel, inod, iel2, inod2, n_face
      integer*8 mm, nn, mm2
      integer*8 vert_1(2), vert_2(2)
c
      n_face = nel
c
      do j=1,npt
        visite(j) = 0
      enddo
c
      mm = 4 ! The first 4 entries of table_N_E_F(*,i) correspond to face and edges
      mm2 = n_edge + n_face
      nn = 0
      do iel=1,nel
        do inod=1,nnodes
          k = table_nod(inod,iel)
          nod_el_p(inod) = k
        enddo
c       P3 element: Vertices
        do inod=1,3
          k = nod_el_p(inod)
          if(visite(k) .eq. 0) then
            visite(k) = iel
            nn = nn + 1
            table_N_E_F(inod+mm,iel) = nn + mm2
          else
            iel2 = visite(k)
            inod2 = 0
            do j=1,3
               j1 = table_nod(j,iel2)
              if (k .eq. j1) inod2 = j
            enddo
            if(inod2 .eq. 0) then
              print*, "list_node_P3: problem with a vertex ", iel, inod
              print*, "k, visite(k) = ", k, visite(k)
              stop
            endif
            table_N_E_F(inod+mm,iel) = table_N_E_F(inod2+mm,iel2)
          endif
        enddo
c       P3 element: nodes on the edges
        do inod=4,6
          k = nod_el_p(inod)
c         Vertices of the edge
          j1 = inod - 3
          vert_1(1) = table_nod(j1,iel)
          j1 = inod - 2
          if (j1 .gt. 3) j1 = j1 - 3 
          vert_1(2) = table_nod(j1,iel)
          if(visite(k) .eq. 0) then
            visite(k) = iel
            m1 = 2*(inod-4)+1
            m2 = 2*(inod-4)+2
            do j1=m1,m2
              nn = nn+1
              table_N_E_F(j1+3+mm,iel) = nn + mm2
            enddo
          else
            iel2 = visite(k)
            inod2 = 0
            do j=4,6
              j1=table_nod(j,iel2)
              if (k .eq. j1) then
                inod2 = j
c               Vertices of the edge
                j1 = inod2 - 3
                vert_2(1) = table_nod(j1,iel2)
                j1 = inod2 - 2
                if (j1 .gt. 3) j1 = j1 - 3 
                vert_2(2) = table_nod(j1,iel2)
              endif
            enddo
            if(inod2 .eq. 0) then
              print*, "list_node_P3: problem with a node ", iel, inod
              stop
            endif
            do j=1,2
c             local numbering along the edges
              if (vert_2(1) .eq. vert_1(1) .and. 
     *          vert_2(2) .eq. vert_1(2)) then
c               The nodes on the edges inod and inod2 are numbered in the same order
c               This is possible only when the elements iel and iel2 have opposite orientations
                m1 = j+2*(inod-4)+3
                m2 = j+2*(inod2-4)+3
                table_N_E_F(m1+mm,iel) = table_N_E_F(m2+mm,iel2)
              elseif (vert_2(1) .eq. vert_1(2) .and. 
     *          vert_2(2) .eq. vert_1(1)) then
c               The nodes on the edges inod and inod2 are numbered in the opposite order
                j1 = 3 - j ! the numbering of the nodes are reversed
                m1 = j1+2*(inod-4)+3
                m2 = j+2*(inod2-4)+3
                table_N_E_F(m1+mm,iel) = table_N_E_F(m2+mm,iel2)
              else
                write(*,*) "list_node_P3: problems: ",
     *            "Check the edge endpoints"
                write(*,*) "inod, table_nod(inod,iel) = ", inod, 
     *            table_nod(inod,iel)
                write(*,*) "inod2, table_nod(inod2,iel2) = ", inod2, 
     *            table_nod(inod2,iel2)
                write(*,*) "iel, iel2 = ", iel, iel2
                write(*,*) "vert_1 = ", vert_1
                write(*,*) "vert_2 = ", vert_2
                write(*,*) "list_node_P3: Aborting..."
                stop
              endif
            enddo
          endif
        enddo
c       Numbering the interior nodes of the triangle
        do j=1,1  ! there is only one interior node for a P3 triangle
          nn = nn+1
          table_N_E_F(j+9+mm,iel) = nn + mm2
        enddo
      enddo
      npt_p3 = nn

      return
      end
