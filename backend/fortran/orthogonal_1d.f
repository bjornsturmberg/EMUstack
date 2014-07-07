C   Calculate the Overlap integral of the prime and adjoint Plane Waves
C
      subroutine orthogonal_1d (nval, nel, npt_P2, 
     *  nb_typ_el, pp, bloch_vec_y, table_nod, 
     *  type_el, x_P2, beta_1, beta_2,
     *  sol_1, sol_2, mat_overlap, overlap_file, PrintAll,
     *  pair_warning, k_0)
c
      implicit none
      integer*8 nval, nel, npt_P2, nb_typ_el

      integer*8 type_el(nel)
      integer*8 table_nod(3,nel)
      double precision x_P2(npt_P2), k_0
      complex*16 sol_1(3+4+4,nval,nel)
      complex*16 sol_2(3+4+4,nval,nel)
      complex*16 pp(nb_typ_el)
      complex*16 beta_1(nval), beta_2(nval)
      complex*16, dimension(nval,nval) :: mat_overlap
      character overlap_file*100
      double precision bloch_vec_y
c     Local variables
      integer*8 nnodes_0, nddl_0
      parameter (nnodes_0 = 3, nddl_0 = 11)
      integer*8 nod_el_p(nnodes_0)
      complex*16 vec_1(nddl_0)
      integer*8 node_P2, node_P3
      parameter (node_P2 = 3, node_P3 = 4)
      double precision mat_Mxx_0(node_P2,node_P2)
      double precision mat_Kxy_0(node_P2,node_P3)
      double precision mat_Kyy_0(node_P3,4)
      double precision mat_Myy_0(4,4)
      double precision mat_Kyx_0(4,node_P2)
c     nddl_0 = node_P2+2*node_P3
      complex*16 mat_el_2(nddl_0,nddl_0)
      complex*16 sol_el_1(nddl_0), sol_el_2(nddl_0)
      double precision xmin, xmax
      integer*8 i, j, j1, typ_e
      integer*8 iel, ival, jval
      integer*8 debug, ui
      double precision r_tmp1
      complex*16 z_tmp1, z_tmp2, z_beta_1
      complex*16 ii
c
C     Mode ordering
      integer*8 skip, PrintAll, pair_warning
      complex*16 betatmp1(1), betatmp2(1)
      complex*16 soltmp1(3+4+4,nel,1)
      complex*16 soltmp2(3+4+4,nel,1)
      integer*8 elcount, nodecount, redo, j2
      double precision val_max_diag, val_max_off
C
C
CCCCCCCCCCCCCCCCCCCCC Start Program CCCCCCCCCCCCCCCCCCCCCCCC
C
c  ii = sqrt(-1)
      ii = dcmplx(0.0d0, 1.0d0)

      ui = 6
      debug = 0
c
c      pi = 3.141592653589793d0
      redo = 0
122   continue                !second rearranged overlap
C
      do jval=1,nval
        do ival=1,nval
          mat_overlap(ival,jval) = 0.0d0
        enddo
      enddo
c
      do iel=1,nel
        typ_e = type_el(iel)
        do j=1,nnodes_0
          j1 = table_nod(j,iel)
          nod_el_p(j) = j1
        enddo
        xmin = x_P2(nod_el_p(1))
        xmax = x_P2(nod_el_p(2))
c       Matrix of overlap integrals betwween the Lagrange polynomials
        call matrix_mxx_1d (xmin, xmax, mat_Mxx_0)
        call matrix_kxy_1d (xmin, xmax, mat_Kxy_0)
        call matrix_kyy_1d (xmin, xmax, mat_Kyy_0)
        call matrix_myy_1d (xmin, xmax, mat_Myy_0)
        call matrix_kyx_1d (xmin, xmax, mat_Kyx_0)
c       Block column 1  -------------   C1
        do j=1,node_P2
          do i=1,node_P2
            mat_el_2(i,j) = pp(typ_e) * mat_Mxx_0(i,j)    ! M_xx
          enddo
          do i=1,node_P3
            mat_el_2(i+node_P2,j) = 0
          enddo
          do i=1,node_P3
            mat_el_2(i+node_P2+node_P3,j) = 0
          enddo
        enddo
c       Block column 2  -------------   C2
        do j=1,node_P3
          do i=1,node_P2
            mat_el_2(i,j+node_P2) = 0
          enddo
        enddo
        do j=1,node_P3
          do i=1,node_P3
            mat_el_2(i+node_P2,j+node_P2) = pp(typ_e) * mat_Myy_0(i,j)  ! M_yy
          enddo
        enddo
        do j=1,node_P3
          do i=1,node_P3
            mat_el_2(i+node_P2+node_P3,j+node_P2) = 0
          enddo
        enddo
c       Block column 3  -------------   C3
        do j=1,node_P3
          do i=1,node_P2
            mat_el_2(i,j+node_P2+node_P3) = - pp(typ_e)  ! M_xz
     *                            * mat_Kxy_0(i,j)
          enddo
        enddo
        do j=1,node_P3
          do i=1,node_P3
            mat_el_2(i+node_P2,j+node_P2+node_P3) = - ii * bloch_vec_y  ! M_yz
     *                           * pp(typ_e) * mat_Myy_0(i,j)
          enddo
        enddo
        do j=1,node_P3
          do i=1,node_P3
            mat_el_2(i+node_P2+node_P3,j+node_P2+node_P3) = 0
          enddo
        enddo
        do ival=1,nval
          do i=1,nddl_0
            z_tmp1 = sol_2(i,ival,iel)
            sol_el_2(i) = z_tmp1
          enddo
          do jval=1,nval
            z_beta_1 = beta_1(jval)
            do i=1,nddl_0
              z_tmp1 = sol_1(i,jval,iel)
              sol_el_1(i) = z_tmp1 * z_beta_1
            enddo
c           Matrix-Vector product
            do i=1,nddl_0
              vec_1(i) = 0.0d0
              do j=1,nddl_0
                z_tmp1 = sol_el_1(j)
                z_tmp2 = mat_el_2(i,j)
                vec_1(i) = vec_1(i) + z_tmp1 * z_tmp2
              enddo
            enddo
c           Scalar product
            z_tmp1 = 0.0d0
            do i=1,nddl_0
              z_tmp1 = vec_1(i) * sol_el_2(i)
              z_tmp1 = z_tmp1 / k_0
              mat_overlap(ival,jval) = mat_overlap(ival,jval)
     *              + z_tmp1
            enddo
          enddo
        enddo
      enddo

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  reorder complex conjugate douplets
      j = 1
      skip = 0
      redo = redo + 1
      if (redo .eq. 2) goto 123

121   continue
      if (j .gt. nval) then 
      goto 122 
C       if all is well - correct orthogonality with its self
      elseif (abs(mat_overlap(j,j)) .gt. 1.0d-7) then
        j = j+1
        goto 121
C       first of a wrongly ordered complex conjugate pair (save values)
      elseif (skip .eq. 0) then
        if (j .eq. nval) then
          pair_warning = pair_warning + 1
        endif
C       find jvals (j and j2) of swaped pair and switch
        do j2 = j+1,nval
          if (abs(mat_overlap(j2,j)) .gt. 1.0d-7) then
          betatmp1(1) = beta_2(j)
          betatmp2(1) = beta_2(j2)
          beta_2(j)  = betatmp2(1)
          beta_2(j2) = betatmp1(1)
ccc            do compcount = 1,3
              do elcount = 1,nel
                do nodecount = 1,nddl_0
              soltmp1(nodecount,elcount,1)
     *          = sol_2(nodecount,j,elcount)
              soltmp2(nodecount,elcount,1)
     *          = sol_2(nodecount,j2,elcount)
              sol_2(nodecount,j,elcount)
     *          = soltmp2(nodecount,elcount,1)
              sol_2(nodecount,j2,elcount)
     *          = soltmp1(nodecount,elcount,1)
               enddo
ccc             enddo
           enddo
          endif  
        enddo
        skip = 1
        j = j+1
        goto 121
C  dont touch second half of pair (as already switched)
      elseif (skip .gt. 0) then
        skip = 0
        j = j+1
        goto 121 
      endif
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
123   continue
      if (PrintAll .eq. 1) then
      open(3,file=overlap_file)
      val_max_diag = 0.0d0
      val_max_off = 0.0d0
      do jval=1,nval
        do ival=1,nval
          r_tmp1 = abs(mat_overlap(ival,jval))
          if (ival .eq. jval) then ! Diagonal
            if (val_max_diag .lt. r_tmp1) then
              val_max_diag = r_tmp1
            endif
          else
            if (val_max_off .lt. r_tmp1) then
            val_max_off = r_tmp1
            endif
          endif
          write(3,12) ival, jval, mat_overlap(ival,jval),
     *      abs(mat_overlap(ival,jval)),
     *      abs(beta_1(ival)-beta_2(jval)), beta_1(ival), beta_2(jval),
     *      abs(mat_overlap(ival,jval) - 
     *      conjg(mat_overlap(jval,ival)))
        enddo
      enddo
      write(3,*) "val_max_diag, val_max_off = ", 
     *      val_max_diag, val_max_off
      close(3)
12    format(2(I4),2(g25.17), 2(g16.8), 8(g18.10))
      endif
C
C
      return
      end
