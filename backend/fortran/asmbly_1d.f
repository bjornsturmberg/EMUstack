

      subroutine asmbly_1d (nel, npt_P2, n_ddl, neq, 
     *  shift, bloch_vec_x, bloch_vec_y, nb_typ_el, 
     *  pp, qq, table_nod, table_ddl, type_el, ineq, 
     *  ip_period_ddl, x_P2, x_ddl, matrix_1, matrix_2)
C   
      implicit none
      integer*8 nel, npt_P2, n_ddl, neq, nb_typ_el
      integer*8 table_nod(3,nel), type_el(nel)
      integer*8 table_ddl(3+4+4,nel), ineq(n_ddl)
      integer*8 ip_period_ddl(n_ddl)
      complex*16 matrix_1(neq,neq), matrix_2(neq,neq)
      complex*16 pp(nb_typ_el), qq(nb_typ_el)
      double precision x_P2(npt_P2), x_ddl(n_ddl)
      complex*16 shift
      double precision bloch_vec_x, bloch_vec_y
C
      integer*8 node_P2, node_P3
      parameter (node_P2 = 3, node_P3 = 4)
      double precision mat_Mxx_0(node_P2,node_P2)
      double precision mat_Kxy_0(node_P2,node_P3)
      double precision mat_Kyy_0(node_P3,4)
      double precision mat_Myy_0(4,4)
      double precision mat_Kyx_0(4,node_P2)
C
      complex*16 mat_el_1(node_P2+2*node_P3,node_P2+2*node_P3)
      complex*16 mat_el_2(node_P2+2*node_P3,node_P2+2*node_P3)
C
      integer allocate_status
      complex*16, dimension(:,:,:), allocatable :: ls_mat_el_1
      complex*16, dimension(:,:,:), allocatable :: ls_mat_el_2
C
      integer*8 nnodes_0, ui, debug, typ_e, nddl_0
      parameter (nnodes_0 = 3, nddl_0 = 11)
      integer*8 ls_nodes(3)
      double precision xmin, xmax
      double precision period_x, delta_x
      complex*16 val_exp_1, val_exp(nddl_0), z_phase_fact
C
      integer*8 dim1, dim2, dim3, info_periodic
      integer*8 i, j, iel, j1
      integer*8 ip, jp, itrial, jtest, eq_ip, eq_jp
      complex*16 z_tmp1, z_tmp2, ii
c
ccccccccccccccccccccccccccccccccccccccc
c
c  ii = sqrt(-1)
      ii = dcmplx(0.0d0, 1.0d0)
      ui = 6
      debug = 0

      do eq_ip=1,neq
        do eq_jp=1,neq
          matrix_1(eq_ip,eq_jp) = 0
          matrix_2(eq_ip,eq_jp) = 0
        enddo
      enddo

      allocate_status = 0


      dim1 = node_P2+2*node_P3
      dim2 = node_P2+2*node_P3
      dim3 = nel

      if (debug .eq. 1) then
        write(*,*) "asmbly_1d: dim1, dim2, dim3 = ", dim1, dim2, dim3
      endif

      allocate(ls_mat_el_1(dim1,dim2,dim3), STAT=allocate_status)
      if (allocate_status /= 0) then
        write(*,*) "The allocation is unsuccessful"
        write(*,*) "allocate_status = ", allocate_status
        write(*,*) "Not enough memory for the array ls_mat_el_1"
        write(*,*) "dim1,dim2,dim3 = ", dim1,dim2,dim3
        write(*,*) "Aborting..."
        stop
      endif

      allocate(ls_mat_el_2(dim1,dim2,dim3), STAT=allocate_status)
      if (allocate_status /= 0) then
        write(*,*) "The allocation is unsuccessful"
        write(*,*) "allocate_status = ", allocate_status
        write(*,*) "Not enough memory for the array ls_mat_el_2"
        write(*,*) "dim1,dim2,dim3 = ", dim1,dim2,dim3
        write(*,*) "Aborting..."
        stop
      endif
c
ccccccccccccccccccccccccccccccccccccccc
c
      period_x = 0
      do iel=1,nel
        j=4
          ip = table_ddl(j,iel)
          j1 = ip_period_ddl(ip)
          if (j1 .ne. 0 .and. j1 .ne. ip) then
            period_x = x_ddl(ip) - x_ddl(j1)
          endif
c        j=7
        j=5
          ip = table_ddl(j,iel)
          j1 = ip_period_ddl(ip)
          if (j1 .ne. 0 .and. j1 .ne. ip) then
            period_x = x_ddl(ip) - x_ddl(j1)
          endif      
      enddo

      val_exp_1 = exp(ii * bloch_vec_x * period_x)
      if (debug .eq. 1) then
        write(ui,*) "asmbly_1d: Period of the unit cell = ", 
     *              period_x
        write(ui,*) "asmbly_1d: val_exp_1 = ",  val_exp_1
      endif
c      npt_P2 = 3 * nel ! Discontinuous P2 polynomial
c      npt_P3 = 2 * nel + nel + 1
c      neq = npt_P2 + 2 * (npt_P3 -1)
c                    !!!!!!!!!!!!!!!! 
      do iel=1,nel   !!!!!!!!!!!!!!!! BEGIN MAIN LOOP
        typ_e = type_el(iel)
        info_periodic = 0  !  info_periodic = 1 if element number iel has a periodic node
        do i=1,3
          ls_nodes(i) = table_nod(i,iel)
        enddo
        do j=1,nddl_0
          val_exp(j) = 1.0d0
        enddo
c       val_exp: Bloch mode phase factor between the origin point and destination point
c       For a pair of periodic points, one is chosen as origin and the other is the destination
        do j=1,nddl_0
          ip = table_ddl(j,iel)
          j1 = ip_period_ddl(ip)
          if (j1 .ne. 0 .and. j1 .ne. ip) then
            info_periodic = 1
            delta_x = x_ddl(ip) - x_ddl(j1)
            val_exp(j) = exp(ii * bloch_vec_x * delta_x)
            if(delta_x /= period_x) then
              write(ui,*) "asmbly_1d: ??? delta_x /= period_x: ",
     *              delta_x, period_x
              write(ui,*) "asmbly_1d: j, ip, j1, iel = ",
     *              j, ip, j1, iel
            endif
          endif
        enddo
c
        xmin = x_P2(ls_nodes(1))
        xmax = x_P2(ls_nodes(2))
        if (debug .eq. 1) then
          if (xmin >= xmax) then
            write(ui,*) "asmbly_1d: ??? xmin >= xmax: ",
     *             xmin, xmax
          endif
        endif
c       Matrix of overlap integrals betwween the Lagrange polynomials
        call matrix_mxx_1d (xmin, xmax, mat_Mxx_0)
        call matrix_kxy_1d (xmin, xmax, mat_Kxy_0)
        call matrix_kyy_1d (xmin, xmax, mat_Kyy_0)
        call matrix_myy_1d (xmin, xmax, mat_Myy_0)
        call matrix_kyx_1d (xmin, xmax, mat_Kyx_0)
c       Block column 1  -------------   C1
        do j=1,node_P2
          do i=1,node_P2
            z_tmp1 = pp(typ_e) * bloch_vec_y**2 - qq(typ_e)
            mat_el_1(i,j) =   z_tmp1 * mat_Mxx_0(i,j)       ! K_xx
            mat_el_2(i,j) = - pp(typ_e) * mat_Mxx_0(i,j)    ! M_xx
          enddo
          do i=1,node_P3
            mat_el_1(i+node_P2,j) = - ii * bloch_vec_y * pp(typ_e)  ! K_yx
     *                            * mat_Kyx_0(i,j)
            mat_el_2(i+node_P2,j) = 0
          enddo
          do i=1,node_P3
            mat_el_1(i+node_P2+node_P3,j) = - pp(typ_e)  ! K_zx
     *                            * mat_Kyx_0(i,j)
            mat_el_2(i+node_P2+node_P3,j) = 0
          enddo
        enddo
c       Block column 2  -------------   C2
        do j=1,node_P3
          do i=1,node_P2
            mat_el_1(i,j+node_P2) = ii * bloch_vec_y * pp(typ_e)  ! K_xy
     *                            * mat_Kxy_0(i,j)
            mat_el_2(i,j+node_P2) = 0
          enddo
        enddo
        do j=1,node_P3
          do i=1,node_P3
            z_tmp1 = pp(typ_e) * mat_Kyy_0(i,j)
            z_tmp2 = qq(typ_e) * mat_Myy_0(i,j)
            mat_el_1(i+node_P2,j+node_P2) = z_tmp1 - z_tmp2   ! K_yy
            mat_el_2(i+node_P2,j+node_P2) = - pp(typ_e) * mat_Myy_0(i,j)  ! M_yy
          enddo
        enddo
        do j=1,node_P3
          do i=1,node_P3
            mat_el_1(i+node_P2+node_P3,j+node_P2) = ii * bloch_vec_y   ! K_zy
     *                           * pp(typ_e) * mat_Myy_0(i,j)
            mat_el_2(i+node_P2+node_P3,j+node_P2) = 0
          enddo
        enddo
c       Block column 3  -------------   C3
        do j=1,node_P3
          do i=1,node_P2
            mat_el_1(i,j+node_P2+node_P3) = 0
            mat_el_2(i,j+node_P2+node_P3) = pp(typ_e)  ! M_xz
     *                            * mat_Kxy_0(i,j)
          enddo
        enddo
        do j=1,node_P3
          do i=1,node_P3
            mat_el_1(i+node_P2,j+node_P2+node_P3) = 0
            mat_el_2(i+node_P2,j+node_P2+node_P3) = ii * bloch_vec_y  ! M_yz
     *                           * pp(typ_e) * mat_Myy_0(i,j)
          enddo
        enddo
        do j=1,node_P3
          do i=1,node_P3
            z_tmp1 = pp(typ_e) * mat_Kyy_0(i,j)
            z_tmp2 = (pp(typ_e) * bloch_vec_y**2 - qq(typ_e))
     *                            * mat_Myy_0(i,j)
            mat_el_1(i+node_P2+node_P3,j+node_P2+node_P3) =  ! K_zz
     *                             z_tmp1 + z_tmp2
            mat_el_2(i+node_P2+node_P3,j+node_P2+node_P3) = 0
          enddo
        enddo
        dim1 = node_P2+2*node_P3
        dim2 = node_P2+2*node_P3
        do j=1,dim1
          do i=1,dim2
            ls_mat_el_1(i,j,iel) = mat_el_1(i,j)
            ls_mat_el_2(i,j,iel) = mat_el_2(i,j)
          enddo
        enddo
        if (info_periodic /= 0) then
          do jtest=1,nddl_0
            do itrial=1,nddl_0
              z_phase_fact = val_exp(jtest) * conjg(val_exp(itrial))
              ls_mat_el_1(itrial,jtest,iel) = z_phase_fact 
     *                              * ls_mat_el_1(itrial,jtest,iel)
              ls_mat_el_2(itrial,jtest,iel) = z_phase_fact 
     *                              * ls_mat_el_2(itrial,jtest,iel)
            enddo
          enddo
        endif
      enddo   !!!!!!!!!!!!!!!! END MAIN LOOP
c             !!!!!!!!!!!!!!!! 
c
ccccccccccccccccccccccccccccccccccccccc
c
c

      do iel=1,nel   !!!!!!!!!!!!!!!! BEGIN MAIN LOOP (2)
        do jtest=1,nddl_0
          jp = table_ddl(jtest,iel)
          eq_jp = ineq(jp)
          if (eq_jp .gt. 0) then
            do itrial=1,nddl_0
              ip = table_ddl(itrial,iel)
              eq_ip = ineq(ip)
              if (eq_ip .gt. 0) then
                z_tmp1 = ls_mat_el_1(itrial,jtest,iel)
                z_tmp2 = ls_mat_el_2(itrial,jtest,iel)
                matrix_1(eq_ip,eq_jp) =  matrix_1(eq_ip,eq_jp) 
     *                              + z_tmp1 - shift*z_tmp2
              matrix_2(eq_ip,eq_jp) =  matrix_2(eq_ip,eq_jp) 
     *                              + z_tmp2
              endif
            enddo
          endif
        enddo
      enddo   !!!!!!!!!!!!!!!! END MAIN LOOP (2)
c
ccccccccccccccccccccccccccccccccccccccc  TEMPORARY
c

      if (.false.) then
      if (debug .eq. 1) then
        do eq_ip=1,neq
          do eq_jp=1,neq
            matrix_1(eq_ip,eq_jp) = 0.0d0
          enddo
            matrix_1(eq_ip,eq_ip) = eq_ip
        enddo
      endif
      if (debug .eq. 1) then
        do eq_ip=1,neq
          do eq_jp=1,neq
            matrix_2(eq_ip,eq_jp) = 0
          enddo
            matrix_2(eq_ip,eq_ip) = 1
        enddo
      endif
      endif
c
ccccccccccccccccccccccccccccccccccccccc
c

      if (.false.) then
      if (debug .eq. 1) then
        write(ui,*)
        write(ui,*) "Writing matrix_1 into a file"
        open (unit=24,file="Output/matrix_1.txt",status="unknown")
        do eq_ip=1,neq
          do eq_jp=1,neq
            write(24,*) eq_ip, eq_jp, matrix_1(eq_ip,eq_jp) 
          enddo
        enddo
        close(24)
      endif
      if (debug .eq. 1) then
        write(ui,*) "Writing matrix_2 into a file"
        open (unit=24,file="Output/matrix_2.txt",status="unknown")
        do eq_ip=1,neq
          do eq_jp=1,neq
            write(24,*) eq_ip, eq_jp, matrix_2(eq_ip,eq_jp) 
          enddo
        enddo
        close(24)
      endif
      endif
c
ccccccccccccccccccccccccccccccccccccccc
c
      if (.false.) then
      if (debug .eq. 1) then
        write(ui,*)
        write(ui,*) "Writing ls_mat_el_1 into a file"
        open (unit=24,file="Output/ls_mat_el_1.txt",status="unknown")
        do iel=1,nel
        write(24,*) "iel = ", iel
        do itrial=1,nddl_0
          do jtest=1,nddl_0
            write(24,*) itrial, jtest, ls_mat_el_1(itrial,jtest,iel) 
          enddo
        enddo
        enddo
        close(24)
      endif
      if (debug .eq. 1) then
        write(ui,*) "Writing ls_mat_el_2 into a file"
        open (unit=24,file="Output/ls_mat_el_2.txt",status="unknown")
        do iel=1,nel
        write(24,*) "iel = ", iel
        do itrial=1,nddl_0
          do jtest=1,nddl_0
            write(24,*) itrial, jtest, ls_mat_el_2(itrial,jtest,iel) 
          enddo
        enddo
        enddo
        close(24)
      endif
      endif
c
ccccccccccccccccccccccccccccccccccccccc
c

      return
      end
