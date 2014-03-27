C	  Calculate the K overlap integral needed for checking
C           the completeness of the PW and BM basis.
C
      subroutine K_overlap (nval, nel, npt, nnodes, 
     *  nb_typ_el, type_el, table_nod, x, sol, pp, qq, 
     *  lambda, freq, K_overlap_mat, neq_PW,
     *  lat_vecs, bloch_vec, beta1, index_pw_inv,
     *  PrintAll, k_0, ordre_ls)
C		 
      implicit none 
C
      double precision lat_vecs(2,2), bloch_vec(2), k_0
C  input output parameters
      integer*8 nval, nel, npt, nnodes, nb_typ_el
      integer*8 type_el(nel)
      integer*8 table_nod(nnodes,nel), nnodes_0
      parameter (nnodes_0 = 6)
      complex*16 x(2,npt), sol(3,nnodes+7,nval,nel)
      complex*16 sol_el_1(2*nnodes_0+10)
      complex*16 pp(nb_typ_el), qq(nb_typ_el)
      integer*8 neq_PW, PrintAll
      complex*16 beta1(nval)
C  local parameters - purely internal
      integer*8 iel, inode, n, i, j, j1, i_eq
      integer*8 l, ltest, typ_e, s, s2, twos, trans
      integer*8 ui, jval, px, py, ordre_ls
      integer*8 jtest, ind_jp, j_eq, E_K
      integer*8 debug, info_curved, n_curved
      integer*8 nod_el_p(nnodes_0)
      double precision xel(2,nnodes_0), r_tmp1
      complex*16 z_val_exp, z_coeff_1
      double precision phi2_list(6), grad2_mat0(2,6), grad2_mat(2,6)
      double precision phi3_list(10), grad3_mat0(2,10), grad3_mat(2,10)
      double precision vec_phi_j(2)
      complex*16 z_vec_phi_i(2,2)
      complex*16 z_tmp1, z_tmp2K, z_tmp2E, z_beta_1
C  variables for quadrature interpolation
      integer*8 nquad, nquad_max, iq
      parameter (nquad_max = 25)
      double precision wq(nquad_max)
      double precision xq(nquad_max), yq(nquad_max)
      double precision xx(2), xx_g(2), ww, det
      double precision mat_B(2,2), mat_T(2,2), d
C
      integer alloc_stat
      complex*16, dimension(:,:), allocatable :: PlaneW_RK
      complex*16, dimension(:,:), allocatable :: PlaneW_RE
      complex*16, dimension(:,:), allocatable :: overlap_K
      complex*16, dimension(:,:), allocatable :: overlap_E
      complex*16, dimension(:,:,:), allocatable :: mat_scal
C
C      integer*8 nval_max, PW_max
C      parameter (nval_max = 2000)
C      parameter (PW_max = 1000)
C      complex*16 PlaneW_RK(PW_max,2)
C      complex*16 PlaneW_RE(PW_max,2)
C      complex*16 overlap_K(PW_max,nval_max)
C      complex*16 overlap_E(PW_max,nval_max)
C      complex*16 mat_scal(PW_max,2*nnodes_0+10,2)
      complex*16 K_overlap_mat(nval,2*neq_PW)
      double precision ZERO, ONE
      parameter ( ZERO = 0.0D0, ONE = 1.0D0)
      complex*16 ii
      integer*8 index_pw_inv(neq_PW)
      double precision r_tmp
      double precision bloch1, bloch2, pi, alpha, beta, norm
      double precision lambda, freq, k, vec_kx, vec_ky
C
CCCCCCCCCCC Start Program CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
	    ui = 6
      debug = 0
C
      if ( nnodes .ne. 6 ) then
        write(ui,*) "K_overlap: problem nnodes = ", nnodes
        write(ui,*) "K_overlap: nnodes should be equal to 6!"
        write(ui,*) "K_overlap: Aborting..."
        stop
      endif
      allocate(PlaneW_RK(neq_PW,2), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "K_overlap: Mem. allocation is unseccesfull"
        write(*,*) "alloc_stat (PlaneW_RK) = ", alloc_stat
        write(*,*) "Aborting..."
        stop
      endif
      allocate(PlaneW_RE(neq_PW,2), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "K_overlap: Mem. allocation is unseccesfull"
        write(*,*) "alloc_stat (PlaneW_RE) = ", alloc_stat
        write(*,*) "Aborting..."
        stop
      endif
      allocate(overlap_K(neq_PW,nval), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "K_overlap: Mem. allocation is unseccesfull"
        write(*,*) "alloc_stat (overlap_K) = ", alloc_stat
        write(*,*) "Aborting..."
        stop
      endif
      allocate(overlap_E(neq_PW,nval), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "K_overlap: Mem. allocation is unseccesfull"
        write(*,*) "alloc_stat (overlap_E) = ", alloc_stat
        write(*,*) "Aborting..."
        stop
      endif
      allocate(mat_scal(neq_PW,2*nnodes_0+10,2), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "K_overlap: Mem. allocation is unseccesfull"
        write(*,*) "alloc_stat (mat_scal) = ", alloc_stat
        write(*,*) "Aborting..."
        stop
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
C
      call quad_triangle(nquad,nquad_max,wq,xq,yq) 
C
C	set up final solution matrix
        do i=1,neq_PW
          do j=1,nval
            overlap_K(i,j) = 0.0D0
            overlap_E(i,j) = 0.0D0			
          enddo
        enddo 
C
c       write(ui,*)  "K_overlap: bloch_vec = ", bloch_vec
C
      ii = dcmplx(0.0d0, 1.0d0)
      d = lat_vecs(1,1)  
      pi = 3.141592653589793d0
      bloch1 = bloch_vec(1)
      bloch2 = bloch_vec(2)
      vec_kx = 2.0d0*pi/d
      vec_ky = 2.0d0*pi/d
C
CCCCCCCCCCCCCCCCC	loop over all elements	CCCCCCCCCCCCCCCC
C
      n_curved = 0
      do iel=1,nel
        typ_e = type_el(iel)
        do j=1,nnodes
          j1 = table_nod(j,iel)
          nod_el_p(j) = j1
          xel(1,j) = x(1,j1)
          xel(2,j) = x(2,j1)
        enddo
        call curved_elem_tri (nnodes, xel, info_curved, r_tmp1)
        if (info_curved .eq. 1) then
          n_curved = n_curved + 1
        endif
cccccccccc
        do i=1,neq_pw
          do j=1,2*nnodes+10
            do E_K = 1,2
              mat_scal(i,j,E_K) = 0.0d0
            enddo
          enddo
        enddo
cccccccccc
        do iq=1,nquad
          xx(1) = xq(iq)
          xx(2) = yq(iq)
          ww = wq(iq)
c         xx   = coordinate on the reference triangle
c         xx_g = coordinate on the actual triangle
c         We will also need the gradients of the P1 element
           call phi2_2d_mat(xx, phi2_list, grad2_mat0)
           call phi3_2d_mat(xx, phi3_list, grad3_mat0)
c         
          if (info_curved .eq. 0) then
c           Rectilinear element
            call jacobian_p1_2d(xx, xel, nnodes, 
     *               xx_g, det, mat_B, mat_T)
            if (det .le. 0 .and. debug .eq. 1 .and. iq .eq. 1) then
              write(*,*) "   !!!"
              write(*,*) "K_overlap: det <= 0: iel, det ", iel, det
            endif
          else
c           Isoparametric element
            call jacobian_p2_2d(xx, xel, nnodes, phi2_list, 
     *               grad2_mat0, xx_g, det, mat_B, mat_T)
          endif
C
           if(abs(det) .lt. 1.0d-10) then
             write(*,*)
             write(*,*) "   ???"
             write(*,*) "K_overlap: det = 0 : iel, det = ", iel, det
             write(*,*) "K_overlap: Aborting..."
             stop
           endif
c          grad_i  = gradient on the actual triangle
c          grad_i  = Transpose(mat_T)*grad_i0
c          Calculation of the matrix-matrix product:
          call DGEMM('Transpose','N', 2, 6, 2, ONE, mat_T, 2,
     *      grad2_mat0, 2, ZERO, grad2_mat, 2)
          call DGEMM('Transpose','N', 2, 10, 2, ONE, mat_T, 2,
     *      grad3_mat0, 2, ZERO, grad3_mat, 2)
          z_coeff_1 = ww * abs(det) * pp(typ_e)
C
C     prefactors of plane waves order s(px,py) 
            s = 1
            do px = -ordre_ls, ordre_ls
              do py = -ordre_ls, ordre_ls
              if (px**2 + py**2 .le. ordre_ls**2) then
              alpha = bloch1 + vec_kx*px	! Bloch vector along x
              beta  = bloch2 + vec_ky*py	! Bloch vector along y
              norm = 1.0d0/DSQRT(alpha**2 + beta**2)	! sqrt term
              r_tmp = alpha*xx_g(1) + beta*xx_g(2)
              z_val_exp = EXP(ii*r_tmp)*norm
              s2 = index_pw_inv(s)
C
C     Plane waves order s
C  	RK    x component
	        PlaneW_RK(s2,1) = alpha*z_val_exp
C           y component
              PlaneW_RK(s2,2) = beta*z_val_exp
C  	RE	x component
              PlaneW_RE(s2,1) = beta*z_val_exp
C           y component
              PlaneW_RE(s2,2) = -alpha*z_val_exp
CCCis it better to lose s2 from mat_scal and just reset to zero after each s loop?
              do i_eq=1,2
C                do i=1,2
C                  do j=1,2
C                    z_vec_phi_i(i,j) = 0.0d0
C                  enddo
C                enddo
                z_vec_phi_i(i_eq,1) = PlaneW_RK(s2,i_eq)
                z_vec_phi_i(i_eq,2) = PlaneW_RE(s2,i_eq)
              enddo
c             Determine the basis vector
                do jtest=1,nnodes_0
                  do j_eq=1,2
                    ind_jp = j_eq + 2*(jtest-1)
c                   Determine the basis vector
                    do i=1,2
                      vec_phi_j(i) = 0.0d0
                    enddo
                    vec_phi_j(j_eq) = phi2_list(jtest)
                    do E_K = 1,2
                      z_tmp1 = z_vec_phi_i(1,E_K)*vec_phi_j(1) + 
     *                            z_vec_phi_i(2,E_K)*vec_phi_j(2)
                      z_tmp1 = (z_coeff_1/k_0) * z_tmp1
                      mat_scal(s2,ind_jp,E_K) = 
     *                  mat_scal(s2,ind_jp,E_K) + z_tmp1
                    enddo
                  enddo
                enddo
                do jtest=1,10
                  j_eq = 3
                  ind_jp = jtest + 2*nnodes_0
c                 Determine the basis vector
                  do i=1,2
                    vec_phi_j(i) = -grad3_mat(i,jtest)
                  enddo
                  do E_K = 1,2
                    z_tmp1 = z_vec_phi_i(1,E_K)*vec_phi_j(1) + 
     *                             z_vec_phi_i(2,E_K)*vec_phi_j(2)
                    z_tmp1 = (z_coeff_1/k_0) * z_tmp1
                    mat_scal(s2,ind_jp,E_K) = 
     *                mat_scal(s2,ind_jp,E_K) + z_tmp1
                  enddo
                enddo
              s = s + 1
            endif
            enddo             ! - plane waves 1
          enddo               ! - plane waves 2
        enddo                 ! - quad points
CCCCCCCCCCCCCC
            do s2=1,neq_PW 
              do jval=1,nval
                z_beta_1 = beta1(jval)
                do i=1,nnodes
                  do j=1,2
c                   The 2 transverse components of the mode jval
                    ind_jp = j + 2*(i-1)
                    z_tmp1 = sol(j,i,jval,iel)
                    sol_el_1(ind_jp) = z_tmp1 * z_beta_1
                  enddo
                enddo
CCCCCCCCCCCCCC
                do i=1,3
c                 The longitudinal component at the vertices (P3 elements)
                  ind_jp = i + 2*nnodes
                  z_tmp1 = sol(3,i,jval,iel)
                  sol_el_1(ind_jp) = z_tmp1 * z_beta_1
                enddo
                do i=nnodes+1,13
c                 The longitudinal component at the edge nodes and interior node (P3 elements)
                  ind_jp = i + 2*nnodes - nnodes + 3
                  z_tmp1 = sol(3,i,jval,iel)
                  sol_el_1(ind_jp) = z_tmp1 * z_beta_1
                enddo
CCCCCCCCCCCCCC
                do j=1,2*nnodes+10
                  z_tmp1 = sol_el_1(j)
                  z_tmp2K = mat_scal(s2,j,1)
                  z_tmp2E = mat_scal(s2,j,2)
                  overlap_K(s2,jval) = overlap_K(s2,jval) 
     *                   + z_tmp1 * z_tmp2K
                  overlap_E(s2,jval) = overlap_E(s2,jval) 
     *                   + z_tmp1 * z_tmp2E
                enddo
              enddo           ! - j eigenvalues
            enddo             ! - plane waves 1 and 2
c              s = s + 1
c            enddo             ! - plane waves 1
c          enddo               ! - plane waves 2
c        enddo                 ! - quad points
      enddo                   ! - elements
C
C  Save as K_mat = | K_E |                
C                  | K_K |
C
      do twos = 1, 2*neq_PW
        do n = 1, nval
          if (twos .le. neq_PW) then
            K_overlap_mat(n,twos) = overlap_E(twos,n)
          else
            K_overlap_mat(n,twos) = overlap_K(twos-neq_PW,n)
          endif
        enddo
      enddo
C
CCCCCCCCCCCCCCCCC	  save results   CCCCCCCCCCCCCCCC
C
      if (PrintAll .eq. 1) then
C
      open (unit=35, file="Normed/K_mat.txt", status='unknown')
      open (unit=33, file="Normed/K_K.txt", status='unknown')
      open (unit=32, file="Normed/K_E.txt", status='unknown')
      write(35,131) lambda, freq     
      write(33,131) lambda, freq   
      write(32,131) lambda, freq
C
      do n = 1, nval
        do twos = 1, 2*neq_PW
          write(35,132)  n,twos,K_overlap_mat(n,twos),
     *         abs(K_overlap_mat(n,twos))
        enddo
      enddo
C
      do n=1,nval
        do s=1,neq_PW
          write(33,132) s, n, overlap_K(s,n), abs(overlap_K(s,n))
          write(32,132) s, n, overlap_E(s,n), abs(overlap_E(s,n))
        enddo
      enddo
C
131    format(2(f12.4))
132    format(2(I4),2(g25.17),g18.10)
C
      close(35)
      close(33)
      close(32) 
      endif
C
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      deallocate(PlaneW_RK)
      deallocate(PlaneW_RE)
      deallocate(overlap_K)
      deallocate(overlap_E)
      deallocate(mat_scal)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      return
      end 
