C  Calculate the Projection integral of the conjugate Plane Waves & Bloch Modes
C
      subroutine J_overlap(nval, nel, npt, nnodes,
     *  type_el, table_nod, x, sol,
     *  lat_vecs, lambda, freq, overlap_J, neq_PW,
     *  bloch_vec, index_pw_inv, PrintAll,
     *  ordre_ls)
C
      implicit none
C  input output parameters
      integer*8 nval, nel, npt, nnodes
      integer*8 type_el(nel)
      integer*8 table_nod(nnodes,nel)
      double precision x(2,npt)
      complex*16 sol(3,nnodes+7,nval,nel)
      double precision lat_vecs(2,2)
      integer*8 neq_PW, PrintAll, ordre_ls
      double precision bloch_vec(2)
C index_pw(neq_PW)
      integer*8 index_pw_inv(neq_PW)
C  local parameters - purely internal
      integer*8 iel, inode, n, i, j, global
      integer*8 ltest, typ_e, s, s2, twos, trans
      integer*8 ui, px, py, nnodes_0
      parameter (nnodes_0 = 6)
      integer*8 nod_el_p(nnodes_0)
      double precision xel(2,nnodes_0)
      double precision phi2_list(6)
      complex*16 vec_phi(2), K_tmp(2), E_tmp(2)
      complex*16 tmp1, tmp2
C  variables for quadrature interpolation
      integer*8 nquad, nquad_max, iq
      parameter (nquad_max = 25)
      double precision wq(nquad_max)
      double precision xq(nquad_max), yq(nquad_max)
      double precision xx(2), xx_g(2), ww, det
      double precision mat_B(2,2), mat_T(2,2), d, dy
C
      integer alloc_stat
      complex*16, dimension(:,:), allocatable :: PlaneW_RK
      complex*16, dimension(:,:), allocatable :: PlaneW_RE
      complex*16, dimension(:,:), allocatable :: overlap_K
      complex*16, dimension(:,:), allocatable :: overlap_E
C
      complex*16 overlap_J(2*neq_PW,nval)
      complex*16 ii, val_exp, coeff_1
      double precision r_tmp, vec_kx, vec_ky, lambda, freq
      double precision bloch1, bloch2, pi, alpha, beta, norm
C
CCCCCCCCCCC Start Program CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      ui = 6
C
      if ( nnodes .ne. 6 ) then
        write(ui,*) "overlap_J: problem nnodes = ", nnodes
        write(ui,*) "overlap_J: nnodes should be equal to 6!"
        write(ui,*) "overlap_J: Aborting..."
        stop
      endif
      allocate(PlaneW_RK(neq_PW,2), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "J_overlap: Mem. allocation is unseccesfull"
        write(*,*) "alloc_stat (PlaneW_RK) = ", alloc_stat
        write(*,*) "Aborting..."
        stop
      endif
      allocate(PlaneW_RE(neq_PW,2), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "J_overlap: Mem. allocation is unseccesfull"
        write(*,*) "alloc_stat (PlaneW_RE) = ", alloc_stat
        write(*,*) "Aborting..."
        stop
      endif
      allocate(overlap_K(neq_PW,nval), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "J_overlap: Mem. allocation is unseccesfull"
        write(*,*) "alloc_stat (overlap_K) = ", alloc_stat
        write(*,*) "Aborting..."
        stop
      endif
      allocate(overlap_E(neq_PW,nval), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "J_overlap: Mem. allocation is unseccesfull"
        write(*,*) "alloc_stat (overlap_E) = ", alloc_stat
        write(*,*) "Aborting..."
        stop
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
C
      call quad_triangle(nquad,nquad_max,wq,xq,yq);
C
C set up final solution matrix
      do j=1,nval
        do i=1,neq_PW
          overlap_K(i,j) = 0.0D0
          overlap_E(i,j) = 0.0D0
        enddo
      enddo
C
      ii = dcmplx(0.0d0, 1.0d0)
      d = lat_vecs(1,1)
      dy = lat_vecs(2,2)
      pi = 3.141592653589793d0
      bloch1 = bloch_vec(1)
      bloch2 = bloch_vec(2)
      vec_kx = 2.0d0*pi/d
      vec_ky = 2.0d0*pi/dy
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCC

C
C
CCCCCCCCCCCCCCCCC	loop over all elements	CCCCCCCCCCCCCCCC
C
      do iel=1,nel
        typ_e = type_el(iel)
        do inode=1,nnodes
          global = table_nod(inode,iel)
          nod_el_p(inode) = global
          xel(1,inode) = x(1,global)
          xel(2,inode) = x(2,global)
        enddo
C
        do iq=1,nquad
          xx(1) = xq(iq)
          xx(2) = yq(iq)
          ww = wq(iq)
          call phi2_2d_mat_J(xx, phi2_list)
          call jacobian_p1_2d(xx, xel, nnodes,
     *               xx_g, det, mat_B, mat_T)
C
C     prefactors of plane waves order s(px,py)
            s = 1
            do px = -ordre_ls, ordre_ls
              do py = -ordre_ls, ordre_ls
              if ((px/d)**2 + (py/dy)**2
     *           .le. (ordre_ls/MAX(d, dy))**2) then
                alpha = bloch1 + vec_kx*px      ! Bloch vector along x
                beta  = bloch2 + vec_ky*py      ! Bloch vector along y
                norm = 1.0d0/SQRT(alpha**2 + beta**2)  ! sqrt term
                r_tmp = alpha*xx_g(1) + beta*xx_g(2)
                val_exp = EXP(-ii*r_tmp)*norm
                s2 = index_pw_inv(s)
C
C     Plane waves order s
C   RK    x component
                PlaneW_RK(s2,1) = alpha*val_exp
C         y component
                PlaneW_RK(s2,2) = beta*val_exp
C   RE    x component
                PlaneW_RE(s2,1) = beta*val_exp
C         y component
                PlaneW_RE(s2,2) = -alpha*val_exp
                s = s + 1
              endif
              enddo
            enddo
C
            coeff_1 = ww * ABS(det)      !!!!!!!!!!!!!!!!   * pp(typ_e)
            do s=1,neq_PW
              do trans=1,2              ! transverse field components
                K_tmp(trans) = PlaneW_RK(s,trans)
                E_tmp(trans) = PlaneW_RE(s,trans)
              enddo
              do n=1,nval
                do trans=1,2
                  vec_phi(trans) = 0.0d0
                enddo
                do ltest=1,nnodes
                  do trans=1,2          ! transverse field components
                    vec_phi(trans) = vec_phi(trans) +
     *                  sol(trans,ltest,n,iel) * phi2_list(ltest)
                  enddo
                enddo
C take dot product of Bloch Lp and PW for each quad pt.
                   tmp1 = vec_phi(1)*K_tmp(1) + vec_phi(2)*K_tmp(2)
                   overlap_K(s,n) = overlap_K(s,n) + coeff_1 * tmp1
C
                   tmp2 = vec_phi(1)*E_tmp(1) + vec_phi(2)*E_tmp(2)
                   overlap_E(s,n) = overlap_E(s,n) + coeff_1 * tmp2
              enddo
            enddo
        enddo
      enddo
C
C  Save as J_mat = | J_E |
C                  | J_K |
      do n = 1, nval
        do twos = 1, 2*neq_PW
          if (twos .le. neq_PW) then
            overlap_J(twos,n) = overlap_E(twos,n)
          else
            overlap_J(twos,n) = overlap_K(twos-neq_PW,n)
          endif
        enddo
      enddo
C
CCCCCCCCCCCCCCCCC	  save results   CCCCCCCCCCCCCCCC
C
      if (PrintAll .eq. 1) then
C
      open (unit=35, file="Matrices/J_mat.txt", status='unknown')
      open (unit=33, file="Matrices/J_K.txt", status='unknown')
      open (unit=32, file="Matrices/J_E.txt", status='unknown')
C      write(34,131) lambda, freq
      write(35,131) lambda, freq
      write(33,131) lambda, freq
      write(32,131) lambda, freq
C
      do n = 1, nval
        do twos = 1, 2*neq_PW
          write(35,132) twos,n,overlap_J(twos,n),
     *         abs(overlap_J(twos,n))
        enddo
      enddo
      do n=1,nval
        do s=1,neq_PW
          write(33,132) s, n, overlap_K(s,n), abs(overlap_K(s,n))
          write(32,132) s, n, overlap_E(s,n), abs(overlap_E(s,n))
        enddo
      enddo
C
      close(35)
      close(34)
      close(33)
      close(32)
C
131   format(2(f12.4))
132   format(2(I4),2(g25.17),g18.10)
      endif
C
C
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      deallocate(PlaneW_RK)
      deallocate(PlaneW_RE)
      deallocate(overlap_K)
      deallocate(overlap_E)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      return
      end
