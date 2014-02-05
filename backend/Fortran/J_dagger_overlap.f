C  Calculate the Projection integral of the Plane Waves & adjoint Bloch modes
C
      subroutine J_dagger_overlap(nval, nel, npt, nnodes, 
     *  type_el, table_nod, x, sol, 
     *  lat_vecs, lambda, freq, overlap_J_dagger, neq_PW,
     *  bloch_vec, index_pw_inv, PrintAll, ordre_ls)
C 
      implicit none 
C  input output parameters
      integer*8 nval, nel, npt, nnodes
      integer*8 type_el(nel)
      integer*8 table_nod(nnodes,nel)
      complex*16 x(2,npt), sol(3,nnodes+7,nval,nel)
      double precision lat_vecs(2,2)
      integer*8 neq_PW, PrintAll
      double precision bloch_vec(2)
C  local parameters - purely internal
      integer*8 iel, inode, n, i, j, global
      integer*8 ltest, typ_e, s, s2, twos, trans
      integer*8 ui, px, py, ordre_ls
      integer*8 nnodes_0, debug
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
      double precision mat_B(2,2), mat_T(2,2), d
C
      integer*8 nval_max, PW_max
      parameter (nval_max = 3000)
      parameter (PW_max = 200)
      complex*16 PlaneW_RK(PW_max,2)
      complex*16 PlaneW_RE(PW_max,2)
      complex*16 overlap_K(PW_max,nval_max)
      complex*16 overlap_E(PW_max,nval_max)
      complex*16 overlap_J_dagger(nval,2*neq_PW)
      integer*8 index_pw_inv(neq_PW)
      complex*16 ii, val_exp, coeff_1
      double precision r_tmp, vec_kx, vec_ky
      double precision bloch1, bloch2, pi, alpha, beta, norm
      double precision lambda, freq
C
CCCCCCCCCCC Start Program CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      ui = 6
      debug = 0
C
      if ( nnodes .ne. 6 ) then
        write(ui,*) "overlap_J_dagger: problem nnodes = ", nnodes
        write(ui,*) "overlap_J_dagger: nnodes should be equal to 6!"
        write(ui,*) "overlap_J_dagger: Aborting..."
        stop
      endif
      if ( nval .gt. nval_max ) then
        write(ui,*) "overlap_J_dagger: nval = ", nval
        write(ui,*) "overlap_J_dagger: > nval_max = ",  nval_max
        write(ui,*) "overlap_J_dagger: Aborting..."
        stop
      endif
      if ( neq_PW .gt. PW_max ) then
        write(ui,*) "overlap_J_dagger: neq_PW = ", neq_PW
        write(ui,*) "overlap_J_dagger: > PW_max = ",  PW_max
        write(ui,*) "overlap_J_dagger: Aborting..."
        stop
      endif
C
      call quad_triangle(nquad,nquad_max,wq,xq,yq);
C        
C set up final solution matrix
        do i=1,neq_PW
          do j=1,nval
            overlap_K(i,j) = 0.0D0
            overlap_E(i,j) = 0.0D0
          enddo
        enddo 
C
      ii = dcmplx(0.0d0, 1.0d0)
      d = lat_vecs(1,1)  
      pi = 3.141592653589793d0
      bloch1 = bloch_vec(1)
      bloch2 = bloch_vec(2)
      vec_kx = 2.0d0*pi/d
      vec_ky = 2.0d0*pi/d
C
CCCCCCCCCCCCCCCCC loop over all elements CCCCCCCCCCCCCCCC
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
              if (px**2 + py**2 .le. ordre_ls**2) then
                alpha = bloch1 + vec_kx*px      ! Bloch vector along x
                beta  = bloch2 + vec_ky*py      ! Bloch vector along y
                norm = 1.0d0/DSQRT(alpha**2 + beta**2) ! sqrt term             
                r_tmp = alpha*xx_g(1) + beta*xx_g(2)
                s2 = index_pw_inv(s)
                val_exp = EXP(ii*r_tmp)*norm
C
C     Plane waves order s (1-neq_PW) not conjugated (-ii) -> ii 
C   RK   x component
                PlaneW_RK(s2,1) = alpha*val_exp
C             y component
                PlaneW_RK(s2,2) = beta*val_exp
C  
C   RE   x component
                PlaneW_RE(s2,1) = beta*val_exp
C             y component
                PlaneW_RE(s2,2) = -alpha*val_exp        
C
                s = s + 1
              endif
              enddo
            enddo
C
            coeff_1 = ww * ABS(det)     ! * pp(typ_e)
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
     *                 sol(trans,ltest,n,iel) * phi2_list(ltest)
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
C
        enddo
      enddo
C
C  Save as J_mat = | J_E |                
C                  | J_K |
C
      do twos = 1, 2*neq_PW
        do n = 1, nval
          if (twos .le. neq_PW) then
            overlap_J_dagger(n,twos) = overlap_E(twos,n)
          else
            overlap_J_dagger(n,twos) = overlap_K(twos-neq_PW,n)
          endif
        enddo
      enddo
C
CCCCCCCCCCCCCCCCC   save results   CCCCCCCCCCCCCCCC
C
      if (PrintAll .eq. 1) then
C
      open (unit=35, file="Matrices/J_dagger_mat.txt",
     * status='unknown')
      open (unit=33, file="Matrices/J_dagger_K.txt", status='unknown')
      open (unit=32, file="Matrices/J_dagger_E.txt", status='unknown')
      write(35,131) lambda, freq     
      write(33,131) lambda, freq   
      write(32,131) lambda, freq
C
      do n = 1, nval
        do twos = 1, 2*neq_PW
          write(35,132)  n,twos,overlap_J_dagger(n,twos),
     *         abs(overlap_J_dagger(n,twos))
        enddo
      enddo
C
      do s=1,neq_PW
        do n=1,nval
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
      return
      end
