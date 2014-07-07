C  Calculate the Projection integral of the Plane Waves & adjoint Bloch modes
C
      subroutine J_dagger_overlap_1d(nval, nel, npt_P2, 
     *  type_el, table_nod, x_P2, sol, period_x, 
     *  lambda, freq, overlap_J_dagger, neq_PW, bloch_vec_x,
     *  bloch_vec_y, index_pw_inv, PrintAll, ordre_ls)

C************************************************************************
      implicit none 
C  input output parameters
      integer*8 nval, nel, npt_P2, ordre_ls
      integer*8 type_el(nel)
      integer*8 table_nod(3,nel)
      double precision x_P2(npt_P2)
      complex*16 sol(3+4+4,nval,nel)
      double precision period_x
      integer*8 neq_PW, PrintAll
      double precision bloch_vec_x, bloch_vec_y
C  local parameters - purely internal
      integer*8 iel, inode, n, i, j, global
      integer*8 ltest, typ_e, s, s2, twos, trans
      integer*8 ui, px
      integer*8 nnodes_0, debug
      parameter (nnodes_0 = 3)
      integer*8 nod_el_p(nnodes_0)
      complex*16 vec_phi(2), K_tmp(2), E_tmp(2)
      complex*16 tmp1, tmp2, ii
C
      integer alloc_stat
      complex*16, dimension(:,:), allocatable :: PlaneW_RK
      complex*16, dimension(:,:), allocatable :: PlaneW_RE
      complex*16, dimension(:,:), allocatable :: overlap_K
      complex*16, dimension(:,:), allocatable :: overlap_E
C
      integer*8 nnode_P2, nnode_P3
      parameter (nnode_P2 = 3, nnode_P3 = 4)
      complex*16 vecP2Exp(nnode_P2), vecP3Exp(nnode_P3)
C
      complex*16, allocatable :: PlaneW_Exp_P2(:,:)
      complex*16, allocatable :: PlaneW_Exp_P3(:,:)
      double precision, allocatable :: tmp_ls_alpha(:)
      double precision xmin, xmax
C
      complex*16 overlap_J_dagger(nval,2*neq_PW)
      integer*8 index_pw_inv(neq_PW)
      double precision r_tmp, vec_kx
      double precision bloch1, pi, alpha, norm
      double precision lambda, freq
C
CCCCCCCCCCC Start Program CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      ui = 6
      debug = 0
C
      allocate(PlaneW_RK(neq_PW,2), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "J_overlap_dagger: Mem. allocation is unseccesfull"
        write(*,*) "alloc_stat (PlaneW_RK) = ", alloc_stat
        write(*,*) "Aborting..."
        stop
      endif
      allocate(PlaneW_RE(neq_PW,2), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "J_overlap_dagger: Mem. allocation is unseccesfull"
        write(*,*) "alloc_stat (PlaneW_RE) = ", alloc_stat
        write(*,*) "Aborting..."
        stop
      endif
      allocate(overlap_K(neq_PW,nval), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "J_overlap_dagger: Mem. allocation is unseccesfull"
        write(*,*) "alloc_stat (overlap_K) = ", alloc_stat
        write(*,*) "Aborting..."
        stop
      endif
      allocate(overlap_E(neq_PW,nval), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "J_overlap_dagger: Mem. allocation is unseccesfull"
        write(*,*) "alloc_stat (overlap_E) = ", alloc_stat
        write(*,*) "Aborting..."
        stop
      endif
      allocate(PlaneW_Exp_P2(neq_PW,nnode_P2), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "J_overlap_1d: Mem. allocation is unseccesfull"
        write(*,*) "alloc_stat (PlaneW_Exp_P2) = ", alloc_stat
        write(*,*) "Aborting..."
        stop
      endif
      allocate(PlaneW_Exp_P3(neq_PW,nnode_P3), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "J_overlap_1d: Mem. allocation is unseccesfull"
        write(*,*) "alloc_stat (PlaneW_Exp_P3) = ", alloc_stat
        write(*,*) "Aborting..."
        stop
      endif
      allocate(tmp_ls_alpha(neq_PW), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "J_overlap_1d: Mem. allocation is unseccesfull"
        write(*,*) "alloc_stat (tmp_ls_alpha) = ", alloc_stat
        write(*,*) "Aborting..."
        stop
      endif
C
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C set up solution matrix
        do i=1,neq_PW
          do j=1,nval
            overlap_K(i,j) = 0.0D0
            overlap_E(i,j) = 0.0D0
          enddo
        enddo 
C
      ii = dcmplx(0.0d0, 1.0d0)
      pi = 3.141592653589793d0
      bloch1 = bloch_vec_x
      vec_kx = 2.0d0*pi/period_x
C
CCCCCCCCCCCCCCCCC loop over all elements CCCCCCCCCCCCCCCC
C
      do iel=1,nel
        typ_e = type_el(iel)
        do inode=1,nnodes_0
          global = table_nod(inode,iel)
          nod_el_p(inode) = global
        enddo
        xmin = x_P2(nod_el_p(1))
        xmax = x_P2(nod_el_p(2))
C     prefactors of plane waves order s(px,py) 
        s = 1
        do px = -ordre_ls, ordre_ls
ccc              do py = -ordre_ls, ordre_ls
ccc              if (px**2 + py**2 .le. ordre_ls**2) then
          alpha = bloch1 + vec_kx*px      ! Bloch vector along x
          s2 = index_pw_inv(s)
          tmp_ls_alpha(s2) = alpha
ccc                beta  = bloch2 + vec_ky*py      ! Bloch vector along y
ccc                r_tmp = alpha*xx_g(1) + beta*xx_g(2)
ccc                val_exp = EXP(ii*r_tmp)*norm
C
C     Plane waves order s (1-neq_PW) not conjugated (-ii) -> ii
          if (alpha .ne. 0.0d0 .or. bloch_vec_y .ne. 0.0d0) then
            norm = 1.0d0/DSQRT(alpha**2 + bloch_vec_y**2) ! sqrt term             
C           RK   x component
              PlaneW_RK(s2,1) = alpha * norm
C           y component
              PlaneW_RK(s2,2) = bloch_vec_y * norm
C           RE   x component
              PlaneW_RE(s2,1) = bloch_vec_y * norm
C           y component
              PlaneW_RE(s2,2) = -alpha * norm
          else
C           RK   x component
              PlaneW_RK(s2,1) = 0
C           y component
              PlaneW_RK(s2,2) = 1
C           RE   x component
              PlaneW_RE(s2,1) = 1
C           y component
              PlaneW_RE(s2,2) = 0
          endif
C         Compuation of the overlap betwween Exp[I alpha x] and the scalar polynomial basis functions
c         vecP2Exp(i) = Integrate[lsP2[[i]] Exp[I alpha x], {x, xmin, xmax}]
c         vecP3Exp(i) = Integrate[lsP3[[i]] Exp[I alpha x], {x, xmin, xmax}]
          call vector_p2_exp_1d(xmin, xmax, alpha, vecP2Exp)
          call vector_p3_exp_1d(xmin, xmax, alpha, vecP3Exp)
          do i=1,nnode_P2
            PlaneW_Exp_P2(s2,i) = vecP2Exp(i)
          enddo
          do i=1,nnode_P3
            PlaneW_Exp_P3(s2,i) = vecP3Exp(i)
          enddo
          s = s + 1
ccc              endif
ccc              enddo
            enddo
cc            coeff_1 = ww * ABS(det)     ! * pp(typ_e)
          do s=1,neq_PW
            do trans=1,2              ! transverse field components
              K_tmp(trans) = PlaneW_RK(s,trans)
              E_tmp(trans) = PlaneW_RE(s,trans)
            enddo
            do n=1,nval
              do trans=1,2
                vec_phi(trans) = 0.0d0
              enddo
ccc              do ltest=1,nnodes
ccc                do trans=1,2          ! transverse field components
              trans = 1 ! x-component of the field
              do ltest=1,3    ! x-component -- P2 FEM
                vec_phi(trans) = vec_phi(trans) +
     *             sol(ltest,n,iel) * PlaneW_Exp_P2(s,ltest)   ! Integral of the Bloch mode Ex * Exp[I alpha x]
              enddo
                trans = 2 ! y-component of the field 
                do ltest=1,4    ! y-component -- P3 FEM
                vec_phi(trans) = vec_phi(trans) +
     *             sol(ltest+3,n,iel) * PlaneW_Exp_P3(s,ltest)   ! Integral of the Bloch mode Ey * Exp[I alpha x]
              enddo
ccc              enddo
C             Take dot product of Bloch Lp and PW for each quad pt.
              tmp1 = vec_phi(1)*K_tmp(1) + vec_phi(2)*K_tmp(2)
              overlap_K(s,n) = overlap_K(s,n) + tmp1
              tmp2 = vec_phi(1)*E_tmp(1) + vec_phi(2)*E_tmp(2)
              overlap_E(s,n) = overlap_E(s,n) + tmp2
            enddo
          enddo
C
ccc        enddo   !  iq
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
