C************************************************************************
C
C  Calculate the Projection integral of the conjugate Plane Waves & Bloch Modes
c
c     sol(1..3,nval, nel)  contains the values of Ex component at P2 interpolation nodes
c     sol(4..7,nval, nel)  contains the values of Ey component at P3 interpolation nodes
c     sol(8..11,nval, nel) contains the values of Ez component at P3 interpolation nodes
C
      subroutine J_overlap_1d(nval, nel, npt_P2, 
     *  type_el, table_nod, x_P2, sol, period_x, 
     *  lambda, freq, overlap_J, neq_PW, bloch_vec_x,
     *  bloch_vec_y, index_pw_inv, PrintAll, ordre_ls)

C************************************************************************
      implicit none 
C  input output parameters
      integer*8 nval, nel, npt_P2
      integer*8 type_el(nel)
      integer*8 table_nod(3,nel)
      complex*16 sol(3+4+4,nval,nel)
      complex*16 overlap_J(2*neq_PW,nval)
      double precision x_P2(npt_P2), period_x
      integer*8 neq_PW, PrintAll, ordre_ls
      double precision bloch_vec_x, bloch_vec_y
      integer*8 index_pw_inv(neq_PW)
C  local parameters - purely internal
      integer*8 iel, inode, n, i, j, global
      integer*8 ltest, typ_e, s, s2, twos, trans
      integer*8 ui, px, nnodes_0
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
      complex*16, allocatable :: PlaneW_Exp_P2(:,:)
      complex*16, allocatable :: PlaneW_Exp_P3(:,:)
C
      double precision, allocatable :: tmp_ls_alpha(:)
      double precision xmin, xmax
      double precision r_tmp, vec_kx, lambda, freq
      double precision bloch1, pi, alpha, norm
C
CCCCCCCCCCC Start Program CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      ui = 6
C
      allocate(PlaneW_RK(neq_PW,2), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "J_overlap_1d: Mem. allocation is unseccesfull"
        write(*,*) "alloc_stat (PlaneW_RK) = ", alloc_stat
        write(*,*) "Aborting..."
        stop
      endif
      allocate(PlaneW_RE(neq_PW,2), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "J_overlap_1d: Mem. allocation is unseccesfull"
        write(*,*) "alloc_stat (PlaneW_RE) = ", alloc_stat
        write(*,*) "Aborting..."
        stop
      endif
      allocate(overlap_K(neq_PW,nval), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "J_overlap_1d: Mem. allocation is unseccesfull"
        write(*,*) "alloc_stat (overlap_K) = ", alloc_stat
        write(*,*) "Aborting..."
        stop
      endif
      allocate(overlap_E(neq_PW,nval), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "J_overlap_1d: Mem. allocation is unseccesfull"
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
      do j=1,nval
        do i=1,neq_PW
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
CCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      do n = 1, nval
        do twos = 1, 2*neq_PW
          if (twos .le. neq_PW) then
            overlap_J(twos,n) = -911
          else
            overlap_J(twos,n) = -911
          endif
        enddo     
      enddo
C
CCCCCCCCCCCCCCCCC	loop over all elements	CCCCCCCCCCCCCCCC
C
      do iel=1,nel
        typ_e = type_el(iel)
        do inode=1,nnodes_0
          global = table_nod(inode,iel)
          nod_el_p(inode) = global
        enddo
        xmin = x_P2(nod_el_p(1))
        xmax = x_P2(nod_el_p(2))
C     prefactors of plane waves order s(px)
        s = 1
        do px = -ordre_ls, ordre_ls
          alpha = bloch1 + vec_kx*px      ! Bloch vector along x
          s2 = index_pw_inv(s)
          tmp_ls_alpha(s2) = alpha
C         Plane waves order s
          if (alpha .ne. 0.0d0 .or. bloch_vec_y .ne. 0.0d0) then
            norm = 1.0d0/SQRT(alpha**2 + bloch_vec_y**2)  ! sqrt term
C           RK    x component
              PlaneW_RK(s2,1) = alpha *norm  ! amplitude only - do not include the exponential function
C           y component
              PlaneW_RK(s2,2) = bloch_vec_y * norm
C           RE    x component
              PlaneW_RE(s2,1) = bloch_vec_y * norm
C           y component
              PlaneW_RE(s2,2) = -alpha * norm
          else
C           RK    x component
              PlaneW_RK(s2,1) = 0
C           y component
              PlaneW_RK(s2,2) = 1
C           RE    x component
              PlaneW_RE(s2,1) = 1
C           y component
              PlaneW_RE(s2,2) = 0
          endif
C         Compuation of the overlap betwween Exp[-I alpha x] and the scalar polynomial basis functions
c         vecP2Exp(i) = Integrate[lsP2[[i]] Exp[-I alpha x], {x, xmin, xmax}]
c         vecP3Exp(i) = Integrate[lsP3[[i]] Exp[-I alpha x], {x, xmin, xmax}]
          call vector_p2_exp_1d(xmin, xmax, -alpha, vecP2Exp)
          call vector_p3_exp_1d(xmin, xmax, -alpha, vecP3Exp)
          do i=1,nnode_P2
            PlaneW_Exp_P2(s2,i) = vecP2Exp(i)
          enddo
          do i=1,nnode_P3
            PlaneW_Exp_P3(s2,i) = vecP3Exp(i)
          enddo
          s = s + 1
        enddo
C
            do s=1,neq_PW 
              do trans=1,2              ! transverse field components
                K_tmp(trans) = PlaneW_RK(s,trans)
                E_tmp(trans) = PlaneW_RE(s,trans)
              enddo
              do n=1,nval
                do trans=1,2
                  vec_phi(trans) = 0.0d0
                enddo
                trans = 1 ! x-component of the field 
                do ltest=1,3    ! x-component -- P2 FEM
                  vec_phi(trans) = vec_phi(trans) 
     *               + sol(ltest,n,iel) * PlaneW_Exp_P2(s,ltest)   ! Integral of the Bloch mode * Exp[-I alpha x]
                enddo
                trans = 2 ! y-component of the field 
                do ltest=1,4    ! y-component -- P3 FEM
                    vec_phi(trans) = vec_phi(trans) 
     *               + sol(ltest+3,n,iel) * PlaneW_Exp_P3(s,ltest)   ! Integral of the Bloch mode * Exp[-I alpha x]
                enddo
C take dot product of Bloch Lp and PW for each quad pt.
                   tmp1 = vec_phi(1)*K_tmp(1) + vec_phi(2)*K_tmp(2)
                   overlap_K(s,n) = overlap_K(s,n) + tmp1
C
                   tmp2 = vec_phi(1)*E_tmp(1) + vec_phi(2)*E_tmp(2)
                   overlap_E(s,n) = overlap_E(s,n) + tmp2
              enddo
            enddo
cc        enddo
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
        open (unit=35, file="Matrices/ls_alpha.txt", status='unknown')
          do s=1,neq_PW
            write(35,*) s, tmp_ls_alpha(s)
          enddo
        close(35)
C
        open (unit=35, file="Matrices/J_mat.txt", status='unknown')
        open (unit=33, file="Matrices/J_K.txt", status='unknown')
        open (unit=32, file="Matrices/J_E.txt", status='unknown')
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
        close(33)
        close(32) 
C
131   format(2(f12.4))
132   format(2(I4),2(g25.17),g18.10)
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
