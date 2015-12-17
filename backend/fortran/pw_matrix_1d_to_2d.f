C
      subroutine pw_matrix_1d_to_2d (neq_PW_1d, neq_PW_2d, nval,
     *  period_x, bloch_vec_x, bloch_vec_y, index_pw_inv_1d,
     *  debug, ordre_ls, overlap_J_1d, overlap_J_2d,
     *  overlap_J_dagger_1d, overlap_J_dagger_2d,
     *  k_0, PrintAll)
C
      implicit none
C  input output parameters
      integer*8 neq_PW_1d, neq_PW_2d, nval
      integer*8 PrintAll, ordre_ls, debug
      double precision bloch_vec_x, bloch_vec_y, k_0
      double precision period_x
      integer*8 index_pw_inv_1d(neq_PW_1d)
      complex*16 overlap_J_1d(2*neq_PW_1d,nval)
      complex*16 overlap_J_2d(2*neq_PW_2d,nval)
      complex*16 overlap_J_dagger_1d(nval,2*neq_PW_1d)
      complex*16 overlap_J_dagger_2d(nval,2*neq_PW_2d)

C  local parameters - purely internal
C
      integer alloc_stat
      integer*8, allocatable :: index_pw_2d(:), index_pw_inv_2d(:)
      complex*16, dimension(:), allocatable :: beta_z_pw_2d
      integer*8, allocatable :: pw_1d_to_2d(:)
      integer*8 ui, ipw, jval
      integer*8 px, py, s, s2, s_1d, s_2d
      double precision vec_kx, vec_ky
c  , d
      double precision bloch1, bloch2, pi, alpha, beta
      complex*16 z_tmp
C
CCCCCCCCCCC Start Program CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      ui = 6
C
      alloc_stat = 0
      allocate(index_pw_inv_2d(neq_PW_2d), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "pw_matrix_1d_to_2d: Mem. allocation ",
     *             "is unseccesfull"
        write(*,*) "alloc_stat (index_pw_inv_2d) = ", alloc_stat
        write(*,*) "Aborting..."
        stop
      endif
      allocate(index_pw_2d(neq_PW_2d), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "pw_matrix_1d_to_2d: Mem. allocation ",
     *             "is unseccesfull"
        write(*,*) "alloc_stat (index_pw_2d) = ", alloc_stat
        write(*,*) "Aborting..."
        stop
      endif
      allocate(beta_z_pw_2d(neq_PW_2d), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "pw_matrix_1d_to_2d: Mem. allocation ",
     *             "is unseccesfull"
        write(*,*) "alloc_stat (beta_z_pw_2d) = ", alloc_stat
        write(*,*) "Aborting..."
        stop
      endif

      allocate(pw_1d_to_2d(neq_PW_1d), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "pw_matrix_1d_to_2d: Mem. allocation ",
     *             "is unseccesfull"
        write(*,*) "alloc_stat (pw_1d_to_2d) = ", alloc_stat
        write(*,*) "Aborting..."
        stop
      endif

C
cc      d = lat_vecs(1,1)
      pi = 3.141592653589793d0
      bloch1 = bloch_vec_x
      bloch2 = bloch_vec_y
      vec_kx = 2.0d0*pi/period_x
      vec_ky = 2.0d0*pi/period_x
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Ordering
      s = 1
      s_1d = 1
      do px = -ordre_ls, ordre_ls
        do py = -ordre_ls, ordre_ls
          if ((px/period_x)**2 + (py/period_x)**2
     *      .le. (ordre_ls/MAX(period_x, period_x))**2) then
            alpha = bloch1 + vec_kx*px  ! Bloch vector along x
            beta  = bloch2 + vec_ky*py  ! Bloch vector along y
            z_tmp = k_0**2 - alpha**2 - beta**2
            beta_z_pw_2d(s) = sqrt(z_tmp)
            if (py == 0) then
              pw_1d_to_2d(px+ordre_ls+1) = s
            endif
            s = s + 1
          endif
        enddo
      enddo
C

      call z_indexx (neq_PW_2d, beta_z_pw_2d, index_pw_2d)
      if (debug .eq. 1) then
        write(ui,*) "index_pw_2d = ", (index_pw_2d(s),s=1,neq_PW_2d)
      endif
C
C       Inverse of index_pw_2d
      do s=1,neq_PW_2d
        s2 = index_pw_2d(s)
        index_pw_inv_2d(s2) = s
      enddo
C
      if (debug .eq. 1) then
        do s=1,neq_PW_2d
          s2 = index_pw_2d(s)
          write(ui,*) beta_z_pw_2d(s2)
        enddo
      endif
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      do jval= 1, nval
        do ipw = 1, 2*neq_PW_2d
          overlap_J_2d(ipw,jval) = 0
        enddo
      enddo
      do jval= 1, nval
        do ipw = 1, neq_PW_1d
          s_1d = index_pw_inv_1d(ipw)
          s2 = pw_1d_to_2d(ipw)
          s_2d = index_pw_inv_2d(s2)
          z_tmp = overlap_J_1d(s_1d,jval)
          overlap_J_2d(s_2d,jval) = z_tmp
          z_tmp = overlap_J_1d(s_1d+neq_PW_1d,jval)
          overlap_J_2d(s_2d+neq_PW_2d,jval) = z_tmp
        enddo
      enddo
c        do ipw = 1, 2*neq_PW_2d
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      do ipw = 1, 2*neq_PW_2d
        do jval= 1, nval
          overlap_J_dagger_2d(jval,ipw) = 0
        enddo
      enddo

      do ipw = 1, neq_PW_1d
        s_1d = index_pw_inv_1d(ipw)
        s2 = pw_1d_to_2d(ipw)
        s_2d = index_pw_inv_2d(s2)
        do jval= 1, nval
          z_tmp = overlap_J_dagger_1d(jval,s_1d)
          overlap_J_dagger_2d(jval,s_2d) = z_tmp
          z_tmp = overlap_J_dagger_1d(jval,s_1d+neq_PW_1d)
          overlap_J_dagger_2d(jval,s_2d+neq_PW_2d) = z_tmp
        enddo
      enddo
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      if (PrintAll .eq. 1) then
      open (unit=35, file="Matrices/J_mat_2d.txt", status='unknown')
      do jval = 1, nval
        do ipw = 1, 2*neq_PW_2d
          write(35,132) ipw,jval,overlap_J_2d(ipw,jval),
     *         abs(overlap_J_2d(ipw,jval))
        enddo
      enddo
      close(35)
      open (unit=34, file="Matrices/J_dagger_mat_2d.txt",
     * status='unknown')
      do jval = 1, nval
        do ipw = 1, 2*neq_PW_2d
          write(34,132)  jval,ipw,overlap_J_dagger_2d(jval,ipw),
     *         abs(overlap_J_dagger_2d(jval,ipw))
        enddo
      enddo
      close(34)
      open (unit=35, file="Matrices/pw_1d_to_2d.txt",
     * status='unknown')
      s = 1
      do px = -ordre_ls, ordre_ls
          write(35,"(4(I7))")  px, s, pw_1d_to_2d(px+ordre_ls+1)
            s = s + 1
      enddo
      close(35)
      endif
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
132   format(2(I4),2(g25.17),g18.10)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
      deallocate(index_pw_2d)
      deallocate(beta_z_pw_2d)
C
      return
      end
