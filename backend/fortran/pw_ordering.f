C
      subroutine pw_ordering (neq_PW, lat_vecs,
     *  bloch_vec, index_pw_inv,
     *  debug, ordre_ls, k_0)
C
      implicit none
C  input output parameters
      integer*8 neq_PW, ordre_ls, debug
      double precision bloch_vec(2), k_0
      double precision lat_vecs(2,2)
      integer*8 index_pw_inv(neq_PW)

C  local parameters - purely internal
C
      integer alloc_stat
      integer*8, dimension(:), allocatable :: index_pw
      complex*16, dimension(:), allocatable :: beta_z_pw
      integer*8 ui
      integer*8 px, py, s, s2
      double precision vec_kx, vec_ky, d, dy
      double precision bloch1, bloch2, pi, alpha, beta
      complex*16 z_tmp
C
CCCCCCCCCCC Start Program CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      ui = 6
C
      alloc_stat = 0
      allocate(index_pw(neq_PW), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "pw_ordering: Mem. allocation is unseccesfull"
        write(*,*) "alloc_stat (index_pw) = ", alloc_stat
        write(*,*) "Aborting..."
        stop
      endif
      allocate(beta_z_pw(neq_PW), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "pw_ordering: Mem. allocation is unseccesfull"
        write(*,*) "alloc_stat (beta_z_pw) = ", alloc_stat
        write(*,*) "Aborting..."
        stop
      endif
C
      d = lat_vecs(1,1)
      dy = lat_vecs(2,2)
      pi = 3.141592653589793d0
      bloch1 = bloch_vec(1)
      bloch2 = bloch_vec(2)
      vec_kx = 2.0d0*pi/d
      vec_ky = 2.0d0*pi/dy
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Ordering
      s = 1
      do px = -ordre_ls, ordre_ls
        do py = -ordre_ls, ordre_ls
            if ((px/d)**2 + (py/dy)**2
     *         .le. (ordre_ls/MAX(d, dy))**2) then
            alpha = bloch1 + vec_kx*px  ! Bloch vector along x
            beta  = bloch2 + vec_ky*py  ! Bloch vector along y
            z_tmp = k_0**2 - alpha**2 - beta**2
            beta_z_pw(s) = sqrt(z_tmp)
            s = s + 1
          endif
        enddo
      enddo
C
      call z_indexx (neq_PW, beta_z_pw, index_pw)
      if (debug .eq. 1) then
        write(ui,*) "pw_ordering.f:"
        write(ui,*) "index_pw = ", (index_pw(s),s=1,neq_PW)
      endif
C
C       Inverse of index_pw
      do s=1,neq_PW
        s2 = index_pw(s)
        index_pw_inv(s2) = s
      enddo
C
      if (debug .eq. 1) then
        write(ui,*) "pw_ordering.f: PW k_z = "
        do s=1,neq_PW
          s2 = index_pw(s)
          write(ui,*) beta_z_pw(s2)
        enddo
      endif
C
      deallocate(index_pw)
      deallocate(beta_z_pw)
C
      return
      end
