C
      subroutine pw_ordering_1d (neq_PW, period_x,
     *  bloch_vec_x, bloch_vec_y, index_pw_inv,
     *  debug, ordre_ls, k_0)
c lat_vecs, bloch_vec, 
C 
      implicit none 
C  input output parameters
      integer*8 Zeroth_Order, Zeroth_Order_inv
      integer*8 neq_PW, PrintAll, ordre_ls, debug
      double precision k_0, bloch_vec_x, bloch_vec_y
c  bloch_vec(2), 
      double precision period_x
c                lat_vecs(2,2)
      integer*8 index_pw_inv(neq_PW)

C  local parameters - purely internal
C
      integer alloc_stat
      integer*8, dimension(:), allocatable :: index_pw
      complex*16, dimension(:), allocatable :: beta_z_pw
      integer*8 ui
      integer*8 px, py, s, s2
      double precision vec_kx, vec_ky
c , d
      double precision bloch1, pi, alpha
c, beta, bloch2
      complex*16 z_tmp
C
CCCCCCCCCCC Start Program CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      ui = 6
C
      alloc_stat = 0
      allocate(index_pw(neq_PW), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "pw_ordering_1d: Mem. allocation is unseccesfull"
        write(*,*) "alloc_stat (index_pw) = ", alloc_stat
        write(*,*) "Aborting..."
        stop
      endif
      allocate(beta_z_pw(neq_PW), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "pw_ordering_1d: Mem. allocation is unseccesfull"
        write(*,*) "alloc_stat (beta_z_pw) = ", alloc_stat
        write(*,*) "Aborting..."
        stop
      endif
C
c      d = period_x
      pi = 3.141592653589793d0
      bloch1 = bloch_vec_x
cc      bloch2 = bloch_vec_y
      vec_kx = 2.0d0*pi/period_x
c      vec_ky = 2.0d0*pi/d
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Ordering
      s = 1
      do px = -ordre_ls, ordre_ls
        alpha = bloch1 + vec_kx*px 
        z_tmp = k_0**2 - alpha**2 - bloch_vec_y**2
        beta_z_pw(s) = sqrt(z_tmp)
        s = s + 1
      enddo


cc        do py = -ordre_ls, ordre_ls
cc          if (px**2 + py**2 .le. ordre_ls**2) then
cc            beta  = bloch2 + vec_ky*py  ! Bloch vector along y
cc          endif
cc        enddo

C      
      call z_indexx (neq_PW, beta_z_pw, index_pw)
      if (debug .eq. 1) then
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
