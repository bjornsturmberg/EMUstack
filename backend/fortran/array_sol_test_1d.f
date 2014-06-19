c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     sol_1(1..3,nval, nel)  contains the values of Ex component at P2 interpolation nodes
c     sol_1(4..7,nval, nel)  contains the values of Ey component at P3 interpolation nodes
c     sol_1(8..11,nval, nel) contains the values of Ez component at P3 interpolation nodes

c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine array_sol_test_1d (nval, nel, n_ddl, 
     *           table_ddl, x_ddl, sol_1, period_x,
     *           bloch_vec_x, bloch_vec_y)
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit none
      integer*8 nval, nel, n_ddl
      complex*16 sol_1(3+4+4,nval,nel)
      integer*8 table_ddl(3+4+4,nel)
      double precision x_ddl(n_ddl), period_x
      double precision bloch_vec_x, bloch_vec_y

c     Local variables
      integer*8 nddl_0
      parameter (nddl_0 = 11)
      integer*8 inod, ip, iel, ival, debug
      integer*8 neq_PW, ordre_min, ordre_max
      double precision xx, vec_kx
      integer px, s, s2
      double precision bloch1, pi, alpha
      complex*16 ii, z_tmp1, z_tmp
      integer alloc_stat
      double precision, allocatable :: ls_alpha(:)
      integer*8, allocatable ::  index_pw_inv(:)  ! (neq_PW)
      integer*8, dimension(:), allocatable :: index_pw
      complex*16, dimension(:), allocatable :: beta_z_pw

c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
C
      ii = dcmplx(0.0d0, 1.0d0)
      pi = 3.141592653589793d0
      bloch1 = bloch_vec_x
      vec_kx = 2.0d0*pi/period_x

      debug = 1
      ordre_min = -(nval-1)/2
      ordre_max = nval + ordre_min - 1
      neq_PW = ordre_max - ordre_min + 1  ! = nval
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
      alloc_stat = 0
      allocate(ls_alpha(nval), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "J_overlap_1d: Mem. allocation is unseccesfull"
        write(*,*) "alloc_stat (ls_alpha) = ", alloc_stat
        write(*,*) "Aborting..."
        stop
      endif
      allocate(index_pw_inv(neq_PW), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "The allocation is unsuccessful"
        write(*,*) "alloc_stat = ", alloc_stat
        write(*,*) "Not enough memory for index_pw_inv"
        write(*,*) "neq_PW = ", neq_PW
        write(*,*) "Aborting..."
        stop
      endif
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
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
C     PW Ordering
      s = 1
      do px = ordre_min, ordre_max
        alpha = bloch1 + vec_kx*px      ! Bloch vector along x
        z_tmp =  - alpha**2 - bloch_vec_y**2
        beta_z_pw(s) = sqrt(z_tmp)
        s = s + 1
      enddo
      call z_indexx (neq_PW, beta_z_pw, index_pw)
C       Inverse of index_pw
      do s=1,neq_PW
        s2 = index_pw(s)
        index_pw_inv(s2) = s
      enddo
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
C     prefactors of plane waves order s(px)
      s = 1
      do px = ordre_min, ordre_max
        alpha = bloch1 + vec_kx*px      ! Bloch vector along x
        s2 = index_pw_inv(s)
        ls_alpha(s2) = alpha
        s = s + 1
      enddo
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do iel=1,nel
        do ival=1,nval
          alpha = ls_alpha(ival)
          do inod=1,3    ! x-component of the field 
            ip = table_ddl(inod,iel)
            xx = x_ddl(ip)
            z_tmp = exp(ii * alpha * xx)
            sol_1(inod,ival,iel) =  z_tmp
          enddo
          do inod=4,7    ! x-component of the field 
            ip = table_ddl(inod,iel)
            xx = x_ddl(ip)
            z_tmp = exp(ii * alpha * xx)
            sol_1(inod,ival,iel) = 2 * z_tmp
          enddo
          do inod=8,11    ! x-component of the field 
            ip = table_ddl(inod,iel)
            xx = x_ddl(ip)
            z_tmp = exp(ii * alpha * xx)
            sol_1(inod,ival,iel) = 3 * z_tmp
          enddo
        enddo
      enddo

      if (.false.) then
      do iel=1,nel
        do ival=1,nval
          alpha = ls_alpha(ival)
          do inod=1,nddl_0
            ip = table_ddl(inod,iel)
            xx = x_ddl(ip)
            z_tmp = exp(ii * alpha * xx)
            sol_1(inod,ival,iel) =  z_tmp
          enddo
        enddo
      enddo
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
      return
      end
