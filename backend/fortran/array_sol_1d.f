
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   sol_0(*,i) : contains the imaginary and real parts of the solution for points such that ineq(i) != 0
c   sol(i) : contains solution for all points
c     sol(1..3,nval, nel)  contains the values of Ex component at P2 interpolation nodes
c     sol(4..7,nval, nel)  contains the values of Ey component at P3 interpolation nodes
c     sol(8..11,nval, nel) contains the values of Ez component at P3 interpolation nodes
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine array_sol_1d (nval, nel, n_ddl, neq,
     *     n_core, bloch_vec_x, index,
     *     table_ddl, type_el, ineq,
     *     ip_period_ddl, x_ddl,
     *     v_cmplx, mode_pol, sol_0, sol, n_k)
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit none
      integer*8 nval, nel, n_ddl, neq, n_k
      integer*8 n_core(2), type_el(nel)
      integer*8 ineq(n_ddl), index(neq)
      integer*8 ip_period_ddl(n_ddl)
      integer*8 table_ddl(3+4+4,nel)
      double precision x_ddl(n_ddl)
      double precision bloch_vec_x
      complex*16 sol_0(neq,nval)
      complex*16 sol(3+4+4,nval,nel)

      complex*16 v_cmplx(nval), mode_pol(4,nval)

c     Local variables
      complex*16, allocatable :: v_tmp(:)
      integer allocate_status
      integer*8 ip, jp, ind_ip
      double precision delta_x

      integer*8 nnodes_0, nddl_0
      parameter (nnodes_0 = 3, nddl_0 = 11)
      complex*16 sol_el(nddl_0)

      integer*8 nod_el_p(nddl_0)

      double precision mode_comp(4)
      double complex val_exp(nddl_0)
      double precision x_min_0, x_max_0

      integer*8 j, k, j1, inod, typ_e, debug, i_sol_max
      integer*8 iel, ival, ival2
      complex*16 ii, z_tmp1, z_tmp2, z_sol_max

      double precision x_min, x_max, x_mid, dx, xx
      complex*16 z_sol_max_L, z_sol_max_R  ! Maximum to left and right of the middle point
      complex*16 z_sol_max_mid  ! Value at (or near) the middle point
c
c  ii = sqrt(-1)
      ii = cmplx(0.0d0, 1.0d0)
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
      debug = 0

c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c

      allocate_status = 0
      allocate(v_tmp(nval), STAT=allocate_status)
      if (allocate_status /= 0) then
        write(*,*) "The allocation is unsuccessful"
        write(*,*) "allocate_status = ", allocate_status
        write(*,*) "Not enough memory for sol_1"
        write(*,*) "nval, nel = ", nval,nel
        write(*,*) "Aborting..."
        stop
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
      x_min = x_ddl(1)
      x_max = x_ddl(1)
      do j=1,n_ddl
        xx = x_ddl(j)
        if(xx .lt. x_min) x_min = xx
        if(xx .gt. x_max) x_max = xx
      enddo
      x_mid = (x_min + x_max) / 2.0d0
      dx = (x_max - x_min)/dble(nel)
      if (debug .eq. 1) then
        write(*,*) "array_sol_1d: x_min = ", x_min
        write(*,*) "array_sol_1d: x_max = ", x_max
        write(*,*) "array_sol_1d: x_mid = ", x_mid
        write(*,*) "array_sol_1d:    dx = ", dx
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c

      do j=1,nval
        j1=index(j)
        v_tmp(j) = v_cmplx(j1)
      enddo
      do j=1,nval
        v_cmplx(j) = v_tmp(j)
      enddo

      do ival=1,nval
        ival2 = index(ival)
        do j=1,4
          mode_pol(j,ival) = 0.0d0
        enddo
        z_sol_max = 0.0d0
        z_sol_max_L = 0.0d0
        z_sol_max_R = 0.0d0
        z_sol_max_mid = 0.0d0
        i_sol_max = 0
        do iel=1,nel
          typ_e = type_el(iel)
          do j=1,4
            mode_comp(j) = 0.0d0
          enddo
          do inod=1,nddl_0
            j = table_ddl(inod,iel)
            nod_el_p(inod) = j
          enddo
c         Periodic boundary condition
          do inod=1,nddl_0
            j = table_ddl(inod,iel)
            k = ip_period_ddl(j)
            if (k .ne. 0) j=k
            nod_el_p(inod) = j
          enddo
          do inod=1,nddl_0
            val_exp(inod) = 1.0d0
          enddo
c         val_exp: Bloch mode phase factor between the origin point and destination point
c         For a pair of periodic points, one is chosen as origin and the other is the destination
          do j=1,nddl_0
            jp = table_ddl(j,iel)
            j1 = ip_period_ddl(jp)
            if (j1 .ne. 0 .and. j1 .ne. jp) then
              delta_x = x_ddl(jp) - x_ddl(j1)
              val_exp(j) = exp(ii * bloch_vec_x * delta_x)
            endif
          enddo
          do inod=1,nddl_0
            ip = table_ddl(inod,iel)
            ind_ip = ineq(ip)
            z_tmp1 = sol_0(ind_ip, ival2)
            z_tmp1 = z_tmp1 * val_exp(inod)
            sol_el(inod) = z_tmp1
          enddo
          do inod=1,nddl_0
              z_tmp2 = sol_el(inod)
              sol(inod,ival,iel) = z_tmp2
              if (abs(z_sol_max) .lt. abs(z_tmp2)) then
                z_sol_max = z_tmp2
                i_sol_max = table_ddl(inod,iel)
              endif
          enddo
cccccccccc NEW ADDITION
          do inod=1,nddl_0
            jp = table_ddl(inod,iel)
            xx = x_ddl(jp)
            z_tmp2 = sol_el(inod)
            if (abs(z_sol_max_L) < abs(z_tmp2) .and.
     *        xx < x_mid) then
              z_sol_max_L = z_tmp2
            endif
            if (abs(z_sol_max_R) < abs(z_tmp2) .and.
     *        xx > x_mid) then
              z_sol_max_R = z_tmp2
            endif
            if (abs(z_sol_max_mid) < abs(z_tmp2) .and.
     *        abs(xx - x_mid) < dx/4.0d0) then
              z_sol_max_mid = z_tmp2
            endif
          enddo
cccccccccc END OF NEW ADDITION
c           Contribution of the element iel to the mode component
          do inod=1,3
            z_tmp2 = abs(sol_el(inod))**2
            mode_comp(1) = mode_comp(1) + z_tmp2   !   X-component
          enddo
          do inod=4,7
            z_tmp2 = abs(sol_el(inod))**2
            mode_comp(2) = mode_comp(2) + z_tmp2   !   Y-component
          enddo
          do inod=8,11
            z_tmp2 = abs(sol_el(inod))**2
            mode_comp(3) = mode_comp(3) + z_tmp2   !   Z-component
          enddo
          j = table_ddl(1,iel)
c          write(*,*) "iel, j = ", iel, j
          x_min_0 = x_ddl(j)
          j = table_ddl(3,iel)
          x_max_0 = x_ddl(j)
c         Avarage values
          do j=1,3
            mode_comp(j) = abs(x_max_0 - x_min_0) * mode_comp(j)
     *                    / dble(nddl_0)
          enddo
c         Add the contribution of the element iel to the mode component
          do j=1,3
            mode_pol(j,ival) = mode_pol(j,ival) + mode_comp(j)
          enddo
          if (typ_e .eq. n_core(1) .or. typ_e .eq. n_core(2)) then
            z_tmp2 = mode_comp(1) + mode_comp(2)
     *        + mode_comp(3)
            mode_pol(4,ival) = mode_pol(4,ival) + z_tmp2
          endif
        enddo
cccccccccccccccccccccccccccccc
c       Total energy and normalization
        z_tmp2 = mode_pol(1,ival) + mode_pol(2,ival)
     *        + mode_pol(3,ival)
        if (abs(z_tmp2) .lt. 1.0d-10) then
          write(*,*) "array_sol_1d: the total energy ",
     *       "is too small : ", z_tmp2
          write(*,*) "array_sol_1d: ival ival2 = ", ival, ival2
          write(*,*) "array_sol_1d: zero eigenvector; aborting..."
          stop
        endif
        do j=1,3
          mode_pol(j,ival) = mode_pol(j,ival) / z_tmp2
        enddo
        j=4
          mode_pol(j,ival) = mode_pol(j,ival) / z_tmp2
c       Check if the eigenvector is nonzero
        if (abs(z_sol_max) .lt. 1.0d-10) then
          z_sol_max = z_tmp2
          write(*,*) "array_sol_1d: z_sol_max is too small"
          write(*,*) "array_sol_1d: z_sol_max = ", z_sol_max
          write(*,*) "ival, ival2, nval = ", ival, ival2, nval
          write(*,*) "array_sol_1d: zero eigenvector; aborting..."
          stop
        endif
        if (debug .eq. 1) then
          write(*,*) "array_sol_1d: "
          write(*,*) "                           ival = ", ival
          write(*,*) "z_sol_max     = ", z_sol_max, abs(z_sol_max)
          write(*,*) "z_sol_max_L   = ", z_sol_max_L, abs(z_sol_max_L)
          write(*,*) "z_sol_max_R   = ", z_sol_max_R, abs(z_sol_max_R)
          write(*,*) "z_sol_max_mid = ", z_sol_max_mid,
     *        abs(z_sol_max_mid)
        endif
c       Normalization so that the maximum field component is 1
        if (abs(z_sol_max_mid/z_sol_max) >= 0.95d0) then
          if (debug .eq. 1) then
            write(*,*) "array_sol_1d: CASE max_mid: n_k = ", n_k
          endif
          do iel=1,nel
            do inod=1,nddl_0
              z_tmp1 = sol(inod,ival,iel) / z_sol_max_mid
              sol(inod,ival,iel) = z_tmp1
            enddo
          enddo
        elseif (n_k == 1) then
          if (debug .eq. 1) then
            write(*,*) "array_sol_1d: CASE max_R: n_k = ", n_k
          endif
          do iel=1,nel
            do inod=1,nddl_0
              z_tmp1 = sol(inod,ival,iel) / z_sol_max_R
              sol(inod,ival,iel) = z_tmp1
            enddo
          enddo
        elseif (n_k == 2) then
          if (debug .eq. 1) then
            write(*,*) "array_sol_1d: CASE max_L: n_k ", n_k
          endif
          do iel=1,nel
            do inod=1,nddl_0
              z_tmp1 = sol(inod,ival,iel) / z_sol_max_L
              sol(inod,ival,iel) = z_tmp1
            enddo
          enddo
        else
          if (debug .eq. 1) then
            write(*,*) "array_sol_1d: CASE 4: n_k", n_k
          endif
          write(*,*)
          write(*,*) "array_sol_1d: no action defined for n_k = ", n_k
          write(*,*) "array_sol_1d: aborting..."
          stop
        endif
      enddo
c
      if (debug .eq. 1) then
        open (unit=24,file="Output/mode_pol.txt",status="unknown")
        do ival=1,nval
            write(24,*) ival, (dble(mode_pol(j,ival)),j=1,4)
        enddo
        close(24)
      endif
c
      return
      end
