c
c     Periodic boundary conditions
c
c
c***********************************************************************
c
      subroutine periodic_cond (i_cond, n_ddl, neq, type_N_E_F, 
     *      ip_period_E_F, ineq, debug)

c
c***********************************************************************
c
c
c   i_cond = 0 => Dirichlet boundary condition
c   i_cond = 1 => Neumann boundary condition
c   i_cond = 2 => Periodic boundary condition
c
c***********************************************************************
c
      implicit none
      integer*8 i_cond, n_ddl, neq
      integer*8 ip_period_E_F(n_ddl), type_N_E_F(2,n_ddl)
      integer*8 ineq(3,n_ddl), debug

      integer*8 i, j, k, i_boundary, i_dim
c
      if (debug .eq. 1) then
        write(*,*) "periodic_cond: i_cond = ", i_cond
      endif
        if(i_cond .eq. 2) then
          if (debug .eq. 1) then
            write(*,*) "periodic_cond: Periodic boundary conditions"
          endif
        else
        write(*,*) "periodic_cond: problem : i_cond !=2 : ", i_cond
        write(*,*) "i_cond should be 2 for Periodic boundary conditions"
        write(*,*) "periodic_cond: Aborting..."
        stop
        endif
c
      do i=1,n_ddl
        i_boundary = type_N_E_F(1,i)
        i_dim = type_N_E_F(2,i)
        j = ip_period_E_F(i)
        if (i_boundary .ne. 0 .and. j .le. 0) then
          write(*,*) "period_cond: ???"
          write(*,*) "i, i_boundary, i_period = ", 
     *         i, i_boundary, j
          write(*,*) "A boundary point should be periodic:"
          write(*,*) "i_period != 0 when i_boundary != 0"
          write(*,*) "period_cond: Aborting..."
          stop
        endif
        if (i_boundary .ne. 0 .and. j .ne. 0) then
          k = ip_period_E_F(j)
          if (k .ne. j) then
            write(*,*) "period_cond:"
            write(*,*) "period_cond: k != j :"
            write(*,*) "i, ip_period_E_F(i) = ", i, ip_period_E_F(i)
            write(*,*) "j, ip_period_E_F(j) = ", j, ip_period_E_F(j)
            write(*,*) "period_cond: Aborting..."
            stop
          endif
        endif
      enddo
c     check ...
      do i=1,n_ddl
        i_boundary = type_N_E_F(1,i)
        i_dim = type_N_E_F(2,i)
        j = ip_period_E_F(i)
        if (i_boundary .eq. 0 .and. j .ne. 0) then
        write(*,*) "period_cond: ???"
        write(*,*) "i, i_boundary, i_period = ", 
     *         i, i_boundary, j
        write(*,*) "An interior point should not be periodic:"
        write(*,*) "i_period = 0 when i_boundary = 0"
        write(*,*) "period_cond: Aborting..."
        stop
        endif
      enddo
c
      if(i_cond .eq. 2) then
c       Periodic boundary condition: all points have a degree of freedom
        neq = 0
        do i=1,n_ddl
          i_boundary = type_N_E_F(1,i)
          i_dim = type_N_E_F(2,i)
          if (i_boundary .eq. 0) then
            if (i_dim .eq. 2) then ! each element is associated to 3 interior Degrees Of Freedom (DOF)
              ineq(1,i) = neq + 1
              ineq(2,i) = neq + 2
              ineq(3,i) = neq + 3
              neq = neq + 3
            elseif (i_dim .eq. 1) then ! each edge is associated to 3 Degrees Of Freedom (DOF)
              ineq(1,i) = neq + 1
              ineq(2,i) = neq + 2
              ineq(3,i) = neq + 3
              neq = neq + 3
            elseif (i_dim .eq. 0) then ! each nodee is associated to 1 Degree Of Freedom (DOF)
              ineq(1,i) = neq + 1
              ineq(2,i) = 0
              ineq(3,i) = 0
              neq = neq + 1
            else
              write(*,*) "bound_cond: i_dim has invalid value : ", i_dim
              write(*,*) "bound_cond: i_cond = ", i_cond
              write(*,*) "bound_cond: i = ", i
              write(*,*) "bound_cond: Aborting..."
              stop
            endif
          endif
        enddo
        do i=1,n_ddl
          i_boundary = type_N_E_F(1,i)
          i_dim = type_N_E_F(2,i)
          if (i_boundary .ne. 0) then
            j = ip_period_E_F(i)
            if(j .eq. i) then
              if (i_dim .eq. 2) then ! each element is associated to 3 interior Degrees Of Freedom (DOF)
                ineq(1,i) = neq + 1
                ineq(2,i) = neq + 2
                ineq(3,i) = neq + 3
                neq = neq + 3
              elseif (i_dim .eq. 1) then ! each edge is associated to 3 Degrees Of Freedom (DOF)
                ineq(1,i) = neq + 1
                ineq(2,i) = neq + 2
                ineq(3,i) = neq + 3
                neq = neq + 3
              elseif (i_dim .eq. 0) then ! each nodee is associated to 1 Degree Of Freedom (DOF)
                ineq(1,i) = neq + 1
                ineq(2,i) = 0
                ineq(3,i) = 0
                neq = neq + 1
              else
                write(*,*) "bound_cond: i_dim has invalid value : ", 
     *            i_dim
                write(*,*) "bound_cond: i_cond = ", i_cond
                write(*,*) "bound_cond: i = ", i
                write(*,*) "bound_cond: Aborting..."
                stop
              endif
            endif
          endif
        enddo
c       set the equation for the "other" (or "destination") periodic node
        do i=1,n_ddl
          j = ip_period_E_F(i)
          i_boundary = type_N_E_F(1,i)
          i_dim = type_N_E_F(2,i)
          if(i_boundary .ne. 0 .and. j .ne. i) then
            do k=1,3
              ineq(k,i) = ineq(k,j)
              if(i_dim .eq. 1 .or. i_dim .eq. 2) then
                if(ineq(k,j) .le. 0 .or. j .le. 0) then
                  write(*,*)
                  write(*,*) "  ???"
                  write(*,*) "period_cond: ineq(j)  <= 0 or j <=0 : "
                  write(*,*) "period_cond: i, j, k, neq(k,j) = ", 
     *              i, j, k, ineq(k,j)
                  write(*,*) "period_cond: Aborting..."
                  stop
                endif
              endif
            enddo
          endif
        enddo
      else
        write(*,*)
        write(*,*) "  ???"
        write(*,*) "period_cond: i_cond has invalid value : ", i_cond
        write(*,*) "period_cond: Aborting..."
        stop
      endif

c
      return
      end
c
