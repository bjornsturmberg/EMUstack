
c
c   type_N_E_F = 0  => interiour ddl (ddl = Degree Of Freedom)
c   type_N_E_F != 0 => boundary ddl
c
c   i_cond = 0 => Dirichlet boundary condition (E-field: electric wall condition)
c   i_cond = 1 => Neumann boundary condition (E-field: magnetic wall condition)
c   i_cond = 2 => Periodic boundary condition
c
c

c     This subroutine set the boundary condition parameters

      subroutine bound_cond (i_cond, n_ddl, neq, type_N_E_F, ineq)

      implicit none
      integer*8 i_cond, n_ddl, neq
      integer*8 ineq(3,n_ddl), type_N_E_F(2,n_ddl)

      integer*8 i, i_boundary, i_dim
c
      if(i_cond .eq. 0) then
c       Dirichlet boundary condition: all points have a degree of freedom
        write(*,*) "bound_cond: Dirichlet boundary condition"
        neq = 0
        do i=1,n_ddl
          i_boundary = type_N_E_F(1,i)
          i_dim = type_N_E_F(2,i)
          if (i_dim .eq. 2) then ! each element is associated to 3 interior Degrees Of Freedom (DOF)
            ineq(1,i) = neq + 1
            ineq(2,i) = neq + 2
            ineq(3,i) = neq + 3
            neq = neq + 3
          elseif (i_dim .eq. 1) then ! each edge is associated to 3 Degrees Of Freedom (DOF)
            if (i_boundary .eq. 0) then
              ineq(1,i) = neq + 1
              ineq(2,i) = neq + 2
              ineq(3,i) = neq + 3
              neq = neq + 3
            else
              ineq(1,i) = 0
              ineq(2,i) = 0
              ineq(3,i) = 0
            endif
          elseif (i_dim .eq. 0) then ! each nodee is associated to 1 Degree Of Freedom (DOF)
            if (i_boundary .eq. 0) then
              ineq(1,i) = neq + 1
              ineq(2,i) = 0
              ineq(3,i) = 0
              neq = neq + 1
            else
              ineq(1,i) = 0
              ineq(2,i) = 0
              ineq(3,i) = 0
            endif
          else
            write(*,*) "bound_cond: i_dim has invalid value : ", i_dim
            write(*,*) "bound_cond: i_cond = ", i_cond
            write(*,*) "bound_cond: i = ", i
            write(*,*) "bound_cond: Aborting..."
            stop
          endif
        enddo
      elseif(i_cond .eq. 1) then
c       Neumann boundary condition: all points have a degree of freedom
        write(*,*) "bound_cond: Neumann boundary condition"
        neq = 0
        do i=1,n_ddl
          i_dim = type_N_E_F(2,i)
c         Each element or edge is associated to 3 Degrees Of Freedom (DOF)
          if (i_dim .eq. 2 .or. i_dim .eq. 1) then
            ineq(1,i) = neq + 1
            ineq(2,i) = neq + 2
            ineq(3,i) = neq + 3
            neq = neq + 3
          elseif (i_dim .eq. 0) then
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
        enddo
      else
        write(*,*) "bound_cond: i_cond has invalid value : ", i_cond
        write(*,*) "bound_cond: Aborting..."
        stop
      endif
c
      return
      end
