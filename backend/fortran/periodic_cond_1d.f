c
c     Periodic boundary conditions
c
c
c***********************************************************************
c
      subroutine periodic_cond_1d (nel, npt, n_ddl, neq, table_nod_P2, 
     *   ineq, x_P2, table_ddl, ip_period_ddl, x_ddl, debug, period_x)

c
c***********************************************************************
c
c
c   Periodic boundary condition
c
c***********************************************************************
c
      implicit none
      integer*8 nel, npt, n_ddl, neq
      integer*8 table_ddl(3+4+4,nel), ip_period_ddl(n_ddl)
      integer*8 ineq(n_ddl), table_nod_P2(3,nel), debug
      double precision x_P2(npt), x_ddl(n_ddl), period_x

      integer*8 i, j, k, k1, i_el, n_ddl_0, n_ddl_0_bak, neq_0, ui
      double precision xmin, xmax, delta_x

c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
      ui = 6
      if (debug .eq. 1) then
        write(ui,*) "periodic_cond_1d: Periodic boundary conditions"
        write(ui,*) "nel, npt, n_ddl, neq = ", nel, npt, n_ddl, neq
      endif
      if (debug .eq. 1) then
        do i=1,n_ddl
          x_ddl(i) = -911
          ineq(i) = -911
          ip_period_ddl(i) = -911
        enddo
      endif
c
      if (debug .eq. 1) then
        write(ui,*) "periodic_cond_1d: writting the file "
        open (unit=24,file="Output/table_nod_P2.txt",status="unknown")
        do i=1,nel
          write(24,*) i, (table_nod_P2(j,i),j=1,3)
        enddo
        close(24)
      endif
c
      n_ddl_0 = 0
      do i_el=1,nel
        j = table_nod_P2(1,i_el)
        xmin = x_P2(j)
        j = table_nod_P2(2,i_el)
        xmax = x_P2(j)
        delta_x = xmax - xmin
        n_ddl_0_bak = n_ddl_0
c       Numbering the P2 node (x-component)
        n_ddl_0  = n_ddl_0 + 1
        table_ddl(1,i_el) = n_ddl_0   ! Discontinuous P2-node
        x_ddl(n_ddl_0) = xmin
        n_ddl_0  = n_ddl_0 + 1
        table_ddl(2,i_el) = n_ddl_0
        x_ddl(n_ddl_0) = xmax
        n_ddl_0  = n_ddl_0 + 1
        table_ddl(3,i_el) = n_ddl_0
        x_ddl(n_ddl_0) = xmin + delta_x / 2.0d0
c       Numbering the P3 node (y-component)
        if(i_el == 1) then
          n_ddl_0  = n_ddl_0 + 1
          table_ddl(4,i_el) = n_ddl_0
          x_ddl(n_ddl_0) = xmin
        else
          if(i_el == 2) then
          table_ddl(4,i_el) = n_ddl_0_bak - 4 - 2
        else
          table_ddl(4,i_el) = n_ddl_0_bak - 3 - 2
        endif
        endif
        n_ddl_0  = n_ddl_0 + 1
        table_ddl(5,i_el) = n_ddl_0
        x_ddl(n_ddl_0) = xmax
        n_ddl_0  = n_ddl_0 + 1
        table_ddl(6,i_el) = n_ddl_0
        x_ddl(n_ddl_0) = xmin + delta_x / 3.0d0
        n_ddl_0  = n_ddl_0 + 1
        table_ddl(7,i_el) = n_ddl_0
        x_ddl(n_ddl_0) = xmin + 2.0d0 * delta_x / 3.0d0
c       Numbering the P3 node (z-component)
        if(i_el == 1) then
          n_ddl_0  = n_ddl_0 + 1
          table_ddl(8,i_el) = n_ddl_0
          x_ddl(n_ddl_0) = xmin
        else
          table_ddl(8,i_el) = n_ddl_0_bak - 2
        endif
        n_ddl_0  = n_ddl_0 + 1
        table_ddl(9,i_el) = n_ddl_0
        x_ddl(n_ddl_0) = xmax
        n_ddl_0  = n_ddl_0 + 1
        table_ddl(10,i_el) = n_ddl_0
        x_ddl(n_ddl_0) = xmin + delta_x / 3.0d0
        n_ddl_0  = n_ddl_0 + 1
        table_ddl(11,i_el) = n_ddl_0
        x_ddl(n_ddl_0) = xmin + 2.0d0 * delta_x / 3.0d0
      enddo
      if (debug .eq. 1) then
        write(ui,*) "periodic_cond_1d: n_ddl_0 = ", n_ddl_0
      endif
      if(n_ddl /= n_ddl_0) then
         write(ui,*)
         write(ui,*) "   ???"
         write(ui,*) "periodic_cond_1d: n_ddl /= n_ddl_0 : ",
     *                n_ddl, n_ddl_0
         write(ui,*) "periodic_cond_1d: Aborting..."
         stop
      endif
c     DDL corresponding to periodic nodes
      do i=1,n_ddl
        ip_period_ddl(i) = 0
      enddo
      i_el = 1
        k = table_ddl(4,i_el)
          ip_period_ddl(k) = k
        k = table_ddl(8,i_el)
          ip_period_ddl(k) = k
      i_el = nel
        k = table_ddl(5,i_el)
        k1 = table_ddl(4,1)
          ip_period_ddl(k) = k1
          period_x = x_ddl(k) - x_ddl(k1)
        k = table_ddl(9,i_el)
        k1 = table_ddl(8,1)
          ip_period_ddl(k) = k1


c     Equation number for the DDL (degree of freedom)
      neq_0 = 0
      do i_el=1,nel
c       Numbering the P2 node (x-component)
        do j=1,3  ! P2-dicontinuous elements
          k = table_ddl(j,i_el)
          neq_0 = neq_0 + 1
          ineq(k) = neq_0
        enddo
c       Numbering the P3 node (y-component)
        j = 1
          k = table_ddl(j+3,i_el)
          neq_0 = neq_0 + 1
          ineq(k) = neq_0
c       The right end P3-node j=2 will be counted by the next element
        do j=3,4
          k = table_ddl(j+3,i_el)
          neq_0 = neq_0 + 1
          ineq(k) = neq_0
        enddo
c       Numbering the P3 node (z-component)
        j = 1
          k = table_ddl(j+3+4,i_el)
          neq_0 = neq_0 + 1
          ineq(k) = neq_0
c       The right end P3-node j=2 will be counted by the next element
        do j=3,4
          k = table_ddl(j+3+4,i_el)
          neq_0 = neq_0 + 1
          ineq(k) = neq_0
        enddo
      enddo
c     Periodic boundary node
c     Note:  No periodic boundary for the discontinuous P2 elements
      i_el = nel
c        j = 4
        j = 2
          k = table_ddl(j+3,i_el)
          k1 = table_ddl(4,1)
          ineq(k) = ineq(k1)
      i_el = nel
c        j = 4
        j = 2
          k = table_ddl(j+3+4,i_el)
          k1 = table_ddl(8,1)
          ineq(k) = ineq(k1)

      if (debug .eq. 1) then
        write(ui,*) "periodic_cond_1d: neq_0 = ", neq_0
      endif

      if(neq /= neq_0) then
         write(ui,*)
         write(ui,*) "   ???"
         write(ui,*) "periodic_cond_1d: neq /= neq_0 : ",
     *                neq, neq_0
         write(ui,*) "periodic_cond_1d: Aborting..."
         stop
      endif
c
      if (debug .eq. 1) then
        write(ui,*) "periodic_cond_1d: writting the file "
        open (unit=24,file="Output/x_ddl.txt",status="unknown")
          write(24,*) "i, x_ddl(i), ineq(i), ip_period_ddl(i)"
        do i=1,n_ddl
          write(24,*) i, x_ddl(i), ineq(i), ip_period_ddl(i)
        enddo
        close(24)
      endif
c
      if (debug .eq. 1) then
        write(ui,*) "periodic_cond_1d: writting the file "
        open (unit=24,file="Output/table_ddl.txt",status="unknown")
      do i_el=1,nel
          write(24,"(I6,' : ',20(I6))") i_el, (table_ddl(i,i_el),i=1,11)
        enddo
        close(24)
      endif
c
      return
      end
c
