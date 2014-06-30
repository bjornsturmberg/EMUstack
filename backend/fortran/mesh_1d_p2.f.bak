
      subroutine mesh_1d_p2(core_radius, nel, mesh_file)

      implicit none
      
      integer*8 i, nel, npt, i_el, allocate_status
C      parameter (nel = 500)        
C      parameter (npt = 2 * nel + 1)
      double precision, allocatable :: ls_x(:)
      integer*8, allocatable :: type_nod(:), type_el(:), table_nod(:,:)
      double precision x, x_1, x_2
      double precision x_min, x_max, delta_x, core_radius
      character mesh_file*100

Cf2py intent(in) nel, mesh_file, core_radius

      npt = 2 * nel + 1

      allocate_status = 0
      allocate(ls_x(npt), STAT=allocate_status)
      if (allocate_status /= 0) then
        write(*,*) "The allocation is unsuccessful"
        write(*,*) "allocate_status = ", allocate_status
        write(*,*) "Not enough memory for ls_x"
        write(*,*) "ls_x = ", ls_x
        write(*,*) "Aborting..."
        stop
      endif
      allocate(type_nod(npt), STAT=allocate_status)
      if (allocate_status /= 0) then
        write(*,*) "The allocation is unsuccessful"
        write(*,*) "allocate_status = ", allocate_status
        write(*,*) "Not enough memory for type_nod"
        write(*,*) "type_nod = ", type_nod
        write(*,*) "Aborting..."
        stop
      endif
      allocate(type_el(nel), STAT=allocate_status)
      if (allocate_status /= 0) then
        write(*,*) "The allocation is unsuccessful"
        write(*,*) "allocate_status = ", allocate_status
        write(*,*) "Not enough memory for type_el"
        write(*,*) "type_el = ", type_el
        write(*,*) "Aborting..."
        stop
      endif
      allocate(table_nod(3,nel), STAT=allocate_status)
      if (allocate_status /= 0) then
        write(*,*) "The allocation is unsuccessful"
        write(*,*) "allocate_status = ", allocate_status
        write(*,*) "Not enough memory for table_nod"
        write(*,*) "table_nod = ", table_nod
        write(*,*) "Aborting..."
        stop
      endif

      x_min = -0.5d0
      x_max =  0.5d0
      delta_x = (x_max - x_min) / dble(nel)
C      core_radius = (x_max - x_min)/10.0d0

c      Coordinate and type of the nodes
      do i_el=1,nel
        x = x_min + (i_el-1) * delta_x
        ls_x(2*i_el-1) = x
        ls_x(2*i_el) = x + delta_x / 2.0d0
        type_nod(2*i_el-1) = 0
        type_nod(2*i_el) = 0
      enddo
c      End-points
      i_el = nel
        x = x_min + i_el * delta_x
        ls_x(2*i_el+1) = x
        type_nod(2*i_el+1) = 1
      i = 1
      type_nod(i) = 1
c
c     Connectivity table
      do i_el=1,nel
        table_nod(1,i_el) = 2*i_el-1
        table_nod(2,i_el) = 2*i_el+1
        table_nod(3,i_el) = 2*i_el  ! Mid-node
c       End-points of the elements
        x_1 = ls_x(2*i_el-1)
        x_2 = ls_x(2*i_el+1)
        if (abs(x_1) <= core_radius .and. 
     *          abs(x_2) <= core_radius) then
          type_el(i_el) = 2
        else
          type_el(i_el) = 1
        endif
      enddo


      open(3,file = mesh_file, status='unknown')
        write(3,*) npt, nel
      do i=1,npt
        write(3,*) i, ls_x(i), type_nod(i)
      enddo
      do i_el=1,nel
        write(3,*) i_el, (table_nod(i,i_el),i=1,3), type_el(i_el)
      enddo
      close(3)

      deallocate(ls_x, type_nod, type_el, table_nod)

      end

