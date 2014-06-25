
c     Construct the FEM mesh
c
c   type_nod = 0  => interiour point
c   type_nod != 0 => boundary point
c
c
      subroutine geometry_1d (nel, npt, nnodes, nb_typ_el,
     *    lx, type_el, table_nod, x, 
     *    mesh_file)
c
      implicit none
      integer*8 nel, npt, nnodes, nb_typ_el
      integer*8 type_nod(npt), type_el(nel)
      integer*8 table_nod(nnodes,nel), nb_typ_el2
      double precision lx, x(npt)
c      complex*16 x(npt)
c      double precision xx, re_tmp, im_tmp

      character mesh_file*100

      integer*8 npt2, nel2, ui
      integer*8 i, j, k, debug
c
      ui = 6
      debug = 0
c
        open (unit=24,file="../backend/fortran/msh/"//mesh_file,
     *     status='old')
        read(24,*) npt2, nel2
c
      if(npt .ne. npt2) then
        write(ui,*) "geometry_1d: npt != npt2 : ",
     *    npt, npt2
      endif
      if(nel .ne. nel2) then
        write(ui,*) "geometry_1d: nel != nel2 : ",
     *    nel, nel2
      endif

c    Coordinate of the FEM points
      do i=1,npt
        read(24,*) k, x(i), type_nod(i)
C        x(i) = x(i) * lx
      enddo
c     Connectivity table
      nb_typ_el2 = 1
      do i=1,nel
        read(24,*) k, (table_nod(j,i),j=1,nnodes), type_el(i)
        j = type_el(i)
        if(nb_typ_el2 .lt. j) nb_typ_el2 = j
        if(j .lt. 0) then
          write(ui,*)
          write(ui,*) "   ???"
          write(ui,*) "geometry_1d: type_el(i) < 0 : ", 
     *    i, type_el(i)
          write(ui,*) "geometry_1d: Aborting..."
          stop
        endif
      enddo
      close(24)

      nb_typ_el = nb_typ_el2

      if (debug .eq. 1) then
        write(*,*) "geometry_1d: nb_typ_el = ", nb_typ_el
      endif

      return
      end

