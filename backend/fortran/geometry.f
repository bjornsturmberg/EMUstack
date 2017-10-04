
c     Construct the FEM mesh
c
c   type_nod = 0  => interiour point
c   type_nod != 0 => boundary point
c
c
      subroutine geometry (nel, npt, nnodes, nb_typ_el,
     *    lx, ly, type_nod, type_el, table_nod, x,
     *    mesh_file)
c
      implicit none
      integer*8 nel, npt, nnodes, nb_typ_el, max_typ_el
      integer*8 type_nod(npt), type_el(nel)
      integer*8 table_nod(nnodes,nel), nb_typ_el2
      double precision lx, ly
      parameter (max_typ_el=10)
      double precision x(2,npt), xx(2)

      character mesh_file*500

      integer*8 npt2, nel2, ui
      integer*8 i, j, k
c
      ui = 6
c
        open (unit=24,file=mesh_file,
     *     status='old')
        read(24,*) npt2, nel2
c
      if(npt .ne. npt2) then
        write(ui,*) "geometry: npt != npt2 : ",
     *    npt, npt2
      endif
      if(nel .ne. nel2) then
        write(ui,*) "geometry: nel != nel2 : ",
     *    nel, nel2
      endif

c    Coordinate of the FEM points
      do i=1,npt
        read(24,*) k, (xx(j),j=1,2), type_nod(i)
        x(1,i) = xx(1)*lx
        x(2,i) = xx(2)*ly
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
          write(ui,*) "geometry: type_el(i) < 0 : ",
     *    i, type_el(i)
          write(ui,*) "geometry: Aborting..."
          stop
        endif
      enddo
      close(24)

      if(nb_typ_el2 .gt. nb_typ_el) then
         write(ui,*)
         write(ui,*) "   ???"
         write(ui,*) "geometry: nb_typ_el2 > nb_typ_el : ",
     *    nb_typ_el2, nb_typ_el
         write(ui,*) "geometry: Aborting..."
         stop
      endif

      return
      end
