c*******************************************************
c
c     conv_gmsh: covert the GMSH mesh format to the FEM mesh format
c
c*******************************************************
c
      subroutine conv_gmsh(geoname)
c
      implicit none
c
      integer i_sym
      integer gmsh_version
      character file0_mesh*500, geoname*500
      character file1_mesh*500, file2_mesh*500
      integer i_mesh(3)
      integer max_ne, max_npt, i_err
      parameter(max_npt=250000, max_ne=120000)
      integer nu_d1(3,max_ne), typ_el_d1(max_ne)
      integer nu_d2(6,max_ne), typ_el_d2(max_ne)
      integer idfn(max_npt)
      integer number_tags, physic_tag
      integer elm_type(max_ne)
      integer tab_ns(max_npt), tab_el(max_ne)
      double precision x(max_npt),  y(max_npt)
c
      integer ne, npt
      integer ne_d1, ne_d2
      integer i, j, k, etc_1(10)
      integer gmsh_type_line, gmsh_type_el
c
      double precision tmp1, tmp2, tmp3
c
      integer debug, ui, namelength2
      double precision time1, time2
      character objet*5
      character file_ui*500
      character*1000 com_line
C
Cf2py intent(in) geoname
c
ccccccccccccccccccccccccc
c
      debug = 0
      ui = 66
      gmsh_version = 2

      namelength2 = len_trim(geoname)
      if (namelength2 .ge. 500) then
        write(*,*) "Name of .geo file is too long extend in ",
     *  "conv_gmsh_py.f"
      endif
C
      file0_mesh = geoname(1:namelength2)//".geo"
      file1_mesh = geoname(1:namelength2)//".msh"
      file2_mesh = geoname(1:namelength2)//".mail"
      file_ui    = geoname(1:namelength2)//".log"
C
      com_line = "gmsh -0 -2  -order 2 -v 0 -o " //
     *    file1_mesh // " " // file0_mesh
C
      call system(com_line)
c
      i_sym = 0
c
      call cpu_time(time1)  ! initial time  in unit = sec./100 (For CPU time)

      gmsh_type_line = 8
      gmsh_type_el = 9

      if(gmsh_version .eq. 2) then
        number_tags = 5 !formerly 6 on windows gmsh 2.5.0
        physic_tag = 4
      else
        number_tags = 5
        physic_tag = 3
      endif
c
c   Initialisation
      i_mesh(1) = -1
      i_mesh(2) = -1
      i_mesh(3) = -1
c
      i_err = 0
      open (unit=24,file=file1_mesh)
        if(gmsh_version .eq. 2) then
          read(24,'(a1)') objet
          read(24,'(a1)') objet
          read(24,'(a1)') objet
        endif
        read(24,'(a5)') objet
        read(24,*) npt
      close(24)
c
      if(max_npt .lt. npt) then
        open (unit=ui,file=file_ui)
        write(*,*) 'CONV_GMSH: ATTENTION: max_npt < npt',
     *              max_npt, npt
        i_err = -1
        close(ui)
        stop
      endif
c
c      write(*,*) 'npt = ', npt
c
      open (unit=24,file=file1_mesh)
        if(gmsh_version .eq. 2) then
          read(24,'(a1)') objet
          read(24,'(a1)') objet
          read(24,'(a1)') objet
        endif
        read(24,'(a5)') objet
        read(24,*) j

      do i=1,npt
        read(24,*) j, x(i), y(i), tmp1
        tab_ns(j) = i
      enddo
c
      read(24,'(a5)') objet
      read(24,'(a5)') objet
      read(24,*) ne
c
      if(max_ne .lt. ne) then
        open (unit=ui,file=file_ui)
        write(*,*) 'CONV_GMSH: ATTENTION: max_ne < ne', max_ne, ne
        i_err = -2
        close(ui)
        stop
      endif
c
      do i=1,ne
        read(24,*) j, elm_type(i)
        tab_el(j) = i
      enddo
c
      close(24)
c
      open (unit=25,file=file1_mesh)
c
c
c     On saute les lignes deja traitees
        if(gmsh_version .eq. 2) then
          read(25,'(a1)') objet
          read(25,'(a1)') objet
          read(25,'(a1)') objet
        endif
        read(25,'(a5)') objet
        read(25,*) j
      do i=1,npt
        read(25,*) j, tmp1, tmp2, tmp3
      enddo
      read(25,'(a5)') objet
      read(25,'(a5)') objet
      read(25,*) j
c
      ne_d1 = 0
      ne_d2 = 0
      do i=1,ne
        if(elm_type(i) .eq. gmsh_type_line) then
          ne_d1 = ne_d1 + 1
          read(25,*) (etc_1(k), k=1,number_tags),
     *      (nu_d1(k,ne_d1), k=1,3)
          do k=1,3
            j = nu_d1(k,ne_d1)
            nu_d1(k,ne_d1) = tab_ns(j)
          enddo
          typ_el_d1(ne_d1) = etc_1(physic_tag)
        elseif(elm_type(i) .eq. gmsh_type_el) then
          ne_d2 = ne_d2 + 1
          read(25,*) (etc_1(k), k=1,number_tags),
     *      (nu_d2(k,ne_d2), k=1,6)
          do k=1,6
            j = nu_d2(k,ne_d2)
            nu_d2(k,ne_d2) = tab_ns(j)
          enddo
          typ_el_d2(ne_d2) = etc_1(physic_tag)
        else
        open (unit=ui,file=file_ui)
        write(*,*) '?? MAIN: elm_type(i), i = ', elm_type(i), i
        close(ui)
          i_err = -5
          stop
        endif
      enddo
c
      close(25)
c
      do i=1,npt
        idfn(i) = 0
      enddo
c
      do i=1,ne_d1
        j = typ_el_d1(i)
        k = nu_d1(1,i) ! tab_ns(nu_d1(1,i))
        idfn(k) = j
        k = nu_d1(2,i) ! tab_ns(nu_d1(2,i))
        idfn(k) = j
        k = nu_d1(3,i)
        idfn(k) = j
      enddo
c
c      call matlab(npt,ne_d2,nu_d2,idfn,typ_el_d2,x,y)
c      stop
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      print*, 'Appel de renumerote'
c
      if(i_sym .ne. 0) then
        call symmetry(npt, ne_d2,
     *      max_ne, max_npt, idfn, nu_d2, typ_el_d2,
     *      x, y, i_sym)
      endif
c
      call renumerote(npt, ne_d2, idfn, nu_d2, x, y, ui)
c
      open (unit=26,file=file2_mesh)
c
      write(26,*) npt,ne_d2
c
      do i=1,npt
        write(26,88) i, x(i), y(i), idfn(i)
      enddo
c
      do i=1,ne_d2
        write(26,*) i, (nu_d2(k,i), k=1,6), typ_el_d2(i)
      enddo
      close(26)
88    format(2x, i7, 2(g25.15), i6)
c
c      call test_border(npt, idfn, x, y, i_err,
c     *                   file1_mesh, ui)
c
c
c###############################################
c	FEM mesh plot by Matlab
c
c      call matlab(npt,ne_d2,nu_d2,idfn,typ_el_d2,x,y)
c
      If(debug .eq. 1) then
        call cpu_time(time2)  ! temps initial (sert pour la durree des calcul)
        open (unit=ui,file=file_ui)
        write(ui,*) "conv_gmsh_m: debug = ", debug
        write(ui,*) "gmsh_version = ", gmsh_version
        write(ui,*) "i_mesh = ", i_mesh
        write(ui,*) " file1_mesh = ", file1_mesh
        write(ui,*) " file2_mesh = ", file2_mesh
        write(ui,*) " i_sym = ", i_sym
        write(ui,*) 'Number of points = ', npt
        write(ui,*) 'Number of elements = ', ne_d2
        write(ui,*) 'The program terminates normally'
        write(ui,*) 'Symmetry code = ', i_sym
        write(ui,*) 'CPU time (sec.)           = ', time2-time1
c     *                                   dble(time2-time1)/100.0
        close(ui)
      EndIf

      continue

      i_mesh(1) = npt
      i_mesh(2) = ne_d2
      i_mesh(3) = gmsh_version

C      stop
      return
      end


c##################################################################################
