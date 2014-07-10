
c   evecs(i) : contains the values of the solution for all points

      subroutine gmsh_plot_field_3d (lambda, h, nval,   
     *     E_H_field, nel, npt, nnodes, 
     *     type_el, nb_typ_el, table_nod, beta,
     *     evecs, vec_coef, x, gmsh_file_pos, extra_name)


      implicit none
      double precision lambda, h
      integer*8 nval, nel, npt, nnodes, E_H_field, nb_typ_el
      integer*8 table_nod(nnodes,nel), type_el(nel)
      double precision x(2,npt)
      complex*16 beta(nval), vec_coef(2*nval)
      complex*16 evecs(3,nnodes+7,nval,nel)
      character(*) gmsh_file_pos
      character(*) extra_name

Cf2py intent(in) lambda, h, nval, E_H_field, nel, npt, nnodes, 
Cf2py intent(in) type_el, nb_typ_el, table_nod
Cf2py intent(in) x, beta, vec_coef, evecs, gmsh_file_pos, extra_name

Cf2py depend(table_nod) nnodes, nel
Cf2py depend(type_el) nel
Cf2py depend(x) npt
Cf2py depend(beta) nval
Cf2py depend(vec_coef) nval
Cf2py depend(evecs) nnodes, nval, nel

c
c     Local variables
      character*20 dir_name
      integer alloc_stat
      complex*16, dimension(:,:,:,:), allocatable :: sol_3d
      integer*8, dimension(:), allocatable :: map_p1, inv_map_p1

      integer*8 n_prism_p1, n_prism_p2
      parameter (n_prism_p1=6, n_prism_p2 = 12)

      double precision xel(3,n_prism_p2), xel_p1(3,n_prism_p1)
      complex*16 sol_el(3,n_prism_p1)
      double precision sol_el_abs2(n_prism_p1)
      integer*8 table_nod_el_p1(n_prism_p1)

      complex*16 P_down, P_up, coef_down, coef_up, coef_t, coef_z
      complex*16 ii, z_tmp1
      double precision hz, dz, zz, r_tmp

      integer*8 nel_3d, npt_3d  ! prism elements
      integer*8 npt_p1  ! Number of vertices of the 2D FEM mesh
      integer*8 npt_h, i_h, npt_3d_p1  ! Resolution: number of points over the thickness h
      integer*8 i1, i, j, iel, inod, ival, debug, ui
      integer*8 namelen_gmsh, namelen_dir, namelen_tchar
      integer*8 namelen_extra

      integer*8 gmsh_type_prism, choice_type
      integer*8 number_tags, physic_tag
      integer*8 list_tag(6)
      integer*8, dimension(:), allocatable :: type_data

      character*100 tchar
      character*1 tE_H

      integer*8 iFrame, nFrame
      double precision pi, phase
      complex*16 exp_phase
      character tval*4, buf*3
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      dir_name = "3d_fields"
C
c  ii = sqrt(-1)
      ii = cmplx(0.0d0, 1.0d0)
c
      ui = 6
      debug = 0

      npt_h = 10 * CEILING((h/lambda + 1)) ! At least 10 nodes per wavelength
      npt_3d = npt * npt_h
      nel_3d = nel * (npt_h - 1)

      if (debug .eq. 1) then
        write(ui,*)
        write(ui,*) "gmsh_plot_field_3d: h/lambda = ", h/lambda
        write(ui,*) "gmsh_plot_field_3d:    npt_h = ", npt_h
        write(ui,*) "gmsh_plot_field_3d:   npt_3d = ", npt_3d
        write(ui,*) "gmsh_plot_field_3d:   nel_3d = ", nel_3d
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if ( nnodes .ne. 6 ) then
        write(ui,*) "gmsh_plot_field_3d: problem nnodes = ", nnodes
        write(ui,*) "nnodes should be equal to 6 !"
        write(ui,*) "gmsh_plot_field_3d: Aborting..."
        stop
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      alloc_stat = 0

      allocate(sol_3d(3,nnodes,nel,npt_h), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*)
        write(*,*) "gmsh_plot_field_3d: ",
     *     "The allocation is unsuccessful"
        write(*,*) "alloc_stat = ", alloc_stat
        write(*,*) "Not enough memory for the array sol_3d"
        write(*,*) "nel, npt_h = ", nel,npt_h
        write(*,*) "gmsh_plot_field_3d: Aborting..."
        stop
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    Initialise sol
      do i_h=1,npt_h
        do iel=1,nel
          do inod=1,nnodes
            do j=1,3
              sol_3d(j,inod,iel,i_h) = 0.0d0
            enddo
          enddo
        enddo
      enddo
c
      do ival=1,nval
        dz = h/dble(npt_h-1)
        do i_h=1,npt_h
          hz = (i_h-1)*dz
          P_down = EXP(ii*beta(ival)*hz)   !  Introduce Propagation in -z
          P_up = EXP(ii*beta(ival)*(h-hz)) !  Introduce Propagation in +z
          coef_down = vec_coef(ival) * P_down
          coef_up = vec_coef(ival+nval) * P_up
          coef_t = coef_up + coef_down
          coef_z = (coef_up - coef_down)/beta(ival) ! Taking into accout the change of variable for Ez
          do iel=1,nel
            do inod=1,nnodes
              do j=1,2
                z_tmp1 = evecs(j,inod,ival,iel) * coef_t
                sol_3d(j,inod,iel,i_h) = sol_3d(j,inod,iel,i_h) + z_tmp1
              enddo
              j=3
                z_tmp1 = evecs(j,inod,ival,iel) * coef_z
                sol_3d(j,inod,iel,i_h) = sol_3d(j,inod,iel,i_h) + z_tmp1
            enddo
          enddo
        enddo
      enddo
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (E_H_field .eq. 1) then
        tE_H = "E"
      elseif(E_H_field .eq. 2) then
        tE_H = "H"
      else
        write(ui,*) "gmsh_plot_field_3d: E_H_field has invalid value: ", 
     *    E_H_field
        write(ui,*) "Aborting..."
        stop
      endif

      namelen_gmsh = len_trim(gmsh_file_pos)
      namelen_dir = len_trim(dir_name)
      namelen_extra = len_trim(extra_name)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      tchar=dir_name(1:namelen_dir)// '/' // 
     *           extra_name(1:namelen_extra) //
     *           tE_H // '_abs2_3d_' //
     *           gmsh_file_pos(1:namelen_gmsh) // '.pos'
      open (unit=26,file=tchar)
        write(26,*) "View.IntervalsType = 3;"
        write(26,*) "View ""|",tE_H,"_t|^2 ", 
     *     " "" {"

      tchar=dir_name(1:namelen_dir)// '/' // 
     *           extra_name(1:namelen_extra) //
     *           tE_H // 'x_re_3d_' //
     *           gmsh_file_pos(1:namelen_gmsh) // '.pos'
      open (unit=27,file=tchar)
        write(27,*) "View.IntervalsType = 3;"
        write(27,*) "View ""Re ",tE_H,"x ", 
     *     " "" {"

      tchar=dir_name(1:namelen_dir)// '/' // 
     *           extra_name(1:namelen_extra) //
     *           tE_H // 'y_re_3d_' //
     *           gmsh_file_pos(1:namelen_gmsh) // '.pos'
      open (unit=28,file=tchar)
        write(28,*) "View.IntervalsType = 3;"
        write(28,*) "View ""Re ",tE_H,"y ", 
     *     " "" {" 

      tchar=dir_name(1:namelen_dir)// '/' // 
     *           extra_name(1:namelen_extra) //
     *           tE_H // 'z_re_3d_' //
     *           gmsh_file_pos(1:namelen_gmsh) // '.pos'
      open (unit=29,file=tchar)
        write(29,*) "View.IntervalsType = 3;"
        write(29,*) "View ""Re ",tE_H,"z ", 
     *     " "" {"

      tchar=dir_name(1:namelen_dir)// '/' // 
     *           extra_name(1:namelen_extra) //
     *           tE_H // 'v_re_3d_' // 
     *           gmsh_file_pos(1:namelen_gmsh) // '.pos'
      open (unit=30,file=tchar)
        write(30,*) "View.IntervalsType = 3;"
        write(30,*) "View.Axes = 2;"
        write(30,*) "View ""Re ",tE_H, " "" {"

      dz = h/dble(npt_h-1)
      do i_h=1,npt_h-1
        zz = - (i_h-1)*dz  ! hz=0 => top interface; hz=-h => bottom interface
        do iel=1,nel
          do inod=1,nnodes
            i1 = table_nod(inod,iel)
            xel(1,inod) = x(1,i1)
            xel(2,inod) = x(2,i1)
            xel(3,inod) = zz
            xel(1,inod+nnodes) = x(1,i1)  ! Prism element
            xel(2,inod+nnodes) = x(2,i1)
            xel(3,inod+nnodes) = zz - dz
          enddo
          do inod=1,3
            i1 = table_nod(inod,iel)
            xel_p1(1,inod) = x(1,i1)
            xel_p1(2,inod) = x(2,i1)
            xel_p1(3,inod) = zz
            xel_p1(1,inod+3) = x(1,i1)  ! Prism element
            xel_p1(2,inod+3) = x(2,i1)
            xel_p1(3,inod+3) = zz - dz
          enddo
          do inod=1,3
            sol_el_abs2(inod) = 0.0
            do j=1,3
              z_tmp1 = sol_3d(j,inod,iel,i_h)
              sol_el(j,inod) = z_tmp1
              sol_el_abs2(inod) = sol_el_abs2(inod) + 
     *           abs(z_tmp1)**2
            enddo
            sol_el_abs2(inod+3) = 0.0
            do j=1,3
              z_tmp1 = sol_3d(j,inod,iel,i_h+1)
              sol_el(j,inod+3) = z_tmp1
              sol_el_abs2(inod+3) = sol_el_abs2(inod+3) + 
     *           abs(z_tmp1)**2
            enddo
          enddo
          write(26,10) xel_p1, sol_el_abs2
          write(27,10) xel_p1, (dble(sol_el(1,j)),j=1,n_prism_p1)
          write(28,10) xel_p1, (dble(sol_el(2,j)),j=1,n_prism_p1)
          write(29,10) xel_p1, (dble(sol_el(3,j)),j=1,n_prism_p1)
          write(30,11) xel_p1, 
     *     ((dble(sol_el(i,j)),i=1,3),j=1,n_prism_p1)
        enddo
      enddo
      write(26,*) "};"
      write(27,*) "};"
      write(28,*) "};"
      write(29,*) "};"
      write(30,*) "};"
      close(26)
      close(27)
      close(28)
      close(29)
      close(30)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c         PART 2: GENERATE THE PRISM FINITE ELEMENT MESH IN GMSH FORMAT
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      alloc_stat = 0

      allocate(map_p1(npt), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*)
        write(*,*) "gmsh_plot_field_3d: ",
     *     "The allocation is unsuccessful"
        write(*,*) "alloc_stat = ", alloc_stat
        write(*,*) "Not enough memory for the array map_p1"
        write(*,*) "npt = ", npt
        write(*,*) "gmsh_plot_field_3d: Aborting..."
        stop
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do inod=1,npt
        map_p1(inod) = 0
      enddo
      npt_p1 = 0
      do iel=1,nel
        do inod=1,3  ! Scan the vertices
          i1 = table_nod(inod,iel)
          if(map_p1(i1) .eq. 0) then
            npt_p1 = npt_p1+1
            map_p1(i1) = npt_p1
          endif
        enddo
      enddo

      if (debug .eq. 1) then
        write(ui,*) "gmsh_plot_field_3d:    npt_p1 = ", npt_p1
        write(ui,*)
      endif

      allocate(inv_map_p1(npt_p1), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*)
        write(*,*) "gmsh_plot_field_3d: ",
     *     "The allocation is unsuccessful"
        write(*,*) "alloc_stat = ", alloc_stat
        write(*,*) "Not enough memory for the array inv_map_p1"
        write(*,*) "npt = ", npt
        write(*,*) "gmsh_plot_field_3d: Aborting..."
        stop
      endif

      allocate(type_data(nb_typ_el), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*)
        write(*,*) "gmsh_plot_field_3d: ",
     *     "The allocation is unsuccessful"
        write(*,*) "alloc_stat = ", alloc_stat
        write(*,*) "Not enough memory for the array type_data"
        write(*,*) "nb_typ_el = ", nb_typ_el
        write(*,*) "gmsh_plot_field_3d: Aborting..."
        stop
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do inod=1,npt
        if(map_p1(inod) .ne. 0) then
          inv_map_p1(map_p1(inod)) = inod
        endif
      enddo
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     elment type: defines the geometrical type
c     6-node prism
      gmsh_type_prism = 6


c      gmsh_version = 2
c     For choice_type = 3, the user provide a map for the types
      choice_type = 3
      number_tags = 6
      physic_tag = 4
      list_tag(2) = gmsh_type_prism
      list_tag(3) = number_tags - 3
      list_tag(5) = 1
      list_tag(6) = 0

      npt_3d_p1 = npt_p1 * npt_h

      if (nb_typ_el .eq. 1) then
        type_data(nb_typ_el) = 10
      else
        do i=1,nb_typ_el
          r_tmp = (dble(i-1)/dble(nb_typ_el-1))
          type_data(i) = int(1.0 + 18.0d0*r_tmp)
        enddo
      endif

      npt_3d_p1 = npt_p1 * npt_h
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      tchar=dir_name(1:namelen_dir)// '/' // 
     *           extra_name(1:namelen_extra) //
     *           '_abs2_3d_' //
     *           gmsh_file_pos(1:namelen_gmsh) // '.msh'
      open (unit=26,file=tchar)
      write(26,'(a11)') "$MeshFormat"
      write(26,'(3(I1,1x))') 2, 0, 8
      write(26,'(a14)') "$EndMeshFormat"
      write(26,'(a6)') "$Nodes"
      write(26,'(I0.1)') npt_3d_p1
      dz = h/dble(npt_h-1)
      do i_h=1,npt_h
        zz = - (i_h-1)*dz  ! hz=0 => top interface; hz=-h => bottom interface
        do inod=1,npt_p1
          i = inv_map_p1(inod)
          i1 = inod + (i_h-1) * npt_p1
          write(26,'(I0.1,3(g25.16))') i1, (dble(x(j,i)),j=1,2), zz
        enddo
      enddo
      write(26,'(a9)') "$EndNodes"
      write(26,'(a9)') "$Elements"
      write(26,'(I0.1)') nel_3d

      do i_h=1,npt_h-1
        do iel=1,nel
          list_tag(1) = iel + (i_h-1) * nel
          if (choice_type .eq. 1) then
            list_tag(physic_tag) = type_el(iel)
            list_tag(physic_tag+1) = type_el(iel)  
C          elseif (choice_type .eq. 2) then ! must change list_tag to double precision
C            typ_e = type_el(iel)           ! and reintroduce n_index, n_eff
C            r_index = int(n_eff(typ_e))
C            list_tag(physic_tag) = r_index
C            list_tag(physic_tag+1) = r_index
          elseif (choice_type .eq. 3) then
            list_tag(physic_tag) = type_data(type_el(iel))
            list_tag(physic_tag+1) = type_data(type_el(iel))
          else
            write(*,*) "gmsh_plot_field_3d: no action is defined ",
     *      "when choice_type = ", choice_type
            write(*,*) "gmsh_plot_field_3d: Aborting..."
            stop
          endif
          do inod=1,3  ! Scan the vertices
            i = table_nod(inod,iel)
            i1 = map_p1(i)
            table_nod_el_p1(inod) = i1 + (i_h-1) * npt_p1
            table_nod_el_p1(inod+3) = i1 + i_h * npt_p1
          enddo
          write(26,'(100(I0.1,2x))') (list_tag(j), j=1,number_tags), 
     *      (table_nod_el_p1(j), j=1,n_prism_p1)
        enddo
      enddo
      write(26,'(a12)') "$EndElements"
      close(26)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c         PART 2: GENERATE ANIMATION FILES
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      pi = 3.141592653589793d0
      nFrame = 50

      tchar=dir_name(1:namelen_dir)// '/' // 
     *           "anim/" 
     *           // 'view_' // tE_H // '_v_3d.geo'
      open (unit=29,file=tchar)

      do iFrame=1,nFrame
        phase = 2.0d0 * pi * (iFrame -1) / dble(nFrame)
c       Exponential time dependence: Exp(-i omega t)
        exp_phase = exp(- ii *  phase)
        if (iFrame .lt. 1000) then
          write(buf,'(i3.3)') iFrame
        else
          buf = "nnn"
        endif
        tval=tE_H//buf
        tchar = 'view_' // tval // '_v_3d.pos'
        namelen_tchar = len_trim(tchar)
        write(29,*) " Include """, tchar(1:namelen_tchar), """;"
        tchar=dir_name(1:namelen_dir)// '/' // 
     *           "anim/" 
     *           // 'view_' // tval // '_v_3d.pos'
        open (unit=30,file=tchar)
        write(30,*) "View.IntervalsType = 3;"
        write(30,*) "View.Axes = 2;"
        write(30,*) "View ""Vector field ",tE_H, " "" {"
        dz = h/dble(npt_h-1)
        do i_h=1,npt_h-1
          zz = - (i_h-1)*dz  ! hz=0 => top interface; hz=-h => bottom interface
          do iel=1,nel
            do inod=1,nnodes
              i1 = table_nod(inod,iel)
              xel(1,inod) = x(1,i1)
              xel(2,inod) = x(2,i1)
              xel(3,inod) = zz
              xel(1,inod+nnodes) = x(1,i1)  ! Prism element
              xel(2,inod+nnodes) = x(2,i1)
              xel(3,inod+nnodes) = zz - dz
            enddo
            do inod=1,3
              i1 = table_nod(inod,iel)
              xel_p1(1,inod) = x(1,i1)
              xel_p1(2,inod) = x(2,i1)
              xel_p1(3,inod) = zz
              xel_p1(1,inod+3) = x(1,i1)  ! Prism element
              xel_p1(2,inod+3) = x(2,i1)
              xel_p1(3,inod+3) = zz - dz
            enddo
            do inod=1,3
              do j=1,3
                z_tmp1 = sol_3d(j,inod,iel,i_h)
                sol_el(j,inod) = z_tmp1 * exp_phase
              enddo
              do j=1,3
                z_tmp1 = sol_3d(j,inod,iel,i_h+1)
                sol_el(j,inod+3) = z_tmp1 * exp_phase
              enddo
            enddo
          write(30,11) xel_p1, 
     *     ((dble(sol_el(i,j)),i=1,3),j=1,n_prism_p1)
          enddo
        enddo
        write(30,*) "};"
        close(30)
      enddo
      write(29,*) "Combine TimeStepsFromAllViews;"
      close(29)

c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     SI : Scalar prism
10    format("SI(",f10.6,17(",",f10.6),"){",
     *     g24.16,5(",",g24.16),"};")

c     VI : Vector prism
11    format("VI(",f10.6,17(",",f10.6),"){",
     *     g24.16,17(",",g24.16),"};")

c
      deallocate(sol_3d, map_p1, inv_map_p1, type_data)
c
      return
      end  subroutine gmsh_plot_field_3d




