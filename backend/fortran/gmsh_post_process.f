c  Write an ASCII file that can be read into GMSH.
c  P1-element is used to represent the  3D vector field.
c  P2-element is used to represent each component of the 3D vector field.
c

      subroutine gmsh_post_process (plot_val, E_H_field, nval,
     *     nel, npt, nnodes, table_nod, type_el, nb_typ_el,
     *     n_eff, x, val_cmplx, sol, visite,
     *     gmsh_file_pos, dir_name,
     *     q_average, plot_real, plot_imag, plot_abs)

      implicit none
      integer*8 nval, nel, npt, nnodes, plot_val, E_H_field
      integer*8 nb_typ_el
      integer*8 table_nod(nnodes,nel), type_el(nel)
      integer*8 visite(npt)
      double precision x(2,npt)
      complex*16 sol(3,nnodes+7,nval,nel), n_eff(nb_typ_el)
      integer alloc_stat
      complex*16, dimension(:,:), allocatable :: sol_avg

      complex*16 val_cmplx(nval)

      integer*8 nnodes_0
      parameter (nnodes_0 = 6)

      double precision xel(3,nnodes_0), xel_p1(3,3)
      complex*16 sol_el(3,nnodes_0), sol_max(4)
      double precision sol_el_abs2(nnodes_0)
      double precision sol_el_abs2_eE(nnodes_0)
C      double precision sol_el_abs2_iD(nnodes_0)
      double precision ls_index(nnodes_0), r_index, zz
C      double precision ls_im_index(nnodes_0), im_index
C      double precision ls_abs_index(nnodes_0), abs_index
      double precision v_im, v_re

      integer*8 i, j, i1, iel, namelen, namelen2, typ_e
      integer*8 q_average, plot_imag, plot_real, plot_abs
      integer*8 debug, ui
      complex*16 z_tmp1
      complex*16 ii
      character*(*) gmsh_file_pos, dir_name
      character*100 tchar
      character tval*4, buf*3
      character*1 tE_H
      integer*8 namelength
c
c  ii = sqrt(-1)
      ii = cmplx(0.0d0, 1.0d0)
c
      ui = 6
      debug = 0
c
      if ( nnodes .ne. 6 ) then
        write(ui,*) "gmsh_post_process: problem nnodes = ", nnodes
        write(ui,*) "gmsh_post_process: nnodes should be equal to 6 !"
        write(ui,*) "gmsh_post_process: Aborting..."
        stop
      endif
C
      alloc_stat = 0
      allocate(sol_avg(3,npt), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "gmsh_post_process: Mem. allocation is unseccesfull"
        write(*,*) "alloc_stat (sol_avg) = ", alloc_stat
        write(*,*) "Aborting..."
        stop
      endif
C

      if (plot_val .eq. 0) return
c
      if (E_H_field .eq. 1) then
        tE_H = "E"
      elseif(E_H_field .eq. 2) then
        tE_H = "H"
      else
        write(ui,*) "gmsh_post_process: E_H_field has invalid value: ",
     *    E_H_field
        write(ui,*) "Aborting..."
        stop
      endif
c
      if (q_average .eq. 1) then
        do i=1,npt
          visite(i) = 0
          do j=1,3
            sol_avg(j,i) = 0.0d0
          enddo
        enddo
        do iel=1,nel
          do i=1,nnodes
            i1 = table_nod(i,iel)
            visite(i1) = visite(i1) + 1
            do j=1,3
              sol_avg(j,i1) = sol_avg(j,i1) +
     *          sol(j,i,plot_val,iel)
            enddo
          enddo
        enddo
        do i=1,npt
          if (visite(i) .eq. 0) then
            write(ui,*) "gmsh_post_process: visite(i) = 0"
            write(ui,*) " i, visite(i) = ", i, visite(i)
            write(ui,*) "gmsh_post_process: Aborting..."
            stop
          endif
          do j=1,3
            sol_avg(j,i) = sol_avg(j,i)/dble(visite(i))
          enddo
        enddo
      endif
c
      v_re = dble(val_cmplx(plot_val))
      v_im = -dble(ii*val_cmplx(plot_val))
c
      if (plot_val .lt. 1000) then
        write(buf,'(i3.3)') plot_val
      else
        buf = "nnn"
      endif
      tval=tE_H//buf

      namelen = len_trim(gmsh_file_pos)
      namelength = len_trim(dir_name)
c
c###############################################
c
      if (plot_real .eq. 1) then

c    All_plots.geo (unit=34)
C      if (plot_val .eq. 1) then
C      tchar = """../../"//dir_name(1:namelength)// "/"
C     *  //"interface_c4.geo"";"
C      namelen2 = len_trim(tchar)
C        write(34,*) "Merge ", tchar(1:namelen2)
C      else
C        write(34,*)
C        write(34,*) "Delete View[0];"
C      endif
C
C
C Convert from gmsh format to pdf
      write(34,*) "Delete View[0];"
      tchar = '../../'//dir_name(1:namelength)// '/'
     *  //gmsh_file_pos(1:namelen) // '_' // tval // '_abs2_eE.pos'
      namelen2 = len_trim(tchar)
      write(34,*) " Include """, tchar(1:namelen2), """;"
      tchar = gmsh_file_pos(1:namelen) // '_' // tval // '_abs2_eE.pdf'
      namelen2 = len_trim(tchar)
      write(34,*) "Print Sprintf(""", tchar(1:namelen2), """);"
C
      write(34,*) "Delete View[0];"
      tchar = '../../'//dir_name(1:namelength)// '/'
     *  //gmsh_file_pos(1:namelen) // '_' // tval // 'v_re.pos'
      namelen2 = len_trim(tchar)
      write(34,*) " Include """, tchar(1:namelen2), """;"
      tchar = gmsh_file_pos(1:namelen) // '_' // tval // 'v_re.pdf'
      namelen2 = len_trim(tchar)
      write(34,*) "Print Sprintf(""", tchar(1:namelen2), """);"
C
      write(34,*) "Delete View[0];"
      tchar = '../../'//dir_name(1:namelength)// '/'
     *  //gmsh_file_pos(1:namelen) // '_' // tval // 'x_re.pos'
      namelen2 = len_trim(tchar)
      write(34,*) " Include """, tchar(1:namelen2), """;"
      tchar = gmsh_file_pos(1:namelen) // '_' // tval // 'x_re.pdf'
      namelen2 = len_trim(tchar)
      write(34,*) "Print Sprintf(""", tchar(1:namelen2), """);"
C
      write(34,*) "Delete View[0];"
      tchar = '../../'//dir_name(1:namelength)// '/'
     *  //gmsh_file_pos(1:namelen) // '_' // tval // 'y_re.pos'
      namelen2 = len_trim(tchar)
      write(34,*) " Include """, tchar(1:namelen2), """;"
      tchar = gmsh_file_pos(1:namelen) // '_' // tval // 'y_re.pdf'
      namelen2 = len_trim(tchar)
      write(34,*) "Print Sprintf(""", tchar(1:namelen2), """);"
C
      write(34,*) "Delete View[0];"
      tchar = '../../'//dir_name(1:namelength)// '/'
     *  //gmsh_file_pos(1:namelen) // '_' // tval // 'z_re.pos'
      namelen2 = len_trim(tchar)
      write(34,*) " Include """, tchar(1:namelen2), """;"
      tchar = gmsh_file_pos(1:namelen) // '_' // tval // 'z_re.pdf'
      namelen2 = len_trim(tchar)
      write(34,*) "Print Sprintf(""", tchar(1:namelen2), """);"
C
      tchar=dir_name(1:namelength)// '/' // gmsh_file_pos(1:namelen)
     *           // '_' // tval // '_abs2.geo'
      open (unit=26,file=tchar)
      tchar = gmsh_file_pos(1:namelen) // '_' // tval // '_abs2.pos'
        namelen2 = len_trim(tchar)
        write(26,*) " Include """, tchar(1:namelen2), """;"
C        write(26,*) "Merge ""interface_c4.geo"";"
      close (unit=26)
C
      tchar=dir_name(1:namelength)// '/' // gmsh_file_pos(1:namelen)
     *           // '_' // tval // '_abs2_eE.geo'
      open (unit=32,file=tchar)
      tchar = gmsh_file_pos(1:namelen) // '_' // tval // '_abs2_eE.pos'
        namelen2 = len_trim(tchar)
        write(32,*) " Include """, tchar(1:namelen2), """;"
C        write(32,*) "Merge ""interface_c4.geo"";"
      close (unit=32)

C      tchar=dir_name(1:namelength)// '/' // gmsh_file_pos(1:namelen)
C     *           // '_' // tval // '_abs2_iE.geo'
C      open (unit=33,file=tchar)
C      tchar = gmsh_file_pos(1:namelen) // '_' // tval // '_abs2_eE.pos'
C        namelen2 = len_trim(tchar)
C        write(33,*) " Include """, tchar(1:namelen2), """;"
C        write(33,*) "Merge ""interface_c4.geo"";"
C      close (unit=33)

      tchar=dir_name(1:namelength)// '/' // gmsh_file_pos(1:namelen)
     *           // '_' // tval // '_abs2.pos'
      open (unit=26,file=tchar)
        write(26,*) "View.ArrowSizeMax = 40;"
        write(26,*) "View.ArrowSizeMin = 0;"
        write(26,*) "General.Color.Background = {255,255,255};"
        write(26,*) "General.Color.BackgroundGradient={255,255,255};"
        write(26,*) "View.AdaptVisualizationGrid =1;"
        write(26,*) "View.MaxRecursionLevel = 2;"
        write(26,*) "View.IntervalsType = 3;"
        write(26,*) "View.Light = 0;"
c        write(26,*) "View ""|E|^2: n = ", "|",tE_H,"|^2",
c        write(26,*) "View ""|E|^2: n = ", "|",tE_H,"|^2: n = ",
        write(26,*) "View ""|",tE_H,"|^2: n = ",
     *     plot_val, ", beta_n =", v_re, "+ I *", v_im, " "" {"
c
      tchar=dir_name(1:namelength)// '/' // gmsh_file_pos(1:namelen)
     *            // '_' // tval // 'x_re.pos'
      open (unit=27,file=tchar)
        write(27,*) "View.ArrowSizeMax = 40;"
        write(27,*) "View.ArrowSizeMin = 0;"
        write(27,*) "General.Color.Background = {255,255,255};"
        write(27,*) "General.Color.BackgroundGradient={255,255,255};"
        write(27,*) "View.AdaptVisualizationGrid =1;"
        write(27,*) "View.MaxRecursionLevel = 2;"
        write(27,*) "View.IntervalsType = 3;"
        write(27,*) "View.Light = 0;"
c        write(27,*) "View ""Re Ex: n = ",
        write(27,*) "View ""Re ",tE_H,"x: n = ",
     *     plot_val, ", beta_n =", v_re, "+ I *", v_im, " "" {"
c
      tchar=dir_name(1:namelength)// '/' // gmsh_file_pos(1:namelen)
     *           // '_' // tval // 'y_re.pos'
      open (unit=28,file=tchar)
        write(28,*) "View.ArrowSizeMax = 40;"
        write(28,*) "View.ArrowSizeMin = 0;"
        write(28,*) "General.Color.Background = {255,255,255};"
        write(28,*) "General.Color.BackgroundGradient={255,255,255};"
        write(28,*) "View.AdaptVisualizationGrid =1;"
        write(28,*) "View.MaxRecursionLevel = 2;"
        write(28,*) "View.IntervalsType = 3;"
        write(28,*) "View.Light = 0;"
        write(28,*) "View ""Re ",tE_H,"y: n = ",
     *     plot_val, ", beta_n =", v_re, "+ I *", v_im, " "" {"
c
      tchar=dir_name(1:namelength)// '/' // gmsh_file_pos(1:namelen)
     *           // '_' // tval // 'z_re.pos'
      open (unit=29,file=tchar)
        write(29,*) "View.ArrowSizeMax = 40;"
        write(29,*) "View.ArrowSizeMin = 0;"
        write(29,*) "General.Color.Background = {255,255,255};"
        write(29,*) "General.Color.BackgroundGradient={255,255,255};"
        write(29,*) "View.AdaptVisualizationGrid =1;"
        write(29,*) "View.MaxRecursionLevel = 2;"
        write(29,*) "View.IntervalsType = 3;"
        write(29,*) "View.Light = 0;"
        write(29,*) "View ""Re ",tE_H,"z: n = ",
     *     plot_val, ", beta_n =", v_re, "+ I *", v_im, " "" {"
c
      tchar=dir_name(1:namelength)// '/' // gmsh_file_pos(1:namelen)
     *           // '_' // tval // 'v_re.pos'
      open (unit=30,file=tchar)
        write(30,*) "View.ArrowSizeMax = 40;"
        write(30,*) "View.ArrowSizeMin = 0;"
        write(30,*) "General.Color.Background = {255,255,255};"
        write(30,*) "General.Color.BackgroundGradient={255,255,255};"
        write(30,*) "View.AdaptVisualizationGrid =1;"
        write(30,*) "View.MaxRecursionLevel = 2;"
        write(30,*) "View.IntervalsType = 3;"
        write(30,*) "View.Light = 0;"
        write(30,*) "View ""Re ",tE_H,": n = ",
     *     plot_val, ", beta_n =", v_re, "+ I *", v_im, " "" {"
c
      tchar=dir_name(1:namelength)// '/' // gmsh_file_pos(1:namelen)
     *           //'_ind.pos'
      open (unit=31,file=tchar)
        write(31,*) "View.ArrowSizeMax = 40;"
        write(31,*) "View.ArrowSizeMin = 0;"
        write(31,*) "General.Color.Background = {255,255,255};"
        write(31,*) "General.Color.BackgroundGradient={255,255,255};"
        write(31,*) "View.AdaptVisualizationGrid =1;"
        write(31,*) "View.MaxRecursionLevel = 2;"
        write(31,*) "View.IntervalsType = 3;"
        write(31,*) "View.Light = 0;"
        write(31,*) "View ""Refrac. index ", " "" {"
C
      tchar=dir_name(1:namelength)// '/' // gmsh_file_pos(1:namelen)
     *           // '_' // tval // '_abs2_eE.pos'
      open (unit=32,file=tchar)
        write(32,*) "View.ArrowSizeMax = 40;"
        write(32,*) "View.ArrowSizeMin = 0;"
        write(32,*) "General.Color.Background = {255,255,255};"
        write(32,*) "General.Color.BackgroundGradient={255,255,255};"
        write(32,*) "View.AdaptVisualizationGrid =1;"
        write(32,*) "View.MaxRecursionLevel = 2;"
        write(32,*) "View.IntervalsType = 3;"
        write(32,*) "View.Light = 0;"
        write(32,*) "View ""|",tE_H,"|^2: n = ",
     *     plot_val, ", beta_n =", v_re, "+ I *", v_im, " "" {"
C
C      tchar=dir_name(1:namelength)// '/' // gmsh_file_pos(1:namelen)
C     *           // '_' // tval // '_abs2_iD.pos'
C      open (unit=33,file=tchar)
C        write(33,*) "View.AdaptVisualizationGrid =1;"
C        write(33,*) "View.IntervalsType = 3;"
C        write(33,*) "View.Light = 0;"
C        write(33,*) "View ""|",tE_H,"|^2: n = ",
C     *     plot_val, ", beta_n =", v_re, "+ I *", v_im, " "" {"
C
        sol_max(4) = 0.0d0
        do iel=1,nel
          typ_e = type_el(iel)
          r_index = real(n_eff(typ_e))
C          im_index = imag(sqrt(eps_eff(typ_e)))
C          abs_index = abs(sqrt(eps_eff(typ_e)))
          zz = 0.0d0
          do i=1,nnodes
            i1 = table_nod(i,iel)
            xel(1,i) = x(1,i1)
            xel(2,i) = x(2,i1)
            xel(3,i) = zz
            ls_index(i) = r_index
C            ls_im_index(i) = im_index
C            ls_abs_index(i) = abs_index
          enddo
          do i=1,3
            i1 = table_nod(i,iel)
            xel_p1(1,i) = x(1,i1)
            xel_p1(2,i) = x(2,i1)
            xel_p1(3,i) = zz
          enddo
          do i=1,nnodes
            sol_el_abs2(i) = 0.0
            sol_el_abs2_eE(i) = 0.0
C            sol_el_abs2_iD(i) = 0.0
            if (q_average .eq. 1) then
              i1 = table_nod(i,iel)
              do j=1,3
                z_tmp1 = sol_avg(j,i1)
                sol_el(j,i) = z_tmp1
                sol_el_abs2(i) = sol_el_abs2(i) + abs(z_tmp1)**2
                sol_el_abs2_eE(i) = sol_el_abs2_eE(i) +
     *               ls_index(i)**2 * abs(z_tmp1)**2
C                sol_el_abs2_iD(i) = sol_el_abs2_iD(i) +
C     *               ls_im_index(i)**2 * abs(z_tmp1)**2
              enddo
            else
              do j=1,3
                z_tmp1 = sol(j,i,plot_val,iel)
                sol_el(j,i) = z_tmp1
                sol_el_abs2(i) = sol_el_abs2(i) + abs(z_tmp1)**2
                sol_el_abs2_eE(i) = sol_el_abs2_eE(i) +
     *               ls_index(i)**2 * abs(z_tmp1)**2
C                sol_el_abs2_iD(i) = sol_el_abs2_iD(i) +
C     *               ls_im_index(i)**2 * abs(z_tmp1)**2
              enddo
            endif
            if (dble(sol_max(4)) .lt. sol_el_abs2(i)) then
              sol_max(1) = sol_el(1,i)
              sol_max(2) = sol_el(2,i)
              sol_max(3) = sol_el(3,i)
              sol_max(4) = sol_el_abs2(i)
            endif
          enddo
          write(26,10) xel, (sol_el_abs2(i),i=1,nnodes)
          write(27,10) xel, (dble(sol_el(1,i)),i=1,nnodes)
          write(28,10) xel, (dble(sol_el(2,i)),i=1,nnodes)
          write(29,10) xel, (dble(sol_el(3,i)),i=1,nnodes)
          write(30,11) xel_p1, ((dble(sol_el(j,i)),j=1,3),i=1,3)
c          write(30,11) xel, ((dble(sol_el(j,i)),j=1,3),i=1,nnodes)
          write(31,10) xel, (ls_index(i),i=1,nnodes)
          write(32,10) xel, (sol_el_abs2_eE(i),i=1,nnodes)
C          write(33,10) xel, (sol_el_abs2_iD(i),i=1,nnodes)
        enddo
        write(26,*) "};"
        write(27,*) "};"
        write(28,*) "};"
        write(29,*) "};"
        write(30,*) "};"
        write(31,*) "};"
        write(32,*) "};"
C        write(33,*) "};"
      close(26)
      close(27)
      close(28)
      close(29)
      close(30)
      close(31)
      close(32)
C      close(33)
      endif
c
c###############################################
c
      if (plot_imag .eq. 1) then
c
      tchar=dir_name(1:namelength)// '/' // gmsh_file_pos(1:namelen)
     *           // '_' // tval // 'x_im.pos'
      open (unit=27,file=tchar)
        write(27,*) "View.ArrowSizeMax = 40;"
        write(27,*) "View.ArrowSizeMin = 0;"
        write(27,*) "General.Color.Background = {255,255,255};"
        write(27,*) "General.Color.BackgroundGradient={255,255,255};"
        write(27,*) "View.AdaptVisualizationGrid =1;"
        write(27,*) "View.MaxRecursionLevel = 2;"
        write(27,*) "View.IntervalsType = 3;"
        write(27,*) "View.Light = 0;"
        write(27,*) "View ""Im ",tE_H,"x: n = ",
     *     plot_val, ", beta_n  =", v_re, "+ I *", v_im, " "" {"
c
      tchar=dir_name(1:namelength)// '/' // gmsh_file_pos(1:namelen)
     *           // '_' // tval // 'y_im.pos'
      open (unit=28,file=tchar)
        write(28,*) "View.ArrowSizeMax = 40;"
        write(28,*) "View.ArrowSizeMin = 0;"
        write(28,*) "General.Color.Background = {255,255,255};"
        write(28,*) "General.Color.BackgroundGradient={255,255,255};"
        write(28,*) "View.AdaptVisualizationGrid =1;"
        write(28,*) "View.MaxRecursionLevel = 2;"
        write(28,*) "View.IntervalsType = 3;"
        write(28,*) "View.Light = 0;"
        write(28,*) "View ""Im ",tE_H,"y: n = ",
     *     plot_val, ", beta_n  =", v_re, "+ I *", v_im, " "" {"
c
      tchar=dir_name(1:namelength)// '/' // gmsh_file_pos(1:namelen)
     *           // '_' // tval // 'z_im.pos'
      open (unit=29,file=tchar)
        write(29,*) "View.ArrowSizeMax = 40;"
        write(29,*) "View.ArrowSizeMin = 0;"
        write(29,*) "General.Color.Background = {255,255,255};"
        write(29,*) "General.Color.BackgroundGradient={255,255,255};"
        write(29,*) "View.AdaptVisualizationGrid =1;"
        write(29,*) "View.MaxRecursionLevel = 2;"
        write(29,*) "View.IntervalsType = 3;"
        write(29,*) "View.Light = 0;"
        write(29,*) "View ""Im ",tE_H,"z: n = ",
     *     plot_val, ", beta_n  =", v_re, "+ I *", v_im, " "" {"
c
      tchar=dir_name(1:namelength)// '/' // gmsh_file_pos(1:namelen)
     *           // '_' // tval // 'v_im.pos'
      open (unit=30,file=tchar)
        write(30,*) "View.ArrowSizeMax = 40;"
        write(30,*) "View.ArrowSizeMin = 0;"
        write(30,*) "General.Color.Background = {255,255,255};"
        write(30,*) "General.Color.BackgroundGradient={255,255,255};"
        write(30,*) "View.AdaptVisualizationGrid =1;"
        write(30,*) "View.MaxRecursionLevel = 2;"
        write(30,*) "View.IntervalsType = 3;"
        write(30,*) "View.Light = 0;"
        write(30,*) "View ""Im ",tE_H,": n = ",
     *     plot_val, ", beta_n  =", v_re, "+ I *", v_im, " "" {"
c
        sol_max(4) = 0.0d0
        do iel=1,nel
          typ_e = type_el(iel)
          zz = 0.0d0
          do i=1,nnodes
            i1 = table_nod(i,iel)
            xel(1,i) = x(1,i1)
            xel(2,i) = x(2,i1)
            xel(3,i) = zz
          enddo
          do i=1,3
            i1 = table_nod(i,iel)
            xel_p1(1,i) = x(1,i1)
            xel_p1(2,i) = x(2,i1)
            xel_p1(3,i) = zz
          enddo
          do i=1,nnodes
            if (q_average .eq. 1) then
              i1 = table_nod(i,iel)
              do j=1,3
                z_tmp1 = sol_avg(j,i1)
                sol_el(j,i) = z_tmp1
              enddo
            else
              do j=1,3
                z_tmp1 = sol(j,i,plot_val,iel)
                sol_el(j,i) = z_tmp1
              enddo
            endif
          enddo
          write(27,10) xel, (imag(sol_el(1,i)),i=1,nnodes)
          write(28,10) xel, (imag(sol_el(2,i)),i=1,nnodes)
          write(29,10) xel, (imag(sol_el(3,i)),i=1,nnodes)
          write(30,11) xel_p1, ((imag(sol_el(j,i)),j=1,3),i=1,3)
        enddo
        write(27,*) "};"
        write(28,*) "};"
        write(29,*) "};"
        write(30,*) "};"
      close(27)
      close(28)
      close(29)
      close(30)
      endif
c
c###############################################
c
c
c###############################################
c
      if (plot_abs .eq. 1) then
c
      tchar=dir_name(1:namelength)// '/' // gmsh_file_pos(1:namelen)
     *           // '_' // tval // 'x_abs.pos'
      open (unit=27,file=tchar)
        write(27,*) "View.ArrowSizeMax = 40;"
        write(27,*) "View.ArrowSizeMin = 0;"
        write(27,*) "General.Color.Background = {255,255,255};"
        write(27,*) "General.Color.BackgroundGradient={255,255,255};"
        write(27,*) "View.AdaptVisualizationGrid =1;"
        write(27,*) "View.MaxRecursionLevel = 2;"
        write(27,*) "View.IntervalsType = 3;"
        write(27,*) "View.Light = 0;"
        write(27,*) "View ""|",tE_H,"x|: n = ",
     *     plot_val, ", beta_n  =", v_re, "+ I *", v_im, " "" {"
c
      tchar=dir_name(1:namelength)// '/' // gmsh_file_pos(1:namelen)
     *           // '_' // tval // 'y_abs.pos'
      open (unit=28,file=tchar)
        write(28,*) "View.ArrowSizeMax = 40;"
        write(28,*) "View.ArrowSizeMin = 0;"
        write(28,*) "General.Color.Background = {255,255,255};"
        write(28,*) "General.Color.BackgroundGradient={255,255,255};"
        write(28,*) "View.AdaptVisualizationGrid =1;"
        write(28,*) "View.MaxRecursionLevel = 2;"
        write(28,*) "View.IntervalsType = 3;"
        write(28,*) "View.Light = 0;"
        write(28,*) "View ""|",tE_H,"y|: n = ",
     *     plot_val, ", beta_n  =", v_re, "+ I *", v_im, " "" {"
c
      tchar=dir_name(1:namelength)// '/' // gmsh_file_pos(1:namelen)
     *           // '_' // tval // 'z_abs.pos'
      open (unit=29,file=tchar)
        write(29,*) "View.ArrowSizeMax = 40;"
        write(29,*) "View.ArrowSizeMin = 0;"
        write(29,*) "General.Color.Background = {255,255,255};"
        write(29,*) "General.Color.BackgroundGradient={255,255,255};"
        write(29,*) "View.AdaptVisualizationGrid =1;"
        write(29,*) "View.MaxRecursionLevel = 2;"
        write(29,*) "View.IntervalsType = 3;"
        write(29,*) "View.Light = 0;"
        write(29,*) "View ""|",tE_H,"z|: n = ",
     *     plot_val, ", beta_n  =", v_re, "+ I *", v_im, " "" {"
c
        sol_max(4) = 0.0d0
        do iel=1,nel
          typ_e = type_el(iel)
          zz = 0.0d0
          do i=1,nnodes
            i1 = table_nod(i,iel)
            xel(1,i) = x(1,i1)
            xel(2,i) = x(2,i1)
            xel(3,i) = zz
          enddo
          do i=1,nnodes
            if (q_average .eq. 1) then
              i1 = table_nod(i,iel)
              do j=1,3
                z_tmp1 = sol_avg(j,i1)
                sol_el(j,i) = z_tmp1
              enddo
            else
              do j=1,3
                z_tmp1 = sol(j,i,plot_val,iel)
                sol_el(j,i) = z_tmp1
              enddo
            endif
          enddo
          write(27,10) xel, (abs(sol_el(1,i)),i=1,nnodes)
          write(28,10) xel, (abs(sol_el(2,i)),i=1,nnodes)
          write(29,10) xel, (abs(sol_el(3,i)),i=1,nnodes)
        enddo
        write(27,*) "};"
        write(28,*) "};"
        write(29,*) "};"
      close(27)
      close(28)
      close(29)
      endif
c
c###############################################
c
c     ST : Scalar triangle
10    format("ST2(",f10.6,17(",",f10.6),"){",
     *     g24.16,5(",",g24.16),"};")

c     VT : Vector triangle
11    format("VT(",f10.6,8(",",f10.6),"){",
     *     g24.16,8(",",g24.16),"};")

c 11    format("VT2(",f10.6,17(",",f10.6),"){",
c     *     g24.16,17(",",g24.16),"};")

      if (debug .eq. 1) then
        write(ui,*)
        write(ui,*) "gmsh_post_process: plot_val = ", plot_val
        write(ui,*) "gmsh_post_process: sol_max = ", sol_max
      endif
C
      deallocate(sol_avg)
C
      return
      end
