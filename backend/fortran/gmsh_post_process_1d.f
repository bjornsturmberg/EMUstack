c  Write an ASCII file that can be read into GMSH.
c  P1-element is used to represent the  3D vector field.
c  P2-element is used to represent each component of the 3D vector field.
c

      subroutine gmsh_post_process_1d (plot_val, E_H_field, nval, 
     *     nel, npt_P2, table_nod, type_el, nb_typ_el, 
     *     n_eff, x_P2, beta, sol_P2,
     *     gmsh_file_pos, dir_name, 
     *     q_average, plot_real, plot_imag, plot_abs)


      implicit none
      integer*8 nval, nel, npt_P2, plot_val, E_H_field
      integer*8 nb_typ_el
      integer*8 table_nod(3,nel), type_el(nel)
      double precision x_P2(npt_P2)
      complex*16 sol_P2(3,3,nval,nel), n_eff(nb_typ_el)
      integer alloc_stat
      complex*16, dimension(:,:), allocatable :: sol_avg

      complex*16, dimension(:,:), allocatable :: tmp_sol_1

      complex*16 beta(nval)


      integer*8, allocatable ::  visite(:)  ! visite(npt_P2)

      integer*8 nnodes_0
      parameter (nnodes_0 = 3)

      double precision xel(3,nnodes_0), xel_p1(3,2)
      complex*16 sol_el(3,nnodes_0)
      double precision sol_el_abs2(nnodes_0), sol_max(4)
      double precision sol_el_abs2_eE(nnodes_0)
C      double precision sol_el_abs2_iD(nnodes_0)
      double precision ls_index(nnodes_0), r_index, yy, zz
C      double precision ls_im_index(nnodes_0), im_index
C      double precision ls_abs_index(nnodes_0), abs_index

      real v_im, v_re

      integer*8 i, j, i1, iel, namelen, namelen2, typ_e
      integer*8 q_average, plot_imag, plot_real, plot_abs
      integer*8 debug, ui
      complex*16 z_tmp1
      complex*16 ii
      character*(*) gmsh_file_pos, dir_name
      character*100 tchar
      character tval*4, buf*3
      character*1 tE_H
      integer*8 namelength, charlength
c
c  ii = sqrt(-1)
      ii = cmplx(0.0d0, 1.0d0)
c
      ui = 6
      debug = 0
c
c      if ( nnodes .ne. 6 ) then
c        write(ui,*) "gmsh_post_process_1d: problem nnodes = ", nnodes
c        write(ui,*) "gmsh_post_process_1d: nnodes should be equal to 6 !"
c        write(ui,*) "gmsh_post_process_1d: Aborting..."
c        stop
c      endif
C
c
c###############################################
c
      alloc_stat = 0
      allocate(sol_avg(3,npt_P2), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "gmsh_post_process_1d: ",
     *             "Mem. allocation is unseccesfull"
        write(*,*) "alloc_stat (sol_avg) = ", alloc_stat
        write(*,*) "npt_P2 = ", npt_P2
        write(*,*) "Aborting..."
        stop
      endif
      allocate(visite(npt_P2), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "gmsh_post_process_1d: ",
     *             "Mem. allocation is unseccesfull"
        write(*,*) "alloc_stat (visite) = ", alloc_stat
        write(*,*) "npt_P2 = ", npt_P2
        write(*,*) "Aborting..."
        stop
      endif

      allocate(tmp_sol_1(3,npt_P2), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "gmsh_post_process_1d: ",
     *             "Mem. allocation is unseccesfull"
        write(*,*) "alloc_stat (tmp_sol_1) = ", alloc_stat
        write(*,*) "npt_P2 = ", npt_P2
        write(*,*) "Aborting..."
        stop
      endif

c
c################# END: TEMPORARY CHANGE ##############################
c
cc      if (plot_val == 3) then
cc        do i=1,npt_P2
cc          do j=1,3
cc            tmp_sol_1(j,i) = -911
cc          enddo
cc        enddo
cc        do iel=1,nel
cc          do i=1,nnodes_0
cc            i1 = table_nod(i,iel)
cc            do j=1,3
cc              tmp_sol_1(j,i1) = sol_P2(j,i,plot_val,iel)
cc            enddo
cc          enddo
cc        enddo
cc       open (unit=335, file="Matrices/mode_1.txt", status='unknown')
cc        do i=1,npt_P2
cc          write(335,"(I8,20(f20.12))") i, x_P2(i), 
cc     *              (tmp_sol_1(j,i),j=1,3)
cc        enddo
cc        close(335)
cc      endif
c
c################# END: TEMPORARY CHANGE ##############################


c
c###############################################
c
C

      if (plot_val .eq. 0) return
c
      if (E_H_field .eq. 1) then
        tE_H = "E"
      elseif(E_H_field .eq. 2) then
        tE_H = "H"
      else
        write(ui,*) "gmsh_post_process_1d: ",
     *             "E_H_field has invalid value: ", 
     *    E_H_field
        write(ui,*) "Aborting..."
        stop
      endif
c
      if (q_average .eq. 1) then
        do i=1,npt_P2
          visite(i) = 0
          do j=1,3
            sol_avg(j,i) = 0.0d0
          enddo
        enddo
        do iel=1,nel
          do i=1,nnodes_0
            i1 = table_nod(i,iel)
            visite(i1) = visite(i1) + 1
            do j=1,3
              sol_avg(j,i1) = sol_avg(j,i1) + 
     *          sol_P2(j,i,plot_val,iel) 
            enddo
          enddo
        enddo
        do i=1,npt_P2
          if (visite(i) .eq. 0) then
            write(ui,*) "gmsh_post_process_1d: visite(i) = 0"
            write(ui,*) " i, visite(i) = ", i, visite(i)
            write(ui,*) "gmsh_post_process_1d: Aborting..."
            stop
          endif
          do j=1,3
            sol_avg(j,i) = sol_avg(j,i)/dble(visite(i))
          enddo
        enddo
      endif
c
      v_re = dble(beta(plot_val))
      v_im = -dble(ii*beta(plot_val))
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
      
      write(34,*) "Delete View[0];"
      
      tchar = '../../'//dir_name(1:namelength)// '/' 
     *  //gmsh_file_pos(1:namelen) // '_' // tval // '_abs2_eE.pos'
      namelen2 = len_trim(tchar)
      write(34,*) " Include """, tchar(1:namelen2), """;"
      tchar = gmsh_file_pos(1:namelen) // '_' // tval // '_abs2_eE.png'
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
          r_index = n_eff(typ_e)
C          im_index = imag(sqrt(eps_eff(typ_e)))
C          abs_index = abs(sqrt(eps_eff(typ_e)))
          yy = 0.0d0
          zz = 0.0d0
          do i=1,nnodes_0
            i1 = table_nod(i,iel)
            xel(1,i) = x_P2(i1)
            xel(2,i) = yy
            xel(3,i) = zz
            ls_index(i) = r_index
          enddo
          do i=1,2
            i1 = table_nod(i,iel)
            xel_p1(1,i) = x_P2(i1)
            xel_p1(2,i) = yy
            xel_p1(3,i) = zz
          enddo
cccccccccccccccccc
c          i = 1  ! Left side end-point
c            i1 = table_nod(i,iel)
c            xel_p1(1,i) = x_P2(i1)
c            xel_p1(2,i) = yy
c            xel_p1(3,i) = zz
c          i = 2  ! Right side end-point
c            i1 = table_nod(i+1,iel)
c            xel_p1(1,i) = x_P2(i1)
c            xel_p1(2,i) = yy
c            xel_p1(3,i) = zz
cccccccccccccccccc          enddo
          do i=1,nnodes_0
            sol_el_abs2(i) = 0.0
            sol_el_abs2_eE(i) = 0.0
            if (q_average .eq. 1) then
              i1 = table_nod(i,iel)
              do j=1,3
                z_tmp1 = sol_avg(j,i1)
                sol_el(j,i) = z_tmp1
                sol_el_abs2(i) = sol_el_abs2(i) + abs(z_tmp1)**2
                sol_el_abs2_eE(i) = sol_el_abs2_eE(i) + 
     *               ls_index(i)**2 * abs(z_tmp1)**2
              enddo
            else
              do j=1,3
                z_tmp1 = sol_P2(j,i,plot_val,iel) 
                sol_el(j,i) = z_tmp1
                sol_el_abs2(i) = sol_el_abs2(i) + abs(z_tmp1)**2
                sol_el_abs2_eE(i) = sol_el_abs2_eE(i) + 
     *               ls_index(i)**2 * abs(z_tmp1)**2
              enddo
            endif
            if (sol_max(4) .lt. sol_el_abs2(i)) then
              sol_max(1) = sol_el(1,i)
              sol_max(2) = sol_el(2,i)
              sol_max(3) = sol_el(3,i)
              sol_max(4) = sol_el_abs2(i)
            endif
          enddo
          write(26,10) xel, (sol_el_abs2(i),i=1,nnodes_0)
          write(27,10) xel, (dble(sol_el(1,i)),i=1,nnodes_0)
          write(28,10) xel, (dble(sol_el(2,i)),i=1,nnodes_0)
          write(29,10) xel, (dble(sol_el(3,i)),i=1,nnodes_0)
          write(30,11) xel_p1, ((dble(sol_el(j,i)),j=1,3),i=1,2)
c          write(30,11) xel, ((dble(sol_el(j,i)),j=1,3),i=1,nnodes_0)
          write(31,10) xel, (ls_index(i),i=1,nnodes_0)
          write(32,10) xel, (sol_el_abs2_eE(i),i=1,nnodes_0)
C          write(33,10) xel, (sol_el_abs2_iD(i),i=1,nnodes_0)
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
          yy = 0.0d0
          zz = 0.0d0
          do i=1,nnodes_0
            i1 = table_nod(i,iel)
            xel(1,i) = x_P2(i1)
            xel(2,i) = yy
            xel(3,i) = zz
          enddo
          do i=1,2
            i1 = table_nod(i,iel)
            xel_p1(1,i) = x_P2(i1)
            xel_p1(2,i) = yy
            xel_p1(3,i) = zz
          enddo
ccccccccccccccccccccccccc
cc          i = 1
cc            i1 = table_nod(i,iel)
cc            xel_p1(1,i) = x_P2(i1)
cc            xel_p1(2,i) = yy
cc            xel_p1(3,i) = zz
cc          i = 2
cc            i1 = table_nod(i+1,iel)
cc            xel_p1(1,i) = x_P2(i1)
cc            xel_p1(2,i) = yy
cc            xel_p1(3,i) = zz
ccccccccccccccccccccccccc          enddo
          do i=1,nnodes_0
            if (q_average .eq. 1) then
              i1 = table_nod(i,iel)
              do j=1,3
                z_tmp1 = sol_avg(j,i1)
                sol_el(j,i) = z_tmp1
              enddo
            else
              do j=1,3
                z_tmp1 = sol_P2(j,i,plot_val,iel) 
                sol_el(j,i) = z_tmp1
              enddo
            endif
          enddo
          write(27,10) xel, (imag(sol_el(1,i)),i=1,nnodes_0)
          write(28,10) xel, (imag(sol_el(2,i)),i=1,nnodes_0)
          write(29,10) xel, (imag(sol_el(3,i)),i=1,nnodes_0)
          write(30,11) xel_p1, ((imag(sol_el(j,i)),j=1,3),i=1,2)
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
          yy = 0.0d0
          zz = 0.0d0
          do i=1,nnodes_0
            i1 = table_nod(i,iel)
            xel(1,i) = x_P2(i1)
            xel(2,i) = yy
            xel(3,i) = zz
          enddo
          do i=1,nnodes_0
            if (q_average .eq. 1) then
              i1 = table_nod(i,iel)
              do j=1,3
                z_tmp1 = sol_avg(j,i1)
                sol_el(j,i) = z_tmp1
              enddo
            else
              do j=1,3
                z_tmp1 = sol_P2(j,i,plot_val,iel) 
                sol_el(j,i) = z_tmp1
              enddo
            endif
          enddo
          write(27,10) xel, (abs(sol_el(1,i)),i=1,nnodes_0)
          write(28,10) xel, (abs(sol_el(2,i)),i=1,nnodes_0)
          write(29,10) xel, (abs(sol_el(3,i)),i=1,nnodes_0)
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
c     SL : Scalar line
10    format("SL2(",f10.6,8(",",f10.6),"){",
     *     g24.16,2(",",g24.16),"};")

c     VT : Vector line
11    format("VL(",f10.6,5(",",f10.6),"){",
     *     g24.16,5(",",g24.16),"};")

c 11    format("VT2(",f10.6,17(",",f10.6),"){",
c     *     g24.16,17(",",g24.16),"};")

      if (debug .eq. 1) then
        write(ui,*)
        write(ui,*) "gmsh_post_process_1d: plot_val = ", plot_val
        write(ui,*) "gmsh_post_process_1d: sol_max = ", sol_max
      endif
C
      deallocate(sol_avg)
C
      return
      end
