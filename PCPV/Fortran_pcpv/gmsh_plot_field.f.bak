c  Write an ASCII file that can be read into GMSH.
c  P1-element is used to represent the  3D vector field.
c  P2-element is used to represent each component of the 3D vector field.
c
      subroutine gmsh_plot_field (plot_val, E_H_field, nval, 
     *     nel, npt, nnodes, table_nod, type_el, eps_eff,
     *     x, Beta, sol, sol_tmp, vec_coef, h, hz, 
     *     gmsh_file_pos, dir_name, nb_typ_el, 
     *       q_average, plot_real, plot_imag, plot_abs)
c
      implicit none
      integer*8 nval, nel, npt, nnodes, plot_val, E_H_field
      integer*8 nb_typ_el
      double precision h, hz
      integer*8 table_nod(nnodes,nel), type_el(nel)
      complex*16 x(2,npt), eps_eff(nb_typ_el)
      complex*16 sol(3,nnodes+7,nval,nel)
      complex*16 sol_tmp(3,nnodes,nel)
      complex*16 vec_coef(2*nval)
      complex*16 Beta(nval)
      character*(*) gmsh_file_pos, dir_name
c
c     Local variables

      integer*8 nnodes_0
      parameter (nnodes_0 = 6)

      double precision xel(3,nnodes_0), xel_p1(3,3)
      complex*16 sol_el(3,nnodes_0)
      double precision sol_el_abs2(nnodes_0), sol_max(4)
      double precision ls_index(nnodes_0), r_index, zz

      complex*16 P_down, P_up, coef_down, coef_up, coef_t, coef_z

c      real v_im, v_re

      integer*8 i, j, i1, iel, ival, namelen, typ_e
      integer*8 q_average, plot_imag, plot_real, plot_abs
      integer*8 debug, ui
      complex*16 z_tmp1
      complex*16 ii
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
      if ( nnodes .ne. 6 ) then
        write(ui,*) "gmsh_plot_field: problem nnodes = ", nnodes
        write(ui,*) "gmsh_plot_field: nnodes should be equal to 6 !"
        write(ui,*) "gmsh_plot_field: Aborting..."
        stop
      endif
c
      if (hz .lt. 0 .or. hz .gt. h ) then
        write(ui,*)
        write(ui,*) "gmsh_plot_field: invalid value for hz"
        write(ui,*) "hz should be in the interval [0,h]"
        write(ui,*) "hz, h = ", hz, h
        write(ui,*) "gmsh_plot_field: Aborting..."
        stop
      endif
c
c      q_average = 0  ! use average value if q_average = 1
c      plot_real = 1  ! plot real part if plot_real = 1
c      plot_imag = 1  ! plot real part if plot_imag = 1
c      plot_abs = 1  ! plot absolute value if plot_abs = 1
c
      if (plot_val .eq. 0) return
c
      if (E_H_field .eq. 1) then
        tE_H = "E"
      elseif(E_H_field .eq. 2) then
        tE_H = "H"
      else
        write(ui,*) "gmsh_plot_field: E_H_field has invalid value: ", 
     *    E_H_field
        write(ui,*) "Aborting..."
        stop
      endif
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

      tchar=dir_name(1:namelength)// '/' // gmsh_file_pos(1:namelen) 
     *           // '_field_' // tval // '_abs2_eD.pos'
      open (unit=26,file=tchar)
        write(26,*) "View.AdaptVisualizationGrid =1;"
        write(26,*) "View.IntervalsType = 3;"
        write(26,*) "View.Light = 0;"
        write(26,*) "View ""|",tE_H,"_t|^2: n = ", plot_val,
     *     "; nval = ", nval, "; hz = ", (-hz), " "" {"
c
      tchar=dir_name(1:namelength)// '/' // gmsh_file_pos(1:namelen) 
     *            // '_field_' // tval // 'x_re.pos'
      open (unit=27,file=tchar)
        write(27,*) "View.AdaptVisualizationGrid =1;"
        write(27,*) "View.IntervalsType = 3;"
        write(27,*) "View.Light = 0;"
c        write(27,*) "View ""Re Ex: n = ", plot_val, 
        write(27,*) "View ""Re ",tE_H,"x: n = ", plot_val, 
     *     "; nval = ", nval, "; hz = ", (-hz), " "" {"
c
      tchar=dir_name(1:namelength)// '/' // gmsh_file_pos(1:namelen)  
     *           // '_field_' // tval // 'y_re.pos'
      open (unit=28,file=tchar)
        write(28,*) "View.AdaptVisualizationGrid =1;"
        write(28,*) "View.IntervalsType = 3;"
        write(28,*) "View.Light = 0;"
        write(28,*) "View ""Re ",tE_H,"y: n = ", plot_val, 
     *     "; nval = ", nval, "; hz = ", (-hz), " "" {"
c
      tchar=dir_name(1:namelength)// '/' // gmsh_file_pos(1:namelen)  
     *           // '_field_' // tval // 'z_re.pos'
      open (unit=29,file=tchar)
        write(29,*) "View.AdaptVisualizationGrid =1;"
        write(29,*) "View.IntervalsType = 3;"
        write(29,*) "View.Light = 0;"
        write(29,*) "View ""Re ",tE_H,"z: n = ", plot_val, 
     *     "; nval = ", nval, "; hz = ", (-hz), " "" {"
c
      tchar=dir_name(1:namelength)// '/' // gmsh_file_pos(1:namelen)  
     *           // '_field_' // tval // 'v_re.pos'
      open (unit=30,file=tchar)
        write(30,*) "View.AdaptVisualizationGrid =1;"
        write(30,*) "View.IntervalsType = 3;"
        write(30,*) "View.Light = 0;"
        write(30,*) "View ""Re ",tE_H,": n = ", plot_val, 
     *     "; nval = ", nval, "; hz = ", (-hz), " "" {"
c
      do iel=1,nel
        do i=1,nnodes
          do j=1,3
            sol_tmp(j,i,iel) = 0.0d0
          enddo
        enddo
      enddo
c     
      do ival=1,nval
        P_down = EXP(ii*Beta(ival)*hz)   !  Introduce Propagation in -z
        P_up = EXP(ii*Beta(ival)*(h-hz)) !  Introduce Propagation in +z
        coef_down = vec_coef(ival) * P_down
        coef_up = vec_coef(ival+nval) * P_up

c        coef_down = vec_coef(ival) * P_up
c        coef_up = vec_coef(ival+nval) * P_down

        coef_t = coef_up + coef_down
        coef_z = (coef_up - coef_down)/Beta(ival) ! Taking into accout the change of variable for Ez
        do iel=1,nel
          do i=1,nnodes
            do j=1,2
              z_tmp1 = sol(j,i,ival,iel) * coef_t
              sol_tmp(j,i,iel) = sol_tmp(j,i,iel) + z_tmp1
            enddo
            j=3
              z_tmp1 = sol(j,i,ival,iel) * coef_z
              sol_tmp(j,i,iel) = sol_tmp(j,i,iel) + z_tmp1
          enddo
        enddo
      enddo
c       
        sol_max(4) = 0.0d0
        do iel=1,nel
          zz = hz
          typ_e = type_el(iel)
          r_index = (eps_eff(typ_e))
          do i=1,nnodes
            i1 = table_nod(i,iel)
            xel(1,i) = x(1,i1)
            xel(2,i) = x(2,i1)
            xel(3,i) = zz
            ls_index(i) = r_index
          enddo
          do i=1,3
            i1 = table_nod(i,iel)
            xel_p1(1,i) = x(1,i1)
            xel_p1(2,i) = x(2,i1)
            xel_p1(3,i) = zz
          enddo
          do i=1,nnodes
            sol_el_abs2(i) = 0.0
            do j=1,2
              z_tmp1 = sol_tmp(j,i,iel)
              sol_el(j,i) = z_tmp1
              sol_el_abs2(i) = sol_el_abs2(i) + 
     *           ls_index(i) * abs(z_tmp1)**2
            enddo
            j=3
              z_tmp1 = sol_tmp(j,i,iel)
              sol_el(j,i) = z_tmp1
            if (sol_max(4) .lt. sol_el_abs2(i)) then
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
      close(31)
      endif
c
c###############################################
c
      if (plot_imag .eq. 1) then
c
      tchar=dir_name(1:namelength)// '/' // gmsh_file_pos(1:namelen)
     *           // '_field_' // tval // 'x_im.pos'
      open (unit=27,file=tchar)
        write(27,*) "View.AdaptVisualizationGrid =1;"
        write(27,*) "View.IntervalsType = 3;"
        write(27,*) "View.Light = 0;"
        write(27,*) "View ""Im ",tE_H,"x: n = ", plot_val, 
     *     "; nval = ", nval, "; hz = ", (-hz), " "" {"
c
      tchar=dir_name(1:namelength)// '/' // gmsh_file_pos(1:namelen)  
     *           // '_field_' // tval // 'y_im.pos'
      open (unit=28,file=tchar)
        write(28,*) "View.AdaptVisualizationGrid =1;"
        write(28,*) "View.IntervalsType = 3;"
        write(28,*) "View.Light = 0;"
        write(28,*) "View ""Im ",tE_H,"y: n = ", plot_val, 
     *     "; nval = ", nval, "; hz = ", (-hz), " "" {"
c
      tchar=dir_name(1:namelength)// '/' // gmsh_file_pos(1:namelen)  
     *           // '_field_' // tval // 'z_im.pos'
      open (unit=29,file=tchar)
        write(29,*) "View.AdaptVisualizationGrid =1;"
        write(29,*) "View.IntervalsType = 3;"
        write(29,*) "View.Light = 0;"
        write(29,*) "View ""Im ",tE_H,"z: n = ", plot_val, 
     *     "; nval = ", nval, "; hz = ", (-hz), " "" {"
c
      tchar=dir_name(1:namelength)// '/' // gmsh_file_pos(1:namelen)  
     *           // '_field_' // tval // 'v_im.pos'
      open (unit=30,file=tchar)
        write(30,*) "View.AdaptVisualizationGrid =1;"
        write(30,*) "View.IntervalsType = 3;"
        write(30,*) "View.Light = 0;"
        write(30,*) "View ""Im ",tE_H,": n = ", plot_val, 
     *     "; nval = ", nval, "; hz = ", (-hz), " "" {"
c
        sol_max(4) = 0.0d0
        do iel=1,nel
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
            do j=1,3
              z_tmp1 = sol_tmp(j,i,iel)
              sol_el(j,i) = z_tmp1
            enddo
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
     *           // '_field_' // tval // 'x_abs.pos'
      open (unit=27,file=tchar)
        write(27,*) "View.AdaptVisualizationGrid =1;"
        write(27,*) "View.IntervalsType = 3;"
        write(27,*) "View.Light = 0;"
        write(27,*) "View ""|",tE_H,"x|: n = ", plot_val, 
     *     "; nval = ", nval, "; hz = ", (-hz), " "" {"
c
      tchar=dir_name(1:namelength)// '/' // gmsh_file_pos(1:namelen)  
     *           // '_field_' // tval // 'y_abs.pos'
      open (unit=28,file=tchar)
        write(28,*) "View.AdaptVisualizationGrid =1;"
        write(28,*) "View.IntervalsType = 3;"
        write(28,*) "View.Light = 0;"
        write(28,*) "View ""|",tE_H,"y|: n = ", plot_val, 
     *     "; nval = ", nval, "; hz = ", (-hz), " "" {"
c
      tchar=dir_name(1:namelength)// '/' // gmsh_file_pos(1:namelen)  
     *           // '_field_' // tval // 'z_abs.pos'
      open (unit=29,file=tchar)
        write(29,*) "View.AdaptVisualizationGrid =1;"
        write(29,*) "View.IntervalsType = 3;"
        write(29,*) "View.Light = 0;"
        write(29,*) "View ""|",tE_H,"z|: n = ", plot_val, 
     *     "; nval = ", nval, "; hz = ", (-hz), " "" {"
c
        sol_max(4) = 0.0d0
        do iel=1,nel
          zz = 0.0d0
          do i=1,nnodes
            i1 = table_nod(i,iel)
            xel(1,i) = x(1,i1)
            xel(2,i) = x(2,i1)
            xel(3,i) = zz
          enddo
          do i=1,nnodes
            do j=1,3
              z_tmp1 = sol_tmp(j,i,iel)
              sol_el(j,i) = z_tmp1
            enddo
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
        write(ui,*) "gmsh_plot_field: plot_val = ", plot_val
        write(ui,*) "gmsh_plot_field: sol_max = ", sol_max
      endif
c
      return
      end
