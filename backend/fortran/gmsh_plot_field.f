c  Write an ASCII file that can be read into GMSH.
c  P1-element is used to represent the  3D vector field.
c  P2-element is used to represent each component of the 3D vector field.
c
      subroutine gmsh_plot_field (nval, 
     *     nel, npt, nnodes, nb_typ_el, table_nod, type_el, 
     *     eps_eff, x, beta, sol, vec_coef, h, hz, 
     *     gmsh_file_pos, plot_real, plot_imag, 
     *     plot_abs, extra_name)
c
      implicit none
      integer*8 nval, nel, npt, nnodes
      integer*8 nb_typ_el
      double precision h, hz
      integer*8 table_nod(nnodes,nel), type_el(nel)
      double precision x(2,npt)
      complex*16 eps_eff(nb_typ_el)
      complex*16 sol(3,nnodes+7,nval,nel)
      complex*16 vec_coef(2*nval)
      complex*16 beta(nval)
      character(*) gmsh_file_pos
      character(*) extra_name
c
c     Local variables

      integer*8 nnodes_0
      parameter (nnodes_0 = 6)

      character*20 dir_name
      integer alloc_stat
      complex*16, dimension(:,:,:), allocatable :: sol_tmp
      double precision xel(3,nnodes_0), xel_p1(3,3)
      complex*16 sol_el(3,nnodes_0), sol_max(4)
      double precision sol_el_abs2(nnodes_0)
      double precision ls_index(nnodes_0), r_index, zz

      complex*16 P_down, P_up, coef_down, coef_up, coef_t, coef_z

c      real v_im, v_re

      integer*8 i, j, i1, iel, ival, namelen, typ_e
      integer*8 plot_imag, plot_real, plot_abs
      integer*8 debug, ui
      complex*16 z_tmp1
      complex*16 ii
      character*100 tchar
      character*1 tE_H
      integer*8 namelength


Cf2py intent(in) nval, nel, npt, nnodes nb_typ_el, table_nod,
Cf2py intent(in) type_el, eps_eff, x, beta, sol, vec_coef, h, hz
Cf2py intent(in) gmsh_file_pos, extra_name
Cf2py intent(in) plot_real, plot_imag, plot_abs

Cf2py depend(table_nod) nnodes, nel
Cf2py depend(type_el) nel
Cf2py depend(x) npt
Cf2py depend(eps_eff) nb_typ_el
Cf2py depend(sol) nnodes, nval, nel
Cf2py depend(vec_coef) nval
Cf2py depend(beta) nval


      dir_name = "in_plane_fields"
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
C
      alloc_stat = 0
      allocate(sol_tmp(3,nnodes,nel), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "gmsh_plot_field: Mem. allocation is unseccesfull"
        write(*,*) "alloc_stat (sol_tmp) = ", alloc_stat
        write(*,*) "Aborting..."
        stop
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
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
c      plot_real = 1  ! plot real part if plot_real = 1
c      plot_imag = 1  ! plot real part if plot_imag = 1
c      plot_abs = 1  ! plot absolute value if plot_abs = 1
c
      namelen = len_trim(gmsh_file_pos)
      namelength = len_trim(dir_name)
c
c###############################################
c
      if (plot_real .eq. 1) then

      tchar=dir_name(1:namelength)// '/' //extra_name//'_BM_abs2_eD_'
     * // gmsh_file_pos(1:namelen) // '.pos'
      open (unit=26,file=tchar)
        write(26,*) "View.AdaptVisualizationGrid =1;"
        write(26,*) "View.IntervalsType = 3;"
        write(26,*) "View.Light = 0;"
        write(26,*) "View ""|",tE_H,"_t|^2:",
     *     "; nval = ", nval, "; hz = ", (-hz), " "" {"
c
      tchar=dir_name(1:namelength)// '/' // extra_name // '_BMx_re_'
     * // gmsh_file_pos(1:namelen) // '.pos'
      open (unit=27,file=tchar)
        write(27,*) "View.AdaptVisualizationGrid =1;"
        write(27,*) "View.IntervalsType = 3;"
        write(27,*) "View.Light = 0;"
c        write(27,*) "View ""Re Ex: n = ", plot_val, 
        write(27,*) "View ""Re ",tE_H,"x:", 
     *     "; nval = ", nval, "; hz = ", (-hz), " "" {"
c
      tchar=dir_name(1:namelength)// '/' // extra_name // '_BMy_re_'
     * // gmsh_file_pos(1:namelen) // '.pos'
      open (unit=28,file=tchar)
        write(28,*) "View.AdaptVisualizationGrid =1;"
        write(28,*) "View.IntervalsType = 3;"
        write(28,*) "View.Light = 0;"
        write(28,*) "View ""Re ",tE_H,"y:", 
     *     "; nval = ", nval, "; hz = ", (-hz), " "" {"
c
      tchar=dir_name(1:namelength)// '/' // extra_name // '_BMz_re_'
     * // gmsh_file_pos(1:namelen) // '.pos'
      open (unit=29,file=tchar)
        write(29,*) "View.AdaptVisualizationGrid =1;"
        write(29,*) "View.IntervalsType = 3;"
        write(29,*) "View.Light = 0;"
        write(29,*) "View ""Re ",tE_H,"z:", 
     *     "; nval = ", nval, "; hz = ", (-hz), " "" {"
c
      tchar=dir_name(1:namelength)// '/' // extra_name // '_BMv_re_'
     * // gmsh_file_pos(1:namelen) // '.pos'
      open (unit=30,file=tchar)
        write(30,*) "View.AdaptVisualizationGrid =1;"
        write(30,*) "View.IntervalsType = 3;"
        write(30,*) "View.Light = 0;"
        write(30,*) "View ""Re ",tE_H,":", 
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
        P_down = EXP(ii*beta(ival)*hz)   !  Introduce Propagation in -z
        P_up = EXP(ii*beta(ival)*(h-hz)) !  Introduce Propagation in +z
        coef_down = vec_coef(ival) * P_down
        coef_up = vec_coef(ival+nval) * P_up

c        coef_down = vec_coef(ival) * P_up
c        coef_up = vec_coef(ival+nval) * P_down

        coef_t = coef_up + coef_down
        coef_z = (coef_up - coef_down)/beta(ival) ! Taking into accout the change of variable for Ez
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
          r_index = real(eps_eff(typ_e))
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
      tchar=dir_name(1:namelength)// '/' // extra_name // '_BMx_im_'
     * // gmsh_file_pos(1:namelen) // '.pos'
      open (unit=27,file=tchar)
        write(27,*) "View.AdaptVisualizationGrid =1;"
        write(27,*) "View.IntervalsType = 3;"
        write(27,*) "View.Light = 0;"
        write(27,*) "View ""Im ",tE_H,"x:", 
     *     "; nval = ", nval, "; hz = ", (-hz), " "" {"
c
      tchar=dir_name(1:namelength)// '/' // extra_name // '_BMy_im_'
     * // gmsh_file_pos(1:namelen) // '.pos'
      open (unit=28,file=tchar)
        write(28,*) "View.AdaptVisualizationGrid =1;"
        write(28,*) "View.IntervalsType = 3;"
        write(28,*) "View.Light = 0;"
        write(28,*) "View ""Im ",tE_H,"y:", 
     *     "; nval = ", nval, "; hz = ", (-hz), " "" {"
c
      tchar=dir_name(1:namelength)// '/' // extra_name // '_BMz_im_'
     * // gmsh_file_pos(1:namelen) // '.pos'
      open (unit=29,file=tchar)
        write(29,*) "View.AdaptVisualizationGrid =1;"
        write(29,*) "View.IntervalsType = 3;"
        write(29,*) "View.Light = 0;"
        write(29,*) "View ""Im ",tE_H,"z:", 
     *     "; nval = ", nval, "; hz = ", (-hz), " "" {"
c
      tchar=dir_name(1:namelength)// '/' // extra_name // '_BMv_im_'
     * // gmsh_file_pos(1:namelen) // '.pos'
      open (unit=30,file=tchar)
        write(30,*) "View.AdaptVisualizationGrid =1;"
        write(30,*) "View.IntervalsType = 3;"
        write(30,*) "View.Light = 0;"
        write(30,*) "View ""Im ",tE_H,":", 
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
      tchar=dir_name(1:namelength)// '/' // extra_name // '_BMx_abs_'
     * // gmsh_file_pos(1:namelen) // '.pos'
      open (unit=27,file=tchar)
        write(27,*) "View.AdaptVisualizationGrid =1;"
        write(27,*) "View.IntervalsType = 3;"
        write(27,*) "View.Light = 0;"
        write(27,*) "View ""|",tE_H,"x|:", 
     *     "; nval = ", nval, "; hz = ", (-hz), " "" {"
c
      tchar=dir_name(1:namelength)// '/' // extra_name // '_BMy_abs_'
     * // gmsh_file_pos(1:namelen) // '.pos'
      open (unit=28,file=tchar)
        write(28,*) "View.AdaptVisualizationGrid =1;"
        write(28,*) "View.IntervalsType = 3;"
        write(28,*) "View.Light = 0;"
        write(28,*) "View ""|",tE_H,"y|:", 
     *     "; nval = ", nval, "; hz = ", (-hz), " "" {"
c
      tchar=dir_name(1:namelength)// '/' // extra_name // '_BMz_abs_'
     * // gmsh_file_pos(1:namelen) // '.pos'
      open (unit=29,file=tchar)
        write(29,*) "View.AdaptVisualizationGrid =1;"
        write(29,*) "View.IntervalsType = 3;"
        write(29,*) "View.Light = 0;"
        write(29,*) "View ""|",tE_H,"z|:", 
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
        write(ui,*) "gmsh_plot_field: sol_max = ", sol_max
      endif
C
      deallocate(sol_tmp)
C
      return
      end
