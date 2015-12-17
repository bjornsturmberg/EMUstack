c  Write an ASCII file that can be read into GMSH.
c  P1-element is used to represent the  3D vector field.
c  P2-element is used to represent each component of the 3D vector field.
c
      subroutine gmsh_plot_PW (
     *     nel, npt, nnodes, neq_PW, bloch_vec,
     *     table_nod, x, lat_vecs, lambda, n_eff_0,
     *     vec_coef_down, vec_coef_up,
     *     index_pw_inv, ordre_ls, hz, gmsh_file_pos,
     *     plot_real, plot_imag, plot_abs, extra_name)
c
      implicit none
      integer*8 nel, npt, nnodes
      integer*8 neq_PW, ordre_ls
      double precision hz, lat_vecs(2,2), bloch_vec(2)
      integer*8 table_nod(nnodes,nel), index_pw_inv(neq_PW)
      complex*16 n_eff_0 ! dielctric constand of the semi-infinite medium
      double precision x(2,npt)
      complex*16 vec_coef_down(2*neq_PW)
      complex*16 vec_coef_up(2*neq_PW)
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
      double precision zz

      integer*8 px, py, s, s2

      double precision d, dy, r_tmp, xx_g(3)
      double precision k_0, k, vec_kx, vec_ky, lambda
      double precision bloch1, bloch2, pi, alpha, beta, norm
      complex*16 chi_TE, chi_TM, gamma, val_exp

      complex*16 P_down, P_up, coef_RE_down, coef_RE_up
      complex*16 coef_RK_down, coef_RK_up
      complex*16 coef_RE_t, coef_RE_z, coef_RK_t, coef_RK_z

      complex*16 PlaneW_RE(3), PlaneW_RK(3)

      integer*8 i, j, i1, iel, namelen
      integer*8 plot_imag, plot_real, plot_abs
      integer*8 debug, ui
      complex*16 z_tmp1
      complex*16 ii
      character*100 tchar
      integer*8 namelength


Cf2py intent(in) nel, npt, nnodes, neq_PW, bloch_vec,
Cf2py intent(in) table_nod, x, lat_vecs, lambda, n_eff_0
Cf2py intent(in) vec_coef_down, vec_coef_up,
Cf2py intent(in) index_pw_inv, ordre_ls, hz,
Cf2py intent(in) gmsh_file_pos, extra_name
Cf2py intent(in) plot_real, plot_imag, plot_abs

Cf2py depend(table_nod) nnodes, nel
Cf2py depend(index_pw_inv) neq_PW
Cf2py depend(x) npt
Cf2py depend(vec_coef_down) neq_PW
Cf2py depend(vec_coef_up) neq_PW


      dir_name = "in_plane_fields"
c
c  ii = sqrt(-1)
      ii = cmplx(0.0d0, 1.0d0)
c
      ui = 6
      debug = 0
c
      d = lat_vecs(1,1)
      dy = lat_vecs(2,2)
      pi = 3.141592653589793d0
      bloch1 = bloch_vec(1)
      bloch2 = bloch_vec(2)
      k_0 = (2.0d0*pi)/lambda
      k = real(n_eff_0*k_0)
      vec_kx = 2.0d0*pi/d
      vec_ky = 2.0d0*pi/dy
      if (debug .eq. 1) then
        write(ui,*) "gmsh_plot_PW: ", lambda, n_eff_0, bloch_vec
      endif
c
      if ( nnodes .ne. 6 ) then
        write(ui,*) "gmsh_plot_PW: problem nnodes = ", nnodes
        write(ui,*) "gmsh_plot_PW: nnodes should be equal to 6 !"
        write(ui,*) "gmsh_plot_PW: Aborting..."
        stop
      endif
C
      alloc_stat = 0
      allocate(sol_tmp(3,nnodes,nel), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "gmsh_plot_PW: Mem. allocation is unseccesfull"
        write(*,*) "alloc_stat (sol_tmp) = ", alloc_stat
        write(*,*) "Aborting..."
        stop
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (hz .lt. 0) then
        write(ui,*)
        write(ui,*) "gmsh_plot_PW: invalid value for hz"
        write(ui,*) "hz should not be negtive (for current version)"
        write(ui,*) "hz = ", hz
        write(ui,*) "gmsh_plot_PW: Aborting..."
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

      tchar=dir_name(1:namelength)// '/' // extra_name //'_PW_abs2_'
     *           // gmsh_file_pos(1:namelen) // '.pos'
      open (unit=26,file=tchar)
        write(26,*) "View.AdaptVisualizationGrid =1;"
        write(26,*) "View.IntervalsType = 3;"
        write(26,*) "View.Light = 0;"
        write(26,*) "View ""|",extra_name(1:2),"_t|^2:",
     *     "; ordre = ", ordre_ls, "; hz = ", hz, " "" {"
c
      tchar=dir_name(1:namelength)// '/' // extra_name // '_PWx_re_'
     * // gmsh_file_pos(1:namelen) // '.pos'
      open (unit=27,file=tchar)
        write(27,*) "View.AdaptVisualizationGrid =1;"
        write(27,*) "View.IntervalsType = 3;"
        write(27,*) "View.Light = 0;"
c        write(27,*) "View ""Re Ex: n = ", plot_val,
        write(27,*) "View ""Re ",extra_name(1:2),"x:",
     *     "; ordre = ", ordre_ls, "; hz = ", hz, " "" {"
c
      tchar=dir_name(1:namelength)// '/' // extra_name // '_PWy_re_'
     * // gmsh_file_pos(1:namelen) // '.pos'
      open (unit=28,file=tchar)
        write(28,*) "View.AdaptVisualizationGrid =1;"
        write(28,*) "View.IntervalsType = 3;"
        write(28,*) "View.Light = 0;"
        write(28,*) "View ""Re ",extra_name(1:2),"y:",
     *     "; ordre = ", ordre_ls, "; hz = ", hz, " "" {"
c
      tchar=dir_name(1:namelength)// '/' // extra_name // '_PWz_re_'
     * // gmsh_file_pos(1:namelen) // '.pos'
      open (unit=29,file=tchar)
        write(29,*) "View.AdaptVisualizationGrid =1;"
        write(29,*) "View.IntervalsType = 3;"
        write(29,*) "View.Light = 0;"
        write(29,*) "View ""Re ",extra_name(1:2),"z:",
     *     "; ordre = ", ordre_ls, "; hz = ", hz, " "" {"
c
      tchar=dir_name(1:namelength)// '/' // extra_name // '_PWv_re_'
     * // gmsh_file_pos(1:namelen) // '.pos'
      open (unit=30,file=tchar)
        write(30,*) "View.AdaptVisualizationGrid =1;"
        write(30,*) "View.IntervalsType = 3;"
        write(30,*) "View.Light = 0;"
        write(30,*) "View ""Re ",extra_name(1:2),":",
     *     "; ordre = ", ordre_ls, "; hz = ", hz, " "" {"
c
      do iel=1,nel
        do i=1,nnodes
          do j=1,3
            sol_tmp(j,i,iel) = 0.0d0
          enddo
        enddo
      enddo
c
      s = 1
      do px = -ordre_ls, ordre_ls
        do py = -ordre_ls, ordre_ls
            if ((px/d)**2 + (py/dy)**2
     *         .le. (ordre_ls/MAX(d, dy))**2) then
            alpha = bloch1 + vec_kx*px  ! Bloch vector along x
            beta  = bloch2 + vec_ky*py  ! Bloch vector along y
            z_tmp1 = k**2 - alpha**2 - beta**2
            gamma = SQRT(z_tmp1)
            chi_TE = SQRT(gamma/k_0)
            chi_TM = SQRT(gamma/(n_eff_0*k))
C            chi_TE = SQRT(gamma/k)
C            chi_TM = SQRT(gamma/k)
            norm = SQRT(alpha**2 + beta**2) !sqrt term
            s2 = index_pw_inv(s)
C
C           Plane waves order s
C  	  RE
            PlaneW_RE(1) = beta/norm   ! x component
            PlaneW_RE(2) = -alpha/norm ! y component
            PlaneW_RE(3) = 0.0d0       ! z component
C  	  RK
            PlaneW_RK(1) = alpha/norm  ! x component
            PlaneW_RK(2) = beta/norm   ! y component
            PlaneW_RK(3) = -norm/gamma ! z component
            s = s + 1

            P_up = EXP(ii*gamma*hz)    !  Introduce Propagation in +z
            P_down = 1.0d0 / P_up      !  Introduce Propagation in -z
c
            coef_RE_down = vec_coef_down(s2) * P_down / chi_TE
            coef_RE_up = vec_coef_up(s2) * P_up / chi_TE
c
            coef_RE_t = coef_RE_up + coef_RE_down  ! Superposition of the counter-porpagating wave
            coef_RE_z = coef_RE_up - coef_RE_down
c
            coef_RK_down = vec_coef_down(s2+neq_PW) * P_down * chi_TM
            coef_RK_up = vec_coef_up(s2+neq_PW) * P_up * chi_TM
c
            coef_RK_t = coef_RK_up + coef_RK_down
            coef_RK_z = coef_RK_up - coef_RK_down

            do iel=1,nel
              do i=1,nnodes
                i1 = table_nod(i,iel)
                xx_g(1) = x(1,i1)
                xx_g(2) = x(2,i1)
                xx_g(3) = hz
                r_tmp = alpha*xx_g(1) + beta*xx_g(2)
                val_exp = EXP(ii*r_tmp)
                do j=1,2 ! Transverse component
                  z_tmp1 = PlaneW_RE(j) * coef_RE_t * val_exp
                  sol_tmp(j,i,iel) = sol_tmp(j,i,iel) + z_tmp1
                  z_tmp1 = PlaneW_RK(j) * coef_RK_t * val_exp
                  sol_tmp(j,i,iel) = sol_tmp(j,i,iel) + z_tmp1
                enddo
                j=3 ! z component
                  z_tmp1 = PlaneW_RE(j) * coef_RE_z * val_exp
                  sol_tmp(j,i,iel) = sol_tmp(j,i,iel) + z_tmp1
                  z_tmp1 = PlaneW_RK(j) * coef_RK_z * val_exp
                  sol_tmp(j,i,iel) = sol_tmp(j,i,iel) + z_tmp1
              enddo
            enddo
          endif
        enddo
      enddo
c
        sol_max(4) = 0.0d0
        do iel=1,nel
          zz = hz
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
            sol_el_abs2(i) = 0.0
            do j=1,2
              z_tmp1 = sol_tmp(j,i,iel)
              sol_el(j,i) = z_tmp1
              sol_el_abs2(i) = sol_el_abs2(i) + abs(z_tmp1)**2
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
      tchar=dir_name(1:namelength)// '/' // extra_name // '_PW_x_im_'
     * // gmsh_file_pos(1:namelen) // '.pos'
      open (unit=27,file=tchar)
        write(27,*) "View.AdaptVisualizationGrid =1;"
        write(27,*) "View.IntervalsType = 3;"
        write(27,*) "View.Light = 0;"
        write(27,*) "View ""Im ",extra_name(1:2),"x:",
     *     "; ordre = ", ordre_ls, "; hz = ", hz, " "" {"
c
      tchar=dir_name(1:namelength)// '/' // extra_name // '_PW_y_im_'
     * // gmsh_file_pos(1:namelen) // '.pos'
      open (unit=28,file=tchar)
        write(28,*) "View.AdaptVisualizationGrid =1;"
        write(28,*) "View.IntervalsType = 3;"
        write(28,*) "View.Light = 0;"
        write(28,*) "View ""Im ",extra_name(1:2),"y:",
     *     "; ordre = ", ordre_ls, "; hz = ", hz, " "" {"
c
      tchar=dir_name(1:namelength)// '/' // extra_name // '_PW_z_im_'
     * // gmsh_file_pos(1:namelen) // '.pos'
      open (unit=29,file=tchar)
        write(29,*) "View.AdaptVisualizationGrid =1;"
        write(29,*) "View.IntervalsType = 3;"
        write(29,*) "View.Light = 0;"
        write(29,*) "View ""Im ",extra_name(1:2),"z:",
     *     "; ordre = ", ordre_ls, "; hz = ", hz, " "" {"
c
      tchar=dir_name(1:namelength)// '/' // extra_name // '_PW_v_im_'
     * // gmsh_file_pos(1:namelen) // '.pos'
      open (unit=30,file=tchar)
        write(30,*) "View.AdaptVisualizationGrid =1;"
        write(30,*) "View.IntervalsType = 3;"
        write(30,*) "View.Light = 0;"
        write(30,*) "View ""Im ",extra_name(1:2),":",
     *     "; ordre = ", ordre_ls, "; hz = ", hz, " "" {"
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
      tchar=dir_name(1:namelength)// '/' // extra_name // '_PW_x_abs_'
     * // gmsh_file_pos(1:namelen) // '.pos'
      open (unit=27,file=tchar)
        write(27,*) "View.AdaptVisualizationGrid =1;"
        write(27,*) "View.IntervalsType = 3;"
        write(27,*) "View.Light = 0;"
        write(27,*) "View ""|",extra_name(1:2),"x|:",
     *     "; ordre = ", ordre_ls, "; hz = ", hz, " "" {"
c
      tchar=dir_name(1:namelength)// '/' // extra_name // '_PW_y_abs_'
     * // gmsh_file_pos(1:namelen) // '.pos'
      open (unit=28,file=tchar)
        write(28,*) "View.AdaptVisualizationGrid =1;"
        write(28,*) "View.IntervalsType = 3;"
        write(28,*) "View.Light = 0;"
        write(28,*) "View ""|",extra_name(1:2),"y|:",
     *     "; ordre = ", ordre_ls, "; hz = ", hz, " "" {"
c
      tchar=dir_name(1:namelength)// '/' // extra_name // '_PW_z_abs_'
     * // gmsh_file_pos(1:namelen) // '.pos'
      open (unit=29,file=tchar)
        write(29,*) "View.AdaptVisualizationGrid =1;"
        write(29,*) "View.IntervalsType = 3;"
        write(29,*) "View.Light = 0;"
        write(29,*) "View ""|",extra_name(1:2),"z|:",
     *     "; ordre = ", ordre_ls, "; hz = ", hz, " "" {"
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
        write(ui,*) "gmsh_plot_PW: sol_max = ", sol_max
      endif
C
      deallocate(sol_tmp)
C
      return
      end
