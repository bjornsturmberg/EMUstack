
c   evecs(i) : contains the values of the solution for all points

      subroutine gmsh_plot_slice (E_H_field, nval, nel, npt, nnodes,
     *     type_el, nb_typ_el, n_eff,
     *     table_nod, x, beta, evecs, vec_coef, h,  lambda,
     *     gmsh_file_pos, Transf_XX, Offset_Z, Offset_X)


      implicit none
      integer*8 nval, nel, npt, nnodes, E_H_field, nb_typ_el
      double precision h,  lambda
      integer*8 table_nod(nnodes,nel), type_el(nel)
      complex*16 x(2,npt), beta(nval), vec_coef(2*nval)
      complex*16 n_eff(nb_typ_el)
c     evecs(3, 1..nnodes,nval, nel)          contains the values of the 3 components at P2 interpolation nodes
c     evecs(3, nnodes+1..nnodes+7,nval, nel) contains the values of Ez component at P3 interpolation nodes (per element: 6 edge-nodes and 1 interior node)
      complex*16 evecs(3,nnodes+7,nval,nel)
      character(*) gmsh_file_pos
      character*30 dir_name
c
c     Local variables
      integer alloc_stat
      complex*16, dimension(:,:,:,:), allocatable :: sol_3d_X
      complex*16, dimension(:,:,:,:), allocatable :: sol_3d_Y
      complex*16, dimension(:,:,:,:), allocatable :: sol_3d_D_1
      complex*16, dimension(:,:,:,:), allocatable :: sol_3d_D_2

      double precision, allocatable :: x_P1_X(:,:), x_P1_Y(:,:)
      double precision, allocatable :: x_P1_D_1(:,:)
      double precision, allocatable :: x_P1_D_2(:,:)

      integer*8, allocatable :: nb_visit_X(:,:), nb_visit_Y(:,:)
      integer*8, allocatable :: nb_visit_D_1(:,:)
      integer*8, allocatable :: nb_visit_D_2(:,:)

      integer*8, allocatable :: info_elem(:)
      integer*8, dimension(:), allocatable :: map_p1, inv_map_p1

      integer*8 nnodes_0, n_quad
      parameter (nnodes_0 = 6, n_quad=4)
      double precision xel_2d(2,nnodes_0)
C      Quadrangle element
      double precision xel(3,n_quad)
      complex*16 sol_el(3,n_quad)
      double precision sol_el_abs2(n_quad)
      double precision Transf_XX, Transf_YY, Transf_ZZ
      double precision Offset_Z, Offset_X

      complex*16 P_down, P_up, coef_down, coef_up, coef_t, coef_z
      complex*16 ii, z_tmp1, r_index
      double precision hz, dz, zz, r_tmp

C      prism elements
      integer*8 nel_3d, npt_3d
C      Number of vertices of the 2D FEM mesh
      integer*8 npt_p1
C      Resolution: number of points over the thickness h
      integer*8 npt_h, i_h, npt_3d_p1
      integer*8 i1, i, j, k, j1, iel, inod, ival, debug, ui
      integer*8 namelen_gmsh, namelen_dir, namelen_tchar

cc      integer*8 gmsh_type_prism, choice_type, typ_e
cc      integer*8 number_tags, physic_tag, list_tag(6)
cc      integer*8, dimension(:), allocatable :: type_data

      character*100 tchar
      character*1 tE_H

      integer*8 iFrame, nFrame
      double precision pi, phase
      complex*16 exp_phase
      character tval*4, buf*3

      integer*8 nx, ny, nnodes_P1
      parameter (nnodes_P1 = 2)
      integer*8 nel_X, nel_Y, nel_D_1, nel_D_2
      double precision x_min, x_max, y_min, y_max, ls_edge(4)
      double precision dx_1, dy_1, dx_2, dy_2
      double precision dx_3, dy_3, dx_4, dy_4
      double precision d, lx, ly, xx, yy, rx, ry

Cf2py intent(in) E_H_field, nval, nel, npt, nnodes
Cf2py intent(in) type_el, nb_typ_el, n_eff, table_nod, x, beta
Cf2py intent(in) evecs, vec_coef, h, lambda
Cf2py intent(in) gmsh_file_pos

Cf2py depend(table_nod) nnodes, nel
Cf2py depend(type_el) nel
Cf2py depend(x) npt
Cf2py depend(n_eff) nb_typ_el
Cf2py depend(evecs) nnodes, nval, nel
Cf2py depend(vec_coef) nval
Cf2py depend(beta) nval

c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c  ii = sqrt(-1)
      ii = cmplx(0.0d0, 1.0d0)
c
      ui = 6
      debug = 0
      dir_name = "fields_vertically/gmsh_BMs"

      npt_h = 10 * (h/lambda + 1) ! At least 10 nodes per wavelength
      npt_3d = npt * npt_h
      nel_3d = nel * (npt_h - 1)

      if (debug .eq. 1) then
        write(ui,*)
        write(ui,*) "gmsh_plot_slice: h/lambda = ", h/lambda
        write(ui,*) "gmsh_plot_slice:    npt_h = ", npt_h
        write(ui,*) "gmsh_plot_slice:   npt_3d = ", npt_3d
        write(ui,*) "gmsh_plot_slice:   nel_3d = ", nel_3d
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if ( nnodes .ne. 6 ) then
        write(ui,*) "gmsh_plot_slice: problem nnodes = ", nnodes
        write(ui,*) "nnodes should be equal to 6 !"
        write(ui,*) "gmsh_plot_slice: Aborting..."
        stop
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      x_min = x(1,1)
      x_max = x(1,1)
      do i=1,npt
        xx = x(1,i)
        if(xx .lt. x_min) x_min = xx
        if(xx .gt. x_max) x_max = xx
      enddo
      y_min = x(2,1)
      y_max = x(2,1)
      do i=1,npt
        yy = x(2,i)
        if(yy .lt. y_min) y_min = yy
        if(yy .gt. y_max) y_max = yy
      enddo

      ls_edge(1) = x_min
      ls_edge(2) = x_max
      ls_edge(3) = y_min
      ls_edge(4) = y_max

      lx = x_max - x_min  !  Length in the x direction
      ly = y_max - y_min  !  Length in the y direction

      rx = sqrt(nel * lx / (2.0d0 * ly))
      ry = rx * ly/lx
c     rx and ry and such that : rx * ry = 2 * nel

      nx = anint(rx)
      ny = anint(ry)

      if (debug .eq. 1) then
        write(ui,*)
        write(ui,*) "gmsh_plot_slice: nx, ny = ", nx, ny
        write(ui,*) "gmsh_plot_slice: rx, ry = ", rx, ry
        write(ui,*) "gmsh_plot_slice: lx, ly = ", lx, ly
        write(ui,*) "gmsh_plot_slice: 2 rx * ry = ", 2 * rx * ry
        write(ui,*) "gmsh_plot_slice: nel = ", nel
        write(ui,*)
      endif


      nel_X = nx
      nel_Y = ny
      nel_D_1 = max(nx ,ny)
      nel_D_2 = max(nx ,ny)


c     Diagonal elements of the coordinate transformation matrix
c     Scaling factor = Transf_XX
C      Transf_XX = 2.0
      Transf_YY = Transf_XX
      Transf_ZZ = Transf_XX

c     Translation of the view along Z-axis
C      Offset_Z = 20.0


c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      alloc_stat = 0

      allocate(sol_3d_X(3,2,nel_X,npt_h),
     *     sol_3d_Y(3,2,nel_Y,npt_h), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*)
        write(*,*) "gmsh_plot_slice: ",
     *     "The allocation is unsuccessful"
        write(*,*) "alloc_stat = ", alloc_stat
        write(*,*) "Not enough memory for the array sol_3d_X"
        write(*,*) "nel_X, nel_Y, npt_h = ", nel_X,nel_Y,npt_h
        write(*,*) "gmsh_plot_slice: Aborting..."
        stop
      endif

      allocate(sol_3d_D_1(3,2,nel_D_1,npt_h),
     *     sol_3d_D_2(3,2,nel_D_2,npt_h), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*)
        write(*,*) "gmsh_plot_slice: ",
     *     "The allocation is unsuccessful"
        write(*,*) "alloc_stat = ", alloc_stat
        write(*,*) "Not enough memory for the array sol_3d_D_1"
        write(*,*) "nel_D_1, nel_D_2, npt_h = ",
     *     nel_D_1,nel_D_2,npt_h
        write(*,*) "gmsh_plot_slice: Aborting..."
        stop
      endif


      allocate(x_P1_X(2,nel_X+1), nb_visit_X(2,nel_X),
     *     x_P1_Y(2,nel_Y+1), nb_visit_Y(2,nel_Y),
     *     STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*)
        write(*,*) "gmsh_plot_slice: ",
     *     "The allocation is unsuccessful"
        write(*,*) "alloc_stat = ", alloc_stat
        write(*,*) "Not enough memory for the array x_P1_X"
        write(*,*) "nel_X, nel_Y = ", nel_X,nel_Y
        write(*,*) "gmsh_plot_slice: Aborting..."
        stop
      endif

      allocate(x_P1_D_1(2,nel_D_1+1), nb_visit_D_1(2,nel_D_1),
     *     x_P1_D_2(2,nel_D_2+1), nb_visit_D_2(2,nel_D_2),
     *     STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*)
        write(*,*) "gmsh_plot_slice: ",
     *     "The allocation is unsuccessful"
        write(*,*) "alloc_stat = ", alloc_stat
        write(*,*) "Not enough memory for the array x_P1_D_1"
        write(*,*) "nel_D_1, nel_D_2 = ", nel_D_1,nel_D_2
        write(*,*) "gmsh_plot_slice: Aborting..."
        stop
      endif

      allocate(info_elem(nel),STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*)
        write(*,*) "gmsh_plot_slice: ",
     *     "The allocation is unsuccessful"
        write(*,*) "alloc_stat = ", alloc_stat
        write(*,*) "Not enough memory for the array info_elem"
        write(*,*) "nel = ", nel
        write(*,*) "gmsh_plot_slice: Aborting..."
        stop
      endif

c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      dx_1 = lx/dble(nel_X)
      dy_1 = 0
      do i=1,nel_X+1
        x_P1_X(1,i) = x_min + (i-1)*dx_1
        x_P1_X(2,i) = (y_min+y_max) / 2.0d0
      enddo

      dx_2 = 0
      dy_2 = ly/dble(nel_Y)
      do i=1,nel_Y+1
        x_P1_Y(1,i) = (x_min+x_max) / 2.0d0
        x_P1_Y(2,i) = y_min + (i-1)*dy_2
      enddo

      dx_3 = lx/dble(nel_D_1)
      dy_3 = ly/dble(nel_D_1)
      do i=1,nel_D_1+1
        x_P1_D_1(1,i) = x_min + (i-1)*dx_3
        x_P1_D_1(2,i) = y_min + (i-1)*dy_3
      enddo

      dx_4 =  lx/dble(nel_D_2)
      dy_4 = -ly/dble(nel_D_2)
      do i=1,nel_D_2+1
        x_P1_D_2(1,i) = x_min + (i-1)*dx_4
        x_P1_D_2(2,i) = y_max + (i-1)*dy_4
      enddo

c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    Initialise sol

      do i_h=1,npt_h
        do iel=1,nel_X
          do inod=1,nnodes_P1
            do j=1,3
              sol_3d_X(j,inod,iel,i_h) = 0.0d0
            enddo
          enddo
        enddo
      enddo

      do i_h=1,npt_h
        do iel=1,nel_Y
          do inod=1,nnodes_P1
            do j=1,3
              sol_3d_Y(j,inod,iel,i_h) = 0.0d0
            enddo
          enddo
        enddo
      enddo

      do i_h=1,npt_h
        do iel=1,nel_D_1
          do inod=1,nnodes_P1
            do j=1,3
              sol_3d_D_1(j,inod,iel,i_h) = 0.0d0
            enddo
          enddo
        enddo
      enddo

      do i_h=1,npt_h
        do iel=1,nel_D_2
          do inod=1,nnodes_P1
            do j=1,3
              sol_3d_D_2(j,inod,iel,i_h) = 0.0d0
            enddo
          enddo
        enddo
      enddo

      do iel=1,nel_X
        do inod=1,nnodes_P1
          nb_visit_X(inod,iel) = 0
        enddo
      enddo

      do iel=1,nel_Y
        do inod=1,nnodes_P1
          nb_visit_Y(inod,iel) = 0
        enddo
      enddo

      do iel=1,nel_D_1
        do inod=1,nnodes_P1
          nb_visit_D_1(inod,iel) = 0
        enddo
      enddo

      do iel=1,nel_D_2
        do inod=1,nnodes_P1
          nb_visit_D_2(inod,iel) = 0
        enddo
      enddo

      do iel=1,nel
        info_elem(iel) = 0
      enddo
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
          ! Taking into account the change of variable for Ez
          coef_z = (coef_up - coef_down)/beta(ival)
          do iel=1,nel
            do j=1,nnodes
              j1 = table_nod(j,iel)
              xel_2d(1,j) = x(1,j1)
              xel_2d(2,j) = x(2,j1)
            enddo
            if (info_elem(iel) >= 0) then
              call slice_interp(nel, nval, iel, ival, i_h,
     *         nnodes, nx, ny, nel_X, nel_Y, nel_D_1, nel_D_2,
     *         npt_h, ls_edge, xel_2d, evecs, coef_t, coef_z,
     *        sol_3d_X, sol_3d_Y, sol_3d_D_1, sol_3d_D_2,
     *        nb_visit_X, nb_visit_Y, nb_visit_D_1, nb_visit_D_2,
     *        info_elem)
            endif
          enddo
        enddo
      enddo
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    Average values
c
      do i_h=1,npt_h
        do iel=1,nel_X
          do inod=1,nnodes_P1
            k = max(1, nb_visit_X(inod,iel))
            do j=1,3
              z_tmp1 = sol_3d_X(j,inod,iel,i_h)
              sol_3d_X(j,inod,iel,i_h) = z_tmp1 / dble(k)
            enddo
          enddo
        enddo
      enddo

      do i_h=1,npt_h
        do iel=1,nel_Y
          do inod=1,nnodes_P1
            k = max(1, nb_visit_Y(inod,iel))
            do j=1,3
              z_tmp1 = sol_3d_Y(j,inod,iel,i_h)
              sol_3d_Y(j,inod,iel,i_h) = z_tmp1 / dble(k)
            enddo
          enddo
        enddo
      enddo

      do i_h=1,npt_h
        do iel=1,nel_D_1
          do inod=1,nnodes_P1
            k = max(1, nb_visit_D_1(inod,iel))
            do j=1,3
              z_tmp1 = sol_3d_D_1(j,inod,iel,i_h)
              sol_3d_D_1(j,inod,iel,i_h) = z_tmp1 / dble(k)
            enddo
          enddo
        enddo
      enddo

      do i_h=1,npt_h
        do iel=1,nel_D_2
          do inod=1,nnodes_P1
            k = max(1, nb_visit_D_2(inod,iel))
            do j=1,3
              z_tmp1 = sol_3d_D_2(j,inod,iel,i_h)
              sol_3d_D_2(j,inod,iel,i_h) = z_tmp1 / dble(k)
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

c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      tchar=dir_name(1:namelen_dir)// '/' //
     *           gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // '_abs2_sl_X.pos'
      open (unit=26,file=tchar)
        write(26,*) "View.IntervalsType = 3;"
        write(26,*) "General.Trackball = 0;"
        write(26,*) "General.RotationX = -90;"
        write(26,*) "General.RotationY = 0;"
        write(26,*) "General.RotationZ = 0;"
        write(26,"('View.TransformXX = ',f10.6,';')") Transf_XX
        write(26,"('View.TransformYY = ',f10.6,';')") Transf_YY
        write(26,"('View.TransformZZ = ',f10.6,';')") Transf_ZZ
        write(26,"('View.OffsetX = ',f10.6,';')") Offset_X
        write(26,"('View.OffsetZ = ',f10.6,';')") Offset_Z
C         write(26,*) "View.Axes = 2;"
        write(26,*) "View ""|",tE_H,"_t|^2 ",
     *     " "" {"

      tchar=dir_name(1:namelen_dir)// '/' //
     *           gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // 'x_re_sl_X.pos'
      open (unit=27,file=tchar)
        write(27,*) "View.IntervalsType = 3;"
        write(27,*) "General.Trackball = 0;"
        write(27,*) "General.RotationX = -90;"
        write(27,*) "General.RotationY = 0;"
        write(27,*) "General.RotationZ = 0;"
        write(27,"('View.TransformXX = ',f10.6,';')") Transf_XX
        write(27,"('View.TransformYY = ',f10.6,';')") Transf_YY
        write(27,"('View.TransformZZ = ',f10.6,';')") Transf_ZZ
        write(27,"('View.OffsetX = ',f10.6,';')") Offset_X
        write(27,"('View.OffsetZ = ',f10.6,';')") Offset_Z
C         write(27,*) "View.Axes = 2;"
        write(27,*) "View ""Re ",tE_H,"x ",
     *     " "" {"

      tchar=dir_name(1:namelen_dir)// '/' //
     *           gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // 'y_re_sl_X.pos'
      open (unit=28,file=tchar)
        write(28,*) "View.IntervalsType = 3;"
        write(28,*) "General.Trackball = 0;"
        write(28,*) "General.RotationX = -90;"
        write(28,*) "General.RotationY = 0;"
        write(28,*) "General.RotationZ = 0;"
        write(28,"('View.TransformXX = ',f10.6,';')") Transf_XX
        write(28,"('View.TransformYY = ',f10.6,';')") Transf_YY
        write(28,"('View.TransformZZ = ',f10.6,';')") Transf_ZZ
        write(28,"('View.OffsetX = ',f10.6,';')") Offset_X
        write(28,"('View.OffsetZ = ',f10.6,';')") Offset_Z
C         write(28,*) "View.Axes = 2;"
        write(28,*) "View ""Re ",tE_H,"y ",
     *     " "" {"

      tchar=dir_name(1:namelen_dir)// '/' //
     *           gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // 'z_re_sl_X.pos'
      open (unit=29,file=tchar)
        write(29,*) "View.IntervalsType = 3;"
        write(29,*) "General.Trackball = 0;"
        write(29,*) "General.RotationX = -90;"
        write(29,*) "General.RotationY = 0;"
        write(29,*) "General.RotationZ = 0;"
        write(29,"('View.TransformXX = ',f10.6,';')") Transf_XX
        write(29,"('View.TransformYY = ',f10.6,';')") Transf_YY
        write(29,"('View.TransformZZ = ',f10.6,';')") Transf_ZZ
        write(29,"('View.OffsetX = ',f10.6,';')") Offset_X
        write(29,"('View.OffsetZ = ',f10.6,';')") Offset_Z
C         write(29,*) "View.Axes = 2;"
        write(29,*) "View ""Re ",tE_H,"z ",
     *     " "" {"

      tchar=dir_name(1:namelen_dir)// '/' //
     *           gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // 'v_re_sl_X.pos'
      open (unit=30,file=tchar)
        write(30,*) "View.IntervalsType = 3;"
        write(30,*) "General.Trackball = 0;"
        write(30,*) "General.RotationX = -90;"
        write(30,*) "General.RotationY = 0;"
        write(30,*) "General.RotationZ = 0;"
        write(30,"('View.TransformXX = ',f10.6,';')") Transf_XX
        write(30,"('View.TransformYY = ',f10.6,';')") Transf_YY
        write(30,"('View.TransformZZ = ',f10.6,';')") Transf_ZZ
        write(30,"('View.OffsetX = ',f10.6,';')") Offset_X
        write(30,"('View.OffsetZ = ',f10.6,';')") Offset_Z
C         write(30,*) "View.Axes = 2;"
        write(30,*) "View ""Re ",tE_H, " "" {"


      dz = h/dble(npt_h-1)
      do i_h=1,npt_h-1
        zz = - (i_h-1)*dz  ! hz=0 => top interface; hz=-h => bottom interface
        do iel=1,nel_X
          do inod=1,2
            xel(1,inod) = x_P1_X(1,iel+inod-1)
            xel(2,inod) = x_P1_X(2,iel+inod-1)
            xel(3,inod) = zz
          enddo
            xel(1,3) = xel(1,2)  ! Quadrangle element
            xel(2,3) = xel(2,2)
            xel(3,3) = zz - dz
            xel(1,4) = xel(1,1)  ! Quadrangle element
            xel(2,4) = xel(2,1)
            xel(3,4) = zz - dz

          do inod=1,2
            sol_el_abs2(inod) = 0.0
            do j=1,3
              z_tmp1 = sol_3d_X(j,inod,iel,i_h)
              sol_el(j,inod) = z_tmp1
              sol_el_abs2(inod) = sol_el_abs2(inod) +
     *           abs(z_tmp1)**2
            enddo
            sol_el_abs2(5-inod) = 0.0
            do j=1,3
              z_tmp1 = sol_3d_X(j,inod,iel,i_h+1)
              sol_el(j,5-inod) = z_tmp1
              sol_el_abs2(5-inod) = sol_el_abs2(5-inod) +
     *           abs(z_tmp1)**2
            enddo
          enddo
          write(26,10) xel, sol_el_abs2
          write(27,10) xel, (dble(sol_el(1,j)),j=1,n_quad)
          write(28,10) xel, (dble(sol_el(2,j)),j=1,n_quad)
          write(29,10) xel, (dble(sol_el(3,j)),j=1,n_quad)
          write(30,11) xel,
     *     ((dble(sol_el(i,j)),i=1,3),j=1,n_quad)
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
      tchar=dir_name(1:namelen_dir)// '/' //
     *           gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // '_abs2_sl_Y.pos'
      open (unit=26,file=tchar)
        write(26,*) "View.IntervalsType = 3;"
        write(26,*) "General.Trackball = 0;"
        write(26,*) "General.RotationX = -90;"
        write(26,*) "General.RotationY = 0;"
        write(26,*) "General.RotationZ = -90;"
        write(26,"('View.TransformXX = ',f10.6,';')") Transf_XX
        write(26,"('View.TransformYY = ',f10.6,';')") Transf_YY
        write(26,"('View.TransformZZ = ',f10.6,';')") Transf_ZZ
        write(26,"('View.OffsetX = ',f10.6,';')") Offset_X
        write(26,"('View.OffsetZ = ',f10.6,';')") Offset_Z
C         write(26,*) "View.Axes = 2;"
        write(26,*) "View ""|",tE_H,"_t|^2 ",
     *     " "" {"

      tchar=dir_name(1:namelen_dir)// '/' //
     *           gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // 'x_re_sl_Y.pos'
      open (unit=27,file=tchar)
        write(27,*) "View.IntervalsType = 3;"
        write(27,*) "General.Trackball = 0;"
        write(27,*) "General.RotationX = -90;"
        write(27,*) "General.RotationY = 0;"
        write(27,*) "General.RotationZ = -90;"
        write(27,"('View.TransformXX = ',f10.6,';')") Transf_XX
        write(27,"('View.TransformYY = ',f10.6,';')") Transf_YY
        write(27,"('View.TransformZZ = ',f10.6,';')") Transf_ZZ
        write(27,"('View.OffsetX = ',f10.6,';')") Offset_X
        write(27,"('View.OffsetZ = ',f10.6,';')") Offset_Z
C         write(27,*) "View.Axes = 2;"
        write(27,*) "View ""Re ",tE_H,"x ",
     *     " "" {"

      tchar=dir_name(1:namelen_dir)// '/' //
     *           gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // 'y_re_sl_Y.pos'
      open (unit=28,file=tchar)
        write(28,*) "View.IntervalsType = 3;"
        write(28,*) "General.Trackball = 0;"
        write(28,*) "General.RotationX = -90;"
        write(28,*) "General.RotationY = 0;"
        write(28,*) "General.RotationZ = -90;"
        write(28,"('View.TransformXX = ',f10.6,';')") Transf_XX
        write(28,"('View.TransformYY = ',f10.6,';')") Transf_YY
        write(28,"('View.TransformZZ = ',f10.6,';')") Transf_ZZ
        write(28,"('View.OffsetX = ',f10.6,';')") Offset_X
        write(28,"('View.OffsetZ = ',f10.6,';')") Offset_Z
C         write(28,*) "View.Axes = 2;"
        write(28,*) "View ""Re ",tE_H,"y ",
     *     " "" {"

      tchar=dir_name(1:namelen_dir)// '/' //
     *           gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // 'z_re_sl_Y.pos'
      open (unit=29,file=tchar)
        write(29,*) "View.IntervalsType = 3;"
        write(29,*) "General.Trackball = 0;"
        write(29,*) "General.RotationX = -90;"
        write(29,*) "General.RotationY = 0;"
        write(29,*) "General.RotationZ = -90;"
        write(29,"('View.TransformXX = ',f10.6,';')") Transf_XX
        write(29,"('View.TransformYY = ',f10.6,';')") Transf_YY
        write(29,"('View.TransformZZ = ',f10.6,';')") Transf_ZZ
        write(29,"('View.OffsetX = ',f10.6,';')") Offset_X
        write(29,"('View.OffsetZ = ',f10.6,';')") Offset_Z
C         write(29,*) "View.Axes = 2;"
        write(29,*) "View ""Re ",tE_H,"z ",
     *     " "" {"

      tchar=dir_name(1:namelen_dir)// '/' //
     *           gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // 'v_re_sl_Y.pos'
      open (unit=30,file=tchar)
        write(30,*) "View.IntervalsType = 3;"
        write(30,*) "General.Trackball = 0;"
        write(30,*) "General.RotationX = -90;"
        write(30,*) "General.RotationY = 0;"
        write(30,*) "General.RotationZ = -90;"
        write(30,"('View.TransformXX = ',f10.6,';')") Transf_XX
        write(30,"('View.TransformYY = ',f10.6,';')") Transf_YY
        write(30,"('View.TransformZZ = ',f10.6,';')") Transf_ZZ
        write(30,"('View.OffsetX = ',f10.6,';')") Offset_X
        write(30,"('View.OffsetZ = ',f10.6,';')") Offset_Z
C         write(30,*) "View.Axes = 2;"
        write(30,*) "View ""Re ",tE_H, " "" {"


      dz = h/dble(npt_h-1)
      do i_h=1,npt_h-1
        zz = - (i_h-1)*dz  ! hz=0 => top interface; hz=-h => bottom interface
        do iel=1,nel_Y
          do inod=1,2
            xel(1,inod) = x_P1_Y(1,iel+inod-1)
            xel(2,inod) = x_P1_Y(2,iel+inod-1)
            xel(3,inod) = zz
          enddo
            xel(1,3) = xel(1,2)  ! Quadrangle element
            xel(2,3) = xel(2,2)
            xel(3,3) = zz - dz
            xel(1,4) = xel(1,1)  ! Quadrangle element
            xel(2,4) = xel(2,1)
            xel(3,4) = zz - dz

          do inod=1,2
            sol_el_abs2(inod) = 0.0
            do j=1,3
              z_tmp1 = sol_3d_Y(j,inod,iel,i_h)
              sol_el(j,inod) = z_tmp1
              sol_el_abs2(inod) = sol_el_abs2(inod) +
     *           abs(z_tmp1)**2
            enddo
            sol_el_abs2(5-inod) = 0.0
            do j=1,3
              z_tmp1 = sol_3d_Y(j,inod,iel,i_h+1)
              sol_el(j,5-inod) = z_tmp1
              sol_el_abs2(5-inod) = sol_el_abs2(5-inod) +
     *           abs(z_tmp1)**2
            enddo
          enddo
          write(26,10) xel, sol_el_abs2
          write(27,10) xel, (dble(sol_el(1,j)),j=1,n_quad)
          write(28,10) xel, (dble(sol_el(2,j)),j=1,n_quad)
          write(29,10) xel, (dble(sol_el(3,j)),j=1,n_quad)
          write(30,11) xel,
     *     ((dble(sol_el(i,j)),i=1,3),j=1,n_quad)
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
      tchar=dir_name(1:namelen_dir)// '/' //
     *           gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // '_abs2_sl_D_1.pos'
      open (unit=26,file=tchar)
        write(26,*) "View.IntervalsType = 3;"
        write(26,*) "General.Trackball = 0;"
        write(26,*) "General.RotationX = -90;"
        write(26,*) "General.RotationY = 0;"
        write(26,*) "General.RotationZ = 0;"
        write(26,"('View.TransformXX = ',f10.6,';')") Transf_XX
        write(26,"('View.TransformYY = ',f10.6,';')") Transf_YY
        write(26,"('View.TransformZZ = ',f10.6,';')") Transf_ZZ
        write(26,"('View.OffsetX = ',f10.6,';')") Offset_X
        write(26,"('View.OffsetZ = ',f10.6,';')") Offset_Z
C         write(26,*) "View.Axes = 2;"
        write(26,*) "View ""|",tE_H,"_t|^2 ",
     *     " "" {"

      tchar=dir_name(1:namelen_dir)// '/' //
     *           gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // 'x_re_sl_D_1.pos'
      open (unit=27,file=tchar)
        write(27,*) "View.IntervalsType = 3;"
        write(27,*) "General.Trackball = 0;"
        write(27,*) "General.RotationX = -90;"
        write(27,*) "General.RotationY = 0;"
        write(27,*) "General.RotationZ = 0;"
        write(27,"('View.TransformXX = ',f10.6,';')") Transf_XX
        write(27,"('View.TransformYY = ',f10.6,';')") Transf_YY
        write(27,"('View.TransformZZ = ',f10.6,';')") Transf_ZZ
        write(27,"('View.OffsetX = ',f10.6,';')") Offset_X
        write(27,"('View.OffsetZ = ',f10.6,';')") Offset_Z
C         write(27,*) "View.Axes = 2;"
        write(27,*) "View ""Re ",tE_H,"x ",
     *     " "" {"

      tchar=dir_name(1:namelen_dir)// '/' //
     *           gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // 'y_re_sl_D_1.pos'
      open (unit=28,file=tchar)
        write(28,*) "View.IntervalsType = 3;"
        write(28,*) "General.Trackball = 0;"
        write(28,*) "General.RotationX = -90;"
        write(28,*) "General.RotationY = 0;"
        write(28,*) "General.RotationZ = 0;"
        write(28,"('View.TransformXX = ',f10.6,';')") Transf_XX
        write(28,"('View.TransformYY = ',f10.6,';')") Transf_YY
        write(28,"('View.TransformZZ = ',f10.6,';')") Transf_ZZ
        write(28,"('View.OffsetX = ',f10.6,';')") Offset_X
        write(28,"('View.OffsetZ = ',f10.6,';')") Offset_Z
C         write(28,*) "View.Axes = 2;"
        write(28,*) "View ""Re ",tE_H,"y ",
     *     " "" {"

      tchar=dir_name(1:namelen_dir)// '/' //
     *           gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // 'z_re_sl_D_1.pos'
      open (unit=29,file=tchar)
        write(29,*) "View.IntervalsType = 3;"
        write(29,*) "General.Trackball = 0;"
        write(29,*) "General.RotationX = -90;"
        write(29,*) "General.RotationY = 0;"
        write(29,*) "General.RotationZ = 0;"
        write(29,"('View.TransformXX = ',f10.6,';')") Transf_XX
        write(29,"('View.TransformYY = ',f10.6,';')") Transf_YY
        write(29,"('View.TransformZZ = ',f10.6,';')") Transf_ZZ
        write(29,"('View.OffsetX = ',f10.6,';')") Offset_X
        write(29,"('View.OffsetZ = ',f10.6,';')") Offset_Z
C         write(29,*) "View.Axes = 2;"
        write(29,*) "View ""Re ",tE_H,"z ",
     *     " "" {"

      tchar=dir_name(1:namelen_dir)// '/' //
     *           gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // 'v_re_sl_D_1.pos'
      open (unit=30,file=tchar)
        write(30,*) "View.IntervalsType = 3;"
        write(30,*) "General.Trackball = 0;"
        write(30,*) "General.RotationX = -90;"
        write(30,*) "General.RotationY = 0;"
        write(30,*) "General.RotationZ = 0;"
        write(30,"('View.TransformXX = ',f10.6,';')") Transf_XX
        write(30,"('View.TransformYY = ',f10.6,';')") Transf_YY
        write(30,"('View.TransformZZ = ',f10.6,';')") Transf_ZZ
        write(30,"('View.OffsetX = ',f10.6,';')") Offset_X
        write(30,"('View.OffsetZ = ',f10.6,';')") Offset_Z
C         write(30,*) "View.Axes = 2;"
        write(30,*) "View ""Re ",tE_H, " "" {"


      dz = h/dble(npt_h-1)
      do i_h=1,npt_h-1
        zz = - (i_h-1)*dz  ! hz=0 => top interface; hz=-h => bottom interface
        do iel=1,nel_D_1
          do inod=1,2
            xel(1,inod) = x_P1_D_1(1,iel+inod-1)
            xel(2,inod) = x_P1_D_1(2,iel+inod-1)
            xel(3,inod) = zz
          enddo
            xel(1,3) = xel(1,2)  ! Quadrangle element
            xel(2,3) = xel(2,2)
            xel(3,3) = zz - dz
            xel(1,4) = xel(1,1)  ! Quadrangle element
            xel(2,4) = xel(2,1)
            xel(3,4) = zz - dz

          do inod=1,2
            sol_el_abs2(inod) = 0.0
            do j=1,3
              z_tmp1 = sol_3d_D_1(j,inod,iel,i_h)
              sol_el(j,inod) = z_tmp1
              sol_el_abs2(inod) = sol_el_abs2(inod) +
     *           abs(z_tmp1)**2
            enddo
            sol_el_abs2(5-inod) = 0.0
            do j=1,3
              z_tmp1 = sol_3d_D_1(j,inod,iel,i_h+1)
              sol_el(j,5-inod) = z_tmp1
              sol_el_abs2(5-inod) = sol_el_abs2(5-inod) +
     *           abs(z_tmp1)**2
            enddo
          enddo
          write(26,10) xel, sol_el_abs2
          write(27,10) xel, (dble(sol_el(1,j)),j=1,n_quad)
          write(28,10) xel, (dble(sol_el(2,j)),j=1,n_quad)
          write(29,10) xel, (dble(sol_el(3,j)),j=1,n_quad)
          write(30,11) xel,
     *     ((dble(sol_el(i,j)),i=1,3),j=1,n_quad)
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
      tchar=dir_name(1:namelen_dir)// '/' //
     *           gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // '_abs2_sl_D_2.pos'
      open (unit=26,file=tchar)
        write(26,*) "View.IntervalsType = 3;"
        write(26,*) "General.Trackball = 0;"
        write(26,*) "General.RotationX = 90;"
        write(26,*) "General.RotationY = 180;"
        write(26,*) "General.RotationZ = 90;"
        write(26,"('View.TransformXX = ',f10.6,';')") Transf_XX
        write(26,"('View.TransformYY = ',f10.6,';')") Transf_YY
        write(26,"('View.TransformZZ = ',f10.6,';')") Transf_ZZ
        write(26,"('View.OffsetX = ',f10.6,';')") Offset_X
        write(26,"('View.OffsetZ = ',f10.6,';')") Offset_Z
C         write(26,*) "View.Axes = 2;"
        write(26,*) "View ""|",tE_H,"_t|^2 ",
     *     " "" {"

      tchar=dir_name(1:namelen_dir)// '/' //
     *           gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // 'x_re_sl_D_2.pos'
      open (unit=27,file=tchar)
        write(27,*) "View.IntervalsType = 3;"
        write(27,*) "General.Trackball = 0;"
        write(27,*) "General.RotationX = 90;"
        write(27,*) "General.RotationY = 180;"
        write(27,*) "General.RotationZ = 90;"
        write(27,"('View.TransformXX = ',f10.6,';')") Transf_XX
        write(27,"('View.TransformYY = ',f10.6,';')") Transf_YY
        write(27,"('View.TransformZZ = ',f10.6,';')") Transf_ZZ
        write(27,"('View.OffsetX = ',f10.6,';')") Offset_X
        write(27,"('View.OffsetZ = ',f10.6,';')") Offset_Z
C         write(27,*) "View.Axes = 2;"
        write(27,*) "View ""Re ",tE_H,"x ",
     *     " "" {"

      tchar=dir_name(1:namelen_dir)// '/' //
     *           gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // 'y_re_sl_D_2.pos'
      open (unit=28,file=tchar)
        write(28,*) "View.IntervalsType = 3;"
        write(28,*) "General.Trackball = 0;"
        write(28,*) "General.RotationX = 90;"
        write(28,*) "General.RotationY = 180;"
        write(28,*) "General.RotationZ = 90;"
        write(28,"('View.TransformXX = ',f10.6,';')") Transf_XX
        write(28,"('View.TransformYY = ',f10.6,';')") Transf_YY
        write(28,"('View.TransformZZ = ',f10.6,';')") Transf_ZZ
        write(28,"('View.OffsetX = ',f10.6,';')") Offset_X
        write(28,"('View.OffsetZ = ',f10.6,';')") Offset_Z
C         write(28,*) "View.Axes = 2;"
        write(28,*) "View ""Re ",tE_H,"y ",
     *     " "" {"

      tchar=dir_name(1:namelen_dir)// '/' //
     *           gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // 'z_re_sl_D_2.pos'
      open (unit=29,file=tchar)
        write(29,*) "View.IntervalsType = 3;"
        write(29,*) "General.Trackball = 0;"
        write(29,*) "General.RotationX = 90;"
        write(29,*) "General.RotationY = 180;"
        write(29,*) "General.RotationZ = 90;"
        write(29,"('View.TransformXX = ',f10.6,';')") Transf_XX
        write(29,"('View.TransformYY = ',f10.6,';')") Transf_YY
        write(29,"('View.TransformZZ = ',f10.6,';')") Transf_ZZ
        write(29,"('View.OffsetX = ',f10.6,';')") Offset_X
        write(29,"('View.OffsetZ = ',f10.6,';')") Offset_Z
C         write(29,*) "View.Axes = 2;"
        write(29,*) "View ""Re ",tE_H,"z ",
     *     " "" {"

      tchar=dir_name(1:namelen_dir)// '/' //
     *           gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // 'v_re_sl_D_2.pos'
      open (unit=30,file=tchar)
        write(30,*) "View.IntervalsType = 3;"
        write(30,*) "General.Trackball = 0;"
        write(30,*) "General.RotationX = 90;"
        write(30,*) "General.RotationY = 180;"
        write(30,*) "General.RotationZ = 90;"
        write(30,"('View.TransformXX = ',f10.6,';')") Transf_XX
        write(30,"('View.TransformYY = ',f10.6,';')") Transf_YY
        write(30,"('View.TransformZZ = ',f10.6,';')") Transf_ZZ
        write(30,"('View.OffsetX = ',f10.6,';')") Offset_X
        write(30,"('View.OffsetZ = ',f10.6,';')") Offset_Z
C         write(30,*) "View.Axes = 2;"
        write(30,*) "View ""Re ",tE_H, " "" {"


      dz = h/dble(npt_h-1)
      do i_h=1,npt_h-1
        zz = - (i_h-1)*dz  ! hz=0 => top interface; hz=-h => bottom interface
        do iel=1,nel_D_2
          do inod=1,2
            xel(1,inod) = x_P1_D_2(1,iel+inod-1)
            xel(2,inod) = x_P1_D_2(2,iel+inod-1)
            xel(3,inod) = zz
          enddo
            xel(1,3) = xel(1,2)  ! Quadrangle element
            xel(2,3) = xel(2,2)
            xel(3,3) = zz - dz
            xel(1,4) = xel(1,1)  ! Quadrangle element
            xel(2,4) = xel(2,1)
            xel(3,4) = zz - dz


          do inod=1,2
            sol_el_abs2(inod) = 0.0
            do j=1,3
              z_tmp1 = sol_3d_D_2(j,inod,iel,i_h)
              sol_el(j,inod) = z_tmp1
              sol_el_abs2(inod) = sol_el_abs2(inod) +
     *           abs(z_tmp1)**2
            enddo
            sol_el_abs2(5-inod) = 0.0
            do j=1,3
              z_tmp1 = sol_3d_D_2(j,inod,iel,i_h+1)
              sol_el(j,5-inod) = z_tmp1
              sol_el_abs2(5-inod) = sol_el_abs2(5-inod) +
     *           abs(z_tmp1)**2
            enddo
          enddo
          write(26,10) xel, sol_el_abs2
          write(27,10) xel, (dble(sol_el(1,j)),j=1,n_quad)
          write(28,10) xel, (dble(sol_el(2,j)),j=1,n_quad)
          write(29,10) xel, (dble(sol_el(3,j)),j=1,n_quad)
          write(30,11) xel,
     *     ((dble(sol_el(i,j)),i=1,3),j=1,n_quad)
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
      tchar = "fields_vertically/to_png-"//
     *   gmsh_file_pos(1:namelen_gmsh) // ".geo"
      open (unit=34,file=tchar)

      tchar = "../"//dir_name(1:namelen_dir)// "/" //
     *           gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // "_abs2_sl_X.pos"
      namelen_tchar = len_trim(tchar)
      write(34,*) " Include """, tchar(1:namelen_tchar), """;"
      tchar= gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // "_abs2_sl_X.png"
      namelen_tchar = len_trim(tchar)
      write(34,*) "Print Sprintf(""", tchar(1:namelen_tchar), """);"

      write(34,*)
      write(34,*) "Delete View[0];"
      tchar = "../"//dir_name(1:namelen_dir)// "/" //
     *           gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // "x_re_sl_X.pos"
      namelen_tchar = len_trim(tchar)
      write(34,*) " Include """, tchar(1:namelen_tchar), """;"
      tchar = gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // "x_re_sl_X.png"
      namelen_tchar = len_trim(tchar)
      write(34,*) "Print Sprintf(""", tchar(1:namelen_tchar), """);"

      write(34,*)
      write(34,*) "Delete View[0];"
      tchar = "../"//dir_name(1:namelen_dir)// "/" //
     *           gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // "y_re_sl_X.pos"
      namelen_tchar = len_trim(tchar)
      write(34,*) " Include """, tchar(1:namelen_tchar), """;"
      tchar = gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // "y_re_sl_X.png"
      namelen_tchar = len_trim(tchar)
      write(34,*) "Print Sprintf(""", tchar(1:namelen_tchar), """);"

      write(34,*)
      write(34,*) "Delete View[0];"
      tchar = "../"//dir_name(1:namelen_dir)// "/" //
     *           gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // "z_re_sl_X.pos"
      namelen_tchar = len_trim(tchar)
      write(34,*) " Include """, tchar(1:namelen_tchar), """;"
      tchar = gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // "z_re_sl_X.png"
      namelen_tchar = len_trim(tchar)
      write(34,*) "Print Sprintf(""", tchar(1:namelen_tchar), """);"

      write(34,*)
      write(34,*) "Delete View[0];"
      tchar = "../"//dir_name(1:namelen_dir)// "/" //
     *           gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // "v_re_sl_X.pos"
      namelen_tchar = len_trim(tchar)
      write(34,*) " Include """, tchar(1:namelen_tchar), """;"
      tchar = gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // "v_re_sl_X.png"
      namelen_tchar = len_trim(tchar)
      write(34,*) "Print Sprintf(""", tchar(1:namelen_tchar), """);"
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      write(34,*)
      write(34,*) "Delete View[0];"
      tchar = "../"//dir_name(1:namelen_dir)// "/" //
     *           gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // "_abs2_sl_Y.pos"
      namelen_tchar = len_trim(tchar)
      write(34,*) " Include """, tchar(1:namelen_tchar), """;"
      tchar= gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // "_abs2_sl_Y.png"
      namelen_tchar = len_trim(tchar)
      write(34,*) "Print Sprintf(""", tchar(1:namelen_tchar), """);"

      write(34,*)
      write(34,*) "Delete View[0];"
      tchar = "../"//dir_name(1:namelen_dir)// "/" //
     *           gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // "x_re_sl_Y.pos"
      namelen_tchar = len_trim(tchar)
      write(34,*) " Include """, tchar(1:namelen_tchar), """;"
      tchar = gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // "x_re_sl_Y.png"
      namelen_tchar = len_trim(tchar)
      write(34,*) "Print Sprintf(""", tchar(1:namelen_tchar), """);"

      write(34,*)
      write(34,*) "Delete View[0];"
      tchar = "../"//dir_name(1:namelen_dir)// "/" //
     *           gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // "y_re_sl_Y.pos"
      namelen_tchar = len_trim(tchar)
      write(34,*) " Include """, tchar(1:namelen_tchar), """;"
      tchar = gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // "y_re_sl_Y.png"
      namelen_tchar = len_trim(tchar)
      write(34,*) "Print Sprintf(""", tchar(1:namelen_tchar), """);"

      write(34,*)
      write(34,*) "Delete View[0];"
      tchar = "../"//dir_name(1:namelen_dir)// "/" //
     *           gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // "z_re_sl_Y.pos"
      namelen_tchar = len_trim(tchar)
      write(34,*) " Include """, tchar(1:namelen_tchar), """;"
      tchar = gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // "z_re_sl_Y.png"
      namelen_tchar = len_trim(tchar)
      write(34,*) "Print Sprintf(""", tchar(1:namelen_tchar), """);"

      write(34,*)
      write(34,*) "Delete View[0];"
      tchar = "../"//dir_name(1:namelen_dir)// "/" //
     *           gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // "v_re_sl_Y.pos"
      namelen_tchar = len_trim(tchar)
      write(34,*) " Include """, tchar(1:namelen_tchar), """;"
      tchar = gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // "v_re_sl_Y.png"
      namelen_tchar = len_trim(tchar)
      write(34,*) "Print Sprintf(""", tchar(1:namelen_tchar), """);"
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      write(34,*)
      write(34,*) "Delete View[0];"
      tchar = "../"//dir_name(1:namelen_dir)// "/" //
     *           gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // "_abs2_sl_D_1.pos"
      namelen_tchar = len_trim(tchar)
      write(34,*) " Include """, tchar(1:namelen_tchar), """;"
      tchar= gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // "_abs2_sl_D_1.png"
      namelen_tchar = len_trim(tchar)
      write(34,*) "Print Sprintf(""", tchar(1:namelen_tchar), """);"

      write(34,*)
      write(34,*) "Delete View[0];"
      tchar = "../"//dir_name(1:namelen_dir)// "/" //
     *           gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // "x_re_sl_D_1.pos"
      namelen_tchar = len_trim(tchar)
      write(34,*) " Include """, tchar(1:namelen_tchar), """;"
      tchar = gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // "x_re_sl_D_1.png"
      namelen_tchar = len_trim(tchar)
      write(34,*) "Print Sprintf(""", tchar(1:namelen_tchar), """);"

      write(34,*)
      write(34,*) "Delete View[0];"
      tchar = "../"//dir_name(1:namelen_dir)// "/" //
     *           gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // "y_re_sl_D_1.pos"
      namelen_tchar = len_trim(tchar)
      write(34,*) " Include """, tchar(1:namelen_tchar), """;"
      tchar = gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // "y_re_sl_D_1.png"
      namelen_tchar = len_trim(tchar)
      write(34,*) "Print Sprintf(""", tchar(1:namelen_tchar), """);"

      write(34,*)
      write(34,*) "Delete View[0];"
      tchar = "../"//dir_name(1:namelen_dir)// "/" //
     *           gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // "z_re_sl_D_1.pos"
      namelen_tchar = len_trim(tchar)
      write(34,*) " Include """, tchar(1:namelen_tchar), """;"
      tchar = gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // "z_re_sl_D_1.png"
      namelen_tchar = len_trim(tchar)
      write(34,*) "Print Sprintf(""", tchar(1:namelen_tchar), """);"

      write(34,*)
      write(34,*) "Delete View[0];"
      tchar = "../"//dir_name(1:namelen_dir)// "/" //
     *           gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // "v_re_sl_D_1.pos"
      namelen_tchar = len_trim(tchar)
      write(34,*) " Include """, tchar(1:namelen_tchar), """;"
      tchar = gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // "v_re_sl_D_1.png"
      namelen_tchar = len_trim(tchar)
      write(34,*) "Print Sprintf(""", tchar(1:namelen_tchar), """);"
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      write(34,*)
      write(34,*) "Delete View[0];"
      tchar = "../"//dir_name(1:namelen_dir)// "/" //
     *           gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // "_abs2_sl_D_2.pos"
      namelen_tchar = len_trim(tchar)
      write(34,*) " Include """, tchar(1:namelen_tchar), """;"
      tchar= gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // "_abs2_sl_D_2.png"
      namelen_tchar = len_trim(tchar)
      write(34,*) "Print Sprintf(""", tchar(1:namelen_tchar), """);"

      write(34,*)
      write(34,*) "Delete View[0];"
      tchar = "../"//dir_name(1:namelen_dir)// "/" //
     *           gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // "x_re_sl_D_2.pos"
      namelen_tchar = len_trim(tchar)
      write(34,*) " Include """, tchar(1:namelen_tchar), """;"
      tchar = gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // "x_re_sl_D_2.png"
      namelen_tchar = len_trim(tchar)
      write(34,*) "Print Sprintf(""", tchar(1:namelen_tchar), """);"

      write(34,*)
      write(34,*) "Delete View[0];"
      tchar = "../"//dir_name(1:namelen_dir)// "/" //
     *           gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // "y_re_sl_D_2.pos"
      namelen_tchar = len_trim(tchar)
      write(34,*) " Include """, tchar(1:namelen_tchar), """;"
      tchar = gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // "y_re_sl_D_2.png"
      namelen_tchar = len_trim(tchar)
      write(34,*) "Print Sprintf(""", tchar(1:namelen_tchar), """);"

      write(34,*)
      write(34,*) "Delete View[0];"
      tchar = "../"//dir_name(1:namelen_dir)// "/" //
     *           gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // "z_re_sl_D_2.pos"
      namelen_tchar = len_trim(tchar)
      write(34,*) " Include """, tchar(1:namelen_tchar), """;"
      tchar = gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // "z_re_sl_D_2.png"
      namelen_tchar = len_trim(tchar)
      write(34,*) "Print Sprintf(""", tchar(1:namelen_tchar), """);"

      write(34,*)
      write(34,*) "Delete View[0];"
      tchar = "../"//dir_name(1:namelen_dir)// "/" //
     *           gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // "v_re_sl_D_2.pos"
      namelen_tchar = len_trim(tchar)
      write(34,*) " Include """, tchar(1:namelen_tchar), """;"
      tchar = gmsh_file_pos(1:namelen_gmsh)
     *           // tE_H // "v_re_sl_D_2.png"
      namelen_tchar = len_trim(tchar)
      write(34,*) "Print Sprintf(""", tchar(1:namelen_tchar), """);"


      close (unit=34)

c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
cc      open (unit=63, file="Output/sol_3d_X.txt",
cc     *         status="unknown")
cc        do i_h=1,npt_h
cc          do iel=1,nel_X
cc            do inod=1,2
cc            write(63,*) i_h, iel, inod,
cc     *                 (sol_3d_X(j,inod,iel,i_h),j=1,3)
cc            enddo
cc          enddo
cc        enddo
cc      close(63)
cc
cc      open (unit=64, file="Output/nb_visit_D_2.txt",
cc     *         status="unknown")
cc        do iel=1,nel_X
cc          write(64,*) iel, (nb_visit_X(inod,iel),inod=1,2)
cc        enddo
cc      close(64)
cc
cc      open (unit=63, file="Output/info_elem.txt",
cc     *         status="unknown")
cc        do iel=1,nel
cc          write(63,*) iel, info_elem(iel)
cc        enddo
cc      close(63)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     SQ : Scalar quadrangle
10    format("SQ(",f10.6,11(",",f10.6),"){",
     *     g24.16,3(",",g24.16),"};")

c     VQ : Vector quadrangle
11    format("VQ(",f10.6,11(",",f10.6),"){",
     *     g24.16,11(",",g24.16),"};")

cc

      return
      end



