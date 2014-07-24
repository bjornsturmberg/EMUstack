
c***********************************************************************
c
c
      subroutine slice_interp (nel, nval, iel, ival, i_h,
     *     nnodes, nx, ny, nel_X, nel_Y, nel_D_1, nel_D_2,
     *     npt_h, ls_edge, xel_2d, evecs, coef_t, coef_z,
     *     sol_3d_X, sol_3d_Y, sol_3d_D_1, sol_3d_D_2,
     *     nb_visit_X, nb_visit_Y, nb_visit_D_1, nb_visit_D_2,
     *     info_elem)

c
      implicit none
      integer*8 iel, nel, ival, nval, i_h
      integer*8 nnodes, nx, ny
      integer*8 nel_X, nel_Y, nel_D_1, nel_D_2, npt_h
      integer*8 info_elem(nel)
      double precision xel_2d(2,6), ls_edge(4)

      double precision x_min, x_max, y_min, y_max

      complex*16 evecs(3,nnodes+7,nval,nel)
      complex*16 coef_t, coef_z

      complex*16 sol_3d_X(3,2,nel_X,npt_h)
      complex*16 sol_3d_Y(3,2,nel_Y,npt_h)
      complex*16 sol_3d_D_1(3,2,nel_D_1,npt_h)
      complex*16 sol_3d_D_2(3,2,nel_D_2,npt_h)

      integer*8 nb_visit_X(2,nel_X), nb_visit_Y(2,nel_Y)
      integer*8 nb_visit_D_1(2,nel_D_1), nb_visit_D_2(2,nel_D_2)

c     Local variables
      integer*8 i, j, nb_intersect
      integer*8 inod, iel_2, inod_2, iel_3, inod_3
      double precision b(2,2), binv(2,2), detb
      double precision x_min2, x_max2, y_min2, y_max2, tol
      integer*8 nx_min2, nx_max2, ny_min2, ny_max2
      double precision x_0, y_0, xx, yy, g1(3), g2(6)
      double precision dx_1, dy_1, dx_2, dy_2
      double precision dx_3, dy_3, dx_4, dy_4
      complex*16 z_tmp1, z_sol(3)

c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c         Contruction of FEM transformation for a linear element
c         construction of b
c         -----------------
          b(1,1) = xel_2d(1,2) - xel_2d(1,1)
          b(1,2) = xel_2d(1,3) - xel_2d(1,1)
          b(2,1) = xel_2d(2,2) - xel_2d(2,1)
          b(2,2) = xel_2d(2,3) - xel_2d(2,1)
cccccccccc
c         inverse of b
c         --------------
          detb = b(1,1) * b(2,2) - b(1,2) * b(2,1)
          binv(1,1) = b(2,2) / detb
          binv(2,2) = b(1,1) / detb
          binv(1,2) = b(1,2) / (-detb)
          binv(2,1) = b(2,1) / (-detb)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      x_min = ls_edge(1)
      x_max = ls_edge(2)
      y_min = ls_edge(3)
      y_max = ls_edge(4)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      tol = 1.0d-5
      x_min2 = xel_2d(1,1)
      x_max2 = xel_2d(1,1)
      do j=1,3
        if(xel_2d(1,j) .lt. x_min2) x_min2 = xel_2d(1,j)
        if(xel_2d(1,j) .gt. x_max2) x_max2 = xel_2d(1,j)
      enddo
      y_min2 = xel_2d(2,1)
      y_max2 = xel_2d(2,1)
      do j=1,3
        if(xel_2d(2,j) .lt. y_min2) y_min2 = xel_2d(2,j)
        if(xel_2d(2,j) .gt. y_max2) y_max2 = xel_2d(2,j)
      enddo
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      nb_intersect = 0
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Check for points of the X-slice which belong to the triangle

      dx_1 = (x_max - x_min)/dble(nel_X)
      dy_1 = 0
      y_0 = (y_min+y_max) / 2.0d0
      if (y_0 > y_min2-tol .and. y_0 < y_max2+tol) then
        nx_min2 = anint((x_min2-x_min)/dx_1) + 1
        nx_max2 = anint((x_max2-x_min)/dx_1) + 1
        if (nx_max2 < 0 .or. nx_max2 > nel_X + 1) then
          write(*,*)
          write(*,*) "slice_interp: problem with nx_max2 = ",
     *        nx_max2
          write(*,*) "slice_interp: nx_max2 > = nel_X + 1 :",
     *        nel_X + 1
          write(*,*) "dx_1, nel_X = ", dx_1, nel_X
          write(*,*) "x_min,  x_max  = ", x_min, x_max
          write(*,*) "x_min2, x_max2 = ", x_min2, x_max2
          write(*,*) "slice_interp: Aborting..."
          write(*,*)
          stop
        endif
        do i=nx_min2,nx_max2
          xx = x_min + (i-1)*dx_1
          yy = (y_min+y_max) / 2.0d0

          g1(1) = 1.0d0 - (binv(1,1)+binv(2,1))*(xx - xel_2d(1,1))
     *        - (binv(1,2)+binv(2,2))*(yy - xel_2d(2,1))
          g1(2) = binv(1,1)*(xx - xel_2d(1,1))
     *        + binv(1,2)*(yy - xel_2d(2,1))
          g1(3) = binv(2,1)*(xx - xel_2d(1,1))
     *        + binv(2,2)*(yy - xel_2d(2,1))
          if(abs(g1(1) - 0.5) .le. (0.5 + tol) .and.
     *      abs(g1(2) - 0.5) .le. (0.5 + tol) .and.
     *      abs(g1(3) - 0.5) .le. (0.5 + tol) ) then
c           Quadratique interpolation
            g2(1) = 2.0d0*g1(1)*(g1(1)-0.5d0)
            g2(2) = 2.0d0*g1(2)*(g1(2)-0.5d0)
            g2(3) = 2.0d0*g1(3)*(g1(3)-0.5d0)
            g2(4) = 4.d0*g1(1)*g1(2)
            g2(5) = 4.d0*g1(2)*g1(3)
            g2(6) = 4.d0*g1(3)*g1(1)

            if (i > 0 .and. i <= nel_X) then
              inod_2 = 1
              iel_2 = i
              if (i_h == 1 .and. ival == 1) then
                nb_intersect = nb_intersect + 1
                nb_visit_X(inod_2,iel_2) = nb_visit_X(inod_2,iel_2)
     *                     + 1
              endif
            elseif (i == nel_X+1) then
              inod_2 = 2
              iel_2 = i-1
              if (i_h == 1 .and. ival == 1) then
                nb_intersect = nb_intersect + 1
                nb_visit_X(inod_2,iel_2) = nb_visit_X(inod_2,iel_2)
     *                     + 1
              endif
            else
              write(*,*)
              write(*,*) "slice_interp: X-slice"
              write(*,*) "slice_interp: problem with i = ", i
              write(*,*) "slice_interp: iel, ival, i_h = ",
     *                                  iel, ival, i_h
              write(*,*) "x_min,  x_max  = ", x_min, x_max
              write(*,*) "x_min2, x_max2 = ", x_min2, x_max2
              write(*,*) "slice_interp: Aborting..."
              write(*,*)
              stop
            endif
            do j=1,3
              z_sol(j) = 0
            enddo
            do inod=1,6
              do j=1,2
                z_tmp1 = g2(inod) * evecs(j,inod,ival,iel) * coef_t
                z_sol(j) = z_sol(j) + z_tmp1
              enddo
              j=3
                z_tmp1 = g2(inod) * evecs(j,inod,ival,iel) * coef_z
                z_sol(j) = z_sol(j) + z_tmp1
            enddo
            do j=1,3
              sol_3d_X(j,inod_2,iel_2,i_h) =
     *                    sol_3d_X(j,inod_2,iel_2,i_h) + z_sol(j)
            enddo
            if (i >= 2 .and. i <= nel_X) then  ! node which belong to two sub-intervals
              inod_2 = 1
              iel_2 = i
              inod_3 = 2
              iel_3 = i-1
              if (i_h == 1 .and. ival == 1) then
                nb_intersect = nb_intersect + 1
                nb_visit_X(inod_3,iel_3) = nb_visit_X(inod_3,iel_3)
     *                     + 1
              endif
              do j=1,3
                sol_3d_X(j,inod_3,iel_3,i_h) =
     *                    sol_3d_X(j,inod_3,iel_3,i_h) + z_sol(j)
              enddo
            endif
          endif
        enddo
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Check for points of the Y-slice which belong to the triangle
      dx_2 = 0
      dy_2 = (y_max - y_min)/dble(nel_Y)
      x_0 = (x_min+x_max) / 2.0d0
      if (x_0 > x_min2-tol .and. x_0 < x_max2+tol) then
        ny_min2 = anint((y_min2-y_min)/dy_2) + 1
        ny_max2 = anint((y_max2-y_min)/dy_2) + 1
        if (ny_max2 < 0 .or. ny_max2 > nel_Y + 1) then
          write(*,*)
          write(*,*) "slice_interp: problem with ny_max2 = ",
     *        ny_max2
          write(*,*) "slice_interp: ny_max2 > = nel_Y + 1 :",
     *        nel_Y + 1
          write(*,*) "dy_2, nel_Y = ", dy_2, nel_Y
          write(*,*) "y_min,  y_max  = ", y_min, y_max
          write(*,*) "y_min2, y_max2 = ", y_min2, y_max2
          write(*,*) "slice_interp: Aborting..."
          write(*,*)
          stop
        endif
        do i=ny_min2,ny_max2
          xx = (x_min+x_max) / 2.0d0
          yy = y_min + (i-1)*dy_2

          g1(1) = 1.0d0 - (binv(1,1)+binv(2,1))*(xx - xel_2d(1,1))
     *        - (binv(1,2)+binv(2,2))*(yy - xel_2d(2,1))
          g1(2) = binv(1,1)*(xx - xel_2d(1,1))
     *        + binv(1,2)*(yy - xel_2d(2,1))
          g1(3) = binv(2,1)*(xx - xel_2d(1,1))
     *        + binv(2,2)*(yy - xel_2d(2,1))
          if(abs(g1(1) - 0.5) .le. (0.5 + tol) .and.
     *      abs(g1(2) - 0.5) .le. (0.5 + tol) .and.
     *      abs(g1(3) - 0.5) .le. (0.5 + tol) ) then
c           Quadratique interpolation
            g2(1) = 2.0d0*g1(1)*(g1(1)-0.5d0)
            g2(2) = 2.0d0*g1(2)*(g1(2)-0.5d0)
            g2(3) = 2.0d0*g1(3)*(g1(3)-0.5d0)
            g2(4) = 4.d0*g1(1)*g1(2)
            g2(5) = 4.d0*g1(2)*g1(3)
            g2(6) = 4.d0*g1(3)*g1(1)

            if (i > 0 .and. i <= nel_Y) then
              inod_2 = 1
              iel_2 = i
              if (i_h == 1 .and. ival == 1) then
                nb_intersect = nb_intersect + 1
                nb_visit_Y(inod_2,iel_2) = nb_visit_Y(inod_2,iel_2)
     *                     + 1
              endif
            elseif (i == nel_Y+1) then
              inod_2 = 2
              iel_2 = i-1
              if (i_h == 1 .and. ival == 1) then
                nb_intersect = nb_intersect + 1
                nb_visit_Y(inod_2,iel_2) = nb_visit_Y(inod_2,iel_2)
     *                     + 1
              endif
            else
              write(*,*)
              write(*,*) "slice_interp: Y-slice"
              write(*,*) "slice_interp: problem with i = ", i
              write(*,*) "slice_interp: iel, ival, i_h = ",
     *                                  iel, ival, i_h
              write(*,*) "y_min,  y_max  = ", y_min, y_max
              write(*,*) "y_min2, y_max2 = ", y_min2, y_max2
              write(*,*) "slice_interp: Aborting..."
              write(*,*)
              stop
            endif
            do j=1,3
              z_sol(j) = 0
            enddo
            do inod=1,6
              do j=1,2
                z_tmp1 = g2(inod) * evecs(j,inod,ival,iel) * coef_t
                z_sol(j) = z_sol(j) + z_tmp1
              enddo
              j=3
                z_tmp1 = g2(inod) * evecs(j,inod,ival,iel) * coef_z
                z_sol(j) = z_sol(j) + z_tmp1
            enddo
            do j=1,3
              sol_3d_Y(j,inod_2,iel_2,i_h) =
     *                    sol_3d_Y(j,inod_2,iel_2,i_h) + z_sol(j)
            enddo
            if (i >= 2 .and. i <= nel_Y) then  ! node which belong to two sub-intervals
              inod_2 = 1
              iel_2 = i
              inod_3 = 2
              iel_3 = i-1
              if (i_h == 1 .and. ival == 1) then
                nb_intersect = nb_intersect + 1
                nb_visit_Y(inod_3,iel_3) = nb_visit_Y(inod_3,iel_3)
     *                     + 1
              endif
              do j=1,3
                sol_3d_Y(j,inod_3,iel_3,i_h) =
     *                    sol_3d_Y(j,inod_3,iel_3,i_h) + z_sol(j)
              enddo
            endif
          endif
        enddo
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Check for points of the D1-slice which belong to the triangle

      dx_3 = (x_max - x_min)/dble(nel_D_1)
      dy_3 = (y_max - y_min)/dble(nel_D_1)
      if (.true.) then
        nx_min2 = anint((x_min2-x_min)/dx_3) + 1
        nx_max2 = anint((x_max2-x_min)/dx_3) + 1
        do i=nx_min2,nx_max2
          xx = x_min + (i-1)*dx_3
          yy = y_min + (i-1)*dy_3

          g1(1) = 1.0d0 - (binv(1,1)+binv(2,1))*(xx - xel_2d(1,1))
     *        - (binv(1,2)+binv(2,2))*(yy - xel_2d(2,1))
          g1(2) = binv(1,1)*(xx - xel_2d(1,1))
     *        + binv(1,2)*(yy - xel_2d(2,1))
          g1(3) = binv(2,1)*(xx - xel_2d(1,1))
     *        + binv(2,2)*(yy - xel_2d(2,1))
          if(abs(g1(1) - 0.5) .le. (0.5 + tol) .and.
     *      abs(g1(2) - 0.5) .le. (0.5 + tol) .and.
     *      abs(g1(3) - 0.5) .le. (0.5 + tol) ) then
c           Quadratique interpolation
            g2(1) = 2.0d0*g1(1)*(g1(1)-0.5d0)
            g2(2) = 2.0d0*g1(2)*(g1(2)-0.5d0)
            g2(3) = 2.0d0*g1(3)*(g1(3)-0.5d0)
            g2(4) = 4.d0*g1(1)*g1(2)
            g2(5) = 4.d0*g1(2)*g1(3)
            g2(6) = 4.d0*g1(3)*g1(1)

            if (i > 0 .and. i <= nel_D_1) then
              inod_2 = 1
              iel_2 = i
              if (i_h == 1 .and. ival == 1) then
                nb_intersect = nb_intersect + 1
                nb_visit_D_1(inod_2,iel_2) = nb_visit_D_1(inod_2,iel_2)
     *                     + 1
              endif
            elseif (i == nel_D_1+1) then
              inod_2 = 2
              iel_2 = i-1
              if (i_h == 1 .and. ival == 1) then
                nb_intersect = nb_intersect + 1
                nb_visit_D_1(inod_2,iel_2) = nb_visit_D_1(inod_2,iel_2)
     *                     + 1
              endif
            else
              write(*,*)
              write(*,*) "slice_interp: D1-slice"
              write(*,*) "slice_interp: problem with i = ", i
              write(*,*) "slice_interp: iel, ival, i_h = ",
     *                                  iel, ival, i_h
              write(*,*) "x_min,  x_max  = ", x_min, x_max
              write(*,*) "x_min2, x_max2 = ", x_min2, x_max2
              write(*,*) "slice_interp: Aborting..."
              write(*,*)
              stop
            endif
            do j=1,3
              z_sol(j) = 0
            enddo
            do inod=1,6
              do j=1,2
                z_tmp1 = g2(inod) * evecs(j,inod,ival,iel) * coef_t
                z_sol(j) = z_sol(j) + z_tmp1
              enddo
              j=3
                z_tmp1 = g2(inod) * evecs(j,inod,ival,iel) * coef_z
                z_sol(j) = z_sol(j) + z_tmp1
            enddo
            do j=1,3
              sol_3d_D_1(j,inod_2,iel_2,i_h) =
     *                    sol_3d_D_1(j,inod_2,iel_2,i_h) + z_sol(j)
            enddo
            if (i >= 2 .and. i <= nel_D_1) then  ! node which belong to two sub-intervals
              inod_2 = 1
              iel_2 = i
              inod_3 = 2
              iel_3 = i-1
              if (i_h == 1 .and. ival == 1) then
                nb_intersect = nb_intersect + 1
                nb_visit_D_1(inod_3,iel_3) = nb_visit_D_1(inod_3,iel_3)
     *                     + 1
              endif
              do j=1,3
                sol_3d_D_1(j,inod_3,iel_3,i_h) =
     *                    sol_3d_D_1(j,inod_3,iel_3,i_h) + z_sol(j)
              enddo
            endif
          endif
        enddo
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Check for points of the D2-slice which belong to the triangle

      dx_4 =  (x_max - x_min)/dble(nel_D_2)
      dy_4 = -(y_max - y_min)/dble(nel_D_2)
      if (.true.) then
        nx_min2 = anint((x_min2-x_min)/dx_4) + 1
        nx_max2 = anint((x_max2-x_min)/dx_4) + 1
        do i=nx_min2,nx_max2
          xx = x_min + (i-1)*dx_4
          yy = y_max + (i-1)*dy_4

          g1(1) = 1.0d0 - (binv(1,1)+binv(2,1))*(xx - xel_2d(1,1))
     *        - (binv(1,2)+binv(2,2))*(yy - xel_2d(2,1))
          g1(2) = binv(1,1)*(xx - xel_2d(1,1))
     *        + binv(1,2)*(yy - xel_2d(2,1))
          g1(3) = binv(2,1)*(xx - xel_2d(1,1))
     *        + binv(2,2)*(yy - xel_2d(2,1))
          if(abs(g1(1) - 0.5) .le. (0.5 + tol) .and.
     *      abs(g1(2) - 0.5) .le. (0.5 + tol) .and.
     *      abs(g1(3) - 0.5) .le. (0.5 + tol) ) then
c           Quadratique interpolation
            g2(1) = 2.0d0*g1(1)*(g1(1)-0.5d0)
            g2(2) = 2.0d0*g1(2)*(g1(2)-0.5d0)
            g2(3) = 2.0d0*g1(3)*(g1(3)-0.5d0)
            g2(4) = 4.d0*g1(1)*g1(2)
            g2(5) = 4.d0*g1(2)*g1(3)
            g2(6) = 4.d0*g1(3)*g1(1)

            if (i > 0 .and. i <= nel_D_2) then
              inod_2 = 1
              iel_2 = i
              if (i_h == 1 .and. ival == 1) then
                nb_intersect = nb_intersect + 1
                nb_visit_D_2(inod_2,iel_2) = nb_visit_D_2(inod_2,iel_2)
     *                     + 1
              endif
            elseif (i == nel_D_2+1) then
              inod_2 = 2
              iel_2 = i-1
              if (i_h == 1 .and. ival == 1) then
                nb_intersect = nb_intersect + 1
                nb_visit_D_2(inod_2,iel_2) = nb_visit_D_2(inod_2,iel_2)
     *                     + 1
              endif
            else
              write(*,*)
              write(*,*) "slice_interp: D2-slice"
              write(*,*) "slice_interp: problem with i = ", i
              write(*,*) "slice_interp: iel, ival, i_h = ",
     *                                  iel, ival, i_h
              write(*,*) "x_min,  x_max  = ", x_min, x_max
              write(*,*) "x_min2, x_max2 = ", x_min2, x_max2
              write(*,*) "slice_interp: Aborting..."
              write(*,*)
              stop
            endif
            do j=1,3
              z_sol(j) = 0
            enddo
            do inod=1,6
              do j=1,2
                z_tmp1 = g2(inod) * evecs(j,inod,ival,iel) * coef_t
                z_sol(j) = z_sol(j) + z_tmp1
              enddo
              j=3
                z_tmp1 = g2(inod) * evecs(j,inod,ival,iel) * coef_z
                z_sol(j) = z_sol(j) + z_tmp1
            enddo
            do j=1,3
              sol_3d_D_2(j,inod_2,iel_2,i_h) =
     *                    sol_3d_D_2(j,inod_2,iel_2,i_h) + z_sol(j)
            enddo
            if (i >= 2 .and. i <= nel_D_2) then  ! node which belong to two sub-intervals
              inod_2 = 1
              iel_2 = i
              inod_3 = 2
              iel_3 = i-1
              if (i_h == 1 .and. ival == 1) then
                nb_intersect = nb_intersect + 1
                nb_visit_D_2(inod_3,iel_3) = nb_visit_D_2(inod_3,iel_3)
     *                     + 1
              endif
              do j=1,3
                sol_3d_D_2(j,inod_3,iel_3,i_h) =
     *                    sol_3d_D_2(j,inod_3,iel_3,i_h) + z_sol(j)
              enddo
            endif
          endif
        enddo
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (i_h == 1 .and. ival == 1) then
        if (nb_intersect >= 1 ) then
          info_elem(iel) = nb_intersect  ! the element contains nb_intersect slice points
        else
          info_elem(iel) = -1 !   the element contains no nb_intersect slice points and will no longer be processed
        endif
      endif
c
      return
      end
