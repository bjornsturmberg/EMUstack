c
c***********************************************************************
c
      subroutine lattice_vec (npt, x, lat_vecs, debug)
c
c***********************************************************************
c
      implicit none
      integer*8 npt
      double precision x(2,npt)
      double precision lat_vecs(2,2)

      double precision x_min, y_min
      double precision x_max, y_max

      double precision x_min_ymin
      double precision x_max_ymin
      double precision x_min_ymax
      double precision x_max_ymax
      double precision x_r, y_r
      double precision lat_vec_1(2), lat_vec_2(2)
      double precision tol
      integer*8 i, n1, n2, debug
c
      x_min = x(1,1)
      x_max = x(1,1)
      do i=1,npt
        x_r = x(1,i)
        if(x_r .lt. x_min) x_min = x_r
        if(x_r .gt. x_max) x_max = x_r
      enddo
      y_min = x(2,1)
      y_max = x(2,1)
      do i=1,npt
        y_r = x(2,i)
        if(y_r .lt. y_min) y_min = y_r
        if(y_r .gt. y_max) y_max = y_r
      enddo
c
      if (debug .eq. 1) then
        write(*,*) "lattice_vec: x_min, x_max = ", x_min, x_max
        write(*,*) "lattice_vec: y_min, y_max = ", y_min, y_max
      endif
c
      tol = 1.0d-6
      n1 = 0
      do i=1,npt
        if (abs(x(2,i)-y_min) .lt. tol) then
          if (n1 .eq. 0) then
            x_min_ymin = x(1,i)
            x_max_ymin = x(1,i)
          endif
          n1 = n1 + 1
          x_r = x(1,i)
          if(x_r .lt. x_min_ymin) x_min_ymin = x_r
          if(x_r .gt. x_max_ymin) x_max_ymin = x_r
        endif
      enddo
c
      n2 = 0
      do i=1,npt
        if (abs(x(2,i)-y_max) .lt. tol) then
          if (n2 .eq. 0) then
            x_min_ymax = x(1,i)
            x_max_ymax = x(1,i)
          endif
          n2 = n2 + 1
          x_r = x(1,i)
          if(x_r .lt. x_min_ymax) x_min_ymax = x_r
          if(x_r .gt. x_max_ymax) x_max_ymax = x_r
        endif
      enddo
c
      lat_vec_1(1) = x_max_ymin - x_min_ymin
      lat_vec_1(2) = 0.0d0
c
      lat_vec_2(1) = x_min_ymax - x_min_ymin
      lat_vec_2(2) = y_max - y_min
c
      do i=1,2
        lat_vecs(i,1) = lat_vec_1(i)
        lat_vecs(i,2) = lat_vec_2(i)
      enddo
c
      if (debug .eq. 1) then
        write(*,*) "lattice_vec: v1 = ", (lat_vecs(i,1),i=1,2)
        write(*,*) "lattice_vec: v2 = ", (lat_vecs(i,2),i=1,2)
      endif
c
      return
      end
