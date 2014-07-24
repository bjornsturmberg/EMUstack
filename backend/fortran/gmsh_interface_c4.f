c*******************************************************
c
c     gmsh_interface_c4: covert the GMSH mesh format to the FEM mesh format
c
c*******************************************************
c
c     nnodes : Number of nodes per element (10-node second order tetrahedron)
c
c*******************************************************
c
      subroutine gmsh_interface_c4 (nel, npt, nnodes, type_el, 
     *  table_nod, x, lat_vecs)

c
      implicit none
      integer*8 nel, npt, nnodes
      integer*8 type_el(nel)
      integer*8 table_nod(nnodes,nel)
      double precision lat_vecs(2,2)
      double precision x(2,npt)
c
      double precision x_min, y_min, x_r
      double precision x_max, y_max, y_r
c
      integer*8 i, j, k, j1, iel, type_cyl
      integer*8 ui, debug

      double precision centre(2,4), rad_cyl(4), r_tmp1

      double precision xx(2)
      double precision period, xx_0(2), ss1, ss2
      double precision xyz(3,5)
c
ccccccccccccccccccccccccccccccccccccc
c
      ui = 6
      debug = 1
      type_cyl = 4 ! Cylinder type
c
      period = lat_vecs(1,1)
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
c     Determination of the cylinder centre (center of the unit cell)
      xx_0(1) = (x_max + x_min)/2.0d0
      xx_0(2) = (y_max + y_min)/2.0d0
      do i=1,2
        if (i .eq. 1) then
          ss1  = -1.0d0  ! Sign
        else
          ss1  = 1.0d0
        endif
        do j=1,2
          if (j .eq. 1) then
            ss2  = -1.0d0  ! Sign
          else
            ss2  = 1.0d0
          endif
          k = i + 2*(j-1)  ! count the centres
          centre(1,k) = xx_0(1) + ss1*period/4.0d0
          centre(2,k) = xx_0(2) + ss2*period/4.0d0
        enddo
      enddo
c
cccccc
c     Determination of the cylinder radius
      do k=1,4
        rad_cyl(k) = 0.0d0
      enddo
      do iel=1,nel
        if (type_el(iel) .eq. type_cyl) then
          do j=1,nnodes
            j1 = table_nod(j,iel)
            do k=1,4
              xx(1) = x(1,j1) - centre(1,k)
              xx(2) = x(2,j1) - centre(2,k)
              r_tmp1 = xx(1)**2 + xx(2)**2
              if (r_tmp1 .gt. rad_cyl(k) .and. 
     *            r_tmp1 .lt. (period/4.0d0)**2) then
                rad_cyl(k) = r_tmp1
              endif
            enddo
          enddo
        endif
      enddo
      do k=1,4
        rad_cyl(k) = sqrt(rad_cyl(k)) ! Cylinder radius
      enddo

c
ccccccccccccccccccccccccccccccccccccc
c

      do k=i,5
        xyz(3,i) = 0.0d0  ! Z-coordinate = 0
      enddo
      open (unit=26,file="Bloch_fields/interface_c4.geo")
          write(26,*) "lc = 0.101;"
          write(26,*) "x_max = ", x_max, ";"
          write(26,*) "y_max = ",  y_max, ";"
          write(26,*) "x_min = ",  x_min, ";"
          write(26,*) "y_min = ",  y_min, ";"
          write(26,*) "zz = 0;"

          k = 1
          xyz(1,1) = centre(1,k)
          xyz(2,1) = centre(2,k)
          xyz(1,2) = centre(1,k) + rad_cyl(k)
          xyz(2,2) = centre(2,k)
          xyz(1,3) = centre(1,k)
          xyz(2,3) = centre(2,k) + rad_cyl(k)
          xyz(1,4) = centre(1,k) - rad_cyl(k)
          xyz(2,4) = centre(2,k)
          xyz(1,5) = centre(1,k)
          xyz(2,5) = centre(2,k) - rad_cyl(k)
          write(26,10) 1, (xyz(i,1),i=1,3)
          write(26,10) 2, (xyz(i,2),i=1,3)
          write(26,10) 3, (xyz(i,3),i=1,3)
          write(26,10) 4, (xyz(i,4),i=1,3)
          write(26,10) 5, (xyz(i,5),i=1,3)

          k = 2
          xyz(1,1) = centre(1,k)
          xyz(2,1) = centre(2,k)
          xyz(1,2) = centre(1,k) + rad_cyl(k)
          xyz(2,2) = centre(2,k)
          xyz(1,3) = centre(1,k)
          xyz(2,3) = centre(2,k) + rad_cyl(k)
          xyz(1,4) = centre(1,k) - rad_cyl(k)
          xyz(2,4) = centre(2,k)
          xyz(1,5) = centre(1,k)
          xyz(2,5) = centre(2,k) - rad_cyl(k)
          write(26,10) 6, (xyz(i,1),i=1,3)
          write(26,10) 7, (xyz(i,2),i=1,3)
          write(26,10) 8, (xyz(i,3),i=1,3)
          write(26,10) 9, (xyz(i,4),i=1,3)
          write(26,10) 10, (xyz(i,5),i=1,3)

          k = 3
          xyz(1,1) = centre(1,k)
          xyz(2,1) = centre(2,k)
          xyz(1,2) = centre(1,k) + rad_cyl(k)
          xyz(2,2) = centre(2,k)
          xyz(1,3) = centre(1,k)
          xyz(2,3) = centre(2,k) + rad_cyl(k)
          xyz(1,4) = centre(1,k) - rad_cyl(k)
          xyz(2,4) = centre(2,k)
          xyz(1,5) = centre(1,k)
          xyz(2,5) = centre(2,k) - rad_cyl(k)
          write(26,10) 11, (xyz(i,1),i=1,3)
          write(26,10) 12, (xyz(i,2),i=1,3)
          write(26,10) 13, (xyz(i,3),i=1,3)
          write(26,10) 14, (xyz(i,4),i=1,3)
          write(26,10) 15, (xyz(i,5),i=1,3)

          k = 4
          xyz(1,1) = centre(1,k)
          xyz(2,1) = centre(2,k)
          xyz(1,2) = centre(1,k) + rad_cyl(k)
          xyz(2,2) = centre(2,k)
          xyz(1,3) = centre(1,k)
          xyz(2,3) = centre(2,k) + rad_cyl(k)
          xyz(1,4) = centre(1,k) - rad_cyl(k)
          xyz(2,4) = centre(2,k)
          xyz(1,5) = centre(1,k)
          xyz(2,5) = centre(2,k) - rad_cyl(k)
          write(26,10) 16, (xyz(i,1),i=1,3)
          write(26,10) 17, (xyz(i,2),i=1,3)
          write(26,10) 18, (xyz(i,3),i=1,3)
          write(26,10) 19, (xyz(i,4),i=1,3)
          write(26,10) 20, (xyz(i,5),i=1,3)

          write(26,*) "Point(21) = {x_min, y_min, zz, lc};"
          write(26,*) "Point(22) = {x_max, y_max, lc};"

          write(26,*) "Circle(1) = {2, 1, 3};"
          write(26,*) "Circle(2) = {3, 1, 4};"
          write(26,*) "Circle(3) = {4, 1, 5};"
          write(26,*) "Circle(4) = {5, 1, 2};"

          write(26,*) "Circle(5) = {7, 6, 8};"
          write(26,*) "Circle(6) = {8, 6, 9};"
          write(26,*) "Circle(7) = {9, 6, 10};"
          write(26,*) "Circle(8) = {10, 6, 7};"

          write(26,*) "Circle(9) = {12, 11, 13};"
          write(26,*) "Circle(10) = {13, 11, 14};"
          write(26,*) "Circle(11) = {14, 11, 15};"
          write(26,*) "Circle(12) = {15, 11, 12};"

          write(26,*) "Circle(13) = {17, 16, 18};"
          write(26,*) "Circle(14) = {18, 16, 19};"
          write(26,*) "Circle(15) = {19, 16, 20};"
          write(26,*) "Circle(16) = {20, 16, 17};"

          write(26,*) "Geometry.LineWidth = 2;"
          write(26,*) "Geometry.Color.Lines = White;"
          write(26,*) "Geometry.Points = 0;"
      close(26)

10    format("Point(",I2,") = {",f10.6,", "f10.6,", "f10.6,", lc};")

      if (debug .eq. 1) then
        write(*,*)
      do k=1,4
        write(*,*) " gmsh_interface_c4: k, rad_cyl(k) = ", 
     *           k, rad_cyl(k)
      enddo
        write(*,*)
      do k=1,4
        write(*,*) " gmsh_interface_c4: k, centre =", 
     *           k, (centre(i,k),i=1,2)
      enddo
        write(*,*)
        write(*,*) " gmsh_interface_c4: period = ", 
     *           period, period/4.0d0
      endif

c
ccccccccccccccccccccccccccccccccccccc
c
      return
      end
