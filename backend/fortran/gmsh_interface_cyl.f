c*******************************************************
c
c     gmsh_interface_cyl: covert the GMSH mesh format to the FEM mesh format
c
c*******************************************************
c
c     nnodes : Number of nodes per element (10-node second order tetrahedron)
c
c*******************************************************
c
      subroutine gmsh_interface_cyl (nel, npt, nnodes, type_el, 
     *  table_nod, x)

c
      implicit none
      integer*8 nel, npt, nnodes
      integer*8 type_el(nel)
      integer*8 table_nod(nnodes,nel)
      double precision x(2,npt)
c
      double precision x_min, y_min, x_r
      double precision x_max, y_max, y_r
c
      integer*8 i, j, j1, iel, type_cyl
      integer*8 ui, debug
      double precision centre(2), rad_cyl, r_tmp1
      double precision xx(2)
c
ccccccccccccccccccccccccccccccccccccc
c
      ui = 6
      debug = 1
      type_cyl = 4 ! Cylinder type
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
      centre(1) = (x_max + x_min)/2.0d0
      centre(2) = (y_max + y_min)/2.0d0
c
cccccc
c     Determination of the cylinder radius
      rad_cyl = 0.0d0
      do iel=1,nel
        if (type_el(iel) .eq. type_cyl) then
        do j=1,nnodes
          j1 = table_nod(j,iel)
          xx(1) = x(1,j1) - centre(1)
          xx(2) = x(2,j1) - centre(2)
          r_tmp1 = xx(1)**2 + xx(2)**2
          if (r_tmp1 .gt. rad_cyl) rad_cyl = r_tmp1
        enddo
        endif
      enddo
      rad_cyl = sqrt(rad_cyl) ! Cylinder radius
c
ccccccccccccccccccccccccccccccccccccc
c
      open (unit=26,file="Bloch_fields/interface_cyl.geo")
          write(26,*) "lc = 0.101;"
          write(26,*) "radius = ",  rad_cyl, ";"
          write(26,*) "x_cent = ",  centre(1), ";"
          write(26,*) "y_cent = ",  centre(2), ";"
          write(26,*) "x_max = ",  x_max, ";"
          write(26,*) "y_max = ",  y_max, ";"
          write(26,*) "x_min = ",  x_min, ";"
          write(26,*) "y_min = ",  y_min, ";"
          write(26,*) "zz = 0;"
          write(26,*) "Point(1) = {x_cent, y_cent, zz, lc};"
          write(26,*) "Point(2) = {x_cent + radius, y_cent, zz, lc};"
          write(26,*) "Point(3) = {x_cent, y_cent + radius, zz, lc};"
          write(26,*) "Point(4) = {x_cent - radius, y_cent, zz, lc};"
          write(26,*) "Point(5) = {x_cent, y_cent - radius, zz, lc};"
          write(26,*) "Point(6) = {x_min, y_min, zz, lc};"
          write(26,*) "Point(7) = {x_max, y_max, lc};"
          write(26,*) "Circle(1) = {2, 1, 3};"
          write(26,*) "Circle(2) = {3, 1, 4};"
          write(26,*) "Circle(3) = {4, 1, 5};"
          write(26,*) "Circle(4) = {5, 1, 2};"
          write(26,*) "Geometry.LineWidth = 2;"
          write(26,*) "Geometry.Color.Lines = White;"
          write(26,*) "Geometry.Points = 0;"
      close(26)
c
ccccccccccccccccccccccccccccccccccccc
c
      return
      end
