
c
c     P2 Lagrange Interpolation nodes for the unit triangle
c
c     unit triangle =  triangle whose vertices are:
c                         (0,0,0), (1,0,0), (0,1,0).
c
      subroutine interp_nod_2d (nnodes, xn)

      integer*8 nnodes
      double precision xn(2,nnodes)
      integer*8 i


      i = 1
      xn(1,i) = 0
      xn(2,i) = 0

      i = 2
      xn(1,i) = 1
      xn(2,i) = 0

      i = 3
      xn(1,i) = 0
      xn(2,i) = 1

      i = 4
      xn(1,i) = 0.5
      xn(2,i) = 0

      i = 5
      xn(1,i) = 0.5
      xn(2,i) = 0.5

      i = 6
      xn(1,i) = 0
      xn(2,i) = 0.5
c
      return
      end
