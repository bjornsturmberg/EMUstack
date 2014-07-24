cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c   phi2_2d_mat evaluates a linear basis function (P3) and its derivative.
c
c     P3 basis function over the unit Triangle
c
c
      subroutine phi2_2d_mat (x, phi, mat_grad)
c
      implicit none
      double precision x(2), phi(6), mat_grad(2,6)
      double precision x0, y0
      integer*8 inode
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      x0 = x(1)
      y0 = x(2)
      inode = 1
        phi(inode) = (-1 + x0 + y0)*(-1 + 2*x0 + 2*y0)
        mat_grad(1,inode) = -3 + 4*x0 + 4*y0
        mat_grad(2,inode) = -3 + 4*x0 + 4*y0
      inode = 2
        phi(inode) = x0*(-1 + 2*x0)
        mat_grad(1,inode) = -1 + 4*x0
        mat_grad(2,inode) = 0
      inode = 3
        phi(inode) = y0*(-1 + 2*y0)
        mat_grad(1,inode) = 0
        mat_grad(2,inode) = -1 + 4*y0
      inode = 4
        phi(inode) = -4*x0*(-1 + x0 + y0)
        mat_grad(1,inode) = -4*(-1 + 2*x0 + y0)
        mat_grad(2,inode) = -4*x0
      inode = 5
        phi(inode) = 4*x0*y0
        mat_grad(1,inode) = 4*y0
        mat_grad(2,inode) = 4*x0
      inode = 6
        phi(inode) = -4*y0*(-1 + x0 + y0)
        mat_grad(1,inode) = -4*y0
        mat_grad(2,inode) = -4*(-1 + x0 + 2*y0)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      return
      end
