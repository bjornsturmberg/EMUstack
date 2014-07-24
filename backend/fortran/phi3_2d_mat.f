cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c   phi3_2d_mat evaluates a linear basis function (P3) and its derivative.
c
c     P3 basis function over the unit Triangle
c
c
      subroutine phi3_2d_mat(x, phi, mat_grad)
c
      implicit none
      double precision x(2), phi(10), mat_grad(2,10)
      double precision x0, y0
      integer*8 inode
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      x0 = x(1)
      y0 = x(2)
      inode = 1
        phi(inode) = -((-1 + x0 + y0)*(-2 + 3*x0 + 3*y0)
     *                *(-1 + 3*x0 + 3*y0))/2.
        mat_grad(1,inode) = (-11 - 27*x0**2 + x0*(36 - 54*y0) 
     *                + 36*y0 - 27*y0**2)/2.
        mat_grad(2,inode) = (-11 - 27*x0**2 + x0*(36 - 54*y0) 
     *                + 36*y0 - 27*y0**2)/2.
      inode = 2
        phi(inode) = (x0*(-2 + 3*x0)*(-1 + 3*x0))/2.
        mat_grad(1,inode) = 1 - 9*x0 + (27*x0**2)/2.
        mat_grad(2,inode) = 0
      inode = 3
        phi(inode) = (y0*(-2 + 3*y0)*(-1 + 3*y0))/2.
        mat_grad(1,inode) = 0
        mat_grad(2,inode) = 1 - 9*y0 + (27*y0**2)/2.
      inode = 4
        phi(inode) = (9*x0*(-1 + x0 + y0)*(-2 + 3*x0 + 3*y0))/2.
        mat_grad(1,inode) = (9*(2 + 9*x0**2 - 5*y0 + 3*y0**2 
     *                + 2*x0*(-5 + 6*y0)))/2.
        mat_grad(2,inode) = (9*x0*(-5 + 6*x0 + 6*y0))/2.
      inode = 5
        phi(inode) = (-9*x0*(-1 + 3*x0)*(-1 + x0 + y0))/2.
        mat_grad(1,inode) =(-9*(1 + 9*x0**2 - y0 + x0*(-8 + 6*y0)))/2.
        mat_grad(2,inode) = (-9*x0*(-1 + 3*x0))/2.
      inode = 6
        phi(inode) = (9*x0*(-1 + 3*x0)*y0)/2.
        mat_grad(1,inode) = (9*(-1 + 6*x0)*y0)/2.
        mat_grad(2,inode) = (9*x0*(-1 + 3*x0))/2.
      inode = 7
        phi(inode) = (9*x0*y0*(-1 + 3*y0))/2.
        mat_grad(1,inode) = (9*y0*(-1 + 3*y0))/2.
        mat_grad(2,inode) = (9*x0*(-1 + 6*y0))/2.
      inode = 8
        phi(inode) = (-9*y0*(-1 + x0 + y0)*(-1 + 3*y0))/2.
        mat_grad(1,inode) = (-9*y0*(-1 + 3*y0))/2.
        mat_grad(2,inode)=(-9*(1 - 8*y0 + 9*y0**2 +x0*(-1 + 6*y0)))/2.
      inode = 9
        phi(inode) = (9*y0*(-1 + x0 + y0)*(-2 + 3*x0 + 3*y0))/2.
        mat_grad(1,inode) = (9*y0*(-5 + 6*x0 + 6*y0))/2.
        mat_grad(2,inode) = (9*(2 + 3*x0**2 - 10*y0 + 9*y0**2 
     *                + x0*(-5 + 12*y0)))/2.
      inode = 10
        phi(inode) = -27*x0*y0*(-1 + x0 + y0)
        mat_grad(1,inode) = -27*y0*(-1 + 2*x0 + y0)
        mat_grad(2,inode) = -27*x0*(-1 + x0 + 2*y0)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      return
      end
