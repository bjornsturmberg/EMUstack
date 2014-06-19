
C************************************************************************
C
      subroutine vector_p2_exp_1d(xmin, xmax, alpha, vecP2Exp)
c
c     vecP2Exp(i) = Integrate[lsP2[[i]] Exp[I alpha x], {x, xmin, xmax}]
c
c	The exact formula for int_x1^x2 (P_i*exp(I alpha x)) is used when 
c	alpha is not too small.
c
c	The exact formula is numerically unstable for small values of alpha; 
c	and an asymptotic approximation based on the Taylor formula of 
c	order 2 is used.
c
C************************************************************************
C
      implicit none

      double precision xmin, xmax, alpha
      complex*16 vecP2Exp(3)

C  Local parameters:

      double precision fact1, alpha_h
      integer*8 dim1, i
      complex*16 ii, z_1, z_2, z_3
c
c  ii = sqrt(-1)
      ii = dcmplx(0.0d0, 1.0d0)

      alpha_h = alpha * (xmax - xmin)
      z_1 = Exp(ii*alpha*xmin)
      z_2 = Exp(ii*alpha*xmax)
      if(abs(alpha_h) > 1.0d-3) then
        vecP2Exp(1) = ii*(-4 + alpha_h**2 - 3*alpha_h*ii)*z_1 
     *              + (alpha_h + 4*ii)*z_2
        vecP2Exp(2) = (alpha_h - 4*ii)*z_1 
     *              - ii*(-4 + alpha_h**2 + 3*alpha_h*ii)*z_2
        vecP2Exp(3) = -4*( 2*ii*(-z_1 + z_2) 
     *              + alpha_h*(z_1 + z_2) )
        fact1 = alpha**3*(xmax - xmin)**2
        dim1 = 3
        do i=1,dim1
          vecP2Exp(i) = vecP2Exp(i) / fact1
        enddo
      else
c       Asymptotic approximation
        z_3 = Exp(ii*alpha*(xmax+xmin)/2.0d0)
        vecP2Exp(1) = (20 + alpha_h**2)*z_1
        vecP2Exp(2) = (20 + alpha_h**2)*z_2
        vecP2Exp(3) = -2*(-40 + alpha_h**2)*z_3
        fact1 = 120/(xmax - xmin)
        dim1 = 3
        do i=1,dim1
          vecP2Exp(i) = vecP2Exp(i) / fact1
        enddo
      endif
C
      end subroutine vector_p2_exp_1d

