
C************************************************************************
C
      subroutine vector_p3_exp_1d(xmin, xmax, alpha, vecP3Exp)
c
c     vecP3Exp(i) = Integrate[lsP3[[i]] Exp[I alpha x], {x, xmin, xmax}]
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
      complex*16 vecP3Exp(4)

C  Local parameters:

      double precision fact1, alpha_h
      integer*8 dim1, i
      complex*16 ii, z_1, z_2, z_3, z_4
c
c  ii = sqrt(-1)
      ii = dcmplx(0.0d0, 1.0d0)

      alpha_h = alpha * (xmax - xmin)
      z_1 = Exp(ii*alpha*xmin)
      z_2 = Exp(ii*alpha*xmax)

      if(abs(alpha_h) > 1.0d-3) then
        vecP3Exp(1) = (-54 + 11*alpha_h**2 - 36*alpha_h*ii 
     *              + 2*alpha_h**3*ii)*z_1 
     *              - 2*(-27 + alpha_h**2 + 9*alpha_h*ii)*z_2
        vecP3Exp(2) = -2*(-27 + alpha_h**2 - 9*alpha_h*ii)*z_1 
     *              + (-54 + 11*alpha_h**2 + 36*alpha_h*ii 
     *              - 2*alpha_h**3*ii)*z_2
        vecP3Exp(3) = 9*(-2*(-9 + alpha_h**2 - 5*alpha_h*ii)*z_1 
     *              + (-18 + alpha_h**2 + 8*alpha_h*ii)*z_2)
        vecP3Exp(4) = -9*((18 - alpha_h**2 + 8*alpha_h*ii)*z_1 
     *              + 2*(-9 + alpha_h**2 + 5*alpha_h*ii)*z_2)
        fact1 = 2*alpha**4*(xmax - xmin)**3
        dim1 = 4
        do i=1,dim1
          vecP3Exp(i) = vecP3Exp(i) / fact1
        enddo
      else
c       Asymptotic approximation
        z_3 = Exp(ii * alpha * (xmax + 2 * xmin) / 3.0d0)
        z_4 = Exp(ii * alpha * (2 * xmax + xmin) / 3.0d0)
        vecP3Exp(1) = (30 - alpha_h**2 + 4*alpha_h*ii) * z_1
        vecP3Exp(2) = (30 - alpha_h**2 - 4*alpha_h*ii) * z_2
        vecP3Exp(3) = (90 + alpha_h**2 - 12*alpha_h*ii) * z_3
        vecP3Exp(4) = (90 + alpha_h**2 + 12*alpha_h*ii) * z_4
        fact1 = 240.0d0/(xmax-xmin)
        dim1 = 4
        do i=1,dim1
          vecP3Exp(i) = vecP3Exp(i) / fact1
        enddo
      endif
C
      end subroutine vector_p3_exp_1d

