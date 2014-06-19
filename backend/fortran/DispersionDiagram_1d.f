C
C     Calculate and save data for a dispersion diagram
C                 
      subroutine DispersionDiagram_1d(lambda, nval, n_conv, 
     *    bloch_vec_x, bloch_vec_y, beta, d_in_nm)

C
      implicit none
C
      integer*8 nval, n_conv
      double precision lambda, d_in_nm
      double precision bloch_vec_x, bloch_vec_y
      complex*16 beta(nval)

      integer*8 ui, i, j, debug
      double precision lambda_nm, pi
C
C
      ui = 6
      debug = 0
      pi = 3.141592653589793d0
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  Rescalling into nm
      lambda_nm = lambda*d_in_nm

C
      if(nval .gt. 1000) write(ui,*) 
     *      "nval>1000 need to make more spaces in DispDiag.f"

C        will need to increase fmt below   \/  to accomodate all Evalues
        write(241,fmt="(2(I6),3(f16.10,2x),1000(2x,f15.10,f11.6,2x))") 
     *       n_conv, nval, lambda_nm, 
     *       bloch_vec_x/(2.0d0*pi), bloch_vec_y/(2.0d0*pi),
     *       (beta(i),i=1,nval)
c

ccccccccccc

C
      return
      end 
