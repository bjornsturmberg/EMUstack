C
C     Calculate and save data for a dispersion diagram
C                 
      subroutine DispersionDiagram(lambda, bloch_vec, shift,
     *    nval, n_conv, beta, mode_pol, d_in_nm)
C
      implicit none
C
      integer*8 nval, ui, i, j, debug, d_in_nm
      integer*8 action, n_conv
      complex*16 beta(nval)
      complex*16 beta_tmp(nval), shift, mode_pol(4,nval)
      double precision lambda, lambda_nm, bloch_vec(2), pi
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

      if(debug .eq. 1) then
        write(ui,*)
        write(ui,*) "lambda, 1/lambda = ", lambda_nm, 1.0d0/lambda_nm
        write(ui,*) (bloch_vec(i)/(2.0d0*pi),i=1,2)
        write(ui,*) "sqrt(shift) = ", sqrt(shift)
        do i=1,nval
          write(ui,"(i4,2(f20.12),4(f18.10))") i, 
     *       beta(i), beta(i)**2
        enddo
        write(ui,*)
        do i=1,nval
          write(ui,"(i4,5(g15.7))") i, 
     *       dble(beta(i)),
     *       (dble(mode_pol(j,i)),j=1,4)
        enddo
      endif

C
      if(nval .gt. 1000) write(ui,*) 
     *      'nval>1000 need to make more spaces in DispDiag.f'
C        will need to increase fmt below   \/  to accomodate all Evalues
        write(241,fmt="(2(I6),3(f16.10,2x),1000(2x,f15.10,f11.6,2x))") 
     *       n_conv, nval, lambda_nm, 
     *       (bloch_vec(i)/(2.0d0*pi),i=1,2),
     *       (beta(i),i=1,nval)
        write(231,fmt="(2(I6),3(f16.10,2x),1000(2x,f15.10,f11.6,2x))") 
     *       n_conv, nval, lambda_nm,
     *       (bloch_vec(i)/(2.0d0*pi),i=1,2),
     *       (dble(mode_pol(4,i)),i=1,nval)
C        do i=1,nval
C          write(231,"(i4,5(g16.8))") i,
C     *       dble(beta(i)),
C     *       (dble(mode_pol(j,i)),j=1,4)
C        enddo

ccccccccccc
c       Collect the modes which have dominant z-component
        do i=1,nval
          j = 3 ! with mode_pol(j,i), the z-component corrresponds to j=3 
          if (abs(mode_pol(j,i)) .gt. 0.5) then
            beta_tmp(i) = beta(i)
          else
            beta_tmp(i) = 0.0d0
          endif
        enddo
C
        write(221,fmt="(2(I6),3(f15.8,2x),1000(2x,f15.10,f11.6,2x))") 
     *       n_conv, nval, lambda_nm,
     *       (bloch_vec(i)/(2.0d0*pi),i=1,2),
     *       (beta_tmp(i),i=1,nval)
C
      if (debug .eq. 1) then
        write(ui,*)
        write(ui,*) "F_z-dominant modes"
        do i=1,nval
          if (beta_tmp(i) .ne. 0.) then
            write(ui,"(i4,5(g15.7))") i, 
     *       dble(beta_tmp(i)),
     *       (dble(mode_pol(j,i)),j=1,4)
          endif
        enddo
      endif
ccccccccccc
c       Collect the modes which are tranverse (have small z-component)
        do i=1,nval
          j = 3 ! with mode_pol(j,i), the z-component corrresponds to j=3 
          if (abs(mode_pol(j,i)) .lt. 0.1) then
            beta_tmp(i) = beta(i)
          else
            beta_tmp(i) = 0.0d0
          endif
        enddo
C
        write(211,fmt="(2(I6),3(f15.8,2x),1000(2x,f15.10,f11.6,2x))") 
     *       n_conv, nval, lambda_nm,
     *       (bloch_vec(i)/(2.0d0*pi),i=1,2),
     *       (beta_tmp(i),i=1,nval)
C
      if (debug .eq. 1) then
        write(ui,*)
        write(ui,*) "F_t-dominant modes (transverse)"
        do i=1,nval
          if (beta_tmp(i) .ne. 0.) then
            write(ui,"(i4,5(g15.7))") i, 
     *       dble(beta_tmp(i)),
     *       (dble(mode_pol(j,i)),j=1,4)
          endif
        enddo
      endif
C
      return
      end 
