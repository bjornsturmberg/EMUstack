C     Calculate the ultimate efficiency as in Lin and Povinelli 09
C                 
      subroutine Ultimate_Efficiency(Lambda_Res, d_in_nm)
C
      implicit none
C
      integer*8 Lambda_Res, n_lambda, i, d_in_nm
      double precision data1(818),data2(818)
      double precision data3(1635),data4(1635)
      double precision Lambda1(818),Lambda2(1635)
      double precision lambda_1, lambda_2
      double precision Limit_values, Ult_Eff
      double precision Tot_Irradiance, Prefactor, Integral
      character*4 buf1
      character*100 tchar, tchar2
C
CCCCCCCCCCC Read in Data CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Total solar irradiance - integral of I(lambda) from 310nm-4000nm
C        done in Mathematica (OtherCode/Silicon_ASTM/ASTMG173.nb)
      Tot_Irradiance = 900.084
      lambda_1 = 310
      lambda_2 = 1127
      Prefactor = (lambda_2-lambda_1)/(lambda_2*Tot_Irradiance)
      Integral = 0.0d0
C
      write(buf1,'(i4.4)') d_in_nm 
C
      if (Lambda_Res .eq. 1) then
        n_lambda = 700
        open (unit=346, file='../Data/irradiance_1.txt',status='old')
        tchar = "Output/d_"//buf1//"p_0000_A_Lambda.txt"
        open (unit=347, file=tchar,status='old')
          do i = 1,n_lambda
            read(346,100) Lambda1(i), data1(i)   ! direct & circumsolar spectrum
            read(347,101) data2(i)               ! calculated absorption spectrum
          enddo
        close(346)
        close(347)
C
        do i=2,(n_lambda-1)
          Integral = Integral + data1(i)*data2(i)*Lambda1(i)
        enddo
        Limit_values = lambda_1*data1(1)*data2(1) +
     *   lambda_2*data1(n_lambda)*data2(n_lambda)
        Ult_Eff = (Integral+Limit_values/2)
     *    *(Prefactor/818)   !where 818 is # actual data points (inclusive zeros)
C
C
      elseif (Lambda_Res .eq. 2) then
        n_lambda = 1400
        open (unit=346, file='../Data/irradiance_2.txt',status='old')
        tchar = "Output/d_"//buf1//"p_0000_A_Lambda.txt"
        open (unit=347, file=tchar,status='old')
          do i = 1,n_lambda
            read(346,100) Lambda2(i), data3(i)  ! direct & circumsolar spectrum
            read(347,101) data4(i)              ! calculated absorption spectrum
          enddo
        close(346)
        close(347)
C
        do i=2,(n_lambda-1)
          Integral = Integral + data3(i)*data4(i)*Lambda2(i)
        enddo
        Limit_values = lambda_1*data3(1)*data4(1) +
     *   lambda_2*data3(n_lambda)*data4(n_lambda)
        Ult_Eff = (Integral+Limit_values/2)
     *    *(Prefactor/1635)   !where 1635 is # actual data points (inclusive zeros)
      endif
C
C
      tchar2 = "Output/d_"//buf1//"_Ultimate_Efficiency.txt"
      open (unit=345,
     * file=tchar2, status='unknown')
        write(345,131) d_in_nm, Ult_Eff 
      close(345)
C
      write(*,*) "Ultimate Efficiency with d=",d_in_nm,"is ",Ult_Eff
C
131   format(I7, G25.17)
100   format(F7.1,F8.5)
101   format(36x,G25.17)
C
      return
      end 
