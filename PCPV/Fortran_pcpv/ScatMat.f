C     Carry out matrix multiplication to calculate Scattering Matrices T, R 
C
      subroutine ScatMat (J, J_dagger, X, neq_PW, nval, Beta, 
     *             T12, R12, T21, R21, PrintAll,
     *             PrintSolution, lx, h_1, h_2, num_h, Checks,
     *             TLambda, RLambda, traLambda, pol, PropModes, 
     *             lambda, d_in_nm,
     *    numberprop_S, freq, Zeroth_Order_inv,
     *    debug, incident, what4incident, out4incident,
     *    title, parallel_1)
C
      implicit none
C
      integer*8 neq_PW, nval, PropModes, d_in_nm, num_h
C     32-bit integers for BLAS and LAPACK
      integer*4 nval_max_32, PW_max_32, h_i
      integer*8 debug, numberprop_S, Zeroth_Order_inv
      integer*8 title, parallel_1
      character buf1*4, buf2*4
      complex*16 Beta(nval)
      double precision h, h_1, h_2, h_int, lambda, freq
      parameter (nval_max_32 = 1000)
      parameter (PW_max_32 = 200)
      complex*16 J(2*neq_PW,nval)
      complex*16 J_dagger(nval,2*neq_PW)
      complex*16 X(2*neq_PW,2*neq_PW)
      complex*16 X_2(2*PW_max_32,2*PW_max_32)
      complex*16 TMP0(2*PW_max_32,nval_max_32)
      complex*16 THETA(nval_max_32,nval_max_32)
      complex*16 THETA_I(nval_max_32,nval_max_32)
      complex*16 PHI(nval_max_32,2*PW_max_32)
      complex*16 PSI(2*PW_max_32,nval_max_32)
      complex*16 TMP1(nval_max_32,2*PW_max_32)
      complex*16 TMP2(nval_max_32,nval_max_32)
C  Fresnel Scattering Coefficients
      complex*16 R12(2*neq_PW,2*neq_PW)
      complex*16 T12(nval,2*neq_PW)
      complex*16 R21(nval,nval)
      complex*16 T21(2*neq_PW,nval)
C   
      complex*16 P(nval_max_32,nval_max_32)
      complex*16 UPSILON(nval_max_32,nval_max_32)
      complex*16 TMP3(nval_max_32,nval_max_32)
      complex*16 TMP3_I(nval_max_32,nval_max_32)
      complex*16 TMP4(2*PW_max_32,nval_max_32)
      complex*16 GAMMA(2*PW_max_32,nval_max_32)
      complex*16 GAMMA_2(nval_max_32,2*PW_max_32)
      complex*16 GAMMA_3(nval_max_32,2*PW_max_32)
C  Transmission and Reflection for each Wavelength
      complex*16 TLambda(2*neq_PW,2*neq_PW)
      complex*16 RLambda(2*neq_PW,2*neq_PW)
C  Single Value Decomposition variables
      complex*16 EVEN(nval_max_32,nval_max_32)
      complex*16 ODD(nval_max_32,nval_max_32)
      double precision SIGMA(nval_max_32,nval_max_32)
      double precision SIGMAe(nval_max_32,nval_max_32)
      double precision SIGMAo(nval_max_32,nval_max_32)
      double precision RWORK(3*nval_max_32)
c     32-but integers for BLAS and LAPACK
      integer*8 LWORK_32, LWMAX_32
      parameter (LWMAX_32 = 50000)
      complex*16 U, WORK(LWMAX_32)
C
C  For Checking Implementation
      complex*16 TMP551(2*PW_max_32,2*PW_max_32)
      complex*16 TMP661(2*PW_max_32,2*PW_max_32)
      complex*16 TMP771(2*PW_max_32,2*PW_max_32)
      complex*16 TMP552(nval_max_32,nval_max_32)
      complex*16 TMP662(nval_max_32,nval_max_32)
      complex*16 TMP772(nval_max_32,nval_max_32)
C      
      integer*8 i, k, PrintAll, Checks, pol
      complex*16 MONE, ZERO, ONE, TWO, ii
C
      double precision total_t, total_r, total_a, tot_0
      integer*8 allincident, incident, what4incident
      integer*8 out4incident, test, traLambda, PrintSolution
      complex*16 detA, detAe, detAo 
      double precision lx

c     32-but integers for BLAS and LAPACK
      integer*4 neq_PW_32, nval_32
      integer*4 INFO_32, NRHS_32, IPIV_32(nval_max_32)
C 
      MONE = -1.0d0
      ZERO = 0.0d0
      ONE = 1.0d0
      TWO = 2.0d0
      ii = dcmplx(0.0d0, 1.0d0)
      neq_PW_32 = neq_PW
      nval_32 = nval
C 
C
C  X^{1/2}
      do k = 1,2*neq_PW
        do i = 1,2*neq_PW
          X_2(i,k) = 0.0d0
        enddo
        X_2(k,k) = SQRT(X(k,k))
      enddo     
C
C  THETA 
      call ZGEMM('N','N', 2*neq_PW_32, nval_32, 2*neq_PW_32, 
     *  ONE, X, 2*neq_PW_32, J, 2*neq_PW_32, ZERO, 
     *  TMP0, 2*PW_max_32)
      call ZGEMM('N','N', nval_32, nval_32, 2*neq_PW_32, ONE, 
     *  J_dagger, nval_32, TMP0, 2*PW_max_32, ZERO, 
     *  THETA, nval_max_32)
      do i = 1, nval
       THETA(i,i) = THETA(i,i) + 1.0d0
      enddo
C      
C  THETA_I
      do i = 1,nval
        do k = 1,nval
          THETA_I(i,k) = 0.0d0       
        enddo
        THETA_I(i,i) = 1.0d0
      enddo
C     Invert 
      NRHS_32 = nval
      call ZGESV( nval_32, NRHS_32, THETA, nval_max_32, 
     *      IPIV_32, THETA_I, nval_max_32, INFO_32 )
      if(INFO_32 .ne. 0) then
        write(*,*) "ScatMat: PROBLEM INFO_32 = ", INFO_32
        write(*,*) "Aborting... pos 1"
        stop
      endif     
C
C  PHI
      call ZGEMM('N','N', nval_32, 2*neq_PW_32, 2*neq_PW_32, ONE, 
     *  J_dagger, nval_32, X_2, 2*PW_max_32, ZERO, 
     *  PHI, nval_max_32)
C
C  PSI
      call ZGEMM('N','N', 2*neq_PW_32, nval_32, 2*neq_PW_32, ONE, 
     *  X_2, 2*PW_max_32, J, 2*neq_PW_32, ZERO, 
     *  PSI, 2*PW_max_32)
C
C
CCCCCCCCCCCCCCC   Scattering Matrices
C
C  T12
      call ZGEMM('N','N', nval_32, 2*neq_PW_32, nval_32, TWO, 
     *  THETA_I, nval_max_32, PHI, nval_max_32, ZERO, 
     *  T12, nval_32)
C  R12
      call ZGEMM('N','N', nval_32, 2*neq_PW_32, nval_32, ONE, 
     *  THETA_I, nval_max_32, PHI, nval_max_32, ZERO, 
     *  TMP1, nval_max_32)
      call ZGEMM('N','N', 2*neq_PW_32, 2*neq_PW_32, nval_32, TWO, 
     *  PSI, 2*PW_max_32, TMP1, nval_max_32, ZERO, 
     *  R12, 2*neq_PW_32)
      do i = 1, 2*neq_PW
        R12(i,i) = R12(i,i) - 1.0d0
      enddo
C  T21
      call ZGEMM('N','N', 2*neq_PW_32, nval_32, nval_32, TWO, 
     *  PSI, 2*PW_max_32, THETA_I, nval_max_32, ZERO, 
     *  T21, 2*neq_PW_32)    
C  R21
      call ZGEMM('N','N', nval_32, nval_32, 2*neq_PW_32, MONE, 
     *  PHI, nval_max_32, PSI, 2*PW_max_32, ZERO, 
     *  TMP2, nval_max_32)
      do i = 1, nval
        TMP2(i,i) = TMP2(i,i) + 1.0d0
      enddo     
      call ZGEMM('N','N', nval_32, nval_32, nval_32, ONE, 
     *  THETA_I, nval_max_32, TMP2, nval_max_32, ZERO, 
     *  R21, nval_32)
C      
CCCCCCCCCCCCCCC   Introduce Propagation in z
C
      do h_i = 1, num_h
        if(h_i .eq. 1) then
          h_int = 0
        else
          h_int = (h_2 - h_1)/(num_h-1)
        endif
        h = h_1 + (h_i-1)*h_int
C
C  P
      do k = 1, nval
        do i = 1,nval
          P(i,k) = 0.0d0
        enddo
        P(k,k) = EXP(ii*Beta(k)*h)
      enddo
C  UPSILON
      call ZGEMM('N','N', nval_32, nval_32, nval_32, ONE, 
     *  R21, nval_32, P, nval_max_32, ZERO, 
     *  UPSILON, nval_max_32)
C  GAMMA
      call ZGEMM('N','N', 2*neq_PW_32, nval_32, nval_32, ONE, 
     *  T21, 2*neq_PW_32, P, nval_max_32, ZERO, 
     *  TMP4, 2*PW_max_32)   
      call ZGEMM('N','N', nval_32, nval_32, nval_32, MONE, 
     *  UPSILON, nval_max_32, UPSILON, nval_max_32, ZERO, 
     *  TMP3, nval_max_32) 
      do i = 1, nval
        TMP3(i,i) = TMP3(i,i) + 1.0d0
      enddo
C     Invert 
      do k = 1,nval
        do i = 1,nval
          TMP3_I(i,k) = 0.0d0
        enddo
        TMP3_I(k,k) = 1.0d0
      enddo
      NRHS_32 = nval
      call ZGESV( nval_32, NRHS_32, TMP3, nval_max_32, 
     *      IPIV_32, TMP3_I, nval_max_32, INFO_32 )
      if(INFO_32 .ne. 0) then
        write(*,*) "ScatMat: PROBLEM INFO_32 = ", INFO_32
        write(*,*) "Aborting... pos 2"
        stop
      endif
      call ZGEMM('N','N', 2*neq_PW_32, nval_32, nval_32, ONE, 
     *  TMP4, 2*PW_max_32, TMP3_I, nval_max_32, ZERO, 
     *  GAMMA, 2*PW_max_32)
C
CCCCCCCCCCCCCCC Reflection, Transmission Matrices for each Lambda
      if (traLambda .eq. 1) then
C  TLambda
        call ZGEMM('N','N', 2*neq_PW_32, 2*neq_PW_32, nval_32, ONE, 
     *    GAMMA, 2*PW_max_32, T12, nval_32, ZERO, 
     *    TLambda, 2*neq_PW_32)       
C  RLambda
        call ZGEMM('N','N', 2*neq_PW_32, nval_32, nval_32, ONE, 
     *    GAMMA, 2*PW_max_32, UPSILON, nval_max_32, ZERO, 
     *    TMP4, 2*PW_max_32)
        call ZGEMM('N','N', 2*neq_PW_32, 2*neq_PW_32, nval_32, ONE, 
     *    TMP4, 2*PW_max_32, T12, nval_32, ZERO, 
     *    RLambda, 2*neq_PW_32)   
        do i = 1, 2*neq_PW
          do k = 1, 2*neq_PW
            RLambda(i,k) = R12(i,k) + RLambda(i,k)
          enddo
        enddo
      endif
C
      if (traLambda .eq. 1) then
      if (debug .eq. 1) then
        write(*,*) "ScatMat: A_and_W_Lambda"
      endif
        call A_and_W_Lambda_sub(TLambda, RLambda, neq_PW, 
     *    numberprop_S, numberprop_S, lambda, pol, 
     *    Zeroth_Order_inv, debug, d_in_nm,
     *    incident, what4incident, out4incident, Checks, h) 
      endif
C
      enddo ! num_h
C      
CCCCCCCCCCCCCCC   Calculating det(I-RPRP)   CCCCCCCCCCCCCCCCC
C
      if (PropModes .ne. 0) then
        do i = 1, nval
          do k = 1, nval
            EVEN(i,k) = UPSILON(i,k)
            ODD(i,k) = UPSILON(i,k)
          enddo
          EVEN(i,i) = EVEN(i,i) - 1.0d0
          ODD(i,i) = ODD(i,i) + 1.0d0
        enddo
      endif
      if (PropModes .eq. 1) then
        call ZGEDET(nval_32, EVEN, nval_max_32, IPIV_32, INFO_32, 
     *         detAe)
        call ZGEDET(nval_32, ODD, nval_max_32, IPIV_32, INFO_32, 
     *         detAo)
        write(565,103) lambda*d_in_nm, detAe
        write(566,103) lambda*d_in_nm, detAo
C  Carry out SVD to find BM composition of megamode
      elseif (PropModes .eq. 2 .and. PrintSolution .ne. 0.0d0) then
        if (3*nval .gt. LWMAX_32) then
          write(*,*) "ScatMat: need to increase SVD workspace"
          write(*,*) "Aborting..."
          stop
        endif
C  Query optimal workspace size
C        LWORK_32 = LWMAX_32 -1
C        CALL ZGESVD('N', 'O', nval_32, nval_32, TMP3, nval_max_32, 
C     *         SIGMA,
C     *         U, 1, U, 1, WORK, LWORK_32, RWORK, INFO_32)
C        LWORK_32 = MIN( LWMAX_32, INT( WORK( 1 ) ) )
        LWORK_32 = LWMAX_32
C  Carry out Single value decomposition
        call ZGESVD('N', 'O', nval_32, nval_32, EVEN, nval_max_32, 
     *         SIGMAe,
     *         U, 1, U, 1, WORK, LWORK_32, RWORK, INFO_32)
        call ZGESVD('N', 'O', nval_32, nval_32, ODD, nval_max_32, 
     *         SIGMAo,
     *         U, 1, U, 1, WORK, LWORK_32, RWORK, INFO_32)

C for new python data structure
C        open (unit=344, file="Matrices/I-RPRP-svmatE.txt",
C     *         status="unknown")
C        open (unit=346, file="Matrices/I-RPRP-svmatO.txt",
C     *         status="unknown")

        open (unit=344, file="I-RPRP-svmatE.txt",
     *         status="unknown")
        open (unit=346, file="I-RPRP-svmatO.txt",
     *         status="unknown")
        do i=1,nval

          write(344,102) i, EVEN(nval,i)
          write(346,102) i, ODD(nval,i)
        enddo 
        close(344)
        close(346)
      endif
C
102   format(I6,2(G25.17))
103   format(F6.1,2(G25.17))
C      
CCCCCCCCCCCCCCC   Coeffecients c+, c-   CCCCCCCCCCCCCCCCC
C
      if (PrintSolution .ne. 0.0d0) then
C
      call ZGEMM('N','N', nval_32, 2*neq_PW_32, nval_32, ONE, 
     *  TMP3_I, nval_max_32, T12, nval_32, ZERO, 
     *  GAMMA_2, nval_max_32)

      call ZGEMM('N','N', nval_32, 2*neq_PW_32, nval_32, ONE, 
     *  UPSILON, nval_max_32, GAMMA_2, nval_max_32, ZERO, 
     *  GAMMA_3, nval_max_32)

C      open (unit=343, file='Matrices/T_c_minus.txt',
      open (unit=343, file='T_c_minus.txt',
     *         status='unknown')
      if (pol .eq. 1) then  !For TE (s) polarised light
      do k=1,neq_PW
        do i=1,nval
          write(343,101) i, k, GAMMA_2(i,k), 
     *     abs(GAMMA_2(i,k))**2
        enddo
      enddo
      elseif (pol .eq. 2) then  !For TM (p) polarised light
        do k=1,neq_PW
          do i=1,nval
            write(343,101) i, k, GAMMA_2(i,neq_PW+k), 
     *       abs(GAMMA_2(i,neq_PW+k))**2
          enddo
        enddo
      endif
      close(343)
c
C      open (unit=343, file='Matrices/T_c_plus.txt',
      open (unit=343, file='T_c_plus.txt',
     *         status='unknown')
      do k=1,neq_PW
        do i=1,nval
          write(343,101) i, k, GAMMA_3(i,k), abs(GAMMA_3(i,k))**2
        enddo
      enddo 
      close(343)
      endif ! PrintSolution .ne. 0.0d0
C      
CCCCCCCCCCCCCCC   Saving Matrices to files   CCCCCCCCCCCCCCCCC
C
C      if (traLambda .eq. 1) then
C      open (unit=341, file='Matrices/TLambda.txt',
C     *         status='unknown')
C      open (unit=342, file='Matrices/RLambda.txt',
C     *         status='unknown')
C     do k=1,2*neq_PW
C        do i=1,2*neq_PW
C          write(341,101) i, k, TLambda(i,k), abs(TLambda(i,k))**2
C          write(342,101) i, k, RLambda(i,k), abs(RLambda(i,k))**2
C        enddo
C      enddo       
C      close(341)      
C      close(342)
C      endif
C
C      if (PrintAll .eq. 1) then
      write(buf1,'(I4.4)') title
      write(buf2,'(I4.4)') parallel_1 

CC      open (unit=345, file='Matrices/R12.txt',   
C      open (unit=345, file="st"//buf1//"_wl"//buf2//"_R12.txt",
C     *         status='unknown')
C      do k=1,2*neq_PW
C        do i=1,2*neq_PW
C          write(345,101) i, k, R12(i,k), abs(R12(i,k))**2
C        enddo
C      enddo 
C      close(345)
CC      open (unit=345, file='Matrices/T12.txt',
C      open (unit=345, file="st"//buf1//"_wl"//buf2//'_T12.txt',
C     *         status='unknown')
C      do k=1,2*neq_PW
C        do i=1,nval
C          write(345,101) i, k, T12(i,k), abs(T12(i,k))**2
C        enddo
C      enddo 
C      close(345)
CC      open (unit=345, file='Matrices/R21.txt',
C      open (unit=345, file="st"//buf1//"_wl"//buf2//'_R21.txt',
C     *         status='unknown')
C      do k=1,nval
C        do i=1,nval
C          write(345,101) i, k, R21(i,k), abs(R21(i,k))**2
C        enddo
C      enddo 
C      close(345)
CC      open (unit=345, file='Matrices/T21.txt',
C      open (unit=345, file="st"//buf1//"_wl"//buf2//'_T21.txt',
C     *         status='unknown')
C      do k=1,nval
C        do i=1,2*neq_PW
C          write(345,101) i, k, T21(i,k), abs(T21(i,k))**2
C        enddo
C      enddo 
C      close(345)
C     P matrix for propagation through finite thickness
      open (unit=345, file="st"//buf1//"_wl"//buf2//"_P.txt",
     *         status='unknown')
      do k=1,nval
        do i=1,nval
          write(345,1000) i, k, P(i,k), abs(P(i,k))**2
        enddo
      enddo 
      close(345)
C     Beta as will be used in P matrix      
      open (unit=345, file="st"//buf1//"_wl"//buf2//'_beta.txt',
     *         status='unknown')
      do k=1,nval
        write(345,1001) Beta(k)
      enddo 
      close(345)

C
101   format(I4,I4,3(G25.17))
1000  format(I4,I4,3(G27.18E3))
1001  format(2(G25.17))
C      endif ! PrintAll .ne. 0.0d0
C
CCCCCCCCCC Check implementation of scattering matrices
C
      if (Checks .eq. 2) then
C
      call ZGEMM('N','N', 2*neq_PW_32, 2*neq_PW_32, nval_32, ONE, 
     *  T21, 2*neq_PW_32, T12, nval_32, ZERO, 
     *  TMP551, 2*PW_max_32)

      do i = 1,2*neq_PW
        do k = 1, 2*neq_PW
            TMP661(i,k) = -(R12(i,k))**2
       enddo
           TMP661(i,i) = TMP661(i,i) + 1.0d0
      enddo

      do i = 1,2*neq_PW
        do k = 1, 2*neq_PW
          TMP771(i,k) = TMP551(i,k) - TMP661(i,k)
        enddo
      enddo

C      open (unit=342, file='Normed/ScatMat_T21T12.txt', 
C     *          status='unknown') 
C      open (unit=344, file='Normed/ScatMat_I_Rsquared1.txt', 
C     *          status='unknown')  
C      open (unit=345, file='Normed/ScatMat_DIFF1.txt', 
C     *          status='unknown')

      open (unit=342, file='ScatMat_T21T12.txt', 
     *          status='unknown') 
      open (unit=344, file='ScatMat_I_Rsquared1.txt', 
     *          status='unknown')  
      open (unit=345, file='ScatMat_DIFF1.txt', 
     *          status='unknown') 
      do i=1,2*neq_PW
        do k=1,2*neq_PW
          write(342,*) TMP551(i,k)        
          write(344,*) TMP661(i,k)        
          write(345,*) TMP771(i,k)
        enddo
      enddo 
      close(342)
      close(344)
      close(345)


      call ZGEMM('N','N', nval_32, nval_32, 2*neq_PW_32, ONE, 
     *  T12, nval_32, T21, 2*neq_PW_32, ZERO, 
     *  TMP552, nval_max_32)

      do i = 1,nval
        do k = 1, nval
            TMP662(i,k) = -(R21(i,k))**2
       enddo
            TMP662(i,i) = TMP662(i,i) + 1.0d0
      enddo

      do i = 1,nval
        do k = 1, nval
          TMP772(i,k) = TMP552(i,k) - TMP662(i,k)
        enddo
      enddo

C      open (unit=342, file='Normed/ScatMat_T12T21.txt', 
C     *          status='unknown') 
C      open (unit=344, file='Normed/ScatMat_I_Rsquared2.txt', 
C     *          status='unknown')  
C      open (unit=345, file='Normed/ScatMat_DIFF2.txt', 
C     *          status='unknown') 

      open (unit=342, file='ScatMat_T12T21.txt', 
     *          status='unknown') 
      open (unit=344, file='ScatMat_I_Rsquared2.txt', 
     *          status='unknown')  
      open (unit=345, file='ScatMat_DIFF2.txt', 
     *          status='unknown')  
      do i=1,nval
        do k=1,nval
          write(342,*) TMP552(i,k)        
          write(344,*) TMP662(i,k)        
          write(345,*) TMP772(i,k)
        enddo
      enddo 
      close(342)
      close(344)
      close(345)
      endif
C
CCCCCCCCCCCCCCCCCCCCCCCCC
C
      return
      end 
