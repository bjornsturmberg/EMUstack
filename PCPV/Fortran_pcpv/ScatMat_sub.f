C     Carry out matrix multiplication to calculate Scattering Matrices T, R 
C
      subroutine ScatMat_sub (J, J_dagger, X, X_b, neq_PW, nval, 
     *             Beta, T12, R12, T21, R21, PrintAll,
     *             PrintSolution, lx, h_1, h_2, num_h, Checks, 
     *             TLambda, RLambda, traLambda, pol, PropModes, 
     *             lambda, d_in_nm,
     *    numberprop_S, numberprop_S_b, freq, Zeroth_Order_inv,
     *    debug, incident, what4incident, out4incident,
     *    title, parallel_1)
C
      implicit none
C
      integer*8 neq_PW, nval, PropModes, d_in_nm, num_h
C     32-but integers for BLAS and LAPACK
      integer*4 nval_max_32, PW_max_32, h_i
      integer*8 debug, numberprop_S, numberprop_S_b
      integer*8 Zeroth_Order_inv
      integer*8 title, parallel_1
      character buf1*4, buf2*4
      complex*16 Beta(nval)
      double precision h, h_1, h_2, h_int, lambda, freq
      parameter (nval_max_32 = 1000)
      parameter (PW_max_32 = 200)
      complex*16 J(2*neq_PW,nval)
      complex*16 J_dagger(nval,2*neq_PW)
      complex*16 X(2*neq_PW,2*neq_PW)
      complex*16 X_b(2*neq_PW,2*neq_PW)
C
      complex*16 X_2(2*PW_max_32,2*PW_max_32)
      complex*16 X_2_b(2*PW_max_32,2*PW_max_32)
      complex*16 TMP0(2*PW_max_32,nval_max_32)
      complex*16 TMP0_b(2*PW_max_32,nval_max_32)
      complex*16 THETA(nval_max_32,nval_max_32)
      complex*16 THETA_b(nval_max_32,nval_max_32)
      complex*16 THETA_I(nval_max_32,nval_max_32)
      complex*16 THETA_b_I(nval_max_32,nval_max_32)
      complex*16 PHI(nval_max_32,2*PW_max_32)
      complex*16 PSI(2*PW_max_32,nval_max_32)
      complex*16 PHI_b(nval_max_32,2*PW_max_32)
      complex*16 PSI_b(2*PW_max_32,nval_max_32)
      complex*16 TMP1(nval_max_32,2*PW_max_32)
      complex*16 TMP2(nval_max_32,nval_max_32)
      complex*16 TMP1_b(nval_max_32,2*PW_max_32)
      complex*16 TMP2_b(nval_max_32,nval_max_32)
C  Fresnel Scattering Coefficients
      complex*16 R12(2*neq_PW,2*neq_PW)
      complex*16 T12(nval,2*neq_PW)
      complex*16 R21(nval,nval)
      complex*16 T21(2*neq_PW,nval)
      complex*16 R23(nval,nval)
      complex*16 T23(2*neq_PW,nval)
C   
      complex*16 P(nval_max_32,nval_max_32)
      complex*16 UPSILON(nval_max_32,nval_max_32)
      complex*16 UPSILON_b(nval_max_32,nval_max_32)
      complex*16 TMP3(nval_max_32,nval_max_32)
      complex*16 TMP3_I(nval_max_32,nval_max_32)
      complex*16 TMP7(2*PW_max_32,nval_max_32)
      complex*16 TMP4(2*PW_max_32,nval_max_32)
      complex*16 TMP5(2*PW_max_32,nval_max_32)
      complex*16 TMP6(2*PW_max_32,nval_max_32)
      complex*16 GAMMA(2*PW_max_32,nval_max_32)
      complex*16 GAMMA_2(nval_max_32,2*PW_max_32)
      complex*16 GAMMA_3(nval_max_32,2*PW_max_32)
C  Transmission and Reflection for each Wavelength
      complex*16 TLambda(2*neq_PW,2*neq_PW)
      complex*16 RLambda(2*neq_PW,2*neq_PW)
C  Single Value Decomposition variables
      complex*16 EVEN(nval_max_32,nval_max_32)
C      double precision SIGMA(nval_max_32,nval_max_32)
      double precision SIGMAe(nval_max_32,nval_max_32)
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
C  X^{1/2}
      do k = 1,2*neq_PW
        do i = 1,2*neq_PW
          X_2(i,k) = 0.0d0
        enddo
        X_2(k,k) = SQRT(X(k,k))
      enddo
C
C  X_b^{1/2}
      do k = 1,2*neq_PW
        do i = 1,2*neq_PW
          X_2_b(i,k) = 0.0d0
        enddo
        X_2_b(k,k) = SQRT(X_b(k,k))
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
        write(*,*) "ScatMat_sub: PROBLEM INFO_32 = ", INFO_32
        write(*,*) "Aborting... pos 1"
        stop
      endif
C
C  THETA_b 
      call ZGEMM('N','N', 2*neq_PW_32, nval_32, 2*neq_PW_32, 
     *  ONE, X_b, 2*neq_PW_32, J, 2*neq_PW_32, ZERO, 
     *  TMP0_b, 2*PW_max_32)
      call ZGEMM('N','N', nval_32, nval_32, 2*neq_PW_32, ONE, 
     *  J_dagger, nval_32, TMP0_b, 2*PW_max_32, ZERO, 
     *  THETA_b, nval_max_32)
      do i = 1, nval
       THETA_b(i,i) = THETA_b(i,i) + 1.0d0
      enddo
C      
C  THETA_b_I
      do i = 1,nval
        do k = 1,nval
          THETA_b_I(i,k) = 0.0d0       
        enddo
        THETA_b_I(i,i) = 1.0d0
      enddo
C     Invert 
      NRHS_32 = nval
      call ZGESV( nval_32, NRHS_32, THETA_b, nval_max_32, 
     *      IPIV_32, THETA_b_I, nval_max_32, INFO_32 )
      if(INFO_32 .ne. 0) then
        write(*,*) "ScatMat_sub: PROBLEM INFO_32 = ", INFO_32
        write(*,*) "Aborting... pos 1"
        stop
      endif    
C
C  PHI
      call ZGEMM('N','N', nval_32, 2*neq_PW_32, 2*neq_PW_32, ONE, 
     *  J_dagger, nval_32, X_2, 2*PW_max_32, ZERO, 
     *  PHI, nval_max_32)
C
C  PHI_b
      call ZGEMM('N','N', nval_32, 2*neq_PW_32, 2*neq_PW_32, ONE, 
     *  J_dagger, nval_32, X_2_b, 2*PW_max_32, ZERO, 
     *  PHI_b, nval_max_32)
C  PSI
      call ZGEMM('N','N', 2*neq_PW_32, nval_32, 2*neq_PW_32, ONE, 
     *  X_2, 2*PW_max_32, J, 2*neq_PW_32, ZERO, 
     *  PSI, 2*PW_max_32)
C
C  PSI_b
      call ZGEMM('N','N', 2*neq_PW_32, nval_32, 2*neq_PW_32, ONE, 
     *  X_2_b, 2*PW_max_32, J, 2*neq_PW_32, ZERO, 
     *  PSI_b, 2*PW_max_32)
C
C
CCCCCCCCCCCCCCC   Scattering Matrices
C
C From top medium to SiNW
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
C
C From SiNW to top medium
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
C From SiNW to bottom medium
C  T23
      call ZGEMM('N','N', 2*neq_PW_32, nval_32, nval_32, TWO, 
     *  PSI_b, 2*PW_max_32, THETA_b_I, nval_max_32, ZERO, 
     *  T23, 2*neq_PW_32)    
C  R23
      call ZGEMM('N','N', nval_32, nval_32, 2*neq_PW_32, MONE, 
     *  PHI_b, nval_max_32, PSI_b, 2*PW_max_32, ZERO, 
     *  TMP2_b, nval_max_32)
      do i = 1, nval
        TMP2_b(i,i) = TMP2_b(i,i) + 1.0d0
      enddo     
      call ZGEMM('N','N', nval_32, nval_32, nval_32, ONE, 
     *  THETA_b_I, nval_max_32, TMP2_b, nval_max_32, ZERO, 
     *  R23, nval_32)
C
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

C  UPSILON_b
      call ZGEMM('N','N', nval_32, nval_32, nval_32, ONE, 
     *  R23, nval_32, P, nval_max_32, ZERO, 
     *  UPSILON_b, nval_max_32)
C
C
       call ZGEMM('N','N', 2*neq_PW_32, nval_32, nval_32, ONE, 
     *  T23, 2*neq_PW_32, P, nval_max_32, ZERO, 
     *  TMP7, 2*PW_max_32)  
C

C  GAMMA
      call ZGEMM('N','N', 2*neq_PW_32, nval_32, nval_32, ONE, 
     *  T21, 2*neq_PW_32, P, nval_max_32, ZERO, 
     *  TMP4, 2*PW_max_32)   
      call ZGEMM('N','N', nval_32, nval_32, nval_32, MONE, 
     *  UPSILON, nval_max_32, UPSILON_b, nval_max_32, ZERO, 
     *  TMP3, nval_max_32)
C TMP3 = I - R21 P R23 P
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
     *  TMP7, 2*PW_max_32, TMP3_I, nval_max_32, ZERO, 
     *  GAMMA, 2*PW_max_32)
C
       call ZGEMM('N','N', 2*neq_PW_32, nval_32, nval_32, ONE, 
     *  TMP4, 2*PW_max_32, UPSILON_b, nval_max_32, ZERO, 
     *  TMP5, 2*PW_max_32)
C
       call ZGEMM('N','N', 2*neq_PW_32, nval_32, nval_32, ONE, 
     *  TMP5, 2*PW_max_32, TMP3_I, nval_max_32, ZERO, 
     *  TMP6, 2*PW_max_32)
C
C
CCCCCCCCCCCCCCC Reflection, Transmission Matrices for each Lambda
      if (traLambda .eq. 1) then
C
C  TLambda
        call ZGEMM('N','N', 2*neq_PW_32, 2*neq_PW_32, nval_32, ONE, 
     *    GAMMA, 2*PW_max_32, T12, nval_32, ZERO, 
     *    TLambda, 2*neq_PW_32)
C     
C  RLambda
        call ZGEMM('N','N', 2*neq_PW_32, 2*neq_PW_32, nval_32, ONE, 
     *    TMP6, 2*PW_max_32, T12, nval_32, ZERO, 
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
        write(*,*) "ScatMat_sub: A_and_W_Lambda_sub"
      endif
      call A_and_W_Lambda_sub(TLambda, RLambda, neq_PW,
     *    numberprop_S, numberprop_S_b, lambda, pol, 
     *    Zeroth_Order_inv, debug, d_in_nm, incident,
     *    what4incident, out4incident, Checks, h)
      endif
C
      enddo ! num_h
C
C      
CCCCCCCCCCCCCCC   Calculating det(I-RPRP)   CCCCCCCCCCCCCCCCC
C
C TMP3 = I - R21 P R23 P
      if (PropModes .eq. 1) then
        call ZGEDET(nval_32, TMP3, nval_max_32, IPIV_32, INFO_32, 
     *         detAe)
        write(565,103) lambda*d_in_nm, detAe
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

        open (unit=344, file="Matrices/I-RPRP-svmatE.txt",
     *         status="unknown")
        do i=1,nval

          write(344,102) i, EVEN(nval,i)
        enddo 
        close(344)
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

      open (unit=343, file='Matrices/T_c_minus.txt',
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
      open (unit=343, file='Matrices/T_c_plus.txt',
     *         status='unknown')
      do k=1,neq_PW
        do i=1,nval
          write(343,101) i, k, GAMMA_3(i,k), abs(GAMMA_3(i,k))**2
        enddo
      enddo 
      close(343)
C      
CCCCCCCCCCCCCCC   Saving Matrices to files   CCCCCCCCCCCCCCCCC
C
      open (unit=341, file='Matrices/TLambda.txt',
     *         status='unknown')
      open (unit=342, file='Matrices/RLambda.txt',
     *         status='unknown')
      do k=1,2*neq_PW
        do i=1,2*neq_PW
          write(341,101) i, k, TLambda(i,k), abs(TLambda(i,k))**2
          write(342,101) i, k, RLambda(i,k), abs(RLambda(i,k))**2
        enddo
      enddo       
      close(341)      
      close(342)
      endif  ! PrintSolution .ne. 0.0d0
C
      write(buf1,'(I4.4)') title
      write(buf2,'(I4.4)') parallel_1
C      open (unit=345, file='Matrices/R12.txt',   
      open (unit=345, file="st"//buf1//"_wl"//buf2//"_R12.txt",
     *         status='unknown')
      do k=1,2*neq_PW
        do i=1,2*neq_PW
          write(345,101) i, k, R12(i,k), abs(R12(i,k))**2
        enddo
      enddo 
      close(345)
C      open (unit=345, file='Matrices/T12.txt',
      open (unit=345, file="st"//buf1//"_wl"//buf2//'_T12.txt',
     *         status='unknown')
      do k=1,2*neq_PW
        do i=1,nval
          write(345,101) i, k, T12(i,k), abs(T12(i,k))**2
        enddo
      enddo 
      close(345)
C      open (unit=345, file='Matrices/R21.txt',
      open (unit=345, file="st"//buf1//"_wl"//buf2//'_R21.txt',
     *         status='unknown')
      do k=1,nval
        do i=1,nval
          write(345,101) i, k, R21(i,k), abs(R21(i,k))**2
        enddo
      enddo 
      close(345)
C      open (unit=345, file='Matrices/T21.txt',
      open (unit=345, file="st"//buf1//"_wl"//buf2//'_T21.txt',
     *         status='unknown')
      do k=1,nval
        do i=1,2*neq_PW
          write(345,101) i, k, T21(i,k), abs(T21(i,k))**2
        enddo
      enddo 
      close(345)
      open (unit=345, file="st"//buf1//"_wl"//buf2//'_R23.txt',
     *         status='unknown')
      do k=1,nval
        do i=1,nval
          write(345,101) i, k, R23(i,k), abs(R23(i,k))**2
        enddo
      enddo 
      close(345)
      open (unit=345, file="st"//buf1//"_wl"//buf2//'_T23.txt',
     *         status='unknown')
      do k=1,nval
        do i=1,2*neq_PW
          write(345,101) i, k, T23(i,k), abs(T23(i,k))**2
        enddo
      enddo 
      close(345)
C
101   format(I4,I4,3(G25.17))
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

      open (unit=342, file='Normed/ScatMat_T21T12.txt', 
     *          status='unknown') 
      open (unit=344, file='Normed/ScatMat_I_Rsquared1.txt', 
     *          status='unknown')  
      open (unit=345, file='Normed/ScatMat_DIFF1.txt', 
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

      open (unit=342, file='Normed/ScatMat_T12T21.txt', 
     *          status='unknown') 
      open (unit=344, file='Normed/ScatMat_I_Rsquared2.txt', 
     *          status='unknown')  
      open (unit=345, file='Normed/ScatMat_DIFF2.txt', 
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
