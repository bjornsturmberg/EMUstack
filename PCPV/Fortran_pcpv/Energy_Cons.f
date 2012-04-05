C	  Calculate the Energy conservation equations 55-58 from M. Byrnes notes.
C       Note: only need to calculate 3 as one Hermitian Conjugate of other.
C
      subroutine Energy_Cons(R12, T12, R21, T21, numberprop_S,
     *             numberprop_N, neq_PW, nval)
C		 
      implicit none 
C
      integer*8 numberprop_S, numberprop_N, neq_PW, nval, n, s
      integer*8 i, j
c     32-but integers for BLAS and LAPACK
      integer*4 nval_max_32, PW_max_32
      parameter (nval_max_32 = 1000)
      parameter (PW_max_32 = 200)
      complex*16 R12(2*neq_PW,2*neq_PW)
      complex*16 T12(nval,2*neq_PW)
      complex*16 R21(nval,nval)
      complex*16 T21(2*neq_PW,nval)
C  Hermitian Conjugate of the above	  
      complex*16 R12H(2*PW_max_32,2*PW_max_32)
      complex*16 T12H(2*PW_max_32,nval_max_32)
      complex*16 R21H(nval_max_32,nval_max_32)
      complex*16 T21H(nval_max_32,2*PW_max_32)
C  Propagating/Evanescent Matrices
      complex*16 I_1(2*PW_max_32,2*PW_max_32)
      complex*16 I_1bar(2*PW_max_32,2*PW_max_32)
      complex*16 I_2(nval_max_32,nval_max_32)
      complex*16 I_2bar(nval_max_32,nval_max_32)
C  Parts of Equation 1
      complex*16 TMP1(2*PW_max_32,2*PW_max_32)
      complex*16 TMP2(2*PW_max_32,nval_max_32)
      complex*16 LHS11(2*PW_max_32,2*PW_max_32)
      complex*16 LHS12(2*PW_max_32,2*PW_max_32)
      complex*16 RHS11(2*PW_max_32,2*PW_max_32)
      complex*16 RHS12(2*PW_max_32,2*PW_max_32)
C  Parts of Equation 2
      complex*16 LHS21(2*PW_max_32,nval_max_32)
      complex*16 LHS22(2*PW_max_32,nval_max_32)
      complex*16 RHS21(2*PW_max_32,nval_max_32)
      complex*16 RHS22(2*PW_max_32,nval_max_32)
C  Parts of Equation 4
      complex*16 TMP3(nval_max_32,nval_max_32)
      complex*16 TMP4(nval_max_32,2*PW_max_32)
      complex*16 LHS41(nval_max_32,nval_max_32)
      complex*16 LHS42(nval_max_32,nval_max_32)
      complex*16 RHS41(nval_max_32,nval_max_32)
      complex*16 RHS42(nval_max_32,nval_max_32)
C
c     32-but integers for BLAS and LAPACK
      integer*4 neq_PW_32, nval_32

C
      complex*16 ZERO, ONE, ii
C
C
CCCCCCCCCCC Start Program CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      ii = dcmplx(0.0d0, 1.0d0)
      ZERO = 0.0d0
      ONE = 1.0d0

      neq_PW_32 = neq_PW
      nval_32 = nval
C
C  I_1 Diagonal with ones for propagating waves in medium 1 (above slab)
C  I_1bar ones for evanescent waves in medium 1
      do i=1,2*neq_PW
        do j=1,2*neq_PW
          I_1(j,i) = 0.0d0
          I_1bar(j,i) = 0.0d0
        enddo
      enddo
      do i=1,neq_PW
        if (i .le. numberprop_S) then
          I_1(i,i) = 1.0d0
          I_1(i+neq_PW,i+neq_PW) = 1.0d0
        else
          I_1bar(i,i) = ii
          I_1bar(i+neq_PW,i+neq_PW) = ii
        endif
      enddo
C
C  I_2 Diagonal with ones for propagating waves in medium 2 (the slab)
C  I_2bar ones for evanescent waves in medium 2
      do i=1,nval
        do j=1,nval
          I_2(j,i) = 0.0d0
          I_2bar(j,i) = 0.0d0
        enddo
        if (i .le. numberprop_N) then
          I_2(i,i) = 1.0d0
        else
          I_2bar(i,i) = ii
        endif
      enddo
C
C  Calculate the Hermitian Conjugate of Scattering Matrices
      do n=1,2*neq_PW
        do s=1,2*neq_PW
          R12H(s,n) = CONJG(R12(n,s))
        enddo
      enddo
      do n=1,nval
        do s=1,2*neq_PW
          T12H(s,n) = CONJG(T12(n,s))
        enddo
      enddo
      do n=1,nval
        do s=1,nval
          R21H(s,n) = CONJG(R21(n,s))
        enddo
      enddo
      do n=1,nval
        do s=1,2*neq_PW
          T21H(n,s) = CONJG(T21(s,n))
        enddo
      enddo
C
C  Carry out the Matrix Multiplication for the first equation
      call ZGEMM('N','N', 2*neq_PW_32, 2*neq_PW_32, 2*neq_PW_32, 
     *  ONE, R12H, 2*PW_max_32, I_1, 2*PW_max_32, ZERO, 
     *  TMP1, 2*PW_max_32)
      call ZGEMM('N','N', 2*neq_PW_32, 2*neq_PW_32, 2*neq_PW_32, 
     *  ONE, TMP1, 2*PW_max_32, R12, 2*neq_PW_32, ZERO, 
     *  LHS11, 2*PW_max_32)
     
      call ZGEMM('N','N', 2*neq_PW_32, nval_32, nval_32, ONE,
     *  T12H, 2*PW_max_32, I_2, nval_max_32, ZERO, 
     *  TMP2, 2*PW_max_32)
      call ZGEMM('N','N', 2*neq_PW_32, 2*neq_PW_32, nval_32, ONE,
     *  TMP2, 2*PW_max_32, T12, nval_32, ZERO, 
     *  LHS12, 2*PW_max_32)

      call ZGEMM('N','N', 2*neq_PW_32, 2*neq_PW_32, 2*neq_PW_32, 
     *  ONE, R12H, 2*PW_max_32, I_1bar, 2*PW_max_32, ZERO, 
     *  RHS11, 2*PW_max_32)

      call ZGEMM('N','N', 2*neq_PW_32, 2*neq_PW_32, 2*neq_PW_32, 
     *  ONE, I_1bar, 2*PW_max_32, R12, 2*neq_PW_32, ZERO, 
     *  RHS12, 2*PW_max_32)
     
      open (unit=342, file='Normed/Energy_Cons_1.txt', 
     *        status='unknown')  
      do n=1,2*neq_PW
        do s=1,2*neq_PW
          write(342,*) LHS11(s,n) + LHS12(s,n) 
     *                 - I_1(s,n) - RHS11(s,n) + RHS12(s,n)
        enddo
      enddo
      close(342)
C
C  Carry out the Matrix Multiplication for the second equation
      call ZGEMM('N','N', 2*neq_PW_32, 2*neq_PW_32, 2*neq_PW_32, 
     *  ONE, R12H, 2*PW_max_32, I_1, 2*PW_max_32, ZERO, 
     *  TMP1, 2*PW_max_32)
      call ZGEMM('N','N', 2*neq_PW_32, nval_32, 2*neq_PW_32, ONE,
     *  TMP1, 2*PW_max_32, T21, 2*neq_PW_32, ZERO, 
     *  LHS21, 2*PW_max_32)
     
      call ZGEMM('N','N', 2*neq_PW_32, nval_32, nval_32, ONE,
     *  T12H, 2*PW_max_32, I_2, nval_max_32, ZERO, 
     *  TMP2, 2*PW_max_32)
      call ZGEMM('N','N', 2*neq_PW_32, nval_32, nval_32, ONE,
     *  TMP2, 2*PW_max_32, R21, nval_32, ZERO, 
     *  LHS22, 2*PW_max_32)

      call ZGEMM('N','N', 2*neq_PW_32, nval_32, nval_32, ONE,
     *  T12H, 2*PW_max_32, I_2bar, nval_max_32, ZERO, 
     *  RHS21, 2*PW_max_32)

      call ZGEMM('N','N', 2*neq_PW_32, nval_32, 2*neq_PW_32, ONE,
     *  I_1bar, 2*PW_max_32, T21, 2*neq_PW_32, ZERO, 
     *  RHS22, 2*PW_max_32)
     
      open (unit=342, file='Normed/Energy_Cons_2.txt', 
     *        status='unknown')  
      do n=1,nval
        do s=1,2*neq_PW
          write(342,*) LHS21(s,n) + LHS22(s,n) 
     *                 - RHS21(s,n) + RHS22(s,n)
        enddo
      enddo
      close(342)
C
C  Carry out the Matrix Multiplication for the fourth equation
      call ZGEMM('N','N', nval_32, nval_32, nval_32, ONE,
     *  R21H, nval_max_32, I_2, nval_max_32, ZERO, 
     *  TMP3, nval_max_32)
      call ZGEMM('N','N', nval_32, nval_32, nval_32, ONE,
     *  TMP3, nval_max_32, R21, nval_32, ZERO, 
     *  LHS41, nval_max_32)
     
      call ZGEMM('N','N', nval_32, 2*neq_PW_32, 2*neq_PW_32, ONE,
     *  T21H, nval_max_32, I_1, 2*PW_max_32, ZERO,
     *  TMP4, nval_max_32)
      call ZGEMM('N','N', nval_32, nval_32, 2*neq_PW_32, ONE,
     *  TMP4, nval_max_32, T21, 2*neq_PW_32, ZERO, 
     *  LHS42, nval_max_32)

      call ZGEMM('N','N', nval_32, nval_32, nval_32, ONE,
     *  R21H, nval_max_32, I_2bar, nval_max_32, ZERO, 
     *  RHS41, nval_max_32)

      call ZGEMM('N','N', nval_32, nval_32, nval_32, ONE,
     *  I_2bar, nval_max_32, R21, nval_32, ZERO, 
     *  RHS42, nval_max_32)
     
      open (unit=342, file='Normed/Energy_Cons_4.txt', 
     *        status='unknown')  
      do n=1,nval
        do s=1,nval
          write(342,*) LHS41(s,n) + LHS42(s,n)
     *                 - I_2(s,n) - RHS41(s,n) + RHS42(s,n)
        enddo
      enddo
      close(342)
C
C
      return
      end 
