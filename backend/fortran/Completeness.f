C	  Carry out matrix multiplication for Completeness Check
C
      subroutine Completeness (nval, neq_PW, 
     *       K_overlap_mat, overlap_J)
C		 
      implicit none 
C
      integer*8 nval, neq_PW
c     32-but integers for BLAS and LAPACK
      integer*8 nval_max_32, PW_max_32
      parameter (nval_max_32 = 1000)
      parameter (PW_max_32 = 200)
      complex*16 K_overlap_mat(nval,2*neq_PW)
      complex*16 overlap_J(2*neq_PW,nval)
      complex*16 JK_sumN(2*PW_max_32,2*PW_max_32)
      complex*16 KJ_sumS(nval_max_32,nval_max_32)
      complex*16 ZERO, ONE
      integer*8 n, s
c     32-but integers for BLAS and LAPACK
      integer*8 neq_PW_32, nval_32

      ZERO = 0.0d0
      ONE = 1.0d0
C
CCCCCCCCCCCCCCCCCCCC Start Program CCCCCCCCCCCCCCCCCCCCCCCCC
C     
      nval_32 = nval
      neq_PW_32 = neq_PW
      call ZGEMM('N','N', nval_32, nval_32, 2*neq_PW_32, 
     *  ONE, K_overlap_mat, nval_32, overlap_J, 2*neq_PW_32, 
     *  ZERO, KJ_sumS, nval_max_32)
     
      call ZGEMM('N','N', 2*neq_PW_32, 2*neq_PW_32, nval_32, 
     *  ONE, overlap_J, 2*neq_PW_32, K_overlap_mat, nval_32, 
     *  ZERO, JK_sumN, 2*PW_max_32)

      open (unit=33, file="Normed/KJ_sumS.txt",
     * status='unknown')
      open (unit=32, file="Normed/JK_sumN.txt", 
     * status='unknown')
      do n=1,nval
        do s=1,nval
          write(33,132) n, s, KJ_sumS(n,s), abs(KJ_sumS(n,s))
        enddo
      enddo
      do n=1,2*neq_PW
        do s=1,2*neq_PW
          write(32,132) n, s, JK_sumN(n,s), abs(JK_sumN(n,s))
        enddo
      enddo
      close(33)
      close(32)
C
132    format(2(I4),2(g25.17),g18.10)      
C
      return
      end 
