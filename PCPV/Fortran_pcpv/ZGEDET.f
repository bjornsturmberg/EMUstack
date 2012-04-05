C     Calculate the determinant of a Matrix with complex entries
C     inputs M-dimension of (square matrix), A
      subroutine ZGEDET( M, A, LDA, IPIV, INFO, detA ) 
C
      implicit none
C  Scalar Arguments
      integer     M           ! dimension of the (square) matrix
      integer     LDA         ! leading dimension of matrix
      integer     INFO        ! info on LPU factorisation see ZGETRF from LAPACK
C  Array Arguments
      integer     IPIV( * )   ! pivot indicies, row i was interchanged with IPIV(i)
      complex*16  A( LDA, * ) ! input matrix to find determinant of
      complex*16  detA        ! output determinant
C  Internal Arguments
      integer     i
C
C  carry out LPU factorisation of matrix
      call ZGETRF( M, M, A, LDA, IPIV, INFO )
      if(INFO .ne. 0) then
        write(*,*) "Failure in Deteriminant"
        write(*,*) "PROBLEM INFO = ", INFO
        write(*,*) "Aborting..."
        stop
      endif
C  detA = product of diagonals, note pivot swaps rows (thus need for *-1)
      detA = 1.0d0
      do i=1,M
         if (IPIV(i) .eq. i) then
            detA = detA * A(i,i)
         else
            detA = -detA * A(i,i)
         endif
      enddo
C
      return
      end 