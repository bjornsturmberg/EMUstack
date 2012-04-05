c
c***********************************************************************
c
c    Compute the product of a matrix in CSC format by a vector
c
c***********************************************************************
c

      subroutine z_mxv_csc (neq, vect1, vect2,
     *     nonz, row_ind, col_ptr, mat)
c
c***********************************************************************
c
      implicit none
      integer*8 neq, nonz
      integer*8 row_ind(nonz), col_ptr(neq+1)
      complex*16 mat(nonz)
      complex*16 vect1(neq), vect2(neq)
c
c     Local variables
      integer*8 i, j, k, col_start, col_end, i_base
c
      do i=1,neq
        vect2(i) = 0.d0
      enddo
c
c     valpr.f has changed the CSC indexing to 0-based indexing
c     so we must add 1 to the CSD row_pointer row_ind
      i_base = 1
c
      do i=1,neq   ! Column index
        col_start = col_ptr(i) + i_base
        col_end = col_ptr(i+1) - 1 + i_base
        do j=col_start,col_end
          k = row_ind(j) + i_base  ! Row number
c          !  mat(j) = value of the matrix entry (k,i)
          vect2(k) = vect2(k) + mat(j)*vect1(i)
        enddo
      enddo
c
      return
      end
