
c
c***********************************************************************
c
      subroutine sort_csr (neq, nonz, max_row_len, col_ind, 
     *  row_ptr, arr, indx, istack)
c
      implicit none
      integer*8 neq, nonz, max_row_len
      integer*8 row_ptr(neq+1), col_ind(nonz)
      integer*8 arr(max_row_len), indx(max_row_len)
      integer*8 istack(max_row_len)
c
      integer*8 row_start, row_end, row_len
      integer*8 i, j, k
c
      do i=1,neq
        row_start = row_ptr(i)
        row_end = row_ptr(i+1) - 1
        row_len = row_end - row_start + 1
cc        print*, "sort_csr: i = ", i
        do j=row_start,row_end
cc          print*, "    ", col_ind(j)
          k = j - row_start + 1
          arr(k) = col_ind(j)
        enddo
        call sort_int (row_len, arr, indx, istack)
        do j=row_start,row_end
          k = j - row_start + 1
cc          print*, "    ", col_ind(j), arr(k), arr(indx(k))
        enddo
        do j=row_start,row_end
          k = j - row_start + 1
          col_ind(j) = arr(indx(k))
        enddo
      enddo

c
c      do i=1,neq
c        row_start = row_ptr(i)
c        row_end = row_ptr(i+1) - 1
c        row_len = row_end - row_start + 1
c        print*, "sort_csr: i = ", i
c        do j=row_start,row_end
c          print*, " ####   ", col_ind(j)
c        enddo
c      enddo
c
c       stop
c
      return
      end

