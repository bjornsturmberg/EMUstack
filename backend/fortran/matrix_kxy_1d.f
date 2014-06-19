C************************************************************************
C
      subroutine matrix_kxy_1d (xmin, xmax, matKxy)
c
c     matKxy(i,j) = Integrate[lsP2[[i]] P3'[[j]], {x, xmin, xmax}]
c
C************************************************************************
C
C
      implicit none

      double precision xmin, xmax, matKxy(3,4)

C  Local parameters:

      double precision fact1
      integer*8 dim1, dim2, i, j

      matKxy(1,1) = -83
      matKxy(1,2) = -7
      matKxy(1,3) = 99
      matKxy(1,4) = -9
      matKxy(2,1) = 7
      matKxy(2,2) = 83
      matKxy(2,3) = 9
      matKxy(2,4) = -99
      matKxy(3,1) = -44
      matKxy(3,2) = 44
      matKxy(3,3) = -108
      matKxy(3,4) = 108

      fact1 = 120
      dim1 = 3
      dim2 = 4
      do j=1,dim2
        do i=1,dim1
          matKxy(i,j) = matKxy(i,j) / fact1
        enddo
      enddo
C
      end subroutine matrix_kxy_1d
