C************************************************************************
C
      subroutine matrix_kyy_1d (xmin, xmax, matKyy)
c
c     matKyy(i,j) - Integrate[P3'[[i]] * P3'[[j]], {x, xmin, xmax}]
c
C************************************************************************
C
C
      implicit none

      double precision xmin, xmax, matKyy(4,4)

C  Local parameters:

      double precision fact1
      integer*8 dim1, dim2, i, j


      matKyy(1,1) = 148
      matKyy(1,2) = -13
      matKyy(1,3) = -189
      matKyy(1,4) = 54
      matKyy(2,1) = -13
      matKyy(2,2) = 148
      matKyy(2,3) = 54
      matKyy(2,4) = -189
      matKyy(3,1) = -189
      matKyy(3,2) = 54
      matKyy(3,3) = 432
      matKyy(3,4) = -297
      matKyy(4,1) = 54
      matKyy(4,2) = -189
      matKyy(4,3) = -297
      matKyy(4,4) = 432

      fact1 = 40 * (xmax-xmin)
      dim1 = 4
      dim2 = 4
      do j=1,dim2
        do i=1,dim1
          matKyy(i,j) = matKyy(i,j) / fact1
        enddo
      enddo
C

      end subroutine matrix_kyy_1d
