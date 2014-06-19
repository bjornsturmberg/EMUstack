C************************************************************************
C
      subroutine matrix_kyx_1d (xmin, xmax, matKyx)
c
c     matKyx(i,j) = Integrate[P3'[[i]] * P2[[j]], {x, xmin, xmax}]
c
C************************************************************************
C
C
      implicit none

      double precision xmin, xmax, matKyx(4,3)

C  Local parameters:

      double precision fact1
      integer*8 dim1, dim2, i, j


      matKyx(1,1) = -83
      matKyx(1,2) = 7
      matKyx(1,3) = -44
      matKyx(2,1) = -7
      matKyx(2,2) = 83
      matKyx(2,3) = 44
      matKyx(3,1) = 99
      matKyx(3,2) = 9
      matKyx(3,3) = -108
      matKyx(4,1) = -9
      matKyx(4,2) = -99
      matKyx(4,3) = 108

      fact1 = 120
      dim1 = 4
      dim2 = 3
      do j=1,dim2
        do i=1,dim1
          matKyx(i,j) = matKyx(i,j) / fact1
        enddo
      enddo
C

      end subroutine matrix_kyx_1d
