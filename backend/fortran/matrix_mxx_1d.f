C************************************************************************
C
      subroutine matrix_mxx_1d (xmin, xmax, matMxx)
c
c     matMxx(i,j) = Integrate[P2[[i]] * P2[[j]], {x, xmin, xmax}]
c
C************************************************************************
C
C
      implicit none

      double precision xmin, xmax, matMxx(3,3)

C  Local parameters:

      double precision fact1
      integer*8 dim1, dim2, i, j


      matMxx(1,1) = 4*(xmax - xmin)
      matMxx(1,2) = -xmax + xmin
      matMxx(1,3) = 2*(xmax - xmin)
      matMxx(2,1) = -xmax + xmin
      matMxx(2,2) = 4*(xmax - xmin)
      matMxx(2,3) = 2*(xmax - xmin)
      matMxx(3,1) = 2*(xmax - xmin)
      matMxx(3,2) = 2*(xmax - xmin)
      matMxx(3,3) = 16*(xmax - xmin)


      fact1 = 30
      dim1 = 3
      dim2 = 3
      do j=1,dim2
        do i=1,dim1
          matMxx(i,j) = matMxx(i,j) / fact1
        enddo
      enddo
C
      end subroutine matrix_mxx_1d
