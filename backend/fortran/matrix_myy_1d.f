C************************************************************************
C
      subroutine matrix_myy_1d (xmin, xmax, matMyy)
c
c      matMyy(j,j) = Integrate[P3[[i]]*P3[[j]], {x, xmin, xmax}]
c
C************************************************************************
C
C
      implicit none

      double precision xmin, xmax, matMyy(4,4)

C  Local parameters:

      double precision fact1
      integer*8 dim1, dim2, i, j


      matMyy(1,1) = 128*(xmax - xmin)
      matMyy(1,2) = 19*(xmax - xmin)
      matMyy(1,3) = 99*(xmax - xmin)
      matMyy(1,4) = -36*(xmax - xmin)
      matMyy(2,1) = 19*(xmax - xmin)
      matMyy(2,2) = 128*(xmax - xmin)
      matMyy(2,3) = -36*(xmax - xmin)
      matMyy(2,4) = 99*(xmax - xmin)
      matMyy(3,1) = 99*(xmax - xmin)
      matMyy(3,2) = -36*(xmax - xmin)
      matMyy(3,3) = 648*(xmax - xmin)
      matMyy(3,4) = -81*(xmax - xmin)
      matMyy(4,1) = -36*(xmax - xmin)
      matMyy(4,2) = 99*(xmax - xmin)
      matMyy(4,3) = -81*(xmax - xmin)
      matMyy(4,4) = 648*(xmax - xmin)



      fact1 = 1680
      dim1 = 4
      dim2 = 4
      do j=1,dim2
        do i=1,dim1
          matMyy(i,j) = matMyy(i,j) / fact1
        enddo
      enddo
C

      end subroutine matrix_myy_1d
