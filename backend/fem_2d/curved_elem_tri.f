c
ccccccccccccccccccccccccccccccccccccccc
c
c     Check if an element (triangle) has a curved face 
c
ccccccccccccccccccccccccccccccccccccccc
c
      subroutine curved_elem_tri (nnodes, xel, info_curved, tmp)
c
ccccccccccccccccccccccccccccccccccccccc
c
      implicit none
      integer*8 nnodes, info_curved
      double precision xel(2,nnodes)

      integer*8 nnd_triangle
      parameter (nnd_triangle = 6)
      double precision xel_triangle(2,nnd_triangle)

      integer*8 i, j, i2
      double precision tol, tmp
c
ccccccccccccccccccccccccccccccccccccccc
c
      if (nnodes .ne. nnd_triangle) then
        write(*,*)
        write(*,*) "   ???"
        write(*,*) "curved_elem_tri: nnodes != nnd_triangle : ",
     *           nnodes, nnd_triangle
        write(*,*) "curved_elem_tri: Aborting..."
        stop
      endif
c     Vertices
      do i=1,3
        do j=1,2
        xel_triangle(j,i) = xel(j,i)
        enddo
      enddo

c     Mid-points
      do i=1,3
         i2 = Mod(i+1,3)
         if(i2 .eq. 0) i2 = 3
        do j=1,2
          xel_triangle(j,i+3) = (xel(j,i)+xel(j,i2))/2.0d0
        enddo
      enddo
c
      tmp = 0.0d0
      do i=1,nnodes
        do j=1,2
        tmp = tmp + (xel_triangle(j,i) - xel(j,i))**2
        enddo
      enddo

      tol = 1.0d-14
      if(abs(tmp) .lt. tol) then
        info_curved = 0
      else
        info_curved = 1
      endif

      info_curved = 0

      return
      end
