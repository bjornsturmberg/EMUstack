c
      subroutine list_face (nel, table_edge_face)
c
c
      implicit none
      integer*8 nel
      integer*8 table_edge_face(14,nel)
      integer*8 i
c
c     Table of connectivity for the face (for 2D FEM, face = triangle element)
c
      do i=1,nel
        table_edge_face(1,i) = i  ! each element is a face
      enddo
c
      return
      end
