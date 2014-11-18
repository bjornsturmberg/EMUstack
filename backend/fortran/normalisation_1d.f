c
      subroutine normalisation_1d (nval, nel,
     *  sol_1, sol_2, mat_overlap)
c
      implicit none
      integer*8 nval, nel
      complex*16 sol_1(3+4+4,nval,nel)
      complex*16 sol_2(3+4+4,nval,nel)
      complex*16, dimension(nval,nval) :: mat_overlap
c     Local variables
      integer*8 i, iel, ival
      integer*8 nddl_0
      parameter (nddl_0 = 11)
      complex*16 z_tmp1, z_tmp2
c
      do iel=1,nel
        do ival=1,nval
          z_tmp1 = sqrt(mat_overlap(ival,ival))
          if (abs(z_tmp1) .gt. 1.0d-8) then
            z_tmp2 =  1.0d0/z_tmp1
C              z_tmp2 =  1.0d0/z_tmp1**2
            do i=1,nddl_0
              sol_1(i,ival,iel) =
     *           sol_1(i,ival,iel) * z_tmp2
              sol_2(i,ival,iel) =
     *           sol_2(i,ival,iel) * z_tmp2
            enddo
          endif
        enddo
      enddo
c
      return
      end

