c
      subroutine normalisation (nval, nel, nnodes,
     *  soln_k1, soln_k2, mat_overlap)
c
      implicit none
      integer*8 nval, nel, nnodes
      complex*16 soln_k1(3,nnodes+7,nval,nel)
      complex*16 soln_k2(3,nnodes+7,nval,nel)
      complex*16 mat_overlap(nval,nval)
c     Local variables
      integer*8 i, j
      integer*8 iel, ival
      complex*16 z_tmp1, z_tmp2
c
      do iel=1,nel
          do ival=1,nval
            z_tmp1 = sqrt(mat_overlap(ival,ival))
            if (abs(z_tmp1) .gt. 1.0d-8) then
c              z_tmp2 =  1.0d0/z_tmp1
              z_tmp2 =  1.0d0/z_tmp1**2
              do i=1,nnodes+7
                do j=1,3
                  soln_k1(j,i,ival,iel) =
     *               soln_k1(j,i,ival,iel)  !  * z_tmp2
                  soln_k2(j,i,ival,iel) =
     *               soln_k2(j,i,ival,iel) * z_tmp2
                enddo
              enddo
          endif
        enddo
      enddo
c
      return
      end

