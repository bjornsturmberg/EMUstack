c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     sol_1(1..3,nval, nel)  contains the values of Ex component at P2 interpolation nodes
c     sol_1(4..7,nval, nel)  contains the values of Ey component at P3 interpolation nodes
c     sol_1(8..11,nval, nel) contains the values of Ez component at P3 interpolation nodes

c     sol_P2(1,1..3,nval, nel)  contains the values of Ex component at P2 interpolation nodes (3 nodes)
c     sol_P2(2,1..3,nval, nel)  contains the values of Ey component at P2 interpolation nodes
c     sol_P2(3,1..3,nval, nel)  contains the values of Ez component at P2 interpolation nodes
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine array_sol_P2_1d (nval, nel, sol_1, sol_P2)
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit none
      integer*8 nval, nel
      complex*16 sol_1(3+4+4,nval,nel)
      complex*16 sol_P2(3,3,nval,nel)

c     Local variables
      integer*8 j, inod, iel, ival, debug
      double precision P3_mid_mode_value(4)
      complex*16 z_tmp1, z_tmp2

c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      debug = 0

c     Values of the P3 basis function at the mid-point
      P3_mid_mode_value(1) = -1.0d0 / 16.0d0
      P3_mid_mode_value(2) = -1.0d0 / 16.0d0
      P3_mid_mode_value(3) =  9.0d0 / 16.0d0
      P3_mid_mode_value(4) =  9.0d0 / 16.0d0

      do iel=1,nel
        do ival=1,nval
          do inod=1,3
            sol_P2(1,inod,ival,iel) = sol_1(inod,ival,iel)  ! x-component
          enddo
          inod = 1
            sol_P2(2,inod,ival,iel) = sol_1(inod+3,ival,iel)    ! y-component
            sol_P2(3,inod,ival,iel) = sol_1(inod+3+4,ival,iel)  ! z-component
          inod = 2
            sol_P2(2,inod,ival,iel) = sol_1(inod+3,ival,iel)    ! y-component
            sol_P2(3,inod,ival,iel) = sol_1(inod+3+4,ival,iel)  ! z-component
c         Interpolated value for the mid-node
c         The initial P3 value of Ey is interpolated to obtain the value at the mid-node
          z_tmp1 = 0
          z_tmp2 = 0
          do j=1,4
            z_tmp1 = z_tmp1 + sol_1(j+3,ival,iel) * P3_mid_mode_value(j)
            z_tmp2 = z_tmp2 + sol_1(j+3+4,ival,iel)*P3_mid_mode_value(j)
          enddo
          inod = 3  ! Mid-node
            sol_P2(2,inod,ival,iel) = z_tmp1   ! y-component
            sol_P2(3,inod,ival,iel) = z_tmp2   ! z-component
        enddo
      enddo
c
      return
      end
