c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c

c     sol_P2(1,1..3,nval, nel)  contains the values of Ex component at P2 interpolation nodes (3 nodes)
c     sol_P2(2,1..3,nval, nel)  contains the values of Ey component at P2 interpolation nodes
c     sol_P2(3,1..3,nval, nel)  contains the values of Ez component at P2 interpolation nodes
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine write_sol_P2_1d (nval, nel, E_H_field,
     *     lambda, beta, sol_P2, mesh_file, dir_name)
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit none
      integer*8 nval, nel, E_H_field
      double precision lambda
      complex*16 sol_P2(3,3,nval,nel), beta(nval)
      character mesh_file*100

c     Local variables

c      integer*8 j, inod, iel, ival, debug
c      double precision P3_mid_mode_value(4)
c      complex*16 z_tmp1, z_tmp2

c      integer*8 nddl_0, nddl_t, dim
c      parameter (nnodes_0 = 6)

      integer*8 i, iel, ival, nnodes_P2
      integer*8 namelen, namelength, namelen2

      character*11 ivalue, jvalue
      character dir_name*100
      character*100 tchar1, tchar2
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
      nnodes_P2 = 3
      namelen = len_trim(mesh_file)
      namelength = len_trim(dir_name)
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      open (unit=63, file="Matrices/beta.txt",
     *         status="unknown")
        do ival=1,nval
          write(63,12) ival, beta(ival)
        enddo
      close(63)

      do ival=1,nval
        jvalue = ivalue(ival)
        namelen2 = len_trim(jvalue)
        tchar1 = dir_name(1:namelength)//"/Sol_P2/r_mode_P2_"
     *     //jvalue(1:namelen2)//".txt"
        tchar2 = dir_name(1:namelength)//"/Sol_P2/i_mode_P2_"
     *     //jvalue(1:namelen2)//".txt"
        open(63, file=tchar1,
     *     form="formatted")
        open(64, file=tchar2,
     *     form="formatted")
        write(63,*) mesh_file(1:namelen)
        write(64,*) mesh_file(1:namelen)
        write(63,"((g25.17),2(I7))") lambda, ival, E_H_field
        write(64,"((g25.17),2(I7))") lambda, ival, E_H_field
        write(63,"(2(g25.17))") beta(ival)
        write(64,"(2(g25.17))") beta(ival)
        do iel=1,nel
          write(63,12) iel, (dble(sol_P2(1,i,ival,iel)),i=1,nnodes_P2)  ! x-component
          write(63,12) iel, (dble(sol_P2(2,i,ival,iel)),i=1,nnodes_P2)  ! y-component
          write(63,12) iel, (dble(sol_P2(3,i,ival,iel)),i=1,nnodes_P2)  ! z-component
c
          write(64,12) iel, (imag(sol_P2(1,i,ival,iel)),i=1,nnodes_P2)  ! x-component
          write(64,12) iel, (imag(sol_P2(2,i,ival,iel)),i=1,nnodes_P2)  ! y-component
          write(64,12) iel, (imag(sol_P2(3,i,ival,iel)),i=1,nnodes_P2)  ! z-component
        enddo
        close(63)
        close(64)
      enddo
c
12    format(I7,13(g25.17) )
c
      return
      end
