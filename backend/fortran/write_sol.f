c

      subroutine write_sol (nval, nel, nnodes, E_H_field,
     *     lambda, beta, sol, mesh_file, dir_name)
c
      implicit none
      integer*8 nval, nel, nnodes, E_H_field
      double precision lambda
      complex*16 sol(3,nnodes+7,nval,nel)
      complex*16 beta(nval)
      character mesh_file*100
c     Local variables
      integer*8 i, iel, ival
      integer*8 namelen, namelength, namelen2

      character*11 ivalue, jvalue
      character dir_name*100
      character*100 tchar1, tchar2
c
      namelen = len_trim(mesh_file)
      namelength = len_trim(dir_name)

      do ival=1,nval
        jvalue = ivalue(ival)
        namelen2 = len_trim(jvalue)
        tchar1 = dir_name(1:namelength)//"/r_mode_"
     *     //jvalue(1:namelen2)//".txt"
        tchar2 = dir_name(1:namelength)//"/i_mode_"
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
          write(63,12) iel, (dble(sol(1,i,ival,iel)),i=1,nnodes+7)  ! x-component
          write(63,12) iel, (dble(sol(2,i,ival,iel)),i=1,nnodes+7)  ! 3-component
          write(63,12) iel, (dble(sol(3,i,ival,iel)),i=1,nnodes+7)  ! z-component
c
          write(64,12) iel, (imag(sol(1,i,ival,iel)),i=1,nnodes+7)  ! x-component
          write(64,12) iel, (imag(sol(2,i,ival,iel)),i=1,nnodes+7)  ! 3-component
          write(64,12) iel, (imag(sol(3,i,ival,iel)),i=1,nnodes+7)  ! z-component
        enddo
        close(63)
        close(64)
      enddo
12    format(I7,13(g25.17) )
c
      return
      end
