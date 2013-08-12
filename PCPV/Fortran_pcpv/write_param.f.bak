c
c     E_H_field = 1 => Electric field formulation (E-Field)
c     E_H_field = 2 => Magnetic field formulation (H-Field)
c
c     i_cond = 0 => Dirichlet boundary condition
c     i_cond = 1 => Neumann boundary condition
c     i_cond = 2 => Periodic boundary condition
c
      subroutine write_param (E_H_field, lambda, npt, nel, i_cond,
     *      nval, nvect, itermax, tol, shift, lx, ly, 
     *      mesh_file, mesh_format, n_conv, nb_typ_el, eps_eff,
     *      bloch_vec, dir_name)
c
c----------------------------------------------------------
c---- write data (lecture des donnees)
c----------------------------------------------------------
c
      implicit none
      integer*8 E_H_field, npt, nel, i_cond
      integer*8 nval, nvect, itermax
      integer*8 n_conv, nb_typ_el
      double precision lambda, tol, lx, ly
      double precision bloch_vec(2)
      complex*16 shift
      complex*16 eps_eff(nb_typ_el)
      integer*8 mesh_format, namelen
      character mesh_file*100
      character dir_name*100
      integer*8 namelength
      character*100 tchar1, tchar2, tchar3
c     Local variables
      integer*8 i
      double precision pi
c
       pi = 3.141592653589793d0
c
      namelength = len_trim(dir_name)
      tchar1 = dir_name(1:namelength)//"/param_out.txt"
C	  
      open(3,file=tchar1)
        write(3,*) E_H_field, "  ! E_H_field"
        write(3,*) lambda, "  ! lambda"
        write(3,*) i_cond, "  ! i_cond"
        write(3,*) lx, "  ! lx"
        write(3,*) ly, "  ! ly"
        write(3,*) nval, "  ! nval"
        write(3,*) nvect, "  ! nvect"
        write(3,*) itermax, "  ! itermax"
        write(3,*) tol, "  ! tol"
        write(3,12) shift, "  ! shift"
        write(3,'(a100)') mesh_file
        write(3,*) mesh_format, "  ! mesh_format"
        write(3,*) npt, nel, "  ! npt, nel"
        write(3,*) n_conv, "  ! n_conv"
      close(3)
c
      tchar2 = dir_name(1:namelength)//"/refrac_index_out.txt"
      open (unit=65,file=tchar2)
        write(65,*) nb_typ_el
        do i=1,nb_typ_el
          write(65,12) sqrt(eps_eff(i))
        enddo
      close(65)
c
      tchar3 = dir_name(1:namelength)//"/bloch_vec_out.txt"
      open(3,file=tchar3)
        write(3,*) bloch_vec(1)/pi, "  ! Kx/pi"
        write(3,*) bloch_vec(2)/pi, "  ! Ky/pi"
      close(3)
c
12    format(6(g25.17) )
c
      return
      end

