c
c     E_H_field = 1 => Electric field formulation (E-Field)
c     E_H_field = 2 => Magnetic field formulation (H-Field)
c
c     i_cond = 0 => Dirichlet boundary condition
c     i_cond = 1 => Neumann boundary condition
c     i_cond = 2 => Periodic boundary condition
c
      subroutine get_param(E_H_field, n_lambda, lambda_1, 
     *      lambda_2, npt, nel, i_cond, nval, nvect, itermax,
     *      tol, lx, ly, mesh_file, mesh_format,
     *      gmsh_file, gmsh_file_pos, log_file, Lambda_Res,
     *      d_in_nm, debug)
c
c----------------------------------------------------------
c---- read data (lecture des donnees)
c----------------------------------------------------------
c
      implicit none
      integer*8 E_H_field, npt, nel, i_cond, d_lat_vec
      integer*8 nval, nvect, itermax, n_lambda, Lambda_Res
      double precision lambda_1, lambda_2, tol, pi, lx, ly
      integer*8 mesh_format, namelen, debug, d_in_nm
      character mesh_file*100, gmsh_file*100, log_file*100
      character gmsh_file_pos*100 ! Post-processing for gmsh
C
      pi = 3.141592653589793d0
      open(3,file='Parameters/param.txt',status='old')
        read(3,*) E_H_field
        read(3,*) n_lambda
        read(3,*) lambda_1, lambda_2
        read(3,*) Lambda_Res    ! sets predefined wavelength resolution
        read(3,*) d_in_nm    ! number of lattice vectors d to calculate
        read(3,*) i_cond
        read(3,*) lx            ! radius (stretching factor) in direction x
        read(3,*) ly            ! radius in direction y
        read(3,*) nval
        read(3,*) itermax
        read(3,*) tol
        read(3,'(a100)') mesh_file
        read(3,*) mesh_format
      close(3)

C
CCCCCCCCCCCCCCC 
C  Convert to correct settings for irradiance_1 and refrac_index_1
      if (Lambda_Res .eq. 1) then
        n_lambda = 818
        lambda_1 = 310
        lambda_2 = 1127
C  Convert to correct settings for irradiance_2 and refrac_index_2
      elseif (Lambda_Res .eq. 2) then
        n_lambda = 1635
        lambda_1 = 310
        lambda_2 = 1127
      endif
C
CCCCCCCCCCCCCCCC
C
      nvect = 2*nval + nval/2 + 2
C
      namelen = len_trim(mesh_file)
      gmsh_file = 'Normed/'//mesh_file(1:namelen)//'.msh'
      gmsh_file_pos = mesh_file(1:namelen)
      log_file = 'Normed/'//mesh_file(1:namelen)//'.log'
      if (debug .eq. 1) then
        write(*,*) "get_param: mesh_file = ", mesh_file
        write(*,*) "get_param: gmsh_file = ", gmsh_file
      endif
c
      open (unit=24,file="../Mesh/"//mesh_file,status='unknown')
        read(24,*) npt, nel
      close(24)
C
      return
      end

