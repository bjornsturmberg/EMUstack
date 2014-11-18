
      subroutine calc_modes_1d(
c     Explicit inputs
     *    lambda, nval, ordre_ls, nb_typ_el,
     *    npt_P2, nel, table_nod, type_el, x_P2,
     *    itermax, debug, mesh_file, n_eff,
     *    bloch_vec_x, bloch_vec_y, shift,
     *    plot_modes, plot_real, plot_imag, plot_abs,
     *    neq_PW, neq_PW_2d, world_1d,
C     Outputs
     *     beta_1, overlap_J, overlap_J_dagger, overlap_J_2d,
     *     overlap_J_dagger_2d, sol_P2)

C************************************************************************
C
C  Program:
C    Finite Element Method for 1D array of planar waveguides
C
C************************************************************************
C
      implicit none

      integer*8 nel, nb_typ_el, world_1d
      integer*8 int_max, cmplx_max
      integer*8 real_max, n_64

      integer*8, dimension(:), allocatable :: a
      complex*16, dimension(:), allocatable :: b
      double precision, dimension(:), allocatable :: c
      integer allocate_status
      integer*8 npt_P2, n_ddl
      integer*8 n_ddl_P2, n_ddl_P3
      integer*8  nnodes_P2, ui
      parameter(nnodes_P2=3)
      integer*8 type_el(nel), table_nod(3,nel)
      integer*8, allocatable :: table_ddl(:,:), ineq(:)
      integer*8, allocatable :: ip_period_ddl(:)

      double precision x_P2(npt_P2)
      double precision, allocatable :: x_ddl(:)

      complex*16 shift
      complex*16, dimension(:,:), allocatable :: matrix_1, matrix_2
      complex*16, dimension(:,:), allocatable :: overlap_L
      complex*16, allocatable :: pp(:), qq(:)
      complex*16, allocatable :: vp(:,:), v(:,:)
      complex*16, allocatable :: mode_pol(:,:)
      complex*16, target :: beta_1(nval), beta_2(nval)
      complex*16, pointer :: beta(:)

      complex*16, target :: sol_1(3+4+4,nval,nel), sol_2(3+4+4,nval,nel)
      complex*16, pointer :: sol(:,:,:)
      complex*16 sol_P2(3,nnodes_P2,nval,nel)
      complex*16 eps_eff(nb_typ_el), n_eff(nb_typ_el)

      complex*16 overlap_J(2*neq_PW, nval)
      complex*16 overlap_J_dagger(nval, 2*neq_PW)
      complex*16 overlap_J_2d(2*neq_PW_2d, nval)
      complex*16 overlap_J_dagger_2d(nval, 2*neq_PW_2d)

      complex*16, allocatable :: beta_PW(:)
      integer*8, allocatable :: index_PW(:)

      integer*8 nval, nvect, n_conv
      double precision k_0, pi, bloch_vec_x, bloch_vec_y
      double precision bloch_vec_x_k, bloch_vec_y_k
      double precision n_eff_0, period_x
      double precision ls_data(10), lambda
c     Normalised wavelength lambda
      integer*8 neq, E_H_field
C  Plane wave parameters
      integer*8 neq_PW, neq_PW_2d, ordre_ls
      integer*8, allocatable :: index_pw_inv(:)
      integer*8 i, n_k
C  Variable used by valpr
      integer*8 itermax, n_core(2)
      integer*8, dimension(:), allocatable :: index
      complex*16 z_beta, z_tmp, z_tmp0
      double precision freq, tol
C  Timing variables
      double precision time1_fact, time2_fact
      double precision time1_postp
      double precision time1_arpack, time2_arpack
      double precision time1_J, time2_J
C  Names and Controls
      character*100 mesh_file
      character*100 overlap_file, dir_name
      character*100 tchar
      integer*8 debug, plot_modes, pair_warning, homogeneous_check
      integer*8 q_average, plot_real, plot_imag, plot_abs

Cf2py intent(in) lambda, nval, ordre_ls, nb_typ_el
Cf2py intent(in) npt_P2, nel, table_nod, type_el, x_P2,
Cf2py intent(in) itermax, debug, mesh_file, n_eff,
Cf2py intent(in) bloch_vec_x, bloch_vec_y, shift,
Cf2py intent(in) plot_modes, plot_real, plot_imag, plot_abs,
Cf2py intent(in) neq_PW, neq_PW_2d,

Cf2py depend(n_eff) nb_typ_el
Cf2py depend(type_el) nel
Cf2py depend(table_nod) nel
Cf2py depend(x_P2) npt_P2

Cf2py intent(out) beta_1, overlap_J, overlap_J_dagger
Cf2py intent(out) overlap_J_2d, overlap_J_dagger_2d
Cf2py intent(out) sol_P2

c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (debug .eq. 1) then
        write(*,*) "WELCOME TO DEBUG MODE"
      endif

      ui = 6
      pi = 3.141592653589793d0

      nvect = 2*nval + nval/2 +3
      E_H_field = 1
c     E_H_field = 1 => Electric field formulation (E-Field)
c     E_H_field = 2 => Magnetic field formulation (H-Field)

      n_64 = 2
      real_max = n_64**1
      int_max  = n_64**1
      cmplx_max = n_64**1

      tol = 0.0 ! ARPACK accuracy (0.0 for machine precision)
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
      allocate_status = 0

      allocate(b(cmplx_max), STAT=allocate_status)
      if (allocate_status /= 0) then
        write(*,*) "The allocation is unsuccessful"
        write(*,*) "allocate_status = ", allocate_status
        write(*,*) "Not enough memory for the complex array b"
        write(*,*) "cmplx_max = ", cmplx_max
        write(*,*) "Aborting..."
        stop
      endif

      allocate(c(real_max), STAT=allocate_status)
      if (allocate_status /= 0) then
        write(*,*) "The allocation is unsuccessful"
        write(*,*) "allocate_status = ", allocate_status
        write(*,*) "Not enough memory for the real array c"
        write(*,*) "real_max = ", real_max
        write(*,*) "Aborting..."
        stop
      endif

      allocate(a(int_max), STAT=allocate_status)
      if (allocate_status /= 0) then
        write(*,*) "The allocation is unsuccessful"
        write(*,*) "allocate_status = ", allocate_status
        write(*,*) "Not enough memory for the integer array a"
        write(*,*) "int_max = ", int_max
        write(*,*) "Aborting..."
        stop
      endif

      allocate(overlap_L(nval,nval), STAT=allocate_status)
      if (allocate_status /= 0) then
        write(*,*) "The allocation is unsuccessful"
        write(*,*) "allocate_status = ", allocate_status
        write(*,*) "Not enough memory for overlap_L"
        write(*,*) "nval = ", nval
        write(*,*) "Aborting..."
        stop
      endif
      allocate(index(nval), STAT=allocate_status)
      if (allocate_status /= 0) then
        write(*,*) "The allocation is unsuccessful"
        write(*,*) "allocate_status = ", allocate_status
        write(*,*) "Not enough memory for index"
        write(*,*) "nval = ", nval
        write(*,*) "Aborting..."
        stop
      endif

c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if( npt_P2 /= (2 * nel + 1)) then
        write(ui,*) "npt_P2 /= (2 * nel + 1) : ",
     *     npt_P2, (2 * nel + 1)
         write(ui,*) "py_calc_modes_1d: Aborting..."
         stop
      endif

      if (debug .eq. 1) then
        write(*,*) "py_calc_modes_1d: nel = ", nel
        write(*,*) "py_calc_modes_1d: npt_P2 = ", npt_P2
      endif

      n_ddl_P2 = 3 * nel ! Discontinuous P2 polynomial
      n_ddl_P3 = 3 * nel + 1
      n_ddl = n_ddl_P2 + 2*n_ddl_P3
      neq = n_ddl_P2 + 2 * (n_ddl_P3 -1)

      if (debug .eq. 1) then
        write(*,*) "py_calc_modes_1d: n_ddl = ", n_ddl
        write(*,*) "py_calc_modes_1d: neq = ", neq
      endif

C      call geometry_1d (nel, npt_P2, nnodes_P2, nb_typ_el,
C     *     1.0, type_el, table_nod,
C     *     x_P2, mesh_file)
C
C      write(ui,*) 'nel', nel
C      write(ui,*) 'npt_P2', npt_P2
C      write(ui,*) 'nb_typ_el', nb_typ_el
C      write(ui,*) 'type_el', type_el
C      write(ui,*) 'table_nod', table_nod
C      write(ui,*) 'x_P2', x_P2

c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c

      allocate(table_ddl(3+4+4,nel), ineq(n_ddl), STAT=allocate_status)
      if (allocate_status /= 0) then
        write(*,*) "The allocation is unsuccessful"
        write(*,*) "allocate_status = ", allocate_status
        write(*,*) "Not enough memory for table_ddl and ineq"
        write(*,*) "nel, n_ddl = ", nel, n_ddl
        write(*,*) "Aborting..."
        stop
      endif

      allocate(matrix_1(neq,neq), matrix_2(neq,neq),
     *     STAT=allocate_status)
      if (allocate_status /= 0) then
        write(*,*) "The allocation is unsuccessful"
        write(*,*) "allocate_status = ", allocate_status
        write(*,*) "Not enough memory for matrix_1 and matrix_2"
        write(*,*) "neq = ", neq
        write(*,*) "Aborting..."
        stop
      endif

      allocate(ip_period_ddl(n_ddl), STAT=allocate_status)
      if (allocate_status /= 0) then
        write(*,*) "The allocation is unsuccessful"
        write(*,*) "allocate_status = ", allocate_status
        write(*,*) "Not enough memory for ip_period_ddl"
        write(*,*) "n_ddl = ", n_ddl
        write(*,*) "Aborting..."
        stop
      endif

      allocate(pp(nb_typ_el), qq(nb_typ_el), STAT=allocate_status)
      if (allocate_status /= 0) then
        write(*,*) "The allocation is unsuccessful"
        write(*,*) "allocate_status = ", allocate_status
        write(*,*) "Not enough memory for pp and qq"
        write(*,*) "nb_typ_el = ", nb_typ_el
        write(*,*) "Aborting..."
        stop
      endif

      allocate(x_ddl(n_ddl), STAT=allocate_status)
      if (allocate_status /= 0) then
        write(*,*) "The allocation is unsuccessful"
        write(*,*) "allocate_status = ", allocate_status
        write(*,*) "Not enough memory for x_ddl"
        write(*,*) "n_ddl = ", n_ddl
        write(*,*) "Aborting..."
        stop
      endif


      allocate(vp(neq,nval), v(neq,nvect), STAT=allocate_status)
      if (allocate_status /= 0) then
        write(*,*) "The allocation is unsuccessful"
        write(*,*) "allocate_status = ", allocate_status
        write(*,*) "Not enough memory for vp, "
        write(*,*) "neq, nval, nvect = ", neq, nval, nvect
        write(*,*) "Aborting..."
        stop
      endif

      allocate(mode_pol(4,nval), STAT=allocate_status)
      if (allocate_status /= 0) then
        write(*,*) "The allocation is unsuccessful"
        write(*,*) "allocate_status = ", allocate_status
        write(*,*) "Not enough memory for mode_pol"
        write(*,*) "nval = ", nval
        write(*,*) "Aborting..."
        stop
      endif

      allocate(index_pw_inv(neq_PW), STAT=allocate_status)
      if (allocate_status /= 0) then
        write(*,*) "The allocation is unsuccessful"
        write(*,*) "allocate_status = ", allocate_status
        write(*,*) "Not enough memory for index_pw_inv"
        write(*,*) "neq_PW = ", neq_PW
        write(*,*) "Aborting..."
        stop
      endif

      allocate(beta_PW(2*nval+1), index_PW(2*nval+1),
     *     STAT=allocate_status)
      if (allocate_status /= 0) then
        write(*,*) "The allocation is unsuccessful"
        write(*,*) "allocate_status = ", allocate_status
        write(*,*) "Not enough memory for beta_PW"
        write(*,*) "neq, nval, nvect = ", nval
        write(*,*) "Aborting..."
        stop
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (debug .eq. 1) then
        write(*,*) "mesh_file = ", mesh_file
      endif

      call periodic_cond_1d (nel, npt_P2, n_ddl, neq, table_nod,
     *      ineq, x_P2, table_ddl, ip_period_ddl, x_ddl, debug,
     *      period_x)
C
      if(debug .eq. 1) then
        write(ui,*) "period_x = ", period_x
      endif
C
      write(ui,*)
      write(ui,*) "---------------------------------------",
     *     "-------"
      write(ui,*) " Wavelength : ", lambda, " (d)"
      write(ui,*) "---------------------------------------",
     *     "-------"
      write(ui,*)
C
c     Calculate effective permittivity
      do i=1,nb_typ_el
        eps_eff(i) = n_eff(i)**2
      enddo

        freq = 1.0d0/lambda
        n_eff_0 = DBLE(n_eff(1))
        k_0 = 2.0d0 * pi * freq

      if(debug .eq. 1) then
        write(ui,*) "lambda = ", lambda
        write(ui,*) "k_0 = ", k_0
        write(ui,*) "n_eff_0 = ", n_eff_0
        write(ui,*) "bloch_vec_x = ", bloch_vec_x
        write(ui,*) "bloch_vec_y = ", bloch_vec_y
      endif
C
C  Index number of the core materials (material with highest Re(eps_eff))
        n_core(1) = 2
        do i=1,nb_typ_el
          if(dble(eps_eff(i)) .gt. dble(eps_eff(n_core(1)))) then
            n_core(1) = i
          endif
        enddo
C  Check that the layer is not in fact homogeneous
        homogeneous_check = 0
        do i=1,nb_typ_el-1
          if(dble(eps_eff(i)) .ne. dble(eps_eff(i+1))) then
            homogeneous_check = 1
          elseif(dimag(eps_eff(i)) .ne. dimag(eps_eff(i+1))) then
            homogeneous_check = 1
          endif
        enddo
        if(homogeneous_check .eq. 0) then
          write(ui,*) "py_calc_modes_1d.f: ",
     *              "FEM routine cannot handle homogeneous layers."
          write(ui,*) "Define layer as object.ThinFilm"
          write(ui,*) "Aborting..."
          stop
        endif
        n_core(2) = n_core(1)
        if(debug .eq. 1) then
          write(ui,*) "n_core = ", n_core
          write(ui,*) "shift = ", shift
          if(E_H_field .eq. 1) then
            write(ui,*) "E-Field formulation"
          else
            write(ui,*) "H-Field formulation"
          endif
        EndIf
C
C
      if(E_H_field .eq. 1) then
        do i=1,nb_typ_el
          qq(i) = eps_eff(i)*k_0**2
          pp(i) = 1.0d0
        enddo
      elseif(E_H_field .eq. 2) then
        do i=1,nb_typ_el
          qq(i) = k_0**2
          pp(i) = 1.0d0/eps_eff(i)
        enddo
      else
        write(ui,*) "py_calc_modes_1d.f: ",
     *              "action indef. avec E_H_field = ",
     *                  E_H_field
        write(ui,*) "Aborting..."
        stop
      endif
C
CCCCCCCCCCCCCCCCCCCC  Loop over Adjoint and Prime  CCCCCCCCCCCCCCCCCCCCCC
C
      do n_k = 1,2
C
      if (n_k .eq. 1) then
        sol => sol_1
        beta => beta_1
        bloch_vec_x_k = bloch_vec_x
        bloch_vec_y_k = bloch_vec_y
      else
        sol => sol_2
        beta => beta_2
        bloch_vec_x_k = -bloch_vec_x
        bloch_vec_y_k = -bloch_vec_y
      endif
C
C     Assemble the coefficient matrix A and the right-hand side F of the
C     finite element equations
      if (debug .eq. 1) then
        write(ui,*) "py_calc_modes_1d.f: Asmbly: call to asmbly_1d"
      endif
C
      call asmbly_1d (nel, npt_P2, n_ddl, neq,
     *  shift, bloch_vec_x_k, bloch_vec_y_k, nb_typ_el,
     *  pp, qq, table_nod, table_ddl, type_el, ineq,
     *  ip_period_ddl, x_P2, x_ddl, matrix_1, matrix_2)


      if (debug .eq. 1) then
        write(ui,*) "py_calc_modes_1d.f: Adjoint(1) / Prime(2)", n_k
      endif

      if (debug .eq. 1) then
        write(ui,*) "py_calc_modes_1d.f: call to valpr_64_1d"
      endif
      call valpr_64_1d (nvect, nval, neq, itermax,
     *  tol, matrix_1, matrix_2, v, beta, vp,
     *  n_conv, ls_data, debug)
C
c
      if (n_conv .ne. nval) then
         write(ui,*) "py_calc_modes_1d.f: convergence problem with "
     *               // "valpr_64_1d"
         write(ui,*) "py_calc_modes_1d.f: n_conv != nval : ",
     *    n_conv, nval
         write(ui,*) "py_calc_modes_1d.f: Aborting..."
         stop
      endif
c
      time1_fact = ls_data(1)
      time2_fact = ls_data(2)
c
      time1_arpack = ls_data(3)
      time2_arpack = ls_data(4)
C
      do i=1,nval
        z_tmp0 = beta(i)
        z_tmp = 1.0d0/z_tmp0+shift
        z_beta = sqrt(z_tmp)
C       Mode classification - we want the forward propagating mode
        if (abs(imag(z_beta)) .lt. 1.0d-8) then
C         re(z_beta) > 0 for forward propagating mode
          if (dble(z_beta) .lt. 0) z_beta = -z_beta
        else
C         im(z_beta) > 0 for forward decaying evanescent mode
          if (imag(z_beta) .lt. 0) z_beta = -z_beta
        endif
C     !  Effective index
C        z_beta = sqrt(z_tmp)/k_0
        beta(i) = z_beta
      enddo
c
      call cpu_time(time1_postp)
C
      call z_indexx (nval, beta, index)
c
C       The eigenvectors will be stored in the array sol
C       The eigenvalues and eigenvectors will be renumbered
C                 using the permutation vector index
      call array_sol_1d (nval, nel, n_ddl, neq,
     *     n_core, bloch_vec_x_k, index,
     *     table_ddl, type_el, ineq,
     *     ip_period_ddl, x_ddl,
     *     beta, mode_pol, vp, sol, n_k)

ccc      if (n_k == 1) then
ccc        call array_sol_1d (nval, nel, n_ddl, neq,
ccc     *     n_core, bloch_vec_x_k, index,
ccc     *     table_ddl, type_el, ineq,
ccc     *     ip_period_ddl, x_ddl,
ccc     *     beta, mode_pol, vp, sol_1)
ccc      else
ccc        call array_sol_1d (nval, nel, n_ddl, neq,
ccc     *     n_core, bloch_vec_x_k, index,
ccc     *     table_ddl, type_el, ineq,
ccc     *     ip_period_ddl, x_ddl,
ccc     *     beta, mode_pol, vp, sol_2)
ccc      endif
C
      if(debug .eq. 1) then
        write(ui,*) 'index = ', (index(i), i=1,nval)
      endif
      if(debug .eq. 1) then
        write(ui,*)
        write(ui,*) "lambda, 1/lambda = ", lambda, 1.0d0/lambda
        write(ui,*) bloch_vec_x_k/(2.0d0*pi), bloch_vec_y_k/(2.0d0*pi)
        write(ui,*) "sqrt(shift) = ", sqrt(shift)
        do i=1,nval
          write(ui,"(i4,2(g22.14),2(g18.10))") i,
     *       beta(i)
        enddo
      endif
C
      enddo
C
CCCCCCCCCCCCCCCCCCCCCCCC  End Prime, Adjoint Loop  CCCCCCCCCCCCCCCCCCCCCC
C
C  Orthogonal integral
      pair_warning = 0
      if (debug .eq. 1) then
        write(ui,*) "py_calc_modes_1d.f: Field product"
      endif
      overlap_file = "Normed/Orthogonal.txt"
      call orthogonal_1d (nval, nel, npt_P2,
     *  nb_typ_el, pp, bloch_vec_y, table_nod,
     *  type_el, x_P2, beta_1, beta_2,
     *  sol_1, sol_2, overlap_L, overlap_file, debug,
     *  pair_warning, k_0)
C
      if (pair_warning .ne. 0 .and. nval .le. 20) then
        write(ui,*) "py_calc_modes_1d.f: Warning found 1 BM
     * of cmplx conj"
        write(ui,*) "pair, increase num_BMs to include the other."
      endif
C
C  Normalisation
      if(debug .eq. 1) then
        write(ui,*) "py_calc_modes_1d.f: Field  Normalisation"
      endif
      call normalisation_1d (nval, nel, sol_1, sol_2, overlap_L)
C
C    Save Original solution
      if (plot_modes .eq. 1) then
        call array_sol_P2_1d (nval, nel, sol_1, sol_P2)
        dir_name = "Bloch_fields"
C        call write_sol_P2_1d (nval, nel, E_H_field,
C     *     lambda, beta_1, sol_P2, mesh_file, dir_name)
C         q_average = 0
C         tchar = "Bloch_fields/PNG/All_plots_png_abs2_eE.geo"
C         open (unit=34,file=tchar)
C         do i=1,nval
C           call gmsh_post_process_1d (i, E_H_field, nval,
C      *       nel, npt_P2, table_nod, type_el, nb_typ_el,
C      *       n_eff, x_P2, beta_1, sol_P2,
C      *       mesh_file, dir_name,
C      *       q_average, plot_real, plot_imag, plot_abs)
C         enddo
C         close (unit=34)
      endif
C
C  Orthonormal integral
      if (debug .eq. 1) then
        write(ui,*) "py_calc_modes_1d.f: Product of normalised field"
        overlap_file = "Normed/Orthogonal_n.txt"
        call cpu_time(time1_J)
        call orthogonal_1d (nval, nel, npt_P2,
     *    nb_typ_el, pp, bloch_vec_y, table_nod,
     *    type_el, x_P2, beta_1, beta_2,
     *    sol_1, sol_2, overlap_L, overlap_file, debug,
     *    pair_warning, k_0)
        call cpu_time(time2_J)
          write(ui,*) "py_calc_modes_1d.f: CPU time for orthogonal :",
     *    (time2_J-time1_J)
      endif
C
C  Plane wave ordering
      if (debug .eq. 1) then
        write(ui,*) "py_calc_modes_1d.f: pw_ordering_1d"
      endif
      call pw_ordering_1d (neq_PW, period_x,
     *  bloch_vec_x, bloch_vec_y, index_pw_inv,
     *  debug, ordre_ls, k_0)
C
C  J_overlap
      if (debug .eq. 1) then
        write(ui,*) "py_calc_modes_1d.f: J_overlap_1d Integral"
      endif
      call J_overlap_1d(nval, nel, npt_P2,
     *  type_el, table_nod, x_P2, sol_1,
     *  period_x, lambda, freq, overlap_J, neq_PW,
     *  bloch_vec_x, bloch_vec_y, index_pw_inv, debug,
     *  ordre_ls)
C
C  J_dagger_overlap
      if (debug .eq. 1) then
        write(ui,*) "py_calc_modes_1d.f: J_dagger_overlap_1d Integral"
      endif
      call J_dagger_overlap_1d(nval, nel, npt_P2,
     *  type_el, table_nod, x_P2, sol_2,
     *  period_x, lambda, freq, overlap_J_dagger, neq_PW,
     *  bloch_vec_x, bloch_vec_y, index_pw_inv, debug, ordre_ls)
C
C  Convert overlap integrals from 1D diffraction order basis to
C  2D diffraction order basis (filling with zeros).
      if (world_1d .eq. 0) then
        if (debug .eq. 1) then
          write(ui,*)"py_calc_modes_1d.f: pw_matrix_1d_to_2d conversion"
        endif
        call pw_matrix_1d_to_2d (neq_PW, neq_PW_2d, nval,
     *    period_x, bloch_vec_x, bloch_vec_y, index_pw_inv,
     *    debug, ordre_ls, overlap_J, overlap_J_2d,
     *    overlap_J_dagger, overlap_J_dagger_2d,
     *    k_0, debug)
      endif

c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    Plane wave diffraction orders
cc      do i=-nval,nval
cc        z_tmp0 = bloch_vec_x + i * 2.0d0 * pi / period_x
cc        z_tmp = (k_0 * n_eff_0)**2 - z_tmp0**2 - bloch_vec_y**2
cc        z_beta = sqrt(z_tmp)
cc        beta_PW(i+nval+1) = z_beta
cc      enddo
cc      call z_indexx (2*nval+1, beta_PW, index_PW)
ccc      if(debug .eq. 1) then
ccc        write(ui,*)
ccc        do i=1,2*nval+1
ccc          j = index_PW(i)
ccc          write(ui,*) "PW", i, beta_PW(j)
ccc        enddo
ccc      endif
cc      if(debug .eq. 1) then
cc        open (unit=24,file="Output/beta_PW.txt",status="unknown")
cc          write(24,*) "k_0, n_eff_0 = ", k_0, n_eff_0
cc          write(24,*) "bloch_vec_x = ", bloch_vec_x
cc          write(24,*) "bloch_vec_y = ", bloch_vec_y
cc          write(24,*) "period_x = ", period_x
cc        do i=1,nval
cc          j = index_PW(i)
cc          write(24,*)  i, beta_PW(j), beta_1(i)
cc        enddo
cc        write(24,*)
cc        do i=nval+1,2*nval+1
cc          j = index_PW(i)
cc          write(24,*)  i, beta_PW(j)
cc        enddo
cc        write(24,*)
cc        do i=1,nval
cc          j = index_PW(i)
cc          write(24,*)  i,
cc     *  sqrt((k_0 * n_eff_0)**2 - beta_PW(j)**2
cc     *  - bloch_vec_y**2),
cc     *  sqrt((k_0 * n_eff_0)**2 - beta_1(i)**2
cc     *  - bloch_vec_y**2)
ccc          write(24,*)  i, beta_PW(j)**2, beta_1(i)**2
cc        enddo
cc        write(24,*)
cc        do i=nval+1,2*nval+1
cc          j = index_PW(i)
cc          write(24,*)  i,
cc     *  sqrt((k_0 * n_eff_0)**2 - beta_PW(j)**2
cc     *  - bloch_vec_y**2)
ccc          write(24,*)  i, beta_PW(j)**2
cc        enddo
cc      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if(debug .eq. 1000) then
        write(ui,*)
        do i=1,nval
          write(ui,*) "b=", i, beta_1(i)
        enddo
      endif
C
C     factorization of the globale matrice
C     -----------------------------------
C
c
c
c      write(*,*) "nel, npt_P2, n_ddl, neq, nnodes_P2 = ",
c     *      nel, npt_P2, n_ddl, neq, nnodes_P2
c      write(*,*) "shift, k_0, bloch_vec_x, nb_typ_el = ",
c     *      shift, k_0, bloch_vec_x, nb_typ_el
c
c      do i=1,10
c        write(*,*) i, x_P2(i), type_nod(i)
c      enddo
c      do i=1,10
c        write(*,"(8(I7))") i, (table_nod(j,i),j=1,nnodes_P2), type_el(i)
c      enddo
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
      deallocate(a,b,c,index,overlap_L)
c
      end subroutine calc_modes_1d
