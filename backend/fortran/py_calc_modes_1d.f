
      subroutine calc_modes_1d(
c     Explicit inputs
     *    lambda, nval, ordre_ls, nb_typ_el,
     *    npt_P2, nel, itermax, debug, mesh_file,
     *    d_in_nm,
     *    plot_modes, plot_real, plot_imag, plot_abs,
     *    neq_PW,
C     Outputs
     *     beta_1, overlap_J, overlap_J_dagger, sol_1, sol_2)



C************************************************************************
C
C  Program:
C    Finite Element Method for 1D array of planear waveguides
C
C************************************************************************
C
      implicit none

      integer*8 nel, nb_typ_el
C  Local parameters:
      integer*8 int_max, cmplx_max, int_used, cmplx_used
      integer*8 real_max, real_used, n_64

      integer*8, dimension(:), allocatable :: a
      complex*16, dimension(:), allocatable :: b
      double precision, dimension(:), allocatable :: c
      integer allocate_status
      integer*8 npt_P2, npt_P3, n_ddl
      integer*8 n_ddl_P2, n_ddl_P3
      integer*8  nnodes_P2, ui
      parameter(nnodes_P2=3)

      integer*8, allocatable :: type_nod(:), type_el(:), table_nod(:, :)
      integer*8, allocatable :: table_ddl(:,:), ineq(:)
      integer*8, allocatable :: ip_period_ddl(:)

      double precision, allocatable :: x_P2(:), x_ddl(:)

      complex*16, dimension(:,:), allocatable :: matrix_1, matrix_2
      complex*16, dimension(:,:), allocatable :: overlap_L
      complex*16, allocatable :: pp(:), qq(:)
      complex*16, allocatable :: vp(:,:), v(:,:)
      complex*16, allocatable :: mode_pol(:,:)
      complex*16, target :: beta_1(nval), beta_2(nval)
      complex*16, pointer :: beta(:)

ccc      complex*16, allocatable :: sol_1(:,:,:), sol_2(:,:,:)
      complex*16, target :: sol_1(3+4+4,nval,nel), sol_2(3+4+4,nval,nel)
      complex*16, pointer :: sol(:,:,:)
c      complex*16, allocatable, target :: 
      complex*16, allocatable :: sol_P2(:,:,:,:)
      complex*16, allocatable :: eps_eff(:), n_eff(:)
      complex*16 overlap_J(2*neq_PW, nval)
      complex*16 overlap_J_dagger(nval, 2*neq_PW)
      complex*16, allocatable :: beta_PW(:)
      integer*8, allocatable :: index_PW(:)

      integer*8 nval, nvect, n_conv
      double precision k_0, pi, lx, bloch_vec_x, bloch_vec_y
      double precision bloch_vec_x_k, bloch_vec_y_k
      double precision n_eff_0, period_x
      double precision ls_data(10)

c     Wavelength lambda is in normalised units of d_in_nm
      double precision lambda, d_in_nm

      integer*8 i_lambda, n_lambda
      double precision lambda_1, lambda_2, d_lambda, d_freq

      complex*16 shift

      character mesh_file*100

      integer*8 neq
c     E_H_field = 1 => Electric field formulation (E-Field)
c     E_H_field = 2 => Magnetic field formulation (H-Field)
      integer*8 E_H_field

C
C  Declare the pointers of the integer super-vector
      integer*8 ip_table_E, ip_table_N_E_F, ip_visite
      integer*8 ip_type_N_E_F, ip_eq
      integer*8 ip_period_N, ip_nperiod_N
      integer*8 ip_period_N_E_F, ip_nperiod_N_E_F
C      integer*8 ip_col_ptr, ip_bandw 
C  Declare the pointers of the real super-vector
      integer*8 jp_x_N_E_F
C      integer*8 jp_matD, jp_matL, jp_matU
C      integer*8 jp_matD2, jp_matL2, jp_matU2
      integer*8 jp_vect1, jp_vect2, jp_workd, jp_resid, jp_vschur
      integer*8 jp_trav, jp_vp
C  Plane wave parameters
      integer*8 neq_PW, nx_PW, ordre_ls
      integer*8, allocatable ::  index_pw_inv(:)


      integer*8 i, j, n_k

c
c
c   , test
c     i_cond = 0 => Dirichlet boundary condition
c     i_cond = 1 => Neumann boundary condition
c     i_cond = 2 => Periodic boundary condition
C     ! Number of nodes per element
c , i_cond
c
C, len_skyl, nsym
c      integer*8 numberprop_N
C  Variable used by valpr
      integer*8 itermax, ltrav
      integer*8 i_base
c      double precision ls_data(10)
c      integer*8 pointer_int(20), pointer_cmplx(20)
C      integer*8 index(2000), n_core(2)
      integer*8, dimension(:), allocatable :: index
      integer*8 n_core(2)
      complex*16 z_beta, z_tmp, z_tmp0
c      integer*8 n_edge, n_face, n_ddl, n_ddl_max
c     variable used by UMFPACK
c      double precision control (20), info_umf (90)
c      integer*8 numeric, status, filenum
C  Renumbering
c      integer*8 ip_row_ptr, ip_bandw_1, ip_adjncy
c      integer*8 len_adj, len_adj_max, len_0_adj_max
c, iout, nonz_1, nonz_2

c    !, mesh_format

      double precision freq, lat_vecs(2,2), tol
c      double precision bloch_vec(2), bloch_vec_k(2)
c
C  Timing variables
      double precision time1, time2
      double precision time1_fact, time2_fact
      double precision time1_asmbl, time2_asmbl
      double precision time1_postp
      double precision time1_arpack, time2_arpack
      double precision time1_J, time2_J
      character*(8) start_date, end_date
      character*(10) start_time, end_time
C  Names and Controls
      character gmsh_file*100, log_file*100
      character gmsh_file_pos*100
      character overlap_file*100, dir_name*100
      character*100 tchar
      integer*8 namelength, debug
      integer*8 plot_modes
      integer*8 pair_warning
      integer*8 q_average, plot_real, plot_imag, plot_abs

c     Declare the pointers of the real super-vector
      integer*8 kp_rhs_re, kp_rhs_im, kp_lhs_re, kp_lhs_im
      integer*8 kp_mat1_re, kp_mat1_im

c     Declare the pointers of for sparse matrix storage
      integer*8 ip_col_ptr, ip_row
      integer*8 jp_mat2
      integer*8 ip_work, ip_work_sort, ip_work_sort2
      integer*8 nonz, nonz_max, max_row_len

      integer*8 ip
      integer i_32

c     new breed of variables to prise out of a, b and c
c      complex*16 x_arr(2,npt)
c      complex*16, target :: sol1(3,nnodes+7,nval,nel)
c      complex*16, target :: sol2(3,nnodes+7,nval,nel)
c      complex*16, pointer :: sol(:,:,:,:)
c      complex*16, target :: beta_1(nval), beta_2(nval)
c      complex*16, pointer :: beta(:)


Cf2py intent(in) lambda, nval, ordre_ls, neq_PW, nb_typ_el
Cf2py intent(in) npt_P2, nel, itermax, debug, mesh_file

Cf2py intent(out) beta_1, overlap_J, overlap_J_dagger
Cf2py intent(out) sol_1, sol_2

c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
      ui = 6
      pi = 3.141592653589793d0

      nvect = 2*nval + nval/2 +3
      tol = 0.0d0
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
      n_64 = 2
      real_max = n_64**1
      int_max  = n_64**1
      cmplx_max = n_64**1
c      3*npt+nel+nnodes*nel 

C      write(*,*) "cmplx_max = ", cmplx_max
C      write(*,*) "real_max = ", real_max
C      write(*,*) "int_max = ", int_max

      tol = 0.0 ! ARPACK accuracy (0.0 for machine precision)
      lx = 1 ! Diameter of unit cell. Default, lx = 1.0.

C     Old inputs now internal to here and commented out by default.
C      mesh_format = 1

      if (debug .eq. 1) then
        write(*,*) "WELCOME TO DEBUG MODE"
      endif
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
      allocate(eps_eff(nb_typ_el), n_eff(nb_typ_el), 
     *           STAT=allocate_status)
      if (allocate_status /= 0) then
        write(*,*) "The allocation is unsuccessful"
        write(*,*) "allocate_status = ", allocate_status
        write(*,*) "Not enough memory for eps_eff, n_eff"
        write(*,*) "nb_typ_el = ", nb_typ_el
        write(*,*) "Aborting..."
        stop
      endif

c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c

CC  tmp shortened tp -4, was -5.
CC  change back if using .mail msh

C     clean mesh_format
      namelength = len_trim(mesh_file)
      gmsh_file = mesh_file(1:namelength-4)//'.msh'
      gmsh_file_pos = mesh_file(1:namelength)
      log_file = mesh_file(1:namelength-4)//'.log'
      if (debug .eq. 1) then
        write(*,*) "mesh_file = ", mesh_file
        write(*,*) "gmsh_file = ", gmsh_file
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
      npt_P3 = 3 * nel + 1
      lx = 1

      if (debug .eq. 1) then
        write(*,*) "py_calc_modes_1d: nel = ", nel
        write(*,*) "py_calc_modes_1d: npt_P2 = ", npt_P2
        write(*,*) "py_calc_modes_1d: npt_P3 = ", npt_P3
      endif

      allocate(type_nod(npt_P2), type_el(nel), x_P2(npt_P2),
     *     table_nod(3, nel), STAT=allocate_status)

      n_ddl_P2 = 3 * nel ! Discontinuous P2 polynomial
      n_ddl_P3 = 3 * nel + 1
      n_ddl = n_ddl_P2 + 2*n_ddl_P3

      neq = n_ddl_P2 + 2 * (n_ddl_P3 -1)

      if (debug .eq. 1) then
        write(*,*) "py_calc_modes_1d: n_ddl = ", n_ddl
        write(*,*) "py_calc_modes_1d: neq = ", neq
      endif

      call geometry_1d (nel, npt_P2, nnodes_P2, nb_typ_el,
     *     lx, type_nod, type_el, table_nod, 
     *     x_P2, mesh_file)

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

      allocate(sol_P2(3,nnodes_P2,nval,nel), STAT=allocate_status)
      if (allocate_status /= 0) then
        write(*,*) "The allocation is unsuccessful"
        write(*,*) "allocate_status = ", allocate_status
        write(*,*) "Not enough memory for sol_P2"
        write(*,*) "nval, nel = ", nval,nel
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
      do i=1,nb_typ_el
        n_eff(i) = 1
      enddo
ccc        n_eff(3) = 3

c     Calculate effective permittivity
      do i=1,nb_typ_el
        eps_eff(i) = n_eff(i)**2
      end do

cc      n_eff_0 = DBLE(n_eff(1))
cc      freq = 1.0d0/lambda
cc      k_0 = 2.0d0 * pi * freq
c      k_0 = 2.0d0*pi*n_eff_0*freq

      bloch_vec_x = 0.1d0 * 1
      bloch_vec_y = 0.2d0 * 1

      E_H_field = 1

      if(debug .eq. 1) then
        write(ui,*) "lambda = ", lambda
        write(ui,*) "k_0 = ", k_0
        write(ui,*) "n_eff_0 = ", n_eff_0
        write(ui,*) "bloch_vec_x = ", bloch_vec_x
        write(ui,*) "bloch_vec_y = ", bloch_vec_y
      endif

c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call periodic_cond_1d (nel, npt_P2, n_ddl, neq, table_nod, 
     *      ineq, x_P2, table_ddl, ip_period_ddl, x_ddl, debug,
     *      period_x)
C
      if(debug .eq. 1) then
        write(ui,*) "period_x = ", period_x
      endif
C
      write(ui,*) 
      write(ui,*) "--------------------------------------------",
     *     "-------"
      write(ui,*) " Wavelength : ", lambda, " (d)"
      write(ui,*) "--------------------------------------------",
     *     "-------"
      write(ui,*) 
C
        freq = 1.0d0/lambda
        n_eff_0 = DBLE(n_eff(1))
        k_0 = 2.0d0 * pi * freq

C  Index number of the core materials (material with highest Re(eps_eff))
        n_core(1) = 3
        do i=3,nb_typ_el
          if(dble(eps_eff(i)) .gt. dble(eps_eff(n_core(1)))) then
            n_core(1) = i
          endif
        enddo
        n_core(2) = n_core(1)
        shift = 1.1d0*Dble(n_eff(n_core(1)))**2 * k_0**2
     *      - bloch_vec_x**2 - bloch_vec_y**2
        If(debug .eq. 1) then
          write(ui,*) "n_core = ", n_core
          write(ui,*) "shift = ", shift
          if(E_H_field .eq. 1) then
            write(ui,*) "E-Field formulation"
          else
            write(ui,*) "H-Field formulation"
          endif
        EndIf

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
     *     beta, mode_pol, vp, sol)

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
     *  nb_typ_el, pp, qq, bloch_vec_y, table_nod, 
     *  type_el, x_P2, beta_1, beta_2,
     *  sol_1, sol_2, overlap_L, overlap_file, debug,
     *  pair_warning, k_0)

      if (pair_warning .ne. 0 .and. nval .le. 20) then
        write(ui,*) "py_calc_modes_1d.f: Warning found 1 BM
     * of cmplx conj"
        write(ui,*) "pair, increase num_BMs to include the other."
      endif

C    Save Original solution
      if (plot_modes .eq. 1) then
        call array_sol_P2_1d (nval, nel, sol_1, sol_P2)
        dir_name = "Bloch_Fields"
C        call write_sol_P2_1d (nval, nel, E_H_field,
C     *     lambda, beta_1, sol_P2, mesh_file, dir_name)
        q_average = 0
        tchar = "Bloch_Fields/PNG/All_plots_png_abs2_eE.geo"
        open (unit=34,file=tchar)
        do i=1,nval
          call gmsh_post_process_1d (i, E_H_field, nval, 
     *       nel, npt_P2, table_nod, type_el, nb_typ_el, 
     *       n_eff, x_P2, beta_1, sol_P2,
     *       gmsh_file_pos, dir_name, 
     *       q_average, plot_real, plot_imag, plot_abs)
        enddo
        close (unit=34)
      endif
C        
C  Normalisation
      if(debug .eq. 1) then
        write(ui,*) "py_calc_modes_1d.f: Field  Normalisation"
      endif
      call normalisation_1d (nval, nel, 
     *  sol_1, sol_2, overlap_L)

C
C  Orthonormal integral
      if (debug .eq. 1) then
        write(ui,*) "py_calc_modes_1d.f: Product of normalised field"
        overlap_file = "Normed/Orthogonal_n.txt"
        call cpu_time(time1_J)
        call orthogonal_1d (nval, nel, npt_P2, 
     *    nb_typ_el, pp, qq, bloch_vec_y, table_nod, 
     *    type_el, x_P2, beta_1, beta_2,
     *    sol_1, sol_2, overlap_L, overlap_file, debug,
     *    pair_warning, k_0)
        call cpu_time(time2_J)
          write(ui,*) "py_calc_modes_1d.f: CPU time for orthogonal :",
     *    (time2_J-time1_J)
      endif

C  Plane wave ordering
      if (debug .eq. 1) then
        write(ui,*) "py_calc_modes_1d.f: pw_ordering_1d"
      endif
      call pw_ordering_1d (neq_PW, period_x,
     *  bloch_vec_x, bloch_vec_y, index_pw_inv,
     *  debug, ordre_ls, k_0)

c      if (debug .eq. 1) then
c        call array_sol_test_1d (nval, nel, n_ddl, 
c     *           table_ddl, x_ddl, sol_1, period_x,
c     *           bloch_vec_x, bloch_vec_y)
c      endif

C  J_overlap
      if (debug .eq. 1) then
        write(ui,*) "py_calc_modes_1d.f: J_overlap_1d Integral"
      endif
      call J_overlap_1d(nval, nel, npt_P2, 
     *  type_el, table_nod, x_P2, sol_1,
     *  period_x, lambda, freq, overlap_J, neq_PW,
     *  bloch_vec_x, bloch_vec_y, index_pw_inv, debug,
     *  ordre_ls)


C  J_dagger_overlap
      if (debug .eq. 1) then
        write(ui,*) "py_calc_modes_1d.f: J_dagger_overlap_1d Integral"
      endif
      call J_dagger_overlap_1d(nval, nel, npt_P2, 
     *  type_el, table_nod, x_P2, sol_2, 
     *  period_x, lambda, freq, overlap_J_dagger, neq_PW,
     *  bloch_vec_x, bloch_vec_y, index_pw_inv, debug, ordre_ls)

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
