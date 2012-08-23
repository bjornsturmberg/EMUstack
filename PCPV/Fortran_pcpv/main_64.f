      program main

C************************************************************************
C
C  Program:
C    Finite Element Method - Scattering Matrix Method (FEM-SMM)
C     for Photonic Crystal Photovoltaics
C
C  Authors:
C    Bjorn Sturmberg & Kokou B. Dossou
C
C************************************************************************
C
      implicit none
C  Local parameters:
      integer*8 int_max, cmplx_max, int_used, cmplx_used
      integer*8 real_max, real_used, n_64
C      parameter (int_max=2**22, cmplx_max=2**26)
C      parameter (real_max=2**21)
      integer*8, dimension(:), allocatable :: a  !   a(int_max)
      complex*16, dimension(:), allocatable :: b  !  b(cmplx_max)
      double precision, dimension(:), allocatable :: c !  c(real_max)
      integer :: allocate_status=0
C
      integer*8 python    !take inputs from python (1) txt files (0)
      parameter (python=1)
C
C  Declare the pointers of the integer super-vector
      integer*8 ip_type_nod, ip_type_el, ip_table_nod
      integer*8 ip_table_E, ip_table_N_E_F, ip_visite
      integer*8 ip_type_N_E_F, ip_eq
      integer*8 ip_period_N, ip_nperiod_N
      integer*8 ip_period_N_E_F, ip_nperiod_N_E_F
      integer*8 ip_index_pw_inv
C      integer*8 ip_col_ptr, ip_bandw 
C  Declare the pointers of the real super-vector
      integer*8 jp_x, jp_x_N_E_F, jp_rhs
C      integer*8 jp_matD, jp_matL, jp_matU
C      integer*8 jp_matD2, jp_matL2, jp_matU2
      integer*8 jp_vect1, jp_vect2, jp_workd, jp_resid, jp_vschur
      integer*8 jp_eigenval_tmp, jp_trav, jp_vp, jp_eigen_pol
      integer*8 jp_overlap_L, jp_overlap_J, jp_overlap_J_dagger
      integer*8 jp_overlap_K, jp_X_mat
      integer*8 jp_sol, jp_sol1, jp_sol2
      integer*8 jp_eigenval, jp_eigenval1, jp_eigenval2
      integer*8 jp_T, jp_R, jp_T12, jp_R12, jp_T21, jp_R21
      integer*8 jp_T_Lambda, jp_R_Lambda
      integer*8 jp_X_mat_b
C  Plane wave parameteres
      integer*8 neq_PW, nx_PW, ny_PW, ordre_ls
      integer*8 Zeroth_Order, Zeroth_Order_inv, max_typ_el, nb_typ_el
      parameter (max_typ_el=10)
      complex*16 pp(max_typ_el),  qq(max_typ_el)
      complex*16 eps_eff(max_typ_el), n_eff(max_typ_el), test
      double precision n_eff_0, n_eff_sub, eps_eff_sub
      integer*8 nel, npt, nnodes, ui, i_cond
C, len_skyl, nsym
      integer*8 neq, debug, E_H_field
      integer*8 npt_p3, numberprop_S, numberprop_N, numberprop_S_b
C  Variable used by valpr
      integer*8 nval, nvect, itermax, ltrav
      integer*8 n_conv, i_base
      double precision ls_data(10)
c      integer*8 pointer_int(20), pointer_cmplx(20)
      integer*8 index(1000), index2(1000), n_core(2)
      complex*16 z_beta, z_tmp, z_tmp0
      integer*8 n_edge, n_face, n_ddl, n_ddl_max, n_k
c     variable used by UMFPACK
      double precision control (20), info_umf (90)
      integer*8 numeric, symbolic, status, sys, filenum
C  Renumbering
c      integer*8 ip_row_ptr, ip_bandw_1, ip_adjncy
c      integer*8 len_adj, len_adj_max, len_0_adj_max
c, iout, nonz_1, nonz_2
      integer*8 i, j, mesh_format
      integer*8 i_lambda, n_lambda, Lambda_Res, debug_i_lambda
      double precision lambda, lambda_1, lambda_2, d_lambda
      double precision d_freq, freq, lat_vecs(2,2), tol, theta, phi
      double precision k_0, pi, lx, ly, bloch_vec(2), bloch_vec0(2)
      complex*16 energ, shift, shift_eps_max
C  Timing variables
      double precision time1, time2
      double precision time1_fact, time2_fact
      double precision time1_asmbl, time2_asmbl
      double precision time1_postp, time2_postp
      double precision time1_arpack, time2_arpack
      double precision time1_J, time2_J
      double precision time_Lambda_end, time_Lambda_start
      character*(8) start_date, end_date
      character*(10) start_time, end_time
C  Names and Controls
      character mesh_file*100, gmsh_file*100, log_file*100
      character gmsh_file_pos*100, get_args*5
      character overlap_file*100, dir_name*100, buf1*4, buf2*4
      character*100 tchar
      integer*8 namelength, PrintAll, PrintOmega, Checks, traLambda
      integer*8 PrintSolution, hole_wire, pol, PrintSupModes
      integer*8 substrate, num_h
      integer*8 Lambda_count, parallel_comp, PropModes
      double precision h_1, h_2, hz
      integer*8 d_in_nm, pair_warning, parallel(2), max_arg, Loss
      complex*16 Complex_refract1(818), Complex_refract2(1635)
      complex*16 Complex_refract3(818), Complex_refract4(1635)
      integer*8 incident, what4incident, out4incident
      integer*8 q_average, plot_real, plot_imag, plot_abs
C  SuperMode Plotting
      complex*16 vec_coef(1000)
      complex*16 vec_coef_down(1000)
      complex*16 vec_coef_up(1000)

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
C  python passing
      integer input
      character*15 get_py
      character*100 get_char
      double precision re_tmp, im_tmp
C
C
CCCCCCCCCCCCCCCCCCCC  Start Program - get parameters  CCCCCCCCCCCCCCCCCCC
C 

      n_64 = 2
      cmplx_max=n_64**26 !n_64**28 on Vayu
      real_max=n_64**21
      int_max=n_64**22

      !write(*,*) "cmplx_max = ", cmplx_max
      !write(*,*) "real_max = ", real_max
      !write(*,*) "int_max = ", int_max

      allocate(b(cmplx_max), STAT=allocate_status)
      if (allocate_status /= 0) then
        write(*,*) "The allocation is unseccesfull"
        write(*,*) "allocate_status = ", allocate_status
        write(*,*) "Not enough memory for the complex array b"
        write(*,*) "cmplx_max = ", cmplx_max
        write(*,*) "Aborting..."
        stop
      endif

      allocate(c(real_max), STAT=allocate_status)
      if (allocate_status /= 0) then
        write(*,*) "The allocation is unseccesfull"
        write(*,*) "allocate_status = ", allocate_status
        write(*,*) "Not enough memory for the real array c"
        write(*,*) "real_max = ", real_max
        write(*,*) "Aborting..."
        stop
      endif

      allocate(a(int_max), STAT=allocate_status)
      if (allocate_status /= 0) then
        write(*,*) "The allocation is unseccesfull"
        write(*,*) "allocate_status = ", allocate_status
        write(*,*) "Not enough memory for the integer array a"
        write(*,*) "int_max = ", int_max
        write(*,*) "Aborting..."
        stop
      endif


      
      ui = 6     !ui = Unite dImpression
      nnodes = 6 ! Number of nodes per element
      pi = 3.141592653589793d0
C      nsym = 1 ! nsym = 0 => symmetric or hermitian matrices
C
C############### using with python for each wavelength #################C
C
      if (python .eq. 1) then
        input = 1
C parallel - int
        CALL GETARG(input, get_args)
        read (get_args,'(I5)') parallel(1)
        parallel(2) = parallel(1)
C        write(ui,*) parallel(1), parallel(2)
        input = input + 1
C lambda - double
        CALL GETARG(input, get_py)
        read (get_py,'(F18.12)') lambda
C        write(ui,*) lambda 
        input = input + 1
C nval - int
        CALL GETARG(input, get_args)
        read (get_args,'(I5)') nval
C        write(ui,*) nval
        input = input + 1
C ordre_ls - int
        CALL GETARG(input,get_args)
        read (get_args,'(I5)') ordre_ls 
C        write(ui,*) ordre_ls
        input = input + 1
C d_in_nm - int
        CALL GETARG(input,get_args)
        read (get_args,'(I5)') d_in_nm
C        write(ui,*) d_in_nm
        input = input + 1
C debug - int
        CALL GETARG(input,get_args)
        read (get_args,'(I5)') debug
C        write(ui,*) debug
        input = input + 1
C mesh_file - char
        CALL GETARG(input, get_char)
        read (get_char,'(a100)') mesh_file 
C        write(ui,*) mesh_file
        input = input + 1
C mesh_format - int
        CALL GETARG(input,get_args)
        read (get_args,'(I5)') mesh_format
C        write(ui,*) mesh_format
        input = input + 1
C
        namelength = len_trim(mesh_file)
C        gmsh_file = 'Normed/'//mesh_file(1:namelength)//'.msh'
        gmsh_file = mesh_file(1:namelength)//'.msh'
        gmsh_file_pos = mesh_file(1:namelength)
C        log_file = 'Normed/'//mesh_file(1:namelength)//'.log'
        log_file = mesh_file(1:namelength)//'.log'
        if (debug .eq. 1) then
          write(*,*) "get_param: mesh_file = ", mesh_file
          write(*,*) "get_param: gmsh_file = ", gmsh_file
        endif    
        open (unit=24,file="../PCPV/Data/"//mesh_file,
     *     status='unknown')
            read(24,*) npt, nel
        close(24)
C nb_typ_el -int
        CALL GETARG(input,get_args)
        read (get_args,'(I5)') nb_typ_el
        if(nb_typ_el .gt. max_typ_el) then
           write(ui,*)
           write(ui,*) "   ???"
           write(ui,*) "MAIN: nb_typ_el > max_typ_el : ",
     *    nb_typ_el, max_typ_el
           write(ui,*) "MAIN: Aborting..."
           stop
        endif
C        write(ui,*) nb_typ_el
        input = input + 1
C n_eff - 2* doubles, per el type
        do i_32 = 1,nb_typ_el
          CALL GETARG(input, get_py)
          read (get_py,'(F18.12)') re_tmp
          CALL GETARG(input+1, get_py)
          read (get_py,'(F18.12)') im_tmp
          n_eff(i_32) = dcmplx(re_tmp,im_tmp)
C          write(ui,*) n_eff(i_32)
          input = input + 2
C eps_eff
          eps_eff(i_32) = n_eff(i_32)**2
C          write(ui,*) eps_eff(i_32)
        enddo
C substrate - int
        CALL GETARG(input,get_args)
        read (get_args,'(I5)') substrate
C        write(ui,*) substrate
        input = input + 1
C theta - double 
        CALL GETARG(input, get_py)
        read (get_py,'(F18.12)') theta
C        write(ui,*) theta
        input = input + 1
C phi - double 
        CALL GETARG(input, get_py)
        read (get_py,'(F18.12)') phi
C        write(ui,*) phi
        input = input + 1
C h_1 - double
        CALL GETARG(input, get_py)
        read (get_py,'(F18.12)') h_1
C        write(ui,*) h_1
        input = input + 1
C h_2 - double
        CALL GETARG(input, get_py)
        read (get_py,'(F18.12)') h_2
C        write(ui,*) h_2
        input = input + 1
C num_h - int
        CALL GETARG(input,get_args)
        read (get_args,'(I5)') num_h
C        write(ui,*) num_h
        input = input + 1
C lx - double 
        CALL GETARG(input, get_py)
        read (get_py,'(F18.12)') lx
C        write(ui,*) lx
        input = input + 1
C ly - double 
        CALL GETARG(input, get_py)
        read (get_py,'(F18.12)') ly
C        write(ui,*) ly
        input = input + 1
C tol - double 
        CALL GETARG(input, get_py)
        read (get_py,'(F18.12)') tol
C        write(ui,*) tol
        input = input + 1
C E_H_field - int
        CALL GETARG(input,get_args)
        read (get_args,'(I5)') E_H_field
C        write(ui,*) E_H_field
        input = input + 1
C i_cond - int
        CALL GETARG(input,get_args)
        read (get_args,'(I5)') i_cond
C        write(ui,*) i_cond
        input = input + 1
C itermax - int
        CALL GETARG(input,get_args)
        read (get_args,'(I5)') itermax
C        write(ui,*) itermax
        input = input + 1
C pol - int
        CALL GETARG(input,get_args)
        read (get_args,'(I5)') pol
C        write(ui,*) pol
        input = input + 1
C traLambda - int
        CALL GETARG(input,get_args)
        read (get_args,'(I5)') traLambda 
C        write(ui,*) traLambda 
        input = input + 1
C PropModes - int
        CALL GETARG(input,get_args)
        read (get_args,'(I5)') PropModes
C        write(ui,*) PropModes
        input = input + 1
C PrintSolution - int
        CALL GETARG(input,get_args)
        read (get_args,'(I5)') PrintSolution
C        write(ui,*) PrintSolution
        input = input + 1
C PrintSupModes - int
        CALL GETARG(input,get_args)
        read (get_args,'(I5)') PrintSupModes
C        write(ui,*) PrintSupModes
        input = input + 1
C PrintOmega - int
        CALL GETARG(input,get_args)
        read (get_args,'(I5)') PrintOmega
C        write(ui,*) PrintOmega
        input = input + 1
C PrintAll - int
        CALL GETARG(input,get_args)
        read (get_args,'(I5)') PrintAll
C        write(ui,*) PrintAll
        input = input + 1
C Checks - int
        CALL GETARG(input,get_args)
        read (get_args,'(I5)') Checks
C        write(ui,*) Checks
        input = input + 1
C q_average - int
        CALL GETARG(input,get_args)
        read (get_args,'(I5)') q_average
C        write(ui,*) q_average
        input = input + 1
C plot_real - int
        CALL GETARG(input,get_args)
        read (get_args,'(I5)') plot_real
C        write(ui,*) plot_real
        input = input + 1
C plot_imag - int
        CALL GETARG(input,get_args)
        read (get_args,'(I5)') plot_imag
C        write(ui,*) plot_imag
        input = input + 1
C plot_abs - int
        CALL GETARG(input,get_args)
        read (get_args,'(I5)') plot_abs
C        write(ui,*) plot_abs
        input = input + 1
C incident - int
        CALL GETARG(input,get_args)
        read (get_args,'(I5)') incident
C        write(ui,*) incident
        input = input + 1
C  what4incident- int
        CALL GETARG(input,get_args)
        read (get_args,'(I5)') what4incident
C        write(ui,*) what4incident
        input = input + 1
C out4incident - int
        CALL GETARG(input,get_args)
        read (get_args,'(I5)') out4incident
C        write(ui,*) out4incident
        input = input + 1
C Loss - int
        CALL GETARG(input,get_args)
        read (get_args,'(I5)') Loss
C        write(ui,*) Loss
        input = input + 1
C
C
        nvect = 2*nval + nval/2 +3
C
C        write(ui,*) "MAIN: Aborting..."
C        stop
C
C############ using as independent program w text files ################C
C
      elseif (python .eq. 0) then
C
      open(3,file='Parameters/controls.txt',status='old')
        read(3,*) debug
        read(3,*) parallel_comp
        read(3,*) Loss
        read(3,*) hole_wire
        read(3,*) substrate
        read(3,*) pol
        read(3,*) traLambda
        read(3,*) PropModes
        read(3,*) PrintSolution
        read(3,*) PrintSupModes
        read(3,*) PrintOmega
        read(3,*) PrintAll
        read(3,*) Checks
      close(3)
C
      open(4,file='Parameters/bloch_vec.txt',status='old')
        read(4,*) theta
        read(4,*) phi
        read(4,*) ordre_ls
        read(4,*) h_1
        read(4,*) h_2
        read(4,*) num_h
      close(4)
C
      open (unit=346, file='Parameters/output.txt', status='old')
        read(346,*) q_average
        read(346,*) plot_real
        read(346,*) plot_imag
        read(346,*) plot_abs
        read(346,*) 
        read(346,*) incident
        read(346,*) what4incident
        read(346,*) out4incident
      close(346)
C
      if (parallel_comp .eq. 1) then
        max_arg = 2
        do i_32=1,max_arg
          CALL GETARG(i_32, get_args)
          read (get_args,'(I5)') parallel(i_32)
        enddo
      else
        parallel(1) = 1       ! n_par_start
        parallel(2) = 0       
      endif
c
      if (debug .eq. 1) then
        write(ui,*) "MAIN: parallel (A) ", parallel(1), parallel(2)
      endif
C
      call get_param(E_H_field, n_lambda, lambda_1, lambda_2, 
     *      npt, nel, i_cond, nval, nvect, itermax, tol,
     *      lx, ly, mesh_file, mesh_format, gmsh_file,
     *      gmsh_file_pos, log_file, Lambda_Res, d_in_nm, debug)
C
      call cmp_ref_index(Lambda_Res, Complex_refract1, Complex_refract2,
     *              Complex_refract3, Complex_refract4, Loss)
C
      else
        write(ui,*) "MAIN: either use python or text files"
        write(ui,*) "MAIN: Aborting..."
        stop
      endif
C
CCCCCCC end of parameter intake options
C
      call cpu_time(time1)  ! initial time  in unit = sec.
      call date_and_time ( start_date, start_time )
C
      if (debug .eq. 2) then
        open(4,file='Parameters/debug_lambda.txt',status='old')
          read(4,*) debug_i_lambda
        close(4)
      else
        debug_i_lambda = 0
      endif
C
      neq_PW = 0
      do nx_PW = -ordre_ls, ordre_ls
        do ny_PW = -ordre_ls, ordre_ls
          if (nx_PW**2 + ny_PW**2 .le. ordre_ls**2) then
            neq_PW = neq_PW + 1
            if (nx_PW .eq. 0 .and. ny_PW .eq. 0) then
              Zeroth_Order = neq_PW       !using neq_PW as counter
            endif
          endif
        enddo
      enddo      
C      
      if (debug .eq. 1) then
        write(ui,*)
        write(ui,*) "start_date = ", start_date
        write(ui,*) "start_time = ", start_time
        write(ui,*) "MAIN: ord_PW = ", ordre_ls
        write(ui,*) "MAIN: neq_PW = ", neq_PW
        write(ui,*)
      endif
C
      pair_warning = 0
C
C####################  Start FEM PRE-PROCESSING  ########################
C
      if (debug .eq. 1) then
        write(ui,*)
        write(ui,*) "lx,ly = ", lx, ly
        write(ui,*) "npt, nel, nnodes = ", npt, nel, nnodes
        write(ui,*) "mesh_file = ", mesh_file
        write(ui,*)
      endif
C
      if ((3*npt+nel+nnodes*nel) .gt. int_max) then
         write(ui,*) "MAIN: (3*npt+nel+nnodes*nel) + npt > int_max : ",
     *    (3*npt+nel+nnodes*nel), int_max
         write(ui,*) "MAIN: increase the size of int_max"
         write(ui,*) "MAIN: Aborting..."
         stop
      endif
      if ((7*npt) .gt. cmplx_max) then
         write(ui,*) "MAIN: (7*npt) > cmplx_max : ",
     *    (7*npt), cmplx_max
         write(ui,*) "MAIN: increase the size of cmplx_max"
         write(ui,*) "MAIN: Aborting..."
         stop
      endif
      ip_type_nod = 1
      ip_type_el = ip_type_nod + npt
      ip_table_nod = ip_type_el + nel ! pointer to FEM connectivity table
      jp_x = 1
C
      if (python .eq. 1) then
      call geometry_py (nel, npt, nnodes, nb_typ_el,
     *     lx, ly, a(ip_type_nod), a(ip_type_el), a(ip_table_nod), 
     *     b(jp_x), mesh_file)
      elseif (python .eq. 0) then
      call geometry (nel, npt, nnodes, nb_typ_el,
     *     lx, ly, a(ip_type_nod), a(ip_type_el), a(ip_table_nod), 
     *     b(jp_x), eps_eff, n_eff, mesh_file)
      else
         write(ui,*) "MAIN: python = 0 or = 1 !!!"
         write(ui,*) "MAIN: Aborting..."
         stop
      endif
C
      if (PrintSupModes + PrintSolution .ge. 1) then
C  Export the mesh to gmsh format
        call mail_to_gmsh (nel, npt, nnodes, a(ip_type_el), 
     *    a(ip_type_nod), a(ip_table_nod), 
     *    nb_typ_el, n_eff, b(jp_x), gmsh_file)
C
C        call gmsh_interface_cyl (nel, npt, nnodes, a(ip_type_el), 
C     *    a(ip_type_nod), a(ip_table_nod), 
C     *    nb_typ_el, b(jp_x))
      endif
C
      call lattice_vec (npt, b(jp_x), lat_vecs)
C
      if (PrintSupModes + PrintSolution .ge. 1) then
        call gmsh_interface_c4 (nel, npt, nnodes, a(ip_type_el), 
     *    a(ip_type_nod), a(ip_table_nod), 
     *    nb_typ_el, b(jp_x), lat_vecs)
      endif
C
C      V = number of vertices
C      E = number of edges
C      F =  number of faces
C      C =  number of cells (3D, tetrahedron)
C
C     From Euler's theorem on 3D graphs: V-E+F-C = 1 - (number of holes)
C     npt = (number of vertices) + (number of mid-edge point) = V + E;
C
      n_face = nel  ! each element is a face
      ip_table_N_E_F = ip_table_nod + nnodes*nel
      call list_face (nel, a(ip_table_N_E_F))

C     n_ddl_max = max(N_Vertices) + max(N_Edge) + max(N_Face)
C     For P2 FEM npt=N_Vertices+N_Edge
C     note: each element has 1 face, 3 edges and 10 P3 nodes
      n_ddl_max = npt + n_face
      ip_visite =  ip_table_N_E_F  + 14*nel 
      ip_table_E = ip_visite + n_ddl_max
C
      call list_edge (nel, npt, nnodes, n_edge, 
     *    a(ip_type_nod), a(ip_table_nod), 
     *    a(ip_table_E), a(ip_table_N_E_F), a(ip_visite))
      call list_node_P3 (nel, npt, nnodes, n_edge, npt_p3, 
     *    a(ip_table_nod), a(ip_table_N_E_F), a(ip_visite))
      n_ddl = n_edge + n_face + npt_p3
C
      if (debug .eq. 1) then
        write(ui,*) "MAIN: npt, nel = ", npt, nel
        write(ui,*) "MAIN: npt_p3 = ", npt_p3
        write(ui,*) "MAIN: n_vertex, n_edge, n_face, nel = ", 
     *    (npt - n_edge), n_edge, n_face, nel
        write(ui,*) "MAIN: 2D case of the Euler characteristic : ",
     *    "V-E+F=1-(number of holes)"
        write(ui,*) "MAIN: Euler characteristic: V - E + F = ", 
     *    (npt - n_edge) - n_edge + n_face !  - nel
      endif
cC
cC-----------
cC
cC     overwriting pointers ip_row_ptr, ..., ip_adjncy
c
      ip_type_N_E_F = ip_table_E + 4*n_edge
C
      jp_x_N_E_F = jp_x + 2*npt
      call type_node_edge_face (nel, npt, nnodes, n_ddl, 
     *      a(ip_type_nod), a(ip_table_nod), a(ip_table_N_E_F), 
     *      a(ip_visite), a(ip_type_N_E_F), 
     *      b(jp_x), b(jp_x_N_E_F))
C
      call get_coord_p3 (nel, npt, nnodes, n_ddl, 
     *      a(ip_table_nod), a(ip_type_nod), a(ip_table_N_E_F), 
     *      a(ip_type_N_E_F), b(jp_x), b(jp_x_N_E_F), a(ip_visite))
C
        ip_period_N = ip_type_N_E_F + 2*n_ddl
        ip_nperiod_N = ip_period_N + npt
        ip_period_N_E_F = ip_nperiod_N + npt
        ip_nperiod_N_E_F = ip_period_N_E_F + n_ddl
        ip_eq = ip_nperiod_N_E_F + n_ddl
C
      if (i_cond .eq. 0 .or. i_cond .eq. 1) then
        call bound_cond (i_cond, n_ddl, neq, a(ip_type_N_E_F), 
     *    a(ip_eq))
      elseif(i_cond .eq. 2) then
        if (debug .eq. 1) then
          write(ui,*) "###### periodic_node"
        endif
        call periodic_node(nel, npt, nnodes, a(ip_type_nod), 
     *      b(jp_x), a(ip_period_N), a(ip_nperiod_N),
     *      a(ip_table_nod), lat_vecs)
        if (debug .eq. 1) then
          write(ui,*) "MAIN: ###### periodic_N_E_F"
        endif
        call periodic_N_E_F (n_ddl, a(ip_type_N_E_F), 
     *      b(jp_x_N_E_F), a(ip_period_N_E_F), 
     *      a(ip_nperiod_N_E_F), lat_vecs)
        call periodic_cond (i_cond, n_ddl, neq, a(ip_type_N_E_F),
     *       a(ip_period_N_E_F), a(ip_eq), debug)
      else
        write(ui,*) "MAIN: i_cond has invalid value : ", i_cond
        write(ui,*) "MAIN: Aborting..."
        stop
      endif
C
      if (debug .eq. 1) then
        write(ui,*) "MAIN: neq, n_ddl = ", neq, n_ddl
      endif
C
C=====calcul du vecteur de localisation des colonnes
C     pour le traitement skyline de la matrice globale
C     Type of sparse storage of the global matrice:
C                                   Symmetric Sparse Skyline format
C     Determine the pointer for the Symmetric Sparse Skyline format
c
c      ip_col_ptr = ip_eq + 3*n_ddl
c      ip_bandw  = ip_col_ptr + neq + 1
c      ip_index_pw_inv = ip_bandw + neq + 1
c      int_used = ip_index_pw_inv + neq_PW 
cC
c      if (int_max .lt. int_used) then
c        write(ui,*)
c        write(ui,*) 'The size of the integer supervector is too small'
c        write(ui,*) 'integer super-vec: int_max  = ', int_max
c        write(ui,*) 'integer super-vec: int_used = ', int_used
c        write(ui,*) 'Aborting...'
c        stop
c      endif
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c     Sparse matrix storage
c
      ip_index_pw_inv = ip_eq + 3*n_ddl
      ip_col_ptr = ip_index_pw_inv + neq_PW 

      call csr_max_length (nel, n_ddl, neq, nnodes, 
     *  a(ip_table_N_E_F), a(ip_eq), a(ip_col_ptr), nonz_max)
c
c      ip = ip_col_ptr + neq + 1 + nonz_max
      ip = ip_col_ptr + neq + 1
      if (ip .gt. int_max) then
         write(ui,*) "main: ip > int_max : ",
     *    ip, int_max
         write(ui,*) "main: nonz_max = ", nonz_max
         write(ui,*) "main: increase the size of int_max"
         write(ui,*) "main: Aborting..."
         stop
      endif
c
      ip_row = ip_col_ptr + neq + 1

      call csr_length (nel, n_ddl, neq, nnodes, a(ip_table_N_E_F), 
     *  a(ip_eq), a(ip_row), a(ip_col_ptr), nonz_max, 
     *  nonz, max_row_len, ip, int_max, debug)

      ip_work = ip_row + nonz
      ip_work_sort = ip_work + 3*n_ddl
      ip_work_sort2 = ip_work_sort + max_row_len

c     sorting csr ...
      call sort_csr (neq, nonz, max_row_len, a(ip_row), 
     *  a(ip_col_ptr), a(ip_work_sort), a(ip_work), 
     *  a(ip_work_sort2))

      if (debug .eq. 1) then
        write(ui,*) "main: nonz_max = ", nonz_max
        write(ui,*) "main: nonz = ", nonz
        write(ui,*) "main: cmplx_max/nonz = ", 
     *    dble(cmplx_max)/dble(nonz)
      endif

      int_used = ip_work_sort2 + max_row_len

      if (int_max .lt. int_used) then
        write(ui,*)
        write(ui,*) 'The size of the integer supervector is too small'
        write(ui,*) 'integer super-vec: int_max  = ', int_max
        write(ui,*) 'integer super-vec: int_used = ', int_used
        write(ui,*) 'Aborting...'
        stop
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
cC
      jp_rhs = jp_x_N_E_F + 3*n_ddl
c     jp_rhs will also be used (in gmsh_post_process) to store a solution
      jp_mat2 = jp_rhs + max(neq, 3*npt)

cC     jp_rhs will also be used (in gmsh_post_process) to store a solution
cC     jp_matD = jp_rhs + max(neq, 3*npt)
c      jp_matD = jp_rhs + max(neq, 3*(nnodes+7)*nel)
c      jp_matL = jp_matD + neq
c      if (nsym .eq. 0) then
c        jp_matU = jp_matL
c      else
c        jp_matU = jp_matL + len_skyl
c      endif
cC
c      jp_matD2 = jp_matU + len_skyl
c      jp_matL2 = jp_matD2 + neq
c      if (nsym .eq. 0) then
c        jp_matU2 = jp_matL2
c      else
c        jp_matU2 = jp_matL2 + len_skyl
c      endif
C
c      jp_vect1 = jp_matU2 + len_skyl

      jp_vect1 = jp_mat2 + nonz
      jp_vect2 = jp_vect1 + neq
      jp_workd = jp_vect2 + neq
      jp_resid = jp_workd + 3*neq

      jp_sol1 = jp_resid + neq
      jp_sol2 = jp_sol1 + 3*(nnodes+7)*nval*nel
      jp_eigenval1 = jp_sol2 + 3*(nnodes+7)*nval*nel
      jp_eigenval2 = jp_eigenval1 + nval + 1
      jp_vschur = jp_eigenval2 + nval + 1     ! Eigenvectors
      jp_eigenval = jp_vschur + neq*nvect
      jp_eigen_pol = jp_eigenval + nval + 1
      jp_eigenval_tmp = jp_eigen_pol + nval*4
      jp_trav = jp_eigenval_tmp + nval + 1

      ltrav = 3*nvect*(nvect+2)
      jp_vp = jp_trav + ltrav
      jp_overlap_L = jp_vp + neq*nval
      jp_overlap_J = jp_overlap_L + nval*nval
      jp_overlap_J_dagger = jp_overlap_J + 2*neq_PW*nval
      jp_overlap_K = jp_overlap_J_dagger + 2*neq_PW*nval
      jp_X_mat = jp_overlap_K + 2*neq_PW*nval
      jp_T12 = jp_X_mat + 2*neq_PW*2*neq_PW  
      jp_T21 = jp_T12 + 2*neq_PW*nval
      jp_R12 = jp_T21 + 2*neq_PW*nval
      jp_R21 = jp_R12 + 4*neq_PW*neq_PW
      jp_T = jp_R21 + nval*nval
      jp_R = jp_T + 2*neq_PW*2*neq_PW
      jp_T_Lambda = jp_R + 2*neq_PW*2*neq_PW
      jp_R_Lambda = jp_T_Lambda + 2*neq_PW*2*neq_PW
 
      if (substrate .eq. 1) then
      jp_X_mat_b = jp_R_Lambda + 2*neq_PW*2*neq_PW
      cmplx_used = jp_X_mat_b + 2*neq_PW*2*neq_PW  
      else
        cmplx_used = jp_R_Lambda + 2*neq_PW*2*neq_PW
      endif
C
      if (cmplx_max .lt. cmplx_used)  then
         write(ui,*) 'The size of the real supervector is too small'
         write(ui,*) 'real super-vec: cmplx_max  = ', cmplx_max
         write(ui,*) 'real super-vec: cmplx_used = ', cmplx_used
         write(ui,*) 'Aborting...'
         stop
      endif
c
      kp_rhs_re = 1
      kp_rhs_im = kp_rhs_re + neq
      kp_lhs_re = kp_rhs_im + neq
      kp_lhs_im = kp_lhs_re + neq
      kp_mat1_re = kp_lhs_im + neq
      kp_mat1_im = kp_mat1_re + nonz
      real_used = kp_mat1_im + nonz

      if (real_max .lt. real_used) then
        write(ui,*)
        write(ui,*) 'The size of the real supervector is too small'
        write(ui,*) '2*nonz  = ', 2*nonz
        write(ui,*) 'real super-vec: real_max  = ', real_max
        write(ui,*) 'real super-vec: real_used = ', real_used
        write(ui,*) 'Aborting...'
        stop
      endif
c
c
c###############################################
c
c       ----------------------------------------------------------------
c       convert from 1-based to 0-based
c       ----------------------------------------------------------------
c
        do 60 j = 1, neq+1
            a(j+ip_col_ptr-1) = a(j+ip_col_ptr-1) - 1
60      continue
        do 70 j = 1, nonz
            a(j+ip_row-1) = a(j+ip_row-1) - 1
70      continue
c
c
c     The CSC indexing, i.e., ip_col_ptr, is 1-based 
c       (but valpr.f will change the CSC indexing to 0-based indexing)
      i_base = 0
c
C
C#####################  End FEM PRE-PROCESSING  #########################
C
      h_1 = h_1/d_in_nm   ! convert h from nm to lattice vector units d
      h_2 = h_2/d_in_nm 
      if (python .eq. 1) then ! lambda from nm - lattice vector units d
C  Avoid directly hitting Wood anomalies
        if (lambda .eq. d_in_nm) then
          lambda = lambda/d_in_nm + 1.0d-15
        elseif (2*lambda .eq. d_in_nm) then
          lambda = lambda/d_in_nm + 1.0d-15
        elseif (3*lambda .eq. d_in_nm) then
          lambda = lambda/d_in_nm + 1.0d-15
        elseif (4*lambda .eq. d_in_nm) then
          lambda = lambda/d_in_nm + 1.0d-15
        elseif (5*lambda .eq. d_in_nm) then
          lambda = lambda/d_in_nm + 1.0d-15
        elseif (6*lambda .eq. d_in_nm) then
          lambda = lambda/d_in_nm + 1.0d-15
        elseif (7*lambda .eq. d_in_nm) then
          lambda = lambda/d_in_nm + 1.0d-15
        else
          lambda = lambda/d_in_nm
        endif
      elseif (python .eq. 0) then
      lambda_1 = lambda_1/d_in_nm
      lambda_2 = lambda_2/d_in_nm

      if(n_lambda .gt. 1) then
        d_lambda = (lambda_2 - lambda_1)/dble(n_lambda-1)
        d_freq   = (1.0d0/lambda_2 - 1.0d0/lambda_1) /
     *             dble(n_lambda-1)
      else
        d_lambda = 0.0d0
        d_freq   = 0.0d0
      endif
      endif
C
      Lambda_count = 1    !  index wavelength loops for A_and_W_Lambda
C
C
      if (python .eq. 1) then
      write(buf2,'(I4.4)') parallel(1)
      if (PrintOmega .eq. 1) then
        open (unit=241, file="p_"//buf2//
     *    "_omega.txt", status='unknown')
        open (unit=231, file="p_"//buf2//
     *    "_omega_pol.txt",  status='unknown')
        open (unit=221, file="p_"//buf2//
     *    "_omega_Fz.txt", status='unknown')
        open (unit=211, file="p_"//buf2//
     *    "_omega_Ft.txt", status='unknown')
      endif
      if (debug .eq. 1) then
        write(ui,*) "MAIN: parallel (B) ", parallel(1), parallel(2)
      endif
      if (parallel(2) .eq. 0) parallel(2)=n_lambda
      if (debug .eq. 1) then
        write(ui,*) "MAIN: parallel (C) ", parallel(1), parallel(2)
      endif
C
      if (traLambda .eq. 1) then
        if (pol .eq. 0) then
          open(643,file="p_"//buf2//
     *    "_T_Lambda.txt",status='unknown')
          open(644,file="p_"//buf2//
     *    "_R_Lambda.txt",status='unknown')
          open(645,file="p_"//buf2//
     *    "_A_Lambda.txt",status='unknown')
          open(660,file="p_"//buf2//
     *    "_T_Lambda_MAT_sp.txt",status='unknown')
          open(661,file="p_"//buf2//
     *    "_R_Lambda_MAT_sp.txt",status='unknown')
        elseif (pol .eq. 5) then
          open(643,file="p_"//buf2//
     *    "_T_Lambda_R.txt",status='unknown')
          open(644,file="p_"//buf2//
     *    "_R_Lambda_R.txt",status='unknown')
          open(645,file="p_"//buf2//
     *    "_A_Lambda_R.txt",status='unknown')
          open(646,file="p_"//buf2//
     *    "_T_Lambda_L.txt",status='unknown')
          open(647,file="p_"//buf2//
     *    "_R_Lambda_L.txt",status='unknown')
          open(648,file="p_"//buf2//
     *    "_A_Lambda_L.txt",status='unknown')
          open(649,file="p_"//buf2//
     *    "_T_Lambda_CD.txt",status='unknown')
          open(650,file="p_"//buf2//
     *    "_R_Lambda_CD.txt",status='unknown')
          open(651,file="p_"//buf2//
     *    "_A_Lambda_CD.txt",status='unknown')
          open(660,file="p_"//buf2//
     *    "_T_Lambda_MAT_lr.txt",status='unknown')
          open(661,file="p_"//buf2//
     *    "_R_Lambda_MAT_lr.txt",status='unknown')
        else
          open(643,file="p_"//buf2//
     *    "_T_Lambda.txt",status='unknown')
          open(644,file="p_"//buf2//
     *    "_R_Lambda.txt",status='unknown')
          open(645,file="p_"//buf2//
     *    "_A_Lambda.txt",status='unknown')
        endif
      endif

      if (PropModes .eq. 1) then
        if (Loss .eq. 0) then
          open(unit=200, file="SuperModes"//buf2//".txt",
     *           status="unknown")   
        endif
        open(unit=565, file="detAe"//buf2//".txt",
     *         status="unknown")
        open(unit=566, file="detAo"//buf2//".txt",
     *         status="unknown")
      endif
C
      elseif (python .eq. 0) then
      write(buf1,'(I4.4)') d_in_nm
      write(buf2,'(I4.4)') parallel(1)
      if (PrintOmega .eq. 1) then
        open (unit=241, file="Output/Dispersion/p_"//buf2//
     *    "_omega.txt", status='unknown')
        open (unit=231, file="Output/Dispersion/p_"//buf2//
     *    "_omega_pol.txt",  status='unknown')
        open (unit=221, file="Output/Dispersion/p_"//buf2//
     *    "_omega_Fz.txt", status='unknown')
        open (unit=211, file="Output/Dispersion/p_"//buf2//
     *    "_omega_Ft.txt", status='unknown')
      endif
      if (debug .eq. 1) then
        write(ui,*) "MAIN: parallel (B) ", parallel(1), parallel(2)
      endif
      if (parallel(2) .eq. 0) parallel(2)=n_lambda
      if (debug .eq. 1) then
        write(ui,*) "MAIN: parallel (C) ", parallel(1), parallel(2)
      endif
C
      if (parallel_comp .ne. 1) write(buf2,'(I4.4)') 0000
C
      if (traLambda .eq. 1) then
        if (pol .eq. 0) then
          open(643,file="Output/d_"//buf1//"p_"//buf2//
     *    "_T_Lambda.txt",status='unknown')
          open(644,file="Output/d_"//buf1//"p_"//buf2//
     *    "_R_Lambda.txt",status='unknown')
          open(645,file="Output/d_"//buf1//"p_"//buf2//
     *    "_A_Lambda.txt",status='unknown')
          open(660,file="p_"//buf2//
     *    "_T_Lambda_MAT_sp.txt",status='unknown')
          open(661,file="p_"//buf2//
     *    "_R_Lambda_MAT_sp.txt",status='unknown')
        elseif (pol .eq. 5) then
          open(643,file="Output/d_"//buf1//"p_"//buf2//
     *    "_T_Lambda_R.txt",status='unknown')
          open(644,file="Output/d_"//buf1//"p_"//buf2//
     *    "_R_Lambda_R.txt",status='unknown')
          open(645,file="Output/d_"//buf1//"p_"//buf2//
     *    "_A_Lambda_R.txt",status='unknown')
          open(646,file="Output/d_"//buf1//"p_"//buf2//
     *    "_T_Lambda_L.txt",status='unknown')
          open(647,file="Output/d_"//buf1//"p_"//buf2//
     *    "_R_Lambda_L.txt",status='unknown')
          open(648,file="Output/d_"//buf1//"p_"//buf2//
     *    "_A_Lambda_L.txt",status='unknown')
          open(649,file="Output/d_"//buf1//"p_"//buf2//
     *    "_T_Lambda_CD.txt",status='unknown')
          open(650,file="Output/d_"//buf1//"p_"//buf2//
     *    "_R_Lambda_CD.txt",status='unknown')
          open(651,file="Output/d_"//buf1//"p_"//buf2//
     *    "_A_Lambda_CD.txt",status='unknown')
        else
          open(643,file="Output/d_"//buf1//"p_"//buf2//
     *    "_T_Lambda.txt",status='unknown')
          open(644,file="Output/d_"//buf1//"p_"//buf2//
     *    "_R_Lambda.txt",status='unknown')
          open(645,file="Output/d_"//buf1//"p_"//buf2//
     *    "_A_Lambda.txt",status='unknown')
        endif
      endif
C
      if (PropModes .eq. 1) then
        if (Loss .eq. 0) then
          open(unit=200, file="Normed/SuperModes"//buf2//".txt",
     *           status="unknown")   
        endif
        open(unit=565, file="Matrices/detAe"//buf2//".txt",
     *         status="unknown")
        open(unit=566, file="Matrices/detAo"//buf2//".txt",
     *         status="unknown")
      endif
      endif !python 1 or 0 (different data structures)
C
C
C#####################  Loop over Wavelengths  ##########################
C
      do i_lambda=parallel(1),parallel(2)
      call cpu_time(time_Lambda_start)
C
      if (python .eq. 1) then
        write(ui,*) 
        write(ui,*) "--------------------------------------------",
     *     "-------"
        write(ui,*) "MAIN: Wavelength Slice", parallel(1)
      elseif (python .eq. 0) then
        write(ui,*) 
        write(ui,*) "--------------------------------------------",
     *     "-------"
        write(ui,*) "MAIN: Wavelength Slice", parallel(1)
        write(ui,*) "MAIN: Slice Progress  ", 
     *     i_lambda-parallel(1)+1,'/', parallel(2)-parallel(1)+1
C
C  For uniformly spaced wavelengths
        lambda = lambda_1 + (i_lambda-1)*d_lambda
C  Avoid directly hitting Wood anomalies
        if (lambda .eq. 1.0d0) lambda = lambda + 1.0d-15
        if (lambda .eq. 0.5d0) lambda = lambda + 1.0d-15
        if (lambda-1 .eq. 0.5d0) lambda = lambda + 1.0d-15
        if (lambda+1 .eq. 0.5d0) lambda = lambda + 1.0d-15
      endif
C
C  For uniformly spaced frequencies
C        freq = 1.0d0/lambda_1 + (i_lambda-1)*d_freq
C        lambda = 1.0d0/freq
C
        n_eff_0 = DBLE(n_eff(1))
        freq = 1.0d0/lambda
        k_0 = 2.0d0*pi*n_eff_0*freq
C
      call Angular_Bloch_Vec(theta, phi, bloch_vec0, k_0)
C
      if (python .eq. 0) then
C  Find refractive index at lambda from cmp_ref_index list
C     - only works for Lambda_Res calculations 
C     - not single lambdas where need to set n const.
      call cmp_ref_ind_lambda (hole_wire, Lambda_Res, eps_eff,
     *    n_eff, Complex_refract1, Complex_refract2,
     *    Complex_refract3, Complex_refract4, 
     *    i_lambda, substrate)
      endif
C
C  Index number of the core materials (material with highest Re(eps_eff))
      if(dble(eps_eff(3)) .gt. dble(eps_eff(4))) then
          n_core(1) = 3
      else
          n_core(1) = 4
      endif
      n_core(2) = n_core(1)
      shift = 1.01d0*Dble(n_eff(n_core(1)))**2 * k_0**2
     *    - bloch_vec0(1)**2 - bloch_vec0(2)**2
      If(debug .eq. 1) then
        write(ui,*) "MAIN: n_core = ", n_core
        if(E_H_field .eq. 1) then
          write(ui,*) "MAIN: E-Field formulation"
        else
          write(ui,*) "MAIN: H-Field formulation"
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
           write(ui,*) "MAIN: action indef. avec E_H_field = ", 
     *                  E_H_field
           write(ui,*) "Aborting..."
           stop
         endif
C
CCCCCCCCCCCCCCCCCCCC  Loop over Prime and Adjoint  CCCCCCCCCCCCCCCCCCCCCC
C
         do n_k = 1,2
C
           if (n_k .eq. 1) then
             dir_name = "Output"
             jp_sol = jp_sol1
             jp_eigenval = jp_eigenval1
             bloch_vec(1) = bloch_vec0(1)
             bloch_vec(2) = bloch_vec0(2)
           else
             dir_name = "Output-"
             jp_sol = jp_sol2
             jp_eigenval = jp_eigenval2
             bloch_vec(1) = -bloch_vec0(1)
             bloch_vec(2) = -bloch_vec0(2)
           endif  
           namelength = len_trim(dir_name)
C
C     Assemble the coefficient matrix A and the right-hand side F of the
C     finite element equations
      if (debug .eq. 1) then
        write(ui,*) "MAIN: Asmbly: call to asmbly"
      endif
      call cpu_time(time1_asmbl)
      call asmbly (i_cond, i_base, nel, npt, n_ddl, neq, nnodes,
     *  shift, bloch_vec, nb_typ_el, pp, qq, a(ip_table_nod), 
     *  a(ip_table_N_E_F), a(ip_type_el), a(ip_eq),
     *   a(ip_period_N), a(ip_period_N_E_F), b(jp_x), b(jp_x_N_E_F), 
     *   nonz, a(ip_row), a(ip_col_ptr), c(kp_mat1_re), 
     *   c(kp_mat1_im), b(jp_mat2), a(ip_work))
      call cpu_time(time2_asmbl)
C
C     factorization of the globale matrice
C     -----------------------------------
C
      if (debug .eq. 1) then
        write(ui,*) "MAIN:        Prime(1) / Adjoint(2)", n_k
c        write(ui,*) "MAIN: factorisation: call to znsy"
      endif
C
      if (debug .eq. 1) then
        write(ui,*) "MAIN: call to valpr"
      endif
      call valpr_64 (i_base, nvect, nval, neq, itermax, ltrav,
     *  tol, nonz, a(ip_row), a(ip_col_ptr), c(kp_mat1_re),
     *  c(kp_mat1_im), b(jp_mat2), b(jp_vect1), b(jp_vect2),
     *  b(jp_workd), b(jp_resid), b(jp_vschur), b(jp_eigenval),
     *  b(jp_trav), b(jp_vp), c(kp_rhs_re), c(kp_rhs_im),
     *  c(kp_lhs_re), c(kp_lhs_im), n_conv, ls_data,
     *  numeric, filenum, status, control, info_umf, debug)
c
      if (n_conv .ne. nval) then
         write(ui,*) "MAIN: convergence problem with valpr_64"
         write(ui,*) "MAIN: n_conv != nval : ",
     *    n_conv, nval
         write(ui,*) "n_core(1), n_eff(n_core(1)) = ",
     *                n_core(1), n_eff(n_core(1))
         write(ui,*) "MAIN: Aborting..."
c         stop
      endif
c
      time1_fact = ls_data(1)
      time2_fact = ls_data(2)
c
      time1_arpack = ls_data(3)
      time2_arpack = ls_data(4)
C
      do i=1,nval
        z_tmp0 = b(jp_eigenval+i-1)
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
C        z_beta = sqrt(z_tmp)/k_0 !  Effective index
        b(jp_eigenval+i-1) = z_beta
      enddo
c
      if (n_conv .ne. nval) then
        do i=1,nval
c          write(ui,"(i4,2(g22.14),2(g18.10))") i,
          write(ui,*) i, b(jp_eigenval+i-1)
        enddo
         write(ui,*) "MAIN: n_conv != nval : ",
     *    n_conv, nval
         write(ui,*) "n_core(1), n_eff(n_core(1)) = ",
     *                n_core(1), n_eff(n_core(1))
         write(ui,*) sqrt(1.01d0*Dble(n_eff(n_core(1)))**2 ),
     *                1.01d0*Dble(n_eff(n_core(1)))**2 
         do i=1,nb_typ_el
            write(ui,*) "i, n_eff(i) = ", i,  n_eff(i)     ! ,  eps_eff(i)
         enddo
      endif
C
      call cpu_time(time1_postp)
C
      call z_indexx (nval, b(jp_eigenval), index)
C
C       The eigenvectors will be stored in the array b(jp_sol)
C       The eigenvalues and eigenvectors will be renumbered  
C                 using the permutation vector index
        call array_sol (i_cond, nval, nel, npt, n_ddl, neq, nnodes, 
     *   n_core, bloch_vec, index, a(ip_table_nod), 
     *   a(ip_table_N_E_F), a(ip_type_el), a(ip_eq), a(ip_period_N), 
     *   a(ip_period_N_E_F), b(jp_x), b(jp_x_N_E_F), b(jp_eigenval), 
     *   b(jp_eigenval_tmp), b(jp_eigen_pol), b(jp_vp), b(jp_sol))
C
      if(debug .eq. 1) then
        write(ui,*) 'index = ', (index(i), i=1,nval)
      endif
      if(debug .eq. 1) then
        write(ui,*)
        write(ui,*) "lambda, 1/lambda = ", lambda, 1.0d0/lambda
        write(ui,*) (bloch_vec(i)/(2.0d0*pi),i=1,2)
        write(ui,*) "sqrt(shift) = ", sqrt(shift)
        do i=1,nval
          write(ui,"(i4,2(g22.14),2(g18.10))") i, 
     *       b(jp_eigenval+i-1)
        enddo
      endif
C
C  Dispersion Diagram
      if (PrintOmega .eq. 1 .and. n_k .eq. 1) then
        call mode_energy (nval, nel, npt, n_ddl, nnodes, 
     *     n_core, a(ip_table_nod), a(ip_type_el), nb_typ_el, eps_eff, 
     *     b(jp_x), b(jp_sol), b(jp_eigenval), b(jp_eigen_pol))
        call DispersionDiagram(lambda, bloch_vec, shift,
     *     nval, n_conv, b(jp_eigenval), b(jp_eigen_pol), d_in_nm)
      endif
C
      enddo
C
CCCCCCCCCCCCCCCCCCCCCCCC  End Prime, Adjoint Loop  CCCCCCCCCCCCCCCCCCCCCC
C

CCCC Hardcore Debuging - Print all arrays + some variables CCCCC
      if (debug .eq. 2 .and. i_lambda .eq. debug_i_lambda) then
        PrintAll = 1
        Checks = 2
        
        open (unit=1111,file="Normed/Debug_data.txt", status='unknown')
        write(1111,*) "i_lambda = ", i_lambda
        write(1111,*) "lambda = ", lambda
        write(1111,*) "eps_eff = ", (eps_eff(i),i=1,nb_typ_el)
        write(1111,*) "shift = ", shift
        write(1111,*) "bloch_vec(1) = ", bloch_vec(1)
        write(1111,*) "bloch_vec(2) = ", bloch_vec(2) 
        write(1111,*) "k_0 = ", k_0
        write(1111,*) 
        do i=1,nval
          write(1111,"(i4,2(g22.14),g18.10)") i, 
     *       b(jp_eigenval1+i-1)
        enddo
      endif
CCCC Hardcore Debuging - End                               CCCCC

      bloch_vec(1) = bloch_vec0(1)
      bloch_vec(2) = bloch_vec0(2)
C  Orthogonal integral
      if (debug .eq. 1) then 
        write(ui,*) "MAIN: Field product"
      endif
      overlap_file = "Normed/Orthogonal.txt"
      call cpu_time(time1_J)
      call orthogonal (nval, nel, npt, nnodes, 
     *  nb_typ_el, pp, qq, a(ip_table_nod), 
     *  a(ip_type_el), b(jp_x), b(jp_eigenval1), b(jp_eigenval2),
     *  b(jp_sol1), b(jp_sol2), b(jp_overlap_L),
     *  overlap_file, PrintAll, d_in_nm, pair_warning, k_0)
      call cpu_time(time2_J)
      if (debug .eq. 1) then
        write(ui,*) "MAIN: CPU time for orthogonal :",
     *  (time2_J-time1_J)
      endif     
C    Save Original solution
      if (PrintSolution .eq. 1) then
        dir_name = "Output/Fields"
        call write_sol (nval, nel, nnodes, E_H_field, lambda,
     *       b(jp_eigenval1), b(jp_sol1), mesh_file, dir_name)
        call write_param (E_H_field, lambda, npt, nel, i_cond,
     *       nval, nvect, itermax, tol, shift, lx, ly, 
     *       mesh_file, mesh_format, n_conv, nb_typ_el, eps_eff,
     *       bloch_vec, dir_name)
      tchar = "Output/FieldsPNG/All_plots_png_abs2_eE.geo"
      open (unit=34,file=tchar)
        do i=1,nval
          call gmsh_post_process (i, E_H_field, nval, nel, npt, 
     *       nnodes, a(ip_table_nod), a(ip_type_el), nb_typ_el,
     *       n_eff, b(jp_x), b(jp_eigenval1), b(jp_sol1),
     *       b(jp_rhs), a(ip_visite), gmsh_file_pos, dir_name, 
     *       q_average, plot_real, plot_imag, plot_abs)
        enddo 
      close (unit=34)
      endif
C        
C  Normalisation
      if(debug .eq. 1) then
        write(ui,*) "MAIN: Field  Normalisation"
      endif 
      call cpu_time(time1_J)
      call normalisation (nval, nel, nnodes, a(ip_table_nod),
     *  b(jp_sol1), b(jp_sol2), b(jp_overlap_L))  
      call cpu_time(time2_J)
      if (debug .eq. 1) then
        write(ui,*) "MAIN: CPU time for normalisation :",
     *  (time2_J-time1_J)
      endif  
C
C  Orthnormal integral
      if (PrintAll .eq. 1) then
        write(ui,*) "MAIN: Product of normalised field"
        overlap_file = "Normed/Orthogonal_n.txt"
        call cpu_time(time1_J)
        call orthogonal (nval, nel, npt, nnodes, 
     *    nb_typ_el, pp, qq, a(ip_table_nod), 
     *    a(ip_type_el), b(jp_x), b(jp_eigenval1), b(jp_eigenval2),
     *    b(jp_sol1), b(jp_sol2), b(jp_overlap_L),
     *    overlap_file, PrintAll, d_in_nm, pair_warning, k_0)
        call cpu_time(time2_J)
          write(ui,*) "MAIN: CPU time for orthogonal :",
     *    (time2_J-time1_J)
      endif
C
C  Plane wave ordering
      call pw_ordering (neq_PW, lat_vecs, bloch_vec, 
     *  a(ip_index_pw_inv), Zeroth_Order, Zeroth_Order_inv, 
     *  debug, ordre_ls, k_0)
C  J_overlap
      if (debug .eq. 1) then
        write(ui,*) "MAIN: J_overlap Integral"
      endif
      call cpu_time(time1_J)
      call J_overlap (nval, nel, npt, nnodes, 
     *  nb_typ_el, a(ip_type_el), a(ip_table_nod), b(jp_x), 
     *  b(jp_sol1), pp, qq, lat_vecs, lambda, freq, n_eff_0,
     *  b(jp_overlap_J), neq_PW, bloch_vec, b(jp_X_mat), numberprop_S,
     *  a(ip_index_pw_inv), PrintAll, debug, ordre_ls, k_0)
      call cpu_time(time2_J)
      if (debug .eq. 1) then
        write(ui,*) "MAIN: CPU time for J_overlap :",
     *  (time2_J-time1_J)
      endif
C
C  J_dagger_overlap
      if (debug .eq. 1) then
        write(ui,*) "MAIN: J_dagger_overlap Integral"
      endif
      call cpu_time(time1_J)
      call J_dagger_overlap (nval, nel, npt, nnodes, 
     *  nb_typ_el, a(ip_type_el), a(ip_table_nod), b(jp_x), 
     *  b(jp_sol2), pp, qq, lat_vecs, lambda, freq,
     *  b(jp_overlap_J_dagger), neq_PW, bloch_vec,
     *  a(ip_index_pw_inv), PrintAll, ordre_ls)
      call cpu_time(time2_J)
      if (debug .eq. 1) then
        write(ui,*) "MAIN: CPU time for J_dagger_overlap :",
     *  (time2_J-time1_J)
      endif
C
C  Overlaps at bottom Substrate
      if (substrate .eq. 1) then
        n_eff_sub = DBLE(n_eff(2))
        eps_eff_sub = DBLE(eps_eff(2))
      call J_overlap_sub ( lat_vecs, lambda, freq, n_eff_sub,
     *  eps_eff_sub, neq_PW, bloch_vec, b(jp_X_mat_b), numberprop_S_b,
     *  a(ip_index_pw_inv), PrintAll, ordre_ls, k_0)
C  Scattering Matrices
      if (debug .eq. 1) then
        write(ui,*) "MAIN: Scattering Matrices substrate"
      endif
C
      call ScatMat_sub ( b(jp_overlap_J), b(jp_overlap_J_dagger),  
     *    b(jp_X_mat), b(jp_X_mat_b), neq_PW, nval, 
     *    b(jp_eigenval1), b(jp_T12), b(jp_R12), b(jp_T21), b(jp_R21),
     *    PrintAll, PrintSolution, 
     *    lx, h_1, h_2, num_h, Checks, b(jp_T_Lambda), 
     *    b(jp_R_Lambda), traLambda, pol, PropModes, lambda, d_in_nm,
     *    numberprop_S, numberprop_S_b, freq, Zeroth_Order_inv,
     *    debug, incident, what4incident, out4incident)
      else  !No Substrate
C  Scattering Matrices
      if (debug .eq. 1) then
        write(ui,*) "MAIN: Scattering Matrices"
      endif
      call ScatMat( b(jp_overlap_J), b(jp_overlap_J_dagger),  
     *    b(jp_X_mat), neq_PW, nval, 
     *    b(jp_eigenval1), b(jp_T12), b(jp_R12), b(jp_T21), b(jp_R21),
     *    PrintAll, PrintSolution, 
     *    lx, h_1, h_2, num_h, Checks, b(jp_T_Lambda), 
     *    b(jp_R_Lambda), traLambda, pol, PropModes, lambda, d_in_nm,
     *    numberprop_S, freq, Zeroth_Order_inv,
     *    debug, incident, what4incident, out4incident)
      endif  !End Substrate Options
C
      if (debug .eq. 2 .and. i_lambda .eq. debug_i_lambda) then
        write(1111,*) 
        write(1111,*) "neq_PW = ", neq_PW
        write(1111,*) "numberprop_S = ", numberprop_S
        if (substrate .eq. 1) then
          write(1111,*) "numberprop_S_b = ", numberprop_S_b
        endif
        write(1111,*) "Zeroth_Order = ", Zeroth_Order
        write(1111,*) "Zeroth_Order_inv = ", Zeroth_Order_inv
      PrintAll = 0
      Checks = 0
      endif
C     
C  Search for number of propagating Bloch Modes
      if (PropModes .eq. 1 .and. Loss .eq. 0) then
      numberprop_N = 0
      do i=1,nval
        test = b(jp_eigenval1 + i - 1)
        if (ABS(IMAG(test)) .lt. 1.0d-5) then
          numberprop_N = numberprop_N + 1
        endif
      enddo
      write(200,*) lambda*d_in_nm, numberprop_N
      endif
C
C  Plot superposition of fields
      if (PrintSupModes .eq. 1) then
C     Coefficent of the modes travelling downward 
      do i=1,nval
        vec_coef(i) = b(jp_T12+i-1)
      enddo
C     Coefficent of the modes travelling upward 
      do i=1,nval
        vec_coef(i+nval) = 0.0d0
      enddo
      dir_name = "Output/Fields"
      hz = 0.0d0  ! hz=0 => top interface; hz=h => bottom interface
      i = 1 ! reference number of the field
      call gmsh_plot_field (i, E_H_field, nval, nel, npt, 
     *     nnodes, a(ip_table_nod), a(ip_type_el), eps_eff, b(jp_x),  
     *     b(jp_eigenval1), b(jp_sol1), b(jp_rhs), 
     *     vec_coef, h_1, hz, gmsh_file_pos, dir_name, nb_typ_el, 
     *       q_average, plot_real, plot_imag, plot_abs)
C
      hz = h_1  ! hz=0 => top interface; hz=h => bottom interface
      i = 2 ! reference number of the field
      call gmsh_plot_field (i, E_H_field, nval, nel, npt, 
     *     nnodes, a(ip_table_nod), a(ip_type_el), eps_eff, b(jp_x),
     *     b(jp_eigenval1), b(jp_sol1), b(jp_rhs), 
     *     vec_coef, h_1, hz, gmsh_file_pos, dir_name, nb_typ_el, 
     *       q_average, plot_real, plot_imag, plot_abs)

      do i=1,2*neq_PW
        vec_coef_down(i) = 0.0d0
      enddo
      vec_coef_down(1) = 1.0d0  ! Incident field
      do i=1,2*neq_PW
        vec_coef_up(i) = b(jp_R12+i-1) ! Reflected field
      enddo

      hz = 0.0d0 ! Upper semi-inifinite medium: Plane wave expansion
      i = 1
      call gmsh_plot_PW (i, E_H_field, 
     *     nel, npt, nnodes, neq_PW, bloch_vec, 
     *  a(ip_table_nod), b(jp_x), lat_vecs, lambda, eps_eff(1),
     *  b(jp_rhs), vec_coef_down, vec_coef_up, 
     *  a(ip_index_pw_inv), ordre_ls, h_1, hz, gmsh_file_pos,
     *  dir_name, q_average, plot_real, plot_imag, plot_abs)

      hz = h_1 ! Upper semi-inifinite medium: Plane wave expansion
      i = 2
      call gmsh_plot_PW (i, E_H_field, 
     *     nel, npt, nnodes, neq_PW, bloch_vec, 
     *  a(ip_table_nod), b(jp_x), lat_vecs, lambda, eps_eff(1),
     *  b(jp_rhs), vec_coef_down, vec_coef_up, 
     *  a(ip_index_pw_inv), ordre_ls, h_1, hz, gmsh_file_pos,
     *  dir_name, q_average, plot_real, plot_imag, plot_abs)
      endif
C
      Lambda_count = Lambda_count + 1
C
      call cpu_time(time_Lambda_end)
      if (python .eq. 0) then
      write(ui,*) "MAIN: Time for Lambda loop:    ", 
     *                time_Lambda_end-time_Lambda_start
      write(ui,*)"-----------------------------------Time-Remaining==",
     *  (parallel(2)-parallel(1)+2-Lambda_count)
     *  *(time_Lambda_end-time_Lambda_start)
      endif
C
      enddo
C
C#########################  End Wavelength Loop  ########################
C
      if (PrintOmega .eq. 1) then
        close(241)
        close(231)
        close(221)
        close(211)
      endif
      if (PropModes .eq. 1) then
        if (Loss .eq. 0) then
          close(200)
        endif
        close(565)
        close(566)
      endif
      if (traLambda .eq. 1) then
        if (pol .eq. 0) then
          close(643)
          close(644)
          close(645)
          close(660)
          close(661)
        elseif (pol .eq. 5) then
          close(643)
          close(644)
          close(645)
          close(646)
          close(647)
          close(648)
          close(649)
          close(650)
          close(651)
          close(660)
          close(661)
        else
          close(643)
          close(644)
          close(645)
        endif
      endif
C
CCCCCCCCCCCCCCCCCCCCC Ultimate Efficiency CCCCCCCCCCCCCCCCCCCCC
C                          as in Lin and Povinelli 09
      if (python .eq. 0) then
      if (parallel_comp .ne. 1 .and. Lambda_Res .ne. 0) then
        if (debug .eq. 1) then
          write(ui,*) "MAIN: Ultimate Efficiency calculation"
        endif
        call Ultimate_Efficiency(Lambda_Res, d_in_nm)
      endif
      endif
C
      if (pair_warning .ne. 0) then
        write(ui,*) "conjugate pair problem", pair_warning,
     *      "times, for d = ", d_in_nm
      endif
C
      call cpu_time(time2_postp)
C
CCCCCCCCCCCCCCCCCCCCC Calculation Checks CCCCCCCCCCCCCCCCCCCCC
C
C  Completeness Check
      if (Checks .eq. 1) then
        write(ui,*) "MAIN: K_overlap Integral"
        call K_overlap(nval, nel, npt, nnodes, 
     *    nb_typ_el, a(ip_type_el), a(ip_table_nod), b(jp_x),   
     *    b(jp_sol2), pp, qq, lambda, freq, b(jp_overlap_K), neq_PW,
     *    lat_vecs, bloch_vec, b(jp_eigenval2), a(ip_index_pw_inv),
     *    PrintAll, k_0, ordre_ls)
        write(ui,*) "MAIN: Completenes Test"
        call Completeness (nval, neq_PW, 
     *    b(jp_overlap_K), b(jp_overlap_J))
C  Search for number of propagating Bloch Modes
      numberprop_N = 0
      do i=1,nval
        test = b(jp_eigenval1 + i - 1)
        if (ABS(IMAG(test)) .lt. 1.0d-5) then
          numberprop_N = numberprop_N + 1
        endif
      enddo
      write(ui,*) "numberprop_N = ", numberprop_N
C  Energy Conservation Check
        write(ui,*) "MAIN: Energy Check"
        call Energy_Cons(b(jp_R12), b(jp_T12), b(jp_R21), b(jp_T21),
     *    numberprop_S, numberprop_N, neq_PW, nval)
      endif 
C
C#########################  End Calculations  ###########################
C
      call date_and_time ( end_date, end_time )
      call cpu_time(time2)
C
        if (python .eq. 0) then
        write(ui,*) 
        write(ui,*) 'Total CPU time (sec.)  = ', (time2-time1)

      open (unit=26,file=log_file)
        write(26,*)
        write(26,*) "Date and time formats = ccyymmdd ; hhmmss.sss"
        write(26,*) "Start date and time   = ", start_date, 
     *    " ; ", start_time
        write(26,*) "End date and time     = ", end_date, 
     *    " ; ", end_time
        write(26,*) "Total CPU time (sec.) = ",  (time2-time1)
        write(26,*) "LU factorisation : CPU time and % Total time = ",  
     *         (time2_fact-time1_fact), 
     *         100*(time2_fact-time1_fact)/(time2-time1),"%"
        write(26,*) "ARPACK : CPU time and % Total time = ",  
     *         (time2_arpack-time1_arpack), 
     *         100*(time2_arpack-time1_arpack)/(time2-time1),"%"
        write(26,*) "Assembly : CPU time and % Total time = ",  
     *         (time2_asmbl-time1_asmbl), 
     *         100*(time2_asmbl-time1_asmbl)/(time2-time1),"%"
        write(26,*) "Post-processsing : CPU time and % Total time = ",  
     *         (time2_postp-time1_postp), 
     *         100*(time2_postp-time1_postp)/(time2-time1),"%"
        write(26,*) "Pre-Assembly : CPU time and % Total time = ",  
     *         (time1_asmbl-time1), 
     *         100*(time1_asmbl-time1)/(time2-time1),"%"
        write(26,*)
        write(26,*) "lambda  = ", lambda
        write(26,*) "npt, nel, nnodes  = ", npt, nel, nnodes
        write(26,*) "neq, i_cond = ", neq, i_cond
        if ( E_H_field .eq. 1) then
          write(26,*) "E_H_field         = ", E_H_field,  
     *                 " (E-Field formulation)"
        elseif ( E_H_field .eq. 2) then
          write(26,*) "E_H_field         = ", E_H_field,  
     *                 " (H-Field formulation)"
       else
          write(ui,*) "MAIN (B): action indef. avec E_H_field = ", 
     *                 E_H_field
          write(ui,*) "Aborting..."
          stop
        endif
        write(26,*) "   bloch_vec = ", bloch_vec
        write(26,*) "bloch_vec/pi = ", (bloch_vec(i)/pi,i=1,2)
        z_tmp = sqrt(shift)/(2.0d0*pi)
        write(26,*) "shift             = ", shift, z_tmp
        write(26,*) "integer super-vector :"
        write(26,*) "int_used, int_max, int_used/int_max         = ", 
     *    int_used , int_max, dble(int_used)/dble(int_max)
        write(26,*) "cmplx super-vector : "
        write(26,*) "cmplx_used, cmplx_max, cmplx_used/cmplx_max = ",
     *     cmplx_used, cmplx_max, dble(cmplx_used)/dble(cmplx_max)
        write(26,*) "Real super-vector : "
        write(26,*) "real_used, real_max, real_max/real_used = ",
     *     real_used, real_max, dble(real_max)/dble(real_used)
         write(26,*)
         write(26,*) "neq_PW = ", neq_PW
         write(26,*) "nval, nvect, n_conv = ", nval, nvect, n_conv
         write(26,*) "nonz, npt*nval, nonz/(npt*nval) = ",
     *   nonz, npt*nval, dble(nonz)/dble(npt*nval)
         write(26,*) "nonz, nonz_max, nonz_max/nonz = ", 
     *   nonz, nonz_max, dble(nonz_max)/dble(nonz)
         write(26,*) "nonz, int_used, int_used/nonz = ", 
     *   nonz, int_used, dble(int_used)/dble(nonz)
c
c         write(26,*) "len_skyl, npt*nval, len_skyl/(npt*nval) = ",
c     *   len_skyl, npt*nval, dble(len_skyl)/dble(npt*nval)
c
        write(26,*) 
        do i=1,nval
          write(26,"(i4,2(g22.14),g18.10)") i, 
     *       b(jp_eigenval1+i-1)
        enddo
        write(26,*)
        write(26,*) "n_core = ", n_core
        write(26,*) "eps_eff = ", (eps_eff(i),i=1,nb_typ_el)
        write(26,*) "n_eff = ", (n_eff(i),i=1,nb_typ_el)
        write(26,*)         
        write(26,*) "conjugate pair problem", pair_warning,
     *      "times, for d = ", d_in_nm
        write(26,*)
        write(26,*) "mesh_file = ", mesh_file
        write(26,*) "gmsh_file = ", gmsh_file
        write(26,*) "log_file  = ", log_file
      close(26)
C
      write(ui,*) "   .      .      ."
      write(ui,*) "   .      .      ."
      write(ui,*) "   .      . (d=",d_in_nm,")"
      write(ui,*) "  and   we're   done!"
      endif
C
      stop
      end 
