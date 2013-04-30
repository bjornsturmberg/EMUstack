      subroutine calc_scat(
c     Explicit inputs
     *    parallel, lambda, ordre_ls, d_in_nm,
     *    debug, 
     *    n_eff,
     *    bloch_vec, h_1, h_2, num_h, lx, ly, 
     *    pol, traLambda, PropModes,
     *    PrintSolution, PrintAll,
     *    Checks,
     *    incident, what4incident, out4incident, title,
     *    neq_PW, Zeroth_Order,
C     New explicit inputs
     *    beta1, sol1, sol2, type_el, table_nod, x_arr, pp, qq,
c     "Optional" inputs (Python guesses these)
     *    nval, npt, nel, nb_typ_el,
c     Outputs
     *    
C     *    J_overlap, J_dagger_overlap,
     *    T12, R12, T21, R21)
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
C  Plane wave parameters
      integer*8 neq_PW, nx_PW, ny_PW, ordre_ls
      integer*8 index_pw_inv(neq_PW)
      integer*8 Zeroth_Order, Zeroth_Order_inv, nb_typ_el
      complex*16 pp(nb_typ_el),  qq(nb_typ_el)
      complex*16 eps_eff(nb_typ_el), n_eff(nb_typ_el)
      double precision n_eff_0
c     i_cond = 0 => Dirichlet boundary condition
c     i_cond = 1 => Neumann boundary condition
c     i_cond = 2 => Periodic boundary condition
      integer*8 nel, npt, nnodes, ui, i_cond
C     ! Number of nodes per element
      parameter(nnodes=6)
      integer*8 type_el(nel), table_nod(nnodes, nel)
C, len_skyl, nsym
      integer*8 neq, debug
      integer*8 npt_p3, numberprop_S, numberprop_N, numberprop_S_b
C  Variable used by valpr
      integer*8 nval, nvect, itermax, ltrav
      integer*8 n_conv, i_base
      double precision ls_data(10)
c      integer*8 pointer_int(20), pointer_cmplx(20)
      integer*8 index(1000), n_core(2)
      integer*8 n_edge, n_face, n_ddl, n_ddl_max, n_k
      integer*8 i, j, mesh_format
c     Wavelength lambda is in normalised units of d_in_nm
      double precision lambda
      double precision freq, lat_vecs(2,2)
      double precision k_0, pi, lx, ly, bloch_vec(2), bloch_vec_k(2)
C  Timing variables
      double precision time1, time2
      double precision time1_J, time2_J
      character*(8) start_date, end_date
      character*(10) start_time, end_time
C  Names and Controls
      character mesh_file*100, gmsh_file*100, log_file*100
      character gmsh_file_pos*100
      character overlap_file*100, dir_name*100, buf1*4, buf2*4
      integer*8 namelength, PrintAll, Checks, traLambda
      integer*8 PrintSolution, pol
      integer*8 num_h
      integer*8 PropModes
C     Thicknesses h_1 and h_2 are in normalised units of d_on_lambda
      double precision h_1, h_2, hz
      integer*8 d_in_nm, parallel, Loss
      integer*8 incident, what4incident, out4incident
      integer*8 title

      integer*8 ip
      integer i_32

c     new breed of variables to prise out of a, b and c
      complex*16 x_arr(2,npt)
      complex*16 sol1(3,nnodes+7,nval,nel)
      complex*16 sol2(3,nnodes+7,nval,nel)
      complex*16 overlap_J(2*neq_PW, nval)
      complex*16 overlap_J_dagger(nval, 2*neq_PW)
      complex*16 overlap_K(nval, 2*neq_PW)
      complex*16 X_mat(2*neq_PW, 2*neq_PW)

      complex*16 beta1(nval), beta2(nval)
c     Fresnel scattering matrices
      complex*16 T12(nval,2*neq_PW)
      complex*16 R12(2*neq_PW,2*neq_PW)
      complex*16 T21(2*neq_PW,nval)
      complex*16 R21(nval,nval)

      complex*16 T_Lambda(2*neq_PW, 2*neq_PW)
      complex*16 R_Lambda(2*neq_PW, 2*neq_PW)
 
Cf2py intent(out) T12, R12, T21, R21

      call cpu_time(time1)

C     !ui = Unite dImpression
      ui = 6
C     ! Number of nodes per element
      pi = 3.141592653589793d0


      write(buf1,'(I4.4)') title
      write(buf2,'(I4.4)') parallel

      n_eff_0 = DBLE(n_eff(1))
      freq = 1.0d0/lambda
      k_0 = 2.0d0*pi*n_eff_0*freq

      call lattice_vec (npt, x_arr, lat_vecs)

C  Plane wave ordering
      call pw_ordering (neq_PW, lat_vecs, bloch_vec, 
     *  index_pw_inv, Zeroth_Order, Zeroth_Order_inv, 
     *  debug, ordre_ls, k_0)
C  J_overlap
      if (debug .eq. 1) then
        write(ui,*) "MAIN: J_overlap Integral"
      endif
      call cpu_time(time1_J)
      call J_overlap (nval, nel, npt, nnodes, 
     *  nb_typ_el, type_el, table_nod, x_arr, 
     *  sol1, pp, qq, lat_vecs, lambda, freq, n_eff_0,
     *  overlap_J, neq_PW, bloch_vec, X_mat, numberprop_S,
     *  index_pw_inv, PrintAll, debug, ordre_ls, k_0)
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
     *  nb_typ_el, type_el, table_nod, x_arr, 
     *  sol2, pp, qq, lat_vecs, lambda, freq,
     *  overlap_J_dagger, neq_PW, bloch_vec,
     *  index_pw_inv, PrintAll, ordre_ls)
      call cpu_time(time2_J)
      if (debug .eq. 1) then
        write(ui,*) "MAIN: CPU time for J_dagger_overlap :",
     *  (time2_J-time1_J)
      endif
C
C  Scattering Matrices
      if (debug .eq. 1) then
        write(ui,*) "MAIN: Scattering Matrices"
      endif
      call ScatMat( overlap_J, overlap_J_dagger,  
     *    X_mat, neq_PW, nval, 
     *    beta1, T12, R12, T21, R21,
     *    PrintAll, PrintSolution, 
     *    lx, h_1, h_2, num_h, Checks, T_Lambda, 
     *    R_Lambda, traLambda, pol, PropModes, lambda, d_in_nm,
     *    numberprop_S, freq, Zeroth_Order_inv,
     *    debug, incident, what4incident, out4incident,
     *    title, parallel)
C
C
C
CCCCCCCCCCCCCCCCCCCCC Calculation Checks CCCCCCCCCCCCCCCCCCCCC
C
C  Completeness Check
      if (Checks .eq. 1) then
        write(ui,*) "MAIN: K_overlap Integral"
        call K_overlap(nval, nel, npt, nnodes, 
     *    nb_typ_el, type_el, table_nod, x_arr,   
     *    sol2, pp, qq, lambda, freq, overlap_K, neq_PW,
     *    lat_vecs, bloch_vec, beta2, index_pw_inv,
     *    PrintAll, k_0, ordre_ls)
        write(ui,*) "MAIN: Completeness Test"
        call Completeness (nval, neq_PW, 
     *    overlap_K, overlap_J)
      write(ui,*) "numberprop_N = ", numberprop_N
C  Energy Conservation Check
        write(ui,*) "MAIN: Energy Check"
        call Energy_Cons(R12, T12, R21, T21,
     *    numberprop_S, numberprop_N, neq_PW, nval)
      endif 
C
C#########################  End Calculations  ###########################
C
      call date_and_time ( end_date, end_time )
      call cpu_time(time2)
C
      if (debug .eq. 1) then
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
        write(26,*)
        write(26,*) "lambda  = ", lambda
c
c         write(26,*) "len_skyl, npt*nval, len_skyl/(npt*nval) = ",
c     *   len_skyl, npt*nval, dble(len_skyl)/dble(npt*nval)
c
        write(26,*) 
        write(26,*) "log_file  = ", log_file
        close(26)
C
        write(ui,*) "  and   we're  done!"
      endif
C
      end subroutine calc_scat
