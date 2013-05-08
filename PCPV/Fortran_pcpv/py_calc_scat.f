      subroutine calc_scat(
c     Explicit inputs
     *    lambda, ordre_ls, 
     *    debug, 
     *    n_eff,
     *    bloch_vec, lx, ly, 
     *    PrintAll,
     *    Checks,
     *    neq_PW, Zeroth_Order,
C     New explicit inputs
     *    sol1, sol2, type_el, table_nod, x_arr, pp, qq,
c     "Optional" inputs (Python guesses these)
     *    nval, npt, nel, nb_typ_el,
c     Outputs
     *    overlap_J, overlap_J_dagger,
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
      complex*16 n_eff(nb_typ_el)
      double precision n_eff_0
      integer*8 nel, npt, nnodes, ui
C     ! Number of nodes per element
      parameter(nnodes=6)
      integer*8 type_el(nel), table_nod(nnodes, nel)
C, len_skyl, nsym
      integer*8 debug
      integer*8 numberprop_S, numberprop_N
      integer*8 nval
c     Wavelength lambda is in normalised units of d_in_nm
      double precision lambda
      double precision freq, lat_vecs(2,2)
      double precision k_0, pi, lx, ly, bloch_vec(2)
C  Timing variables
      double precision time1, time2
      double precision time1_J, time2_J
      character*(8) start_date, end_date
      character*(10) start_time, end_time
C  Names and Controls
      character log_file*100
      integer*8 PrintAll, Checks

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

Cf2py intent(out) overlap_J, overlap_J_dagger
Cf2py intent(out) T12, R12, T21, R21

      call cpu_time(time1)

C     !ui = Unite dImpression
      ui = 6
C     ! Number of nodes per element
      pi = 3.141592653589793d0

      if (n_eff(1) .ne. 1) then
        write(ui,*) "The code seems to assume that n_eff(1) = 1,"
        write(ui,*) "but today, n_eff(1) = ", n_eff(1)
        write(ui,*) "Aborting..."
        stop
      endif
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
     *    T12, R12, T21, R21,
     *    Checks, 
     *    debug)
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
