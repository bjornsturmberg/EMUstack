c     ------------------------------------------------------------------
c
c                    sous-routine VALPR_64.
c                                 ------
c     gere l'utilisation de ARPACK
c
c     ------------------------------------------------------------------
c
c     Parametres d'entree:
c     --------------------
c
c     nvect (I)              : dimension de l'espace de Krylov
c     nval (I)               : nombre de valeurs propres desirees.
c     neq (I)                : nombre d'equations
c     workd (DP)             : matrice de travail pour dnaupd, 
c                              taille  3 neq 
c     resid   (DP)           : vecteur de trvail pour dnaupd,
c                              taille neq
c     v (DP)                 : matrice des vecteurs de Schur,
c                              taille neq*nvect
c     asup,ainf,adiag (DP)   : stockage de la matrice tangente.
c     msup,minf,mdiag (DP)   : stockage de la matrice masse.
c     kld (I)                : vecteur de localisation des debuts de
c                              colonnes
c     vect1,vect2,vect3 (DP) : vecteurs de travail,
c                              taille neq
c     long (I)               : longueur des super-vecteurs pour
c                              les matrices.
c     ddot (DP)              : fonction appelee pour calculer
c                              le produit scalaire
c     trav(DP)               : espace de travail pour dnaupd, 
c                              taille ltrav>= 3*nvect^2+6*nvect
c     action                 : 0: redemarrage avec un vecteur aleatoire
c                              1: lire une vecteur initial                                 c
c     Parametres de sortie:
c     ---------------------
c
c     reel (DP)              : parties reelles, taille nvect
c     imag (DP)              : parties imaginaires, taille nvect
c     
c     ------------------------------------------------------------------
c
      subroutine valpr_64 (i_base, nvect, nval, neq, itermax, ltrav, 
     *  tol, nonz, row_ind, col_ptr, mat1_re, mat1_im, mat2,
     *  vect1, vect2, workd, resid, v, d, trav, vp,   
     *  rhs_re, rhs_im, lhs_re, lhs_im, n_conv, ls_data,
     *  numeric, filenum, status, control, info_umf, debug)
c
c     ------------------------------------------------------------------
c
      implicit none
c
      integer*8 neq, nonz, n_conv, i_base
      integer*8 row_ind(nonz), col_ptr(neq+1)
      complex*16 mat2(nonz)
      double precision mat1_re(nonz), mat1_im(nonz)
      double precision rhs_re(neq), rhs_im(neq)
      double precision lhs_re(neq), lhs_im(neq)
      double precision ls_data(10)
c
      double precision time1_fact, time2_fact
      double precision time1_arpack, time2_arpack
c
      double precision control (20), info_umf (90)
      integer*8 numeric, symbolic, status, sys, filenum
c
      integer*8 itermax, nvect, nval, i, j, ltrav
      integer*8 compteur
      complex*16 resid(neq), v(neq,nvect), workd(3*neq)
      complex*16 vect1(neq), vect2(neq), trav(ltrav)
      complex*16 d(nval+1), shift2, vp(neq,nval)
c
      double precision tol
c
c      integer*8 max_nvect
c      parameter(max_nvect=3000) ! previously 1500
c
      integer alloc_stat
      complex*16, dimension(:), allocatable :: workev !  (3*max_nvect), 
      double precision, dimension(:), allocatable :: rwork  !  (max_nvect)
      logical, dimension(:), allocatable :: select  !  (max_nvect)


c     Local variables
c     32-bit integers for ARPACK
      integer*4 neq_32, nval_32, nvect_32
      integer*4 ido_32, info_32, ierr_32, iparam_32(11)
      integer*4 ipntr_32(14), ltrav_32
c
      logical rvec
      character bmat*1, which*2

c      data bmat/'G'/
c      data which/'SR'/
c      data which/'SM'/
      data bmat/'I'/
      data which/'LM'/
c
      integer*8 ui, debug
c      common/imp/ui, debug
c
c     ------------------------------------------------------------------
c
      ui = 6
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      alloc_stat = 0
      allocate(workev(3*nvect), rwork(nvect), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "VALPR_64: Mem. allocation is unseccesfull"
        write(*,*) "for the arrays workev, rwork"
        write(*,*) "alloc_stat, nvect = ", alloc_stat, nvect
        write(*,*) "Aborting..."
        stop
      endif

      allocate(select(nvect), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "VALPR_64: Mem. allocation is unseccesfull"
        write(*,*) "for the array select"
        write(*,*) "alloc_stat, nvect = ", alloc_stat, nvect
        write(*,*) "Aborting..."
        stop
      endif

c
      shift2 = (0.0d0,0.0d0)
c
      do i=1,nval
        d(i) = 0.0d0
      enddo
c
c       ##################################################################
c
      if (i_base .ne. 0) then
        write(ui,*) "valpr_64: i_base != 0 : ", i_base
        write(ui,*) "valpr_64: UMFPACK requires 0-based indexing"
        write(ui,*) "valpr_64l: Aborting..."
        stop
      endif
c
c       ----------------------------------------------------------------
c       factor the matrix and save to a file
c       ----------------------------------------------------------------
c
      if (debug .eq. 1) then
        write(ui,*) "valpr_64: factorisation (UMFPACK)"
      endif
      call cpu_time(time1_fact)
      ls_data(1) = time1_fact
c
c       set default parameters
        call umf4zdef (control)

c
c     umfpack * report status (print level = control(1)) :
c     print level = 0 or less : No output, even when an error occurs.
c     print level = 1 (default value) : then error messages are printed, 
c                      and nothing is printed if the status is UMFPACK OK.
c     print level = 2 or more : then the status is always printed.
c
c       print control parameters.  set control (1) to 1 to print
c       error messages only
        control (1) = 1
        call umf4zpcon (control)

c       pre-order and symbolic analysis
        call umf4zsym (neq, neq, col_ptr, row_ind, 
     *         mat1_re, mat1_im, symbolic, control, info_umf)
c
c       print statistics computed so far
c       call umf4zpinf (control, info_umf) could also be done.
      if (debug .eq. 1) then
        write(ui,80) info_umf (1), info_umf (16),
     $      (info_umf (21) * info_umf (4)) / 2**20,
     $      (info_umf (22) * info_umf (4)) / 2**20,
     $      info_umf (23), info_umf (24), info_umf (25)
80      format ('  symbolic analysis:',/,
     $      '   status:  ', f5.0, /,
     $      '   time:    ', e10.2, ' (sec)'/,
     $      '   estimates (upper bound) for numeric LU:', /,
     $      '   size of LU:    ', f12.2, ' (MB)', /,
     $      '   memory needed: ', f12.2, ' (MB)', /,
     $      '   flop count:    ', e12.2, /
     $      '   nnz (L):       ', f12.0, /
     $      '   nnz (U):       ', f12.0)
      endif

c       check umf4zsym error condition
        if (info_umf (1) .lt. 0) then
            write(ui,*) 'Error occurred in umf4zsym: ', info_umf (1)
            stop
        endif

c       numeric factorization
        call umf4znum (col_ptr, row_ind, mat1_re, 
     *         mat1_im, symbolic, numeric, control, info_umf)

c       print statistics for the numeric factorization
c       call umf4zpinf (control, info_umf) could also be done.
      if (debug .eq. 1) then
        write(ui,90) info_umf (1), info_umf (66),
     $      (info_umf (41) * info_umf (4)) / 2**20,
     $      (info_umf (42) * info_umf (4)) / 2**20,
     $      info_umf (43), info_umf (44), info_umf (45)
90      format ('  numeric factorization:',/,
     $      '   status:  ', f5.0, /,
     $      '   time:    ', e10.2, /,
     $      '   actual numeric LU statistics:', /,
     $      '   size of LU:    ', f12.2, ' (MB)', /,
     $      '   memory needed: ', f12.2, ' (MB)', /,
     $      '   flop count:    ', e12.2, /
     $      '   nnz (L):       ', f12.0, /
     $      '   nnz (U):       ', f12.0)
      endif

c       check umf4znum error condition
        if (info_umf (1) .lt. 0) then
            write(ui,*) 'Error occurred in umf4znum: ', info_umf (1)
            stop
        endif

c       save the symbolic analysis to the file s42.umf
c       note that this is not needed until another matrix is
c       factorized, below.
c	filenum = 42
c        call umf4zssym (symbolic, filenum, status)
c        if (status .lt. 0) then
c            write(ui,*) 'Error occurred in umf4zssym: ', status
c            stop
c        endif
c
c       save the LU factors to the file n0.umf
c        call umf4zsnum (numeric, filenum, status)
c        if (status .lt. 0) then
c            write(ui,*) 'Error occurred in umf4zsnum: ', status
c            stop
c        endif

c       free the symbolic analysis
        call umf4zfsym (symbolic)

cc       free the numeric factorization
c        call umf4zfnum (numeric)
c
      call cpu_time(time2_fact)
      ls_data(2) = time2_fact
      if (debug .eq. 1) then
        write(ui,*) "valpr_64: factorisation completed"
        write(ui,*) "LU factorisation : CPU time = ",  
     *         (time2_fact-time1_fact)
c , 
c     *         100*(time2_fact-time1_fact)/(time2-time1),"%"
      endif
c
c       No LU factors (symbolic or numeric) are in memory at this point.
c
cc       ----------------------------------------------------------------
cc       load the LU factors back in, and solve the system
cc       ----------------------------------------------------------------
c
cc       At this point the program could terminate and load the LU
cC       factors (numeric) from the n0.umf file, and solve the
cc       system (see below).  Note that the symbolic object is not
cc       required.
c
cc       load the numeric factorization back in (filename: n0.umf)
c        call umf4zlnum (numeric, filenum, status)
c        if (status .lt. 0) then
c            write(ui,*) 'Error occurred in umf4zlnum: ', status
c            stop
c        endif
c
c       ##################################################################
c
c     On commence le travail avec znaupd
c     ----------------------------------
c
c      ------------------------------------------------------------ 
c     | Le parametre IDO_32 est utilise pour la communication.        |
c     | A l'etape initiale il doit valoir 0.                       |
c     | Le choix INFO_32=0 correspond a la construction par           |
c     | Arpack d'un vecteur initial a composantes aleatoires       |
c     | iparam_32(1)=1 signifie que Arpack calcule les translations   |
c     | a partir de la matrice projetee et conformement au critere |
c     | "which".                                                   |
c      ------------------------------------------------------------

      neq_32 = neq
      nval_32 = nval
      nvect_32 = nvect
      ltrav_32 = ltrav
      ido_32 = 0 
      iparam_32(1) = 1
      iparam_32(3) = itermax
c      iparam_32(7) = 3
      iparam_32(7) = 1
      info_32 = 0
ccccccccccccccccccc
c
      compteur = 0
 
c      ----------------------------------------------------
c     | Boucle principale en mode de communication inverse | 
c      ----------------------------------------------------
c
      call cpu_time(time1_arpack)
      ls_data(3) = time1_arpack
c
20    continue
c
      call znaupd (ido_32, bmat, neq_32, which, nval_32, tol, 
     *             resid, nvect_32, v, neq_32, iparam_32, 
     *             ipntr_32, workd, trav, ltrav_32, rwork, info_32)
c
      compteur = compteur + 1

c      if (ido_32.eq.-1) then
         if (ido_32 .eq. -1 .or. ido_32 .eq. 1) then

c      ------------------------------------------------------
c     | On execute  y <--- OP*x = inv[A-SIGMA*M]*M*x         |
c     | pour obtenir un vecteur de depart dans l'image de OP |
c     | x = workd(ipntr_32(1)) et y = workd(ipntr_32(2))           |
c      ------------------------------------------------------

         call zcopy(neq_32, workd(ipntr_32(1)), 1,vect1, 1)
         call z_mxv_csc (neq, vect1, vect2, nonz, row_ind, 
     *     col_ptr, mat2)
c
         do i=1,neq
           rhs_re(i) = dble(vect2(i))
           rhs_im(i) = imag(vect2(i))
         enddo
c
c       solve Ax=b, without iterative refinement
        sys = 0
        call umf4zsol (sys, lhs_re, lhs_im, rhs_re, rhs_im, 
     *     numeric, control, info_umf)
        if (info_umf (1) .lt. 0) then
            write(ui,*) 'Error occurred in umf4zsol: ', info_umf (1)
            stop
        endif
        do i=1,neq
          vect2 (i) = dcmplx (lhs_re (i), lhs_im (i))
        enddo
c
         call zcopy(neq_32, vect2, 1, workd(ipntr_32(2)), 1)
         go to 20
c            
         else if (ido_32.eq.2) then
c
         write(ui,*) 'VALPR_64: ATTENTION ido_32 = ', ido_32
         write(ui,*) 'check the results...'

c           ----------------------------------------------
c          | On execute  y <--- M*x                       |
c          | x = workd(ipntr_32(1))  et  y = workd(ipntr_32(2)) |
c           ----------------------------------------------

            call zcopy(neq_32, workd(ipntr_32(1)), 1, vect1, 1)
            call z_mxv_csc (neq, vect1, vect2, nonz, row_ind, 
     *        col_ptr, mat2)
c
         do i=1,neq
           rhs_re(i) = dble(vect2(i))
           rhs_im(i) = imag(vect2(i))
         enddo
c
c       solve Ax=b, without iterative refinement
        sys = 0
        call umf4zsol (sys, lhs_re, lhs_im, rhs_re, rhs_im, 
     *     numeric, control, info_umf)
        if (info_umf (1) .lt. 0) then
            write(ui,*) 'Error occurred in umf4zsol: ', info_umf (1)
            stop
        endif
        do i=1,neq
          vect2 (i) = dcmplx (lhs_re (i), lhs_im (i))
        enddo
            call zcopy(neq_32, vect2,1, workd(ipntr_32(2)), 1)
            go to 20
           
         end if 

c      ---------------------------------------------------
c     | Either we have convergence, or there is an error. |
c      ---------------------------------------------------

      n_conv = iparam_32(5)

      if (info_32 .gt. 0) then
        write(ui,*)
        write(ui,*) "VALPR_64: info_32 != 0 : ", info_32
        write(ui,*) "VALPR_64: iparam_32(5) = ", iparam_32(5), nval_32
        write(ui,*) "VALPR_64: number of converged values = ", 
     *                iparam_32(5)
        write(ui,*)
      endif

      if (info_32.lt.0) then

c      ---------------------------------------------------
c     | Error message, check the documentation in DNAUPD. |

c  INFO    Integer.  (INPUT/OUTPUT)
c          If INFO .EQ. 0, a randomly initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          Error flag on output.
c          =  0: Normal exit.
c          =  1: Maximum number of iterations taken.
c                All possible eigenvalues of OP has been found. IPARAM(5)  
c                returns the number of wanted converged Ritz values.
c          =  2: No longer an informational error. Deprecated starting
c                with release 2 of ARPACK.
c          =  3: No shifts could be applied during a cycle of the 
c                Implicitly restarted Arnoldi iteration. One possibility 
c                is to increase the size of NCV relative to NEV. 
c                See remark 4 below.
c          = -1: N must be positive.
c          = -2: NEV must be positive.
c          = -3: NCV-NEV >= 2 and less than or equal to N.
c          = -4: The maximum number of Arnoldi update iteration 
c                must be greater than zero.
c          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
c          = -6: BMAT must be one of 'I' or 'G'.
c          = -7: Length of private work array is not sufficient.
c          = -8: Error return from LAPACK eigenvalue calculation;
c          = -9: Starting vector is zero.
c          = -10: IPARAM(7) must be 1,2,3,4.
c          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatable.
c          = -12: IPARAM(1) must be equal to 0 or 1.
c          = -9999: Could not build an Arnoldi factorization.
c                   IPARAM(5) returns the size of the current Arnoldi
c                   factorization.

c      ---------------------------------------------------

         write(ui,*) 'VALPR_64:'
         write(ui,*) ' Error with _naupd, info_32 = ', info_32
         write(ui,*) ' Check the documentation in _naupd.'
         write(ui,*) 'Aborting...'
         stop

      else 

c      -------------------------------------
c     | Ici on recupere les valeurs propres |
c      -------------------------------------

         rvec = .true.

         call zneupd (rvec, 'A', select, d, v, neq_32, shift2, 
     *                workev, bmat, neq_32, which, nval_32, tol, 
     *                resid, nvect_32, v, neq_32, iparam_32, ipntr_32, 
     *                workd, trav, ltrav_32, rwork, ierr_32)        
c      ------------------------------------------------------------
c     | La partie reelle d'une valeur propre se trouve dans la     |
c     | premiere colonne du tableau D, la partie imaginaire est    |
c     | dans la seconde.                                           |
c     | Les vecteurs propres sont dans les premieres nval_32 colonnes |
c     | du tableau V, lorsque demande (rvec=.true.). Sinon, on y   |
c     | trouve une base orthogonale de l'espace propre.            |
c      ------------------------------------------------------------

         if (ierr_32.ne.0) then

c      -----------------------------------------------------
c     | Error condition: Check the documentation of DNEUPD. |
c      -----------------------------------------------------

            write(ui,*) 'VALPR_64:' 
            write(ui,*) ' Error with _neupd, info_32 = ', ierr_32
            write(ui,*) ' Check the documentation of _neupd. '
            write(ui,*) 'Aborting...'
            stop

         else
           do i = 1, nval
             do j = 1, neq
               vp(j,i) = v(j,i)
             enddo
           enddo
         endif
      endif
c
c       free the numeric factorization
        call umf4zfnum (numeric)
c
      call cpu_time(time2_arpack)
      ls_data(4) = time2_arpack
c
c      if (debug .eq. 1) then
c        do i=1,nval
c          write (*,*) "i, d(i) = ", i, d(i)
c        enddo
c      endif
c
      deallocate(workev, rwork, STAT=alloc_stat)
      deallocate(select)
c
      return 
      end
