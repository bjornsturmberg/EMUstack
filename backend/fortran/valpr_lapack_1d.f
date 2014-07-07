c     ------------------------------------------------------------------
c
c                    sous-routine valpr_64_1d.
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
      subroutine valpr_64_1d (nvect, nval, neq, itermax, 
     *  tol, mat1, mat2, v, d, vp,   
     *  n_conv, ls_data, debug)

c     rhs_re, rhs_im, lhs_re, lhs_im, 
c     vect1, vect2, workd, resid, ltrav, 
c     i_base, trav, 
c     nonz, row_ind, col_ptr, 
c     mat1_re, mat1_im, 
c     numeric, filenum, status, control, info_umf, 

c
c     ------------------------------------------------------------------
c
      implicit none
c
      integer*8 neq, n_conv
      complex*16 mat1(neq,neq), mat2(neq,neq)
      complex*16 v(neq,nvect), shift2
      complex*16 d(nval+1), vp(neq,nval)
      double precision ls_data(10)

c
c      integer*8 row_ind(nonz), col_ptr(neq+1)
c     nonz, i_base
c      double precision mat1_re(nonz), mat1_im(nonz)
c      double precision rhs_re(neq), rhs_im(neq)
c      double precision lhs_re(neq), lhs_im(neq)

c
      double precision time1_fact, time2_fact
      double precision time1_arpack, time2_arpack
c
c      double precision control (20), info_umf (90)
c      integer*8 numeric, symbolic, status, sys, filenum
c
      integer*8 itermax, nvect, nval, i, j, col1, col2
      integer*8 compteur


c      complex*16 energ
c
      complex*16 z_tmp0, z_tmp
      double precision tol
c
c      integer*8 max_nvect
c      parameter(max_nvect=3000) ! previously 1500
c
      integer alloc_stat
      complex*16, dimension(:), allocatable :: workev !  (3*max_nvect), 
      double precision, dimension(:), allocatable :: rwork  !  (max_nvect)
      logical , dimension(:), allocatable :: select  !  (max_nvect)


c     Local variables
      integer*8 ltrav
      complex*16, allocatable ::  vect1(:), vect2(:), resid(:)
      complex*16, allocatable ::  trav(:), workd(:)


c     32-b1t integers for ARPACK
      integer*4 neq_32, nval_32, nvect_32
      integer*4 ido_32, info_32, ierr_32, iparam_32(11)
      integer*4 ipntr_32(14), ltrav_32
c     32-b1t integers for LAPACK
      integer*4, dimension(:), allocatable :: IPIV_32   !  (nval_max_32)
      integer*4 INFO_32_Lapk, LDA_32, LDB_32
      complex*16 ZERO, ONE

 

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
      ZERO = 0.0d0
      ONE = 1.0d0

      ltrav = 3*nvect*(nvect+2)

      if (debug .eq. 1) then
        write(ui,*) "valpr_64_1d: nvect, nval = ", nvect, nval
      endif

c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      allocate(IPIV_32(neq), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "valpr_64_1d: Mem. allocation is unseccesfull"
        write(*,*) "alloc_stat (IPIV_32) = ", alloc_stat
        write(*,*) "neq = ", neq
        write(*,*) "Aborting..."
        stop
      endif

      allocate(vect1(neq), vect2(neq), resid(neq), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "valpr_64_1d: Mem. allocation is unseccesfull"
        write(*,*) "alloc_stat (vect1, vect2,...) = ", alloc_stat
        write(*,*) "neq = ", neq
        write(*,*) "Aborting..."
        stop
      endif

      allocate(trav(ltrav), workd(3*neq), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "valpr_64_1d: Mem. allocation is unseccesfull"
        write(*,*) "alloc_stat (trav,workd) = ", alloc_stat
        write(*,*) "neq, ltrav = ", neq, ltrav
        write(*,*) "Aborting..."
        stop
      endif
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      alloc_stat = 0
      allocate(workev(3*nvect), rwork(nvect), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "valpr_64_1d: Mem. allocation is unseccesfull"
        write(*,*) "for the arrays workev, rwork"
        write(*,*) "alloc_stat, nvect = ", alloc_stat, nvect
        write(*,*) "Aborting..."
        stop
      endif

      allocate(select(nvect), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "valpr_64_1d: Mem. allocation is unseccesfull"
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
c      if (i_base .ne. 0) then
c        write(ui,*) "valpr_64_1d: i_base != 0 : ", i_base
c        write(ui,*) "valpr_64_1d: UMFPACK requires 0-based indexing"
c        write(ui,*) "valpr_64_1d: Aborting..."
c        stop
c      endif
c
c       ----------------------------------------------------------------
c       factor the matrix and save to a file
c       ----------------------------------------------------------------
c
c
c
ccccccccccccccccccccccccccccccccccccccccccc
c       set default parameters
c        call umf4zdef (control)
c
c
c     umfpack * report status (print level = control(1)) :
c     print level = 0 or less : No output, even when an error occurs.
c     print level = 1 (default value) : then error messages are printed, 
c                      and nothing is printed if the status is UMFPACK OK.
c     print level = 2 or more : then the status is always printed.
c
c       print control parameters.  set control (1) to 1 to print
c       error messages only
c        control (1) = 1
c        call umf4zpcon (control)
c
c       pre-order and symbolic analysis
c        call umf4zsym (neq, neq, col_ptr, row_ind, 
c     *         mat1_re, mat1_im, symbolic, control, info_umf)
c
c       print statistics computed so far
c       call umf4zpinf (control, info_umf) could also be done.
c      if (debug .eq. 1) then
c        write(ui,80) info_umf (1), info_umf (16),
c     $      (info_umf (21) * info_umf (4)) / 2**20,
c     $      (info_umf (22) * info_umf (4)) / 2**20,
c     $      info_umf (23), info_umf (24), info_umf (25)
c80      format ('  symbolic analysis:',/,
c     $      '   status:  ', f5.0, /,
c     $      '   time:    ', e10.2, ' (sec)'/,
c     $      '   estimates (upper bound) for numeric LU:', /,
c     $      '   size of LU:    ', f12.2, ' (MB)', /,
c     $      '   memory needed: ', f12.2, ' (MB)', /,
c     $      '   flop count:    ', e12.2, /
c     $      '   nnz (L):       ', f12.0, /
c     $      '   nnz (U):       ', f12.0)
c      endif
c
cc       check umf4zsym error condition
c        if (info_umf (1) .lt. 0) then
c            write(ui,*) 'Error occurred in umf4zsym: ', info_umf (1)
c            stop
c        endif
c
cc       numeric factorization
c        call umf4znum (col_ptr, row_ind, mat1_re, 
c     *         mat1_im, symbolic, numeric, control, info_umf)
c
cc       print statistics for the numeric factorization
cc       call umf4zpinf (control, info_umf) could also be done.
c      if (debug .eq. 1) then
c        write(ui,90) info_umf (1), info_umf (66),
c     $      (info_umf (41) * info_umf (4)) / 2**20,
c     $      (info_umf (42) * info_umf (4)) / 2**20,
c     $      info_umf (43), info_umf (44), info_umf (45)
c90      format ('  numeric factorization:',/,
c     $      '   status:  ', f5.0, /,
c     $      '   time:    ', e10.2, /,
c     $      '   actual numeric LU statistics:', /,
c     $      '   size of LU:    ', f12.2, ' (MB)', /,
c     $      '   memory needed: ', f12.2, ' (MB)', /,
c     $      '   flop count:    ', e12.2, /
c     $      '   nnz (L):       ', f12.0, /
c     $      '   nnz (U):       ', f12.0)
c      endif
c
cc       check umf4znum error condition
c        if (info_umf (1) .lt. 0) then
c            write(ui,*) 'Error occurred in umf4znum: ', info_umf (1)
c            stop
c        endif
ccccccccccccccccccccccccccccccccccccccccccc
c
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
c
c       free the symbolic analysis
c        call umf4zfsym (symbolic)
c
cc       free the numeric factorization
c        call umf4zfnum (numeric)
c
ccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
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

c
ccccccccccccccccccccccccccccccccccccccccccc
c
      if (debug .eq. 1) then
        write(ui,*) "valpr_64_1d: factorisation (LAPACK)"
      endif

c
c      if (debug .eq. 1) then
c        write(ui,*) "valpr_64_1d: Writing mat1 into a file"
c        open (unit=24,file="Output/mat1_a.txt",status="unknown")
c        do i=1,neq
c          do j=1,neq
c            write(24,*) i, j, mat1(i,j) 
c          enddo
c        enddo
c        close(24)
c      endif
c
c      if (debug .eq. 1) then
c        write(ui,*) "valpr_64_1d: Writing mat2 into a file"
c        open (unit=24,file="Output/mat2_a.txt",status="unknown")
c        do i=1,neq
c          do j=1,neq
c            write(24,*) i, j, mat2(i,j) 
c          enddo
c        enddo
c        close(24)
c      endif

c       LU Factorization
c     Compute the LU factorization of A.
c
      LDA_32 = neq
c      neq_32 = neq
      call cpu_time(time1_fact)
      ls_data(1) = time1_fact
      CALL ZGETRF( neq_32, neq_32, mat1, LDA_32, IPIV_32, INFO_32_Lapk)
      call cpu_time(time2_fact)

c      if (debug .eq. 1) then
c        write(ui,*) "valpr_64_1d: Writing mat1 (fact.) into a file"
c        open (unit=24,file="Output/mat1_b.txt",status="unknown")
c        do i=1,neq
c          do j=1,neq
c            write(24,*) i, j, mat1(i,j) 
c          enddo
c        enddo
c        close(24)
c      endif

      ls_data(2) = time2_fact
       if(INFO_32_Lapk .ne. 0) then
         write(*,*) "valpr_64_1d: PROBLEM INFO_32_Lapk = ", INFO_32_Lapk
         write(*,*) "Aborting... pos 1"
         stop
       endif
      if (debug .eq. 1) then
        write(ui,*) "valpr_64_1d: factorisation completed"
        write(ui,*) "LU factorisation : CPU time = ",  
     *         (time2_fact-time1_fact)
      endif

c      if (debug .eq. 1) then
c        write(ui,*) "valpr_64_1d: "
c        do i=1,10
c        write(ui,*) "(B) valpr_64_1d: mi, mat1(i,i) = ", i, mat1(5,i)
c        enddo
c        write(ui,*)
c        do i=1,10
c        write(ui,*) "(B) valpr_64_1d: mi, mat2(i,i) = ", i, mat2(5,i)
c        enddo
c      endif

c
c
ccccccccccccccccccccccccccccccccccccccccccc
c

      if (debug .eq. 1) then
        write(ui,*) "valpr_64_1d: call to ARPACK"
      endif
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
c         call z_mxv_csc (neq, vect1, vect2, nonz, row_ind, 
c     *     col_ptr, mat2)

c      if (debug .eq. 1) then
c           write(ui,*) "valpr_64_1d:", compteur
c         do i=1,10
c           write(ui,*) "(A): vect1: ", vect1(i)
c         enddo
c      endif

         call ZGEMV('N', neq_32, neq_32, ONE, mat2, LDA_32, 
     *     vect1, 1, ZERO, vect2, 1)

c      if (debug .eq. 1) then
c           write(ui,*) "valpr_64_1d:", compteur
c         do i=1,10
c           write(ui,*) "(B): vect2: ", vect2(i)
c         enddo
c      endif

c
c        Solve the system mat1*X = B, overwriting B with X.
c
         LDB_32 = neq_32
         CALL ZGETRS( 'N', neq_32, 1, mat1, LDA_32, IPIV_32, vect2, 
     *     LDB_32, INFO_32_Lapk)
         if(INFO_32_Lapk .ne. 0) then
           write(*,*) "valpr_64_1d: PROBLEM INFO_32_Lapk = ", 
     *              INFO_32_Lapk
           write(*,*) "Aborting... pos 2"
           stop
         endif

c      if (debug .eq. 1) then
c           write(ui,*) "valpr_64_1d:"
c         do i=1,10
c           write(ui,*) "Ratio (C): ", compteur, vect2(i)
c         enddo
c         stop
c      endif

c        write(ui,*) "LU factorisation : CPU time = ",  
c     *         (time2_fact-time1_fact)
c
c
c         do i=1,neq
c           rhs_re(i) = dble(vect2(i))
c           rhs_im(i) = imag(vect2(i))
c         enddo
c
c       solve Ax=b, without iterative refinement
c        sys = 0
c        call umf4zsol (sys, lhs_re, lhs_im, rhs_re, rhs_im, 
c     *     numeric, control, info_umf)
c        if (info_umf (1) .lt. 0) then
c            write(ui,*) 'Error occurred in umf4zsol: ', info_umf (1)
c            stop
c        endif
c        do i=1,neq
c          vect2 (i) = dcmplx (lhs_re (i), lhs_im (i))
c        enddo
c
         call zcopy(neq_32, vect2, 1, workd(ipntr_32(2)), 1)
         go to 20
c            
         else if (ido_32.eq.2) then
c
         write(ui,*) 'valpr_64_1d: ATTENTION ido_32 = ', ido_32
         write(ui,*) 'check the results...'

c           ----------------------------------------------
c          | On execute  y <--- M*x                       |
c          | x = workd(ipntr_32(1))  et  y = workd(ipntr_32(2)) |
c           ----------------------------------------------

            call zcopy(neq_32, workd(ipntr_32(1)), 1, vect1, 1)
         call ZGEMV('N', neq_32, neq_32, ONE, mat2, LDA_32, 
     *     vect1, 1, ZERO, vect2, 1)

c            call z_mxv_csc (neq, vect1, vect2, nonz, row_ind, 
c     *        col_ptr, mat2)
c
c
c        Solve the system mat1*X = B, overwriting B with X.
c
         LDB_32 = neq_32
         CALL ZGETRS( 'N', neq_32, 1, mat1, LDA_32, IPIV_32, vect2, 
     *     LDB_32, INFO_32 )
         if(INFO_32 .ne. 0) then
           write(*,*) "valpr_64_1d: PROBLEM INFO_32 = ", INFO_32
           write(*,*) "Aborting... pos 1"
           stop
         endif
c
c         do i=1,neq
c           rhs_re(i) = dble(vect2(i))
c           rhs_im(i) = imag(vect2(i))
c         enddo
cc
cc       solve Ax=b, without iterative refinement
c        sys = 0
c        call umf4zsol (sys, lhs_re, lhs_im, rhs_re, rhs_im, 
c     *     numeric, control, info_umf)
c        if (info_umf (1) .lt. 0) then
c            write(ui,*) 'Error occurred in umf4zsol: ', info_umf (1)
c            stop
c        endif
c        do i=1,neq
c          vect2 (i) = dcmplx (lhs_re (i), lhs_im (i))
c        enddo
            call zcopy(neq_32, vect2,1, workd(ipntr_32(2)), 1)
            go to 20
           
         end if 

c      ---------------------------------------------------
c     | Either we have convergence, or there is an error. |
c      ---------------------------------------------------

      n_conv = iparam_32(5)

      if (info_32 .gt. 0) then
        write(ui,*)
        write(ui,*) "valpr_64_1d: info_32 != 0 : ", info_32
        write(ui,*) "valpr_64_1d: iparam_32(5) = ", iparam_32(5), 
     *                 nval_32
        write(ui,*) "valpr_64_1d: number of converged values = ", 
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

         write(ui,*) 'valpr_64_1d:'
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

            write(ui,*) 'valpr_64_1d:' 
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
c        call umf4zfnum (numeric)
c
      call cpu_time(time2_arpack)
      ls_data(4) = time2_arpack
c
cc      if (debug .eq. 1) then
cc       do i=1,nval
cc          write (*,*) "i, d(i) = ", i, d(i)
cc        enddo
cc      endif
c
      deallocate(workev, rwork, STAT=alloc_stat)
      deallocate(select)
c
      return 
      end
