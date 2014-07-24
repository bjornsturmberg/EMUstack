CCC A collection of subroutines from ARPACK 2.4 1996. It includes;
C znaupd
C zneupd
C znaup2
C zstatn
C zngets
C zgetv0
C znaitr
C zneigh
C zsortc
C znapps

c\BeginDoc
c
c\Name: znaupd
c
c\Description: 
c  Reverse communication interface for the Implicitly Restarted Arnoldi
c  iteration. This is intended to be used to find a few eigenpairs of a 
c  complex linear operator OP with respect to a semi-inner product defined 
c  by a hermitian positive semi-definite real matrix B. B may be the identity 
c  matrix.  NOTE: if both OP and B are real, then dsaupd or dnaupd should
c  be used.
c
c
c  The computed approximate eigenvalues are called Ritz values and
c  the corresponding approximate eigenvectors are called Ritz vectors.
c
c  znaupd is usually called iteratively to solve one of the 
c  following problems:
c
c  Mode 1:  A*x = lambda*x.
c           ===> OP = A  and  B = I.
c
c  Mode 2:  A*x = lambda*M*x, M hermitian positive definite
c           ===> OP = inv[M]*A  and  B = M.
c           ===> (If M can be factored see remark 3 below)
c
c  Mode 3:  A*x = lambda*M*x, M hermitian semi-definite
c           ===> OP =  inv[A - sigma*M]*M   and  B = M. 
c           ===> shift-and-invert mode 
c           If OP*x = amu*x, then lambda = sigma + 1/amu.
c  
c
c  NOTE: The action of w <- inv[A - sigma*M]*v or w <- inv[M]*v
c        should be accomplished either by a direct method
c        using a sparse matrix factorization and solving
c
c           [A - sigma*M]*w = v  or M*w = v,
c
c        or through an iterative method for solving these
c        systems.  If an iterative method is used, the
c        convergence test must be more stringent than
c        the accuracy requirements for the eigenvalue
c        approximations.
c
c\Usage:
c  call znaupd
c     ( IDO, BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM,
c       IPNTR, WORKD, WORKL, LWORKL, RWORK, INFO )
c
c\Arguments
c  IDO     Integer.  (INPUT/OUTPUT)
c          Reverse communication flag.  IDO must be zero on the first 
c          call to znaupd.  IDO will be set internally to
c          indicate the type of operation to be performed.  Control is
c          then given back to the calling routine which has the
c          responsibility to carry out the requested operation and call
c          znaupd with the result.  The operand is given in
c          WORKD(IPNTR(1)), the result must be put in WORKD(IPNTR(2)).
c          -------------------------------------------------------------
c          IDO =  0: first call to the reverse communication interface
c          IDO = -1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c                    This is for the initialization phase to force the
c                    starting vector into the range of OP.
c          IDO =  1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c                    In mode 3, the vector B * X is already
c                    available in WORKD(ipntr(3)).  It does not
c                    need to be recomputed in forming OP * X.
c          IDO =  2: compute  Y = M * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c          IDO =  3: compute and return the shifts in the first 
c                    NP locations of WORKL.
c          IDO = 99: done
c          -------------------------------------------------------------
c          After the initialization phase, when the routine is used in 
c          the "shift-and-invert" mode, the vector M * X is already 
c          available and does not need to be recomputed in forming OP*X.
c             
c  BMAT    Character*1.  (INPUT)
c          BMAT specifies the type of the matrix B that defines the
c          semi-inner product for the operator OP.
c          BMAT = 'I' -> standard eigenvalue problem A*x = lambda*x
c          BMAT = 'G' -> generalized eigenvalue problem A*x = lambda*M*x
c
c  N       Integer.  (INPUT)
c          Dimension of the eigenproblem.
c
c  WHICH   Character*2.  (INPUT)
c          'LM' -> want the NEV eigenvalues of largest magnitude.
c          'SM' -> want the NEV eigenvalues of smallest magnitude.
c          'LR' -> want the NEV eigenvalues of largest real part.
c          'SR' -> want the NEV eigenvalues of smallest real part.
c          'LI' -> want the NEV eigenvalues of largest imaginary part.
c          'SI' -> want the NEV eigenvalues of smallest imaginary part.
c
c  NEV     Integer.  (INPUT)
c          Number of eigenvalues of OP to be computed. 0 < NEV < N-1.
c
c  TOL     Double precision  scalar.  (INPUT)
c          Stopping criteria: the relative accuracy of the Ritz value 
c          is considered acceptable if BOUNDS(I) .LE. TOL*ABS(RITZ(I))
c          where ABS(RITZ(I)) is the magnitude when RITZ(I) is complex.
c          DEFAULT = dlamch('EPS')  (machine precision as computed
c                    by the LAPACK auxiliary subroutine dlamch).
c
c  RESID   Complex*16 array of length N.  (INPUT/OUTPUT)
c          On INPUT: 
c          If INFO .EQ. 0, a random initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          On OUTPUT:
c          RESID contains the final residual vector.
c
c  NCV     Integer.  (INPUT)
c          Number of columns of the matrix V. NCV must satisfy the two
c          inequalities 1 <= NCV-NEV and NCV <= N.
c          This will indicate how many Arnoldi vectors are generated 
c          at each iteration.  After the startup phase in which NEV 
c          Arnoldi vectors are generated, the algorithm generates 
c          approximately NCV-NEV Arnoldi vectors at each subsequent update 
c          iteration. Most of the cost in generating each Arnoldi vector is 
c          in the matrix-vector operation OP*x. (See remark 4 below.)
c
c  V       Complex*16 array N by NCV.  (OUTPUT)
c          Contains the final set of Arnoldi basis vectors. 
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling program.
c
c  IPARAM  Integer array of length 11.  (INPUT/OUTPUT)
c          IPARAM(1) = ISHIFT: method for selecting the implicit shifts.
c          The shifts selected at each iteration are used to filter out
c          the components of the unwanted eigenvector.
c          -------------------------------------------------------------
c          ISHIFT = 0: the shifts are to be provided by the user via
c                      reverse communication.  The NCV eigenvalues of 
c                      the Hessenberg matrix H are returned in the part
c                      of WORKL array corresponding to RITZ.
c          ISHIFT = 1: exact shifts with respect to the current
c                      Hessenberg matrix H.  This is equivalent to 
c                      restarting the iteration from the beginning 
c                      after updating the starting vector with a linear
c                      combination of Ritz vectors associated with the 
c                      "wanted" eigenvalues.
c          ISHIFT = 2: other choice of internal shift to be defined.
c          -------------------------------------------------------------
c
c          IPARAM(2) = No longer referenced 
c
c          IPARAM(3) = MXITER
c          On INPUT:  maximum number of Arnoldi update iterations allowed. 
c          On OUTPUT: actual number of Arnoldi update iterations taken. 
c
c          IPARAM(4) = NB: blocksize to be used in the recurrence.
c          The code currently works only for NB = 1.
c
c          IPARAM(5) = NCONV: number of "converged" Ritz values.
c          This represents the number of Ritz values that satisfy
c          the convergence criterion.
c
c          IPARAM(6) = IUPD
c          No longer referenced. Implicit restarting is ALWAYS used.  
c
c          IPARAM(7) = MODE
c          On INPUT determines what type of eigenproblem is being solved.
c          Must be 1,2,3; See under \Description of znaupd for the 
c          four modes available.
c
c          IPARAM(8) = NP
c          When ido = 3 and the user provides shifts through reverse
c          communication (IPARAM(1)=0), _naupd returns NP, the number
c          of shifts the user is to provide. 0 < NP < NCV-NEV.
c
c          IPARAM(9) = NUMOP, IPARAM(10) = NUMOPB, IPARAM(11) = NUMREO,
c          OUTPUT: NUMOP  = total number of OP*x operations,
c                  NUMOPB = total number of B*x operations if BMAT='G',
c                  NUMREO = total number of steps of re-orthogonalization.
c
c  IPNTR   Integer array of length 14.  (OUTPUT)
c          Pointer to mark the starting locations in the WORKD and WORKL
c          arrays for matrices/vectors used by the Arnoldi iteration.
c          -------------------------------------------------------------
c          IPNTR(1): pointer to the current operand vector X in WORKD.
c          IPNTR(2): pointer to the current result vector Y in WORKD.
c          IPNTR(3): pointer to the vector B * X in WORKD when used in 
c                    the shift-and-invert mode.
c          IPNTR(4): pointer to the next available location in WORKL
c                    that is untouched by the program.
c          IPNTR(5): pointer to the NCV by NCV upper Hessenberg
c                    matrix H in WORKL.
c          IPNTR(6): pointer to the  ritz value array  RITZ
c          IPNTR(7): pointer to the (projected) ritz vector array Q
c          IPNTR(8): pointer to the error BOUNDS array in WORKL.
c          IPNTR(14): pointer to the NP shifts in WORKL. See Remark 5 below.
c
c          Note: IPNTR(9:13) is only referenced by zneupd. See Remark 2 below.
c
c          IPNTR(9): pointer to the NCV RITZ values of the 
c                    original system.
c          IPNTR(10): Not Used
c          IPNTR(11): pointer to the NCV corresponding error bounds.
c          IPNTR(12): pointer to the NCV by NCV upper triangular
c                     Schur matrix for H.
c          IPNTR(13): pointer to the NCV by NCV matrix of eigenvectors
c                     of the upper Hessenberg matrix H. Only referenced by
c                     zneupd if RVEC = .TRUE. See Remark 2 below.
c
c          -------------------------------------------------------------
c          
c  WORKD   Complex*16 work array of length 3*N.  (REVERSE COMMUNICATION)
c          Distributed array to be used in the basic Arnoldi iteration
c          for reverse communication.  The user should not use WORKD 
c          as temporary workspace during the iteration !!!!!!!!!!
c          See Data Distribution Note below.  
c
c  WORKL   Complex*16 work array of length LWORKL.  (OUTPUT/WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.  See Data Distribution Note below.
c
c  LWORKL  Integer.  (INPUT)
c          LWORKL must be at least 3*NCV**2 + 5*NCV.
c
c  RWORK   Double precision  work array of length NCV (WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.
c
c
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
c          = -3: NCV-NEV >= 1 and less than or equal to N.
c          = -4: The maximum number of Arnoldi update iteration 
c                must be greater than zero.
c          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
c          = -6: BMAT must be one of 'I' or 'G'.
c          = -7: Length of private work array is not sufficient.
c          = -8: Error return from LAPACK eigenvalue calculation;
c          = -9: Starting vector is zero.
c          = -10: IPARAM(7) must be 1,2,3.
c          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
c          = -12: IPARAM(1) must be equal to 0 or 1.
c          = -9999: Could not build an Arnoldi factorization.
c                   User input error highly likely.  Please
c                   check actual array dimensions and layout.
c                   IPARAM(5) returns the size of the current Arnoldi
c                   factorization.
c
c\Remarks
c  1. The computed Ritz values are approximate eigenvalues of OP. The
c     selection of WHICH should be made with this in mind when using
c     Mode = 3.  When operating in Mode = 3 setting WHICH = 'LM' will
c     compute the NEV eigenvalues of the original problem that are
c     closest to the shift SIGMA . After convergence, approximate eigenvalues 
c     of the original problem may be obtained with the ARPACK subroutine zneupd.
c
c  2. If a basis for the invariant subspace corresponding to the converged Ritz 
c     values is needed, the user must call zneupd immediately following 
c     completion of znaupd. This is new starting with release 2 of ARPACK.
c
c  3. If M can be factored into a Cholesky factorization M = LL`
c     then Mode = 2 should not be selected.  Instead one should use
c     Mode = 1 with  OP = inv(L)*A*inv(L`).  Appropriate triangular 
c     linear systems should be solved with L and L` rather
c     than computing inverses.  After convergence, an approximate
c     eigenvector z of the original problem is recovered by solving
c     L`z = x  where x is a Ritz vector of OP.
c
c  4. At present there is no a-priori analysis to guide the selection
c     of NCV relative to NEV.  The only formal requirement is that NCV > NEV + 1.
c     However, it is recommended that NCV .ge. 2*NEV.  If many problems of
c     the same type are to be solved, one should experiment with increasing
c     NCV while keeping NEV fixed for a given test problem.  This will 
c     usually decrease the required number of OP*x operations but it
c     also increases the work and storage required to maintain the orthogonal
c     basis vectors.  The optimal "cross-over" with respect to CPU time
c     is problem dependent and must be determined empirically. 
c     See Chapter 8 of Reference 2 for further information.
c
c  5. When IPARAM(1) = 0, and IDO = 3, the user needs to provide the
c     NP = IPARAM(8) complex shifts in locations
c     WORKL(IPNTR(14)), WORKL(IPNTR(14)+1), ... , WORKL(IPNTR(14)+NP).
c     Eigenvalues of the current upper Hessenberg matrix are located in
c     WORKL(IPNTR(6)) through WORKL(IPNTR(6)+NCV-1). They are ordered
c     according to the order defined by WHICH.  The associated Ritz estimates
c     are located in WORKL(IPNTR(8)), WORKL(IPNTR(8)+1), ... ,
c     WORKL(IPNTR(8)+NCV-1).
c
c-----------------------------------------------------------------------
c
c\Data Distribution Note: 
c
c  Fortran-D syntax:
c  ================
c  Complex*16 resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
c  decompose  d1(n), d2(n,ncv)
c  align      resid(i) with d1(i)
c  align      v(i,j)   with d2(i,j)
c  align      workd(i) with d1(i)     range (1:n)
c  align      workd(i) with d1(i-n)   range (n+1:2*n)
c  align      workd(i) with d1(i-2*n) range (2*n+1:3*n)
c  distribute d1(block), d2(block,:)
c  replicated workl(lworkl)
c
c  Cray MPP syntax:
c  ===============
c  Complex*16 resid(n), v(ldv,ncv), workd(n,3), workl(lworkl)
c  shared     resid(block), v(block,:), workd(block,:)
c  replicated workl(lworkl)
c  
c  CM2/CM5 syntax:
c  ==============
c  
c-----------------------------------------------------------------------
c
c     include   'ex-nonsym.doc'
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  Complex*16
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
c     Restarted Arnoldi Iteration", Rice University Technical Report
c     TR95-13, Department of Computational and Applied Mathematics.
c  3. B.N. Parlett & Y. Saad, "_Complex_ Shift and Invert Strategies for
c     Double precision Matrices", Linear Algebra and its Applications, vol 88/89,
c     pp 575-595, (1987).
c
c\Routines called:
c     znaup2  ARPACK routine that implements the Implicitly Restarted
c             Arnoldi Iteration.
c     zstatn  ARPACK routine that initializes the timing variables.
c     ivout   ARPACK utility routine that prints integers.
c     zvout   ARPACK utility routine that prints vectors.
c     second  ARPACK utility routine for timing.
c     dlamch  LAPACK routine that determines machine constants.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics 
c     Rice University           
c     Houston, Texas 
c 
c\SCCS Information: @(#)
c FILE: naupd.F   SID: 2.9   DATE OF SID: 07/21/02   RELEASE: 2
c
c\Remarks
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine znaupd
     &   ( ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, 
     &     ipntr, workd, workl, lworkl, rwork, info )
c
c     %----------------------------------------------------%
c     | Include files for debugging and timing information |
c     %----------------------------------------------------%
c
      include   'debug.h'
      include   'stat.h'
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character  bmat*1, which*2
      integer    ido, info, ldv, lworkl, n, ncv, nev
      Double precision 
     &           tol
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    iparam(11), ipntr(14)
      Complex*16
     &           resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
      Double precision  
     &           rwork(ncv)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Complex*16
     &           one, zero
      parameter (one = (1.0D+0, 0.0D+0), zero = (0.0D+0, 0.0D+0))
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer    bounds, ierr, ih, iq, ishift, iupd, iw, 
     &           ldh, ldq, levec, mode, msglvl, mxiter, nb,
     &           nev0, next, np, ritz, j
      save       bounds, ih, iq, ishift, iupd, iw,
     &           ldh, ldq, levec, mode, msglvl, mxiter, nb,
     &           nev0, next, np, ritz
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   znaup2, zvout, ivout, second, zstatn
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision 
     &           dlamch
      external   dlamch
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c 
      if (ido .eq. 0) then
c 
c        %-------------------------------%
c        | Initialize timing statistics  |
c        | & message level for debugging |
c        %-------------------------------%
c
         call zstatn
         call second (t0)
         msglvl = mcaupd
c
c        %----------------%
c        | Error checking |
c        %----------------%
c
         ierr   = 0
         ishift = iparam(1)
c         levec  = iparam(2)
         mxiter = iparam(3)
c         nb     = iparam(4)
         nb     = 1
c
c        %--------------------------------------------%
c        | Revision 2 performs only implicit restart. |
c        %--------------------------------------------%
c
         iupd   = 1
         mode   = iparam(7)
c
         if (n .le. 0) then
             ierr = -1
         else if (nev .le. 0) then
             ierr = -2
         else if (ncv .le. nev .or.  ncv .gt. n) then
             ierr = -3
         else if (mxiter .le. 0) then
             ierr = -4
         else if (which .ne. 'LM' .and.
     &       which .ne. 'SM' .and.
     &       which .ne. 'LR' .and.
     &       which .ne. 'SR' .and.
     &       which .ne. 'LI' .and.
     &       which .ne. 'SI') then
            ierr = -5
         else if (bmat .ne. 'I' .and. bmat .ne. 'G') then
            ierr = -6
         else if (lworkl .lt. 3*ncv**2 + 5*ncv) then
            ierr = -7
         else if (mode .lt. 1 .or. mode .gt. 3) then
                                                ierr = -10
         else if (mode .eq. 1 .and. bmat .eq. 'G') then
                                                ierr = -11
         end if
c 
c        %------------%
c        | Error Exit |
c        %------------%
c
         if (ierr .ne. 0) then
            info = ierr
            ido  = 99
            go to 9000
         end if
c 
c        %------------------------%
c        | Set default parameters |
c        %------------------------%
c
         if (nb .le. 0)				nb = 1
         if (tol .le. 0.0D+0 )			tol = dlamch('EpsMach')
         if (ishift .ne. 0  .and.  
     &       ishift .ne. 1  .and.
     &       ishift .ne. 2) 			ishift = 1
c
c        %----------------------------------------------%
c        | NP is the number of additional steps to      |
c        | extend the length NEV Lanczos factorization. |
c        | NEV0 is the local variable designating the   |
c        | size of the invariant subspace desired.      |
c        %----------------------------------------------%
c
         np     = ncv - nev
         nev0   = nev 
c 
c        %-----------------------------%
c        | Zero out internal workspace |
c        %-----------------------------%
c
         do 10 j = 1, 3*ncv**2 + 5*ncv
            workl(j) = zero
  10     continue
c 
c        %-------------------------------------------------------------%
c        | Pointer into WORKL for address of H, RITZ, BOUNDS, Q        |
c        | etc... and the remaining workspace.                         |
c        | Also update pointer to be used on output.                   |
c        | Memory is laid out as follows:                              |
c        | workl(1:ncv*ncv) := generated Hessenberg matrix             |
c        | workl(ncv*ncv+1:ncv*ncv+ncv) := the ritz values             |
c        | workl(ncv*ncv+ncv+1:ncv*ncv+2*ncv)   := error bounds        |
c        | workl(ncv*ncv+2*ncv+1:2*ncv*ncv+2*ncv) := rotation matrix Q |
c        | workl(2*ncv*ncv+2*ncv+1:3*ncv*ncv+5*ncv) := workspace       |
c        | The final workspace is needed by subroutine zneigh called   |
c        | by znaup2. Subroutine zneigh calls LAPACK routines for      |
c        | calculating eigenvalues and the last row of the eigenvector |
c        | matrix.                                                     |
c        %-------------------------------------------------------------%
c
         ldh    = ncv
         ldq    = ncv
         ih     = 1
         ritz   = ih     + ldh*ncv
         bounds = ritz   + ncv
         iq     = bounds + ncv
         iw     = iq     + ldq*ncv
         next   = iw     + ncv**2 + 3*ncv
c
         ipntr(4) = next
         ipntr(5) = ih
         ipntr(6) = ritz
         ipntr(7) = iq
         ipntr(8) = bounds
         ipntr(14) = iw
      end if
c
c     %-------------------------------------------------------%
c     | Carry out the Implicitly restarted Arnoldi Iteration. |
c     %-------------------------------------------------------%
c
      call znaup2 
     &   ( ido, bmat, n, which, nev0, np, tol, resid, mode, iupd,
     &     ishift, mxiter, v, ldv, workl(ih), ldh, workl(ritz), 
     &     workl(bounds), workl(iq), ldq, workl(iw), 
     &     ipntr, workd, rwork, info )
c 
c     %--------------------------------------------------%
c     | ido .ne. 99 implies use of reverse communication |
c     | to compute operations involving OP.              |
c     %--------------------------------------------------%
c
      if (ido .eq. 3) iparam(8) = np
      if (ido .ne. 99) go to 9000
c 
      iparam(3) = mxiter
      iparam(5) = np
      iparam(9) = nopx
      iparam(10) = nbx
      iparam(11) = nrorth
c
c     %------------------------------------%
c     | Exit if there was an informational |
c     | error within znaup2.               |
c     %------------------------------------%
c
      if (info .lt. 0) go to 9000
      if (info .eq. 2) info = 3
c
      if (msglvl .gt. 0) then
         call ivout (logfil, 1, mxiter, ndigit,
     &               '_naupd: Number of update iterations taken')
         call ivout (logfil, 1, np, ndigit,
     &               '_naupd: Number of wanted "converged" Ritz values')
         call zvout (logfil, np, workl(ritz), ndigit, 
     &               '_naupd: The final Ritz values')
         call zvout (logfil, np, workl(bounds), ndigit, 
     &               '_naupd: Associated Ritz estimates')
      end if
c
      call second (t1)
      tcaupd = t1 - t0
c
      if (msglvl .gt. 0) then
c
c        %--------------------------------------------------------%
c        | Version Number & Version Date are defined in version.h |
c        %--------------------------------------------------------%
c
         write (6,1000)
         write (6,1100) mxiter, nopx, nbx, nrorth, nitref, nrstrt,
     &                  tmvopx, tmvbx, tcaupd, tcaup2, tcaitr, titref,
     &                  tgetv0, tceigh, tcgets, tcapps, tcconv, trvec
 1000    format (//,
     &      5x, '=============================================',/
     &      5x, '= Complex implicit Arnoldi update code      =',/
     &      5x, '= Version Number: ', ' 2.3', 21x, ' =',/
     &      5x, '= Version Date:   ', ' 07/31/96', 16x,   ' =',/
     &      5x, '=============================================',/
     &      5x, '= Summary of timing statistics              =',/
     &      5x, '=============================================',//)
 1100    format (
     &      5x, 'Total number update iterations             = ', i5,/
     &      5x, 'Total number of OP*x operations            = ', i5,/
     &      5x, 'Total number of B*x operations             = ', i5,/
     &      5x, 'Total number of reorthogonalization steps  = ', i5,/
     &      5x, 'Total number of iterative refinement steps = ', i5,/
     &      5x, 'Total number of restart steps              = ', i5,/
     &      5x, 'Total time in user OP*x operation          = ', f12.6,/
     &      5x, 'Total time in user B*x operation           = ', f12.6,/
     &      5x, 'Total time in Arnoldi update routine       = ', f12.6,/
     &      5x, 'Total time in naup2 routine                = ', f12.6,/
     &      5x, 'Total time in basic Arnoldi iteration loop = ', f12.6,/
     &      5x, 'Total time in reorthogonalization phase    = ', f12.6,/
     &      5x, 'Total time in (re)start vector generation  = ', f12.6,/
     &      5x, 'Total time in Hessenberg eig. subproblem   = ', f12.6,/
     &      5x, 'Total time in getting the shifts           = ', f12.6,/
     &      5x, 'Total time in applying the shifts          = ', f12.6,/
     &      5x, 'Total time in convergence testing          = ', f12.6,/
     &      5x, 'Total time in computing final Ritz vectors = ', f12.6/)
      end if
c
 9000 continue
c
      return
c
c     %---------------%
c     | End of znaupd |
c     %---------------%
c
      end



c\BeginDoc
c 
c\Name: zneupd 
c 
c\Description: 
c  This subroutine returns the converged approximations to eigenvalues 
c  of A*z = lambda*B*z and (optionally): 
c 
c      (1) The corresponding approximate eigenvectors; 
c 
c      (2) An orthonormal basis for the associated approximate 
c          invariant subspace; 
c 
c      (3) Both.  
c
c  There is negligible additional cost to obtain eigenvectors.  An orthonormal 
c  basis is always computed.  There is an additional storage cost of n*nev
c  if both are requested (in this case a separate array Z must be supplied). 
c
c  The approximate eigenvalues and eigenvectors of  A*z = lambda*B*z
c  are derived from approximate eigenvalues and eigenvectors of
c  of the linear operator OP prescribed by the MODE selection in the
c  call to ZNAUPD.  ZNAUPD must be called before this routine is called.
c  These approximate eigenvalues and vectors are commonly called Ritz
c  values and Ritz vectors respectively.  They are referred to as such 
c  in the comments that follow.   The computed orthonormal basis for the 
c  invariant subspace corresponding to these Ritz values is referred to as a 
c  Schur basis. 
c 
c  The definition of OP as well as other terms and the relation of computed
c  Ritz values and vectors of OP with respect to the given problem
c  A*z = lambda*B*z may be found in the header of ZNAUPD.  For a brief 
c  description, see definitions of IPARAM(7), MODE and WHICH in the
c  documentation of ZNAUPD.
c
c\Usage:
c  call zneupd 
c     ( RVEC, HOWMNY, SELECT, D, Z, LDZ, SIGMA, WORKEV, BMAT, 
c       N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR, WORKD, 
c       WORKL, LWORKL, RWORK, INFO )
c
c\Arguments:
c  RVEC    LOGICAL  (INPUT)
c          Specifies whether a basis for the invariant subspace corresponding
c          to the converged Ritz value approximations for the eigenproblem 
c          A*z = lambda*B*z is computed.
c
c             RVEC = .FALSE.     Compute Ritz values only.
c
c             RVEC = .TRUE.      Compute Ritz vectors or Schur vectors.
c                                See Remarks below.
c
c  HOWMNY  Character*1  (INPUT)
c          Specifies the form of the basis for the invariant subspace 
c          corresponding to the converged Ritz values that is to be computed.
c
c          = 'A': Compute NEV Ritz vectors;
c          = 'P': Compute NEV Schur vectors;
c          = 'S': compute some of the Ritz vectors, specified
c                 by the logical array SELECT.
c
c  SELECT  Logical array of dimension NCV.  (INPUT)
c          If HOWMNY = 'S', SELECT specifies the Ritz vectors to be
c          computed. To select the  Ritz vector corresponding to a
c          Ritz value D(j), SELECT(j) must be set to .TRUE.. 
c          If HOWMNY = 'A' or 'P', SELECT need not be initialized 
c          but it is used as internal workspace.
c
c  D       Complex*16 array of dimension NEV+1.  (OUTPUT)
c          On exit, D contains the  Ritz  approximations 
c          to the eigenvalues lambda for A*z = lambda*B*z.
c
c  Z       Complex*16 N by NEV array    (OUTPUT)
c          On exit, if RVEC = .TRUE. and HOWMNY = 'A', then the columns of 
c          Z represents approximate eigenvectors (Ritz vectors) corresponding 
c          to the NCONV=IPARAM(5) Ritz values for eigensystem
c          A*z = lambda*B*z.
c
c          If RVEC = .FALSE. or HOWMNY = 'P', then Z is NOT REFERENCED.
c
c          NOTE: If if RVEC = .TRUE. and a Schur basis is not required, 
c          the array Z may be set equal to first NEV+1 columns of the Arnoldi 
c          basis array V computed by ZNAUPD.  In this case the Arnoldi basis 
c          will be destroyed and overwritten with the eigenvector basis.
c
c  LDZ     Integer.  (INPUT)
c          The leading dimension of the array Z.  If Ritz vectors are
c          desired, then  LDZ .ge.  max( 1, N ) is required.  
c          In any case,  LDZ .ge. 1 is required.
c
c  SIGMA   Complex*16  (INPUT)
c          If IPARAM(7) = 3 then SIGMA represents the shift. 
c          Not referenced if IPARAM(7) = 1 or 2.
c
c  WORKEV  Complex*16 work array of dimension 2*NCV.  (WORKSPACE)
c
c  **** The remaining arguments MUST be the same as for the   ****
c  **** call to ZNAUPD that was just completed.               ****
c
c  NOTE: The remaining arguments 
c
c           BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR, 
c           WORKD, WORKL, LWORKL, RWORK, INFO 
c
c         must be passed directly to ZNEUPD following the last call 
c         to ZNAUPD.  These arguments MUST NOT BE MODIFIED between
c         the the last call to ZNAUPD and the call to ZNEUPD.
c
c  Three of these parameters (V, WORKL and INFO) are also output parameters:
c
c  V       Complex*16 N by NCV array.  (INPUT/OUTPUT)
c
c          Upon INPUT: the NCV columns of V contain the Arnoldi basis
c                      vectors for OP as constructed by ZNAUPD .
c
c          Upon OUTPUT: If RVEC = .TRUE. the first NCONV=IPARAM(5) columns
c                       contain approximate Schur vectors that span the
c                       desired invariant subspace.
c
c          NOTE: If the array Z has been set equal to first NEV+1 columns
c          of the array V and RVEC=.TRUE. and HOWMNY= 'A', then the
c          Arnoldi basis held by V has been overwritten by the desired
c          Ritz vectors.  If a separate array Z has been passed then
c          the first NCONV=IPARAM(5) columns of V will contain approximate
c          Schur vectors that span the desired invariant subspace.
c
c  WORKL   Double precision work array of length LWORKL.  (OUTPUT/WORKSPACE)
c          WORKL(1:ncv*ncv+2*ncv) contains information obtained in
c          znaupd.  They are not changed by zneupd.
c          WORKL(ncv*ncv+2*ncv+1:3*ncv*ncv+4*ncv) holds the
c          untransformed Ritz values, the untransformed error estimates of 
c          the Ritz values, the upper triangular matrix for H, and the
c          associated matrix representation of the invariant subspace for H.
c
c          Note: IPNTR(9:13) contains the pointer into WORKL for addresses
c          of the above information computed by zneupd.
c          -------------------------------------------------------------
c          IPNTR(9):  pointer to the NCV RITZ values of the
c                     original system.
c          IPNTR(10): Not used
c          IPNTR(11): pointer to the NCV corresponding error estimates.
c          IPNTR(12): pointer to the NCV by NCV upper triangular
c                     Schur matrix for H.
c          IPNTR(13): pointer to the NCV by NCV matrix of eigenvectors
c                     of the upper Hessenberg matrix H. Only referenced by
c                     zneupd if RVEC = .TRUE. See Remark 2 below.
c          -------------------------------------------------------------
c
c  INFO    Integer.  (OUTPUT)
c          Error flag on output.
c          =  0: Normal exit.
c
c          =  1: The Schur form computed by LAPACK routine csheqr
c                could not be reordered by LAPACK routine ztrsen.
c                Re-enter subroutine zneupd with IPARAM(5)=NCV and
c                increase the size of the array D to have
c                dimension at least dimension NCV and allocate at least NCV
c                columns for Z. NOTE: Not necessary if Z and V share
c                the same space. Please notify the authors if this error
c                occurs.
c
c          = -1: N must be positive.
c          = -2: NEV must be positive.
c          = -3: NCV-NEV >= 1 and less than or equal to N.
c          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
c          = -6: BMAT must be one of 'I' or 'G'.
c          = -7: Length of private work WORKL array is not sufficient.
c          = -8: Error return from LAPACK eigenvalue calculation.
c                This should never happened.
c          = -9: Error return from calculation of eigenvectors.
c                Informational error from LAPACK routine ztrevc.
c          = -10: IPARAM(7) must be 1,2,3
c          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
c          = -12: HOWMNY = 'S' not yet implemented
c          = -13: HOWMNY must be one of 'A' or 'P' if RVEC = .true.
c          = -14: ZNAUPD did not find any eigenvalues to sufficient
c                 accuracy.
c          = -15: ZNEUPD got a different count of the number of converged
c                 Ritz values than ZNAUPD got.  This indicates the user
c                 probably made an error in passing data from ZNAUPD to
c                 ZNEUPD or that the data was modified before entering
c                 ZNEUPD
c
c\BeginLib
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
c     Restarted Arnoldi Iteration", Rice University Technical Report
c     TR95-13, Department of Computational and Applied Mathematics.
c  3. B. Nour-Omid, B. N. Parlett, T. Ericsson and P. S. Jensen,
c     "How to Implement the Spectral Transformation", Math Comp.,
c     Vol. 48, No. 178, April, 1987 pp. 664-673. 
c
c\Routines called:
c     ivout   ARPACK utility routine that prints integers.
c     zmout   ARPACK utility routine that prints matrices
c     zvout   ARPACK utility routine that prints vectors.
c     zgeqr2  LAPACK routine that computes the QR factorization of 
c             a matrix.
c     zlacpy  LAPACK matrix copy routine.
c     zlahqr  LAPACK routine that computes the Schur form of a
c             upper Hessenberg matrix.
c     zlaset  LAPACK matrix initialization routine.
c     ztrevc  LAPACK routine to compute the eigenvectors of a matrix
c             in upper triangular form.
c     ztrsen  LAPACK routine that re-orders the Schur form.
c     zunm2r  LAPACK routine that applies an orthogonal matrix in 
c             factored form.
c     dlamch  LAPACK routine that determines machine constants.
c     ztrmm   Level 3 BLAS matrix times an upper triangular matrix.
c     zgeru   Level 2 BLAS rank one update to a matrix.
c     zcopy   Level 1 BLAS that copies one vector to another .
c     zscal   Level 1 BLAS that scales a vector.
c     zdscal  Level 1 BLAS that scales a complex vector by a real number.
c     dznrm2  Level 1 BLAS that computes the norm of a complex vector.
c
c\Remarks
c
c  1. Currently only HOWMNY = 'A' and 'P' are implemented. 
c
c  2. Schur vectors are an orthogonal representation for the basis of
c     Ritz vectors. Thus, their numerical properties are often superior.
c     If RVEC = .true. then the relationship
c             A * V(:,1:IPARAM(5)) = V(:,1:IPARAM(5)) * T, and
c       transpose( V(:,1:IPARAM(5)) ) * V(:,1:IPARAM(5)) = I
c     are approximately satisfied.
c     Here T is the leading submatrix of order IPARAM(5) of the 
c     upper triangular matrix stored workl(ipntr(12)). 
c
c\Authors
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Chao Yang                    Houston, Texas 
c     Dept. of Computational & 
c     Applied Mathematics 
c     Rice University 
c     Houston, Texas
c
c\SCCS Information: @(#)
c FILE: neupd.F   SID: 2.8   DATE OF SID: 07/21/02   RELEASE: 2
c
c\EndLib
c
c-----------------------------------------------------------------------
      subroutine zneupd(rvec , howmny, select, d     ,
     &                   z    , ldz   , sigma , workev,
     &                   bmat , n     , which , nev   ,
     &                   tol  , resid , ncv   , v     ,
     &                   ldv  , iparam, ipntr , workd ,
     &                   workl, lworkl, rwork , info  )
c
c     %----------------------------------------------------%
c     | Include files for debugging and timing information |
c     %----------------------------------------------------%
c
      include   'debug.h'
      include   'stat.h'
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character  bmat, howmny, which*2
      logical    rvec
      integer    info, ldz, ldv, lworkl, n, ncv, nev
      Complex*16     
     &           sigma
      Double precision 
     &           tol
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    iparam(11), ipntr(14)
      logical    select(ncv)
      Double precision
     &           rwork(ncv)
      Complex*16
     &           d(nev)     , resid(n)     , v(ldv,ncv),
     &           z(ldz, nev), 
     &           workd(3*n) , workl(lworkl), workev(2*ncv)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Complex*16
     &           one, zero
      parameter  (one = (1.0D+0, 0.0D+0), zero = (0.0D+0, 0.0D+0))
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character  type*6
      integer    bounds, ierr  , ih    , ihbds, iheig , nconv ,
     &           invsub, iuptri, iwev  , j    , ldh   , ldq   ,
     &           mode  , msglvl, ritz  , wr   , k     , irz   ,
     &           ibd   , outncv, iq    , np   , numcnv, jj    ,
     &           ishift
      Complex*16
     &           rnorm, temp, vl(1)
      Double precision
     &           conds, sep, rtemp, eps23
      logical    reord
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   zcopy , zgeru, zgeqr2, zlacpy, zmout,
     &           zunm2r, ztrmm, zvout, ivout,
     &           zlahqr
c  
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision
     &           dznrm2, dlamch, dlapy2
      external   dznrm2, dlamch, dlapy2
c
      Complex*16
     &           zdotc
      external   zdotc
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c 
c     %------------------------%
c     | Set default parameters |
c     %------------------------%
c
      msglvl = mceupd
      mode = iparam(7)
      nconv = iparam(5)
      info = 0
c
c
c     %---------------------------------%
c     | Get machine dependent constant. |
c     %---------------------------------%
c
      eps23 = dlamch('Epsilon-Machine')
      eps23 = eps23**(2.0D+0 / 3.0D+0)
c
c     %-------------------------------%
c     | Quick return                  |
c     | Check for incompatible input  |
c     %-------------------------------%
c
      ierr = 0
c
      if (nconv .le. 0) then
         ierr = -14
      else if (n .le. 0) then
         ierr = -1
      else if (nev .le. 0) then
         ierr = -2
      else if (ncv .le. nev .or.  ncv .gt. n) then
         ierr = -3
      else if (which .ne. 'LM' .and.
     &        which .ne. 'SM' .and.
     &        which .ne. 'LR' .and.
     &        which .ne. 'SR' .and.
     &        which .ne. 'LI' .and.
     &        which .ne. 'SI') then
         ierr = -5
      else if (bmat .ne. 'I' .and. bmat .ne. 'G') then
         ierr = -6
      else if (lworkl .lt. 3*ncv**2 + 4*ncv) then
         ierr = -7
      else if ( (howmny .ne. 'A' .and.
     &           howmny .ne. 'P' .and.
     &           howmny .ne. 'S') .and. rvec ) then
         ierr = -13
      else if (howmny .eq. 'S' ) then
         ierr = -12
      end if
c     
      if (mode .eq. 1 .or. mode .eq. 2) then
         type = 'REGULR'
      else if (mode .eq. 3 ) then
         type = 'SHIFTI'
      else 
                                              ierr = -10
      end if
      if (mode .eq. 1 .and. bmat .eq. 'G')    ierr = -11
c
c     %------------%
c     | Error Exit |
c     %------------%
c
      if (ierr .ne. 0) then
         info = ierr
         go to 9000
      end if
c 
c     %--------------------------------------------------------%
c     | Pointer into WORKL for address of H, RITZ, WORKEV, Q   |
c     | etc... and the remaining workspace.                    |
c     | Also update pointer to be used on output.              |
c     | Memory is laid out as follows:                         |
c     | workl(1:ncv*ncv) := generated Hessenberg matrix        |
c     | workl(ncv*ncv+1:ncv*ncv+ncv) := ritz values            |
c     | workl(ncv*ncv+ncv+1:ncv*ncv+2*ncv) := error bounds     |
c     %--------------------------------------------------------%
c
c     %-----------------------------------------------------------%
c     | The following is used and set by ZNEUPD.                 |
c     | workl(ncv*ncv+2*ncv+1:ncv*ncv+3*ncv) := The untransformed |
c     |                                      Ritz values.         |
c     | workl(ncv*ncv+3*ncv+1:ncv*ncv+4*ncv) := The untransformed |
c     |                                      error bounds of      |
c     |                                      the Ritz values      |
c     | workl(ncv*ncv+4*ncv+1:2*ncv*ncv+4*ncv) := Holds the upper |
c     |                                      triangular matrix    |
c     |                                      for H.               |
c     | workl(2*ncv*ncv+4*ncv+1: 3*ncv*ncv+4*ncv) := Holds the    |
c     |                                      associated matrix    |
c     |                                      representation of    |
c     |                                      the invariant        |
c     |                                      subspace for H.      |
c     | GRAND total of NCV * ( 3 * NCV + 4 ) locations.           |
c     %-----------------------------------------------------------%
c     
      ih     = ipntr(5)
      ritz   = ipntr(6)
      iq     = ipntr(7)
      bounds = ipntr(8)
      ldh    = ncv
      ldq    = ncv
      iheig  = bounds + ldh
      ihbds  = iheig  + ldh
      iuptri = ihbds  + ldh
      invsub = iuptri + ldh*ncv
      ipntr(9)  = iheig
      ipntr(11) = ihbds
      ipntr(12) = iuptri
      ipntr(13) = invsub
      wr = 1
      iwev = wr + ncv
c
c     %-----------------------------------------%
c     | irz points to the Ritz values computed  |
c     |     by _neigh before exiting _naup2.    |
c     | ibd points to the Ritz estimates        |
c     |     computed by _neigh before exiting   |
c     |     _naup2.                             |
c     %-----------------------------------------%
c
      irz = ipntr(14) + ncv*ncv
      ibd = irz + ncv
c
c     %------------------------------------%
c     | RNORM is B-norm of the RESID(1:N). |
c     %------------------------------------%
c
      rnorm = workl(ih+2)
      workl(ih+2) = zero
c
      if (msglvl .gt. 2) then
         call zvout(logfil, ncv, workl(irz), ndigit,
     &   '_neupd: Ritz values passed in from _NAUPD.')
         call zvout(logfil, ncv, workl(ibd), ndigit,
     &   '_neupd: Ritz estimates passed in from _NAUPD.')
      end if
c
      if (rvec) then
c
         reord = .false.
c
c        %---------------------------------------------------%
c        | Use the temporary bounds array to store indices   |
c        | These will be used to mark the select array later |
c        %---------------------------------------------------%
c
         do 10 j = 1,ncv
            workl(bounds+j-1) = j
            select(j) = .false.
   10    continue
c
c        %-------------------------------------%
c        | Select the wanted Ritz values.      |
c        | Sort the Ritz values so that the    |
c        | wanted ones appear at the tailing   |
c        | NEV positions of workl(irr) and     |
c        | workl(iri).  Move the corresponding |
c        | error estimates in workl(ibd)       |
c        | accordingly.                        |
c        %-------------------------------------%
c
         np     = ncv - nev
         ishift = 0
         call zngets(ishift, which     , nev          ,
     &                np    , workl(irz), workl(bounds))
c
         if (msglvl .gt. 2) then
            call zvout (logfil, ncv, workl(irz), ndigit,
     &      '_neupd: Ritz values after calling _NGETS.')
            call zvout (logfil, ncv, workl(bounds), ndigit,
     &      '_neupd: Ritz value indices after calling _NGETS.')
         end if
c
c        %-----------------------------------------------------%
c        | Record indices of the converged wanted Ritz values  |
c        | Mark the select array for possible reordering       |
c        %-----------------------------------------------------%
c
         numcnv = 0
         do 11 j = 1,ncv
            rtemp = max(eps23,
     &                 dlapy2 ( dble(workl(irz+ncv-j)),
     &                          dimag(workl(irz+ncv-j)) ))
            jj = workl(bounds + ncv - j)
            if (numcnv .lt. nconv .and.
     &          dlapy2( dble(workl(ibd+jj-1)),
     &          dimag(workl(ibd+jj-1)) )
     &          .le. tol*rtemp) then
               select(jj) = .true.
               numcnv = numcnv + 1
               if (jj .gt. nev) reord = .true.
            endif
   11    continue
c
c        %-----------------------------------------------------------%
c        | Check the count (numcnv) of converged Ritz values with    |
c        | the number (nconv) reported by dnaupd.  If these two      |
c        | are different then there has probably been an error       |
c        | caused by incorrect passing of the dnaupd data.           |
c        %-----------------------------------------------------------%
c
         if (msglvl .gt. 2) then
             call ivout(logfil, 1, numcnv, ndigit,
     &            '_neupd: Number of specified eigenvalues')
             call ivout(logfil, 1, nconv, ndigit,
     &            '_neupd: Number of "converged" eigenvalues')
         end if
c
         if (numcnv .ne. nconv) then
            info = -15
            go to 9000
         end if
c
c        %-------------------------------------------------------%
c        | Call LAPACK routine zlahqr to compute the Schur form |
c        | of the upper Hessenberg matrix returned by ZNAUPD.   |
c        | Make a copy of the upper Hessenberg matrix.           |
c        | Initialize the Schur vector matrix Q to the identity. |
c        %-------------------------------------------------------%
c
         call zcopy(ldh*ncv, workl(ih), 1, workl(iuptri), 1)
         call zlaset('All', ncv, ncv          , 
     &                zero , one, workl(invsub),
     &                ldq)
         call zlahqr(.true., .true.       , ncv          , 
     &                1     , ncv          , workl(iuptri),
     &                ldh   , workl(iheig) , 1            ,
     &                ncv   , workl(invsub), ldq          ,
     &                ierr)
         call zcopy(ncv         , workl(invsub+ncv-1), ldq,
     &               workl(ihbds), 1)
c
         if (ierr .ne. 0) then
            info = -8
            go to 9000
         end if
c
         if (msglvl .gt. 1) then
            call zvout (logfil, ncv, workl(iheig), ndigit,
     &           '_neupd: Eigenvalues of H')
            call zvout (logfil, ncv, workl(ihbds), ndigit,
     &           '_neupd: Last row of the Schur vector matrix')
            if (msglvl .gt. 3) then
               call zmout (logfil       , ncv, ncv   , 
     &                     workl(iuptri), ldh, ndigit,
     &              '_neupd: The upper triangular matrix ')
            end if
         end if
c
         if (reord) then
c
c           %-----------------------------------------------%
c           | Reorder the computed upper triangular matrix. |
c           %-----------------------------------------------%
c
            call ztrsen('None'       , 'V'          , select      ,
     &                   ncv          , workl(iuptri), ldh         ,
     &                   workl(invsub), ldq          , workl(iheig),
     &                   nconv        , conds        , sep         , 
     &                   workev       , ncv          , ierr)
c
            if (ierr .eq. 1) then
               info = 1
               go to 9000
            end if
c
            if (msglvl .gt. 2) then
                call zvout (logfil, ncv, workl(iheig), ndigit,
     &           '_neupd: Eigenvalues of H--reordered')
                if (msglvl .gt. 3) then
                   call zmout(logfil       , ncv, ncv   ,
     &                         workl(iuptri), ldq, ndigit,
     &              '_neupd: Triangular matrix after re-ordering')
                end if
            end if
c
         end if
c
c        %---------------------------------------------%
c        | Copy the last row of the Schur basis matrix |
c        | to workl(ihbds).  This vector will be used  |
c        | to compute the Ritz estimates of converged  |
c        | Ritz values.                                |
c        %---------------------------------------------%
c
         call zcopy(ncv         , workl(invsub+ncv-1), ldq,
     &               workl(ihbds), 1)
c 
c        %--------------------------------------------%
c        | Place the computed eigenvalues of H into D |
c        | if a spectral transformation was not used. |
c        %--------------------------------------------%
c
         if (type .eq. 'REGULR') then
            call zcopy(nconv, workl(iheig), 1, d, 1)
         end if
c
c        %----------------------------------------------------------%
c        | Compute the QR factorization of the matrix representing  |
c        | the wanted invariant subspace located in the first NCONV |
c        | columns of workl(invsub,ldq).                            |
c        %----------------------------------------------------------%
c
         call zgeqr2(ncv , nconv , workl(invsub),
     &                ldq , workev, workev(ncv+1),
     &                ierr)
c
c        %--------------------------------------------------------%
c        | * Postmultiply V by Q using zunm2r.                    |
c        | * Copy the first NCONV columns of VQ into Z.           |
c        | * Postmultiply Z by R.                                 |
c        | The N by NCONV matrix Z is now a matrix representation |
c        | of the approximate invariant subspace associated with  |
c        | the Ritz values in workl(iheig). The first NCONV       | 
c        | columns of V are now approximate Schur vectors         |
c        | associated with the upper triangular matrix of order   |
c        | NCONV in workl(iuptri).                                |
c        %--------------------------------------------------------%
c
         call zunm2r('Right', 'Notranspose', n            ,
     &                ncv    , nconv        , workl(invsub),
     &                ldq    , workev       , v            ,
     &                ldv    , workd(n+1)   , ierr)
         call zlacpy('All', n, nconv, v, ldv, z, ldz)
c
         do 20 j=1, nconv
c
c           %---------------------------------------------------%
c           | Perform both a column and row scaling if the      |
c           | diagonal element of workl(invsub,ldq) is negative |
c           | I'm lazy and don't take advantage of the upper    |
c           | triangular form of workl(iuptri,ldq).             |
c           | Note that since Q is orthogonal, R is a diagonal  |
c           | matrix consisting of plus or minus ones.          |
c           %---------------------------------------------------%
c
            if ( dble( workl(invsub+(j-1)*ldq+j-1) ) .lt. 
     &                  dble(zero) ) then
               call zscal(nconv, -one, workl(iuptri+j-1), ldq)
               call zscal(nconv, -one, workl(iuptri+(j-1)*ldq), 1)
            end if
c
 20      continue
c
         if (howmny .eq. 'A') then
c
c           %--------------------------------------------%
c           | Compute the NCONV wanted eigenvectors of T |
c           | located in workl(iuptri,ldq).              |
c           %--------------------------------------------%
c
            do 30 j=1, ncv
               if (j .le. nconv) then
                  select(j) = .true.
               else
                  select(j) = .false.
               end if
 30         continue
c
            call ztrevc('Right', 'Select'     , select       ,
     &                   ncv    , workl(iuptri), ldq          ,
     &                   vl     , 1            , workl(invsub),
     &                   ldq    , ncv          , outncv       ,
     &                   workev , rwork        , ierr)
c
            if (ierr .ne. 0) then
                info = -9
                go to 9000
            end if
c
c           %------------------------------------------------%
c           | Scale the returning eigenvectors so that their |
c           | Euclidean norms are all one. LAPACK subroutine |
c           | ztrevc returns each eigenvector normalized so  |
c           | that the element of largest magnitude has      |
c           | magnitude 1.                                   |
c           %------------------------------------------------%
c
            do 40 j=1, nconv
                  rtemp = dznrm2(ncv, workl(invsub+(j-1)*ldq), 1)
                  rtemp = dble(one) / rtemp
                  call zdscal ( ncv, rtemp,
     &                 workl(invsub+(j-1)*ldq), 1 )
c
c                 %------------------------------------------%
c                 | Ritz estimates can be obtained by taking |
c                 | the inner product of the last row of the |
c                 | Schur basis of H with eigenvectors of T. |
c                 | Note that the eigenvector matrix of T is |
c                 | upper triangular, thus the length of the |
c                 | inner product can be set to j.           |
c                 %------------------------------------------%
c 
                  workev(j) = zdotc(j, workl(ihbds), 1,
     &                        workl(invsub+(j-1)*ldq), 1)
 40         continue
c
            if (msglvl .gt. 2) then
               call zcopy(nconv, workl(invsub+ncv-1), ldq,
     &                    workl(ihbds), 1)
               call zvout (logfil, nconv, workl(ihbds), ndigit,
     &            '_neupd: Last row of the eigenvector matrix for T')
               if (msglvl .gt. 3) then
                  call zmout(logfil       , ncv, ncv   ,
     &                        workl(invsub), ldq, ndigit,
     &               '_neupd: The eigenvector matrix for T')
               end if
            end if
c
c           %---------------------------------------%
c           | Copy Ritz estimates into workl(ihbds) |
c           %---------------------------------------%
c 
            call zcopy(nconv, workev, 1, workl(ihbds), 1)
c
c           %----------------------------------------------%
c           | The eigenvector matrix Q of T is triangular. |
c           | Form Z*Q.                                    |
c           %----------------------------------------------%
c
            call ztrmm('Right'   , 'Upper'      , 'No transpose',
     &                  'Non-unit', n            , nconv         ,
     &                  one       , workl(invsub), ldq           ,
     &                  z         , ldz)
         end if 
c
      else
c
c        %--------------------------------------------------%
c        | An approximate invariant subspace is not needed. |
c        | Place the Ritz values computed ZNAUPD into D.    |
c        %--------------------------------------------------%
c
         call zcopy(nconv, workl(ritz), 1, d, 1)
         call zcopy(nconv, workl(ritz), 1, workl(iheig), 1)
         call zcopy(nconv, workl(bounds), 1, workl(ihbds), 1)
c
      end if
c
c     %------------------------------------------------%
c     | Transform the Ritz values and possibly vectors |
c     | and corresponding error bounds of OP to those  |
c     | of A*x = lambda*B*x.                           |
c     %------------------------------------------------%
c
      if (type .eq. 'REGULR') then
c
         if (rvec) 
     &      call zscal(ncv, rnorm, workl(ihbds), 1)
c      
      else
c     
c        %---------------------------------------%
c        |   A spectral transformation was used. |
c        | * Determine the Ritz estimates of the |
c        |   Ritz values in the original system. |
c        %---------------------------------------%
c
         if (rvec) 
     &      call zscal(ncv, rnorm, workl(ihbds), 1)
c    
         do 50 k=1, ncv
            temp = workl(iheig+k-1)
            workl(ihbds+k-1) = workl(ihbds+k-1) / temp / temp
  50     continue
c  
      end if
c
c     %-----------------------------------------------------------%
c     | *  Transform the Ritz values back to the original system. |
c     |    For TYPE = 'SHIFTI' the transformation is              |
c     |             lambda = 1/theta + sigma                      |
c     | NOTES:                                                    |
c     | *The Ritz vectors are not affected by the transformation. |
c     %-----------------------------------------------------------%
c    
      if (type .eq. 'SHIFTI') then
         do 60 k=1, nconv
            d(k) = one / workl(iheig+k-1) + sigma
  60     continue
      end if
c
      if (type .ne. 'REGULR' .and. msglvl .gt. 1) then
         call zvout (logfil, nconv, d, ndigit,
     &     '_neupd: Untransformed Ritz values.')
         call zvout (logfil, nconv, workl(ihbds), ndigit,
     &     '_neupd: Ritz estimates of the untransformed Ritz values.')
      else if ( msglvl .gt. 1) then
         call zvout (logfil, nconv, d, ndigit,
     &     '_neupd: Converged Ritz values.')
         call zvout (logfil, nconv, workl(ihbds), ndigit,
     &     '_neupd: Associated Ritz estimates.')
      end if
c
c     %-------------------------------------------------%
c     | Eigenvector Purification step. Formally perform |
c     | one of inverse subspace iteration. Only used    |
c     | for MODE = 3. See reference 3.                  |
c     %-------------------------------------------------%
c
      if (rvec .and. howmny .eq. 'A' .and. type .eq. 'SHIFTI') then
c
c        %------------------------------------------------%
c        | Purify the computed Ritz vectors by adding a   |
c        | little bit of the residual vector:             |
c        |                      T                         |
c        |          resid(:)*( e    s ) / theta           |
c        |                      NCV                       |
c        | where H s = s theta.                           |
c        %------------------------------------------------%
c
         do 100 j=1, nconv
            if (workl(iheig+j-1) .ne. zero) then
               workev(j) =  workl(invsub+(j-1)*ldq+ncv-1) /
     &                      workl(iheig+j-1)
            endif
 100     continue

c        %---------------------------------------%
c        | Perform a rank one update to Z and    |
c        | purify all the Ritz vectors together. |
c        %---------------------------------------%
c
         call zgeru (n, nconv, one, resid, 1, workev, 1, z, ldz)
c
      end if
c
 9000 continue
c
      return
c     
c     %---------------%
c     | End of zneupd|
c     %---------------%
c
      end


c\BeginDoc
c
c\Name: znaup2 
c
c\Description: 
c  Intermediate level interface called by znaupd .
c
c\Usage:
c  call znaup2 
c     ( IDO, BMAT, N, WHICH, NEV, NP, TOL, RESID, MODE, IUPD,
c       ISHIFT, MXITER, V, LDV, H, LDH, RITZ, BOUNDS, 
c       Q, LDQ, WORKL, IPNTR, WORKD, RWORK, INFO )
c
c\Arguments
c
c  IDO, BMAT, N, WHICH, NEV, TOL, RESID: same as defined in znaupd .
c  MODE, ISHIFT, MXITER: see the definition of IPARAM in znaupd .
c
c  NP      Integer.  (INPUT/OUTPUT)
c          Contains the number of implicit shifts to apply during
c          each Arnoldi iteration.
c          If ISHIFT=1, NP is adjusted dynamically at each iteration
c          to accelerate convergence and prevent stagnation.
c          This is also roughly equal to the number of matrix-vector
c          products (involving the operator OP) per Arnoldi iteration.
c          The logic for adjusting is contained within the current
c          subroutine.
c          If ISHIFT=0, NP is the number of shifts the user needs
c          to provide via reverse comunication. 0 < NP < NCV-NEV.
c          NP may be less than NCV-NEV since a leading block of the current
c          upper Hessenberg matrix has split off and contains "unwanted"
c          Ritz values.
c          Upon termination of the IRA iteration, NP contains the number
c          of "converged" wanted Ritz values.
c
c  IUPD    Integer.  (INPUT)
c          IUPD .EQ. 0: use explicit restart instead implicit update.
c          IUPD .NE. 0: use implicit update.
c
c  V       Complex*16  N by (NEV+NP) array.  (INPUT/OUTPUT)
c          The Arnoldi basis vectors are returned in the first NEV 
c          columns of V.
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling 
c          program.
c
c  H       Complex*16  (NEV+NP) by (NEV+NP) array.  (OUTPUT)
c          H is used to store the generated upper Hessenberg matrix
c
c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling 
c          program.
c
c  RITZ    Complex*16  array of length NEV+NP.  (OUTPUT)
c          RITZ(1:NEV)  contains the computed Ritz values of OP.
c
c  BOUNDS  Complex*16  array of length NEV+NP.  (OUTPUT)
c          BOUNDS(1:NEV) contain the error bounds corresponding to 
c          the computed Ritz values.
c          
c  Q       Complex*16  (NEV+NP) by (NEV+NP) array.  (WORKSPACE)
c          Private (replicated) work array used to accumulate the
c          rotation in the shift application step.
c
c  LDQ     Integer.  (INPUT)
c          Leading dimension of Q exactly as declared in the calling
c          program.
c
c  WORKL   Complex*16  work array of length at least 
c          (NEV+NP)**2 + 3*(NEV+NP).  (WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.  It is used in shifts calculation, shifts
c          application and convergence checking.
c
c
c  IPNTR   Integer array of length 3.  (OUTPUT)
c          Pointer to mark the starting locations in the WORKD for 
c          vectors used by the Arnoldi iteration.
c          -------------------------------------------------------------
c          IPNTR(1): pointer to the current operand vector X.
c          IPNTR(2): pointer to the current result vector Y.
c          IPNTR(3): pointer to the vector B * X when used in the 
c                    shift-and-invert mode.  X is the current operand.
c          -------------------------------------------------------------
c          
c  WORKD   Complex*16  work array of length 3*N.  (WORKSPACE)
c          Distributed array to be used in the basic Arnoldi iteration
c          for reverse communication.  The user should not use WORKD
c          as temporary workspace during the iteration !!!!!!!!!!
c          See Data Distribution Note in ZNAUPD .
c
c  RWORK   Double precision    work array of length  NEV+NP ( WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.
c
c  INFO    Integer.  (INPUT/OUTPUT)
c          If INFO .EQ. 0, a randomly initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          Error flag on output.
c          =     0: Normal return.
c          =     1: Maximum number of iterations taken.
c                   All possible eigenvalues of OP has been found.  
c                   NP returns the number of converged Ritz values.
c          =     2: No shifts could be applied.
c          =    -8: Error return from LAPACK eigenvalue calculation;
c                   This should never happen.
c          =    -9: Starting vector is zero.
c          = -9999: Could not build an Arnoldi factorization.
c                   Size that was built in returned in NP.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  Complex*16 
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
c     Restarted Arnoldi Iteration", Rice University Technical Report
c     TR95-13, Department of Computational and Applied Mathematics.
c
c\Routines called:
c     zgetv0   ARPACK initial vector generation routine. 
c     znaitr   ARPACK Arnoldi factorization routine.
c     znapps   ARPACK application of implicit shifts routine.
c     zneigh   ARPACK compute Ritz values and error bounds routine. 
c     zngets   ARPACK reorder Ritz values and error bounds routine.
c     zsortc   ARPACK sorting routine.
c     ivout   ARPACK utility routine that prints integers.
c     second  ARPACK utility routine for timing.
c     zmout    ARPACK utility routine that prints matrices
c     zvout    ARPACK utility routine that prints vectors.
c     dvout    ARPACK utility routine that prints vectors.
c     dlamch   LAPACK routine that determines machine constants.
c     dlapy2   LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     zcopy    Level 1 BLAS that copies one vector to another .
c     zdotc    Level 1 BLAS that computes the scalar product of two vectors. 
c     zswap    Level 1 BLAS that swaps two vectors.
c     dznrm2   Level 1 BLAS that computes the norm of a vector.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice Universitya
c     Chao Yang                    Houston, Texas
c     Dept. of Computational &
c     Applied Mathematics 
c     Rice University           
c     Houston, Texas 
c 
c\SCCS Information: @(#)
c FILE: naup2.F   SID: 2.6   DATE OF SID: 06/01/00   RELEASE: 2
c
c\Remarks
c     1. None
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine znaup2 
     &   ( ido, bmat, n, which, nev, np, tol, resid, mode, iupd, 
     &     ishift, mxiter, v, ldv, h, ldh, ritz, bounds, 
     &     q, ldq, workl, ipntr, workd, rwork, info )
c
c     %----------------------------------------------------%
c     | Include files for debugging and timing information |
c     %----------------------------------------------------%
c
      include   'debug.h'
      include   'stat.h'
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character  bmat*1, which*2
      integer    ido, info, ishift, iupd, mode, ldh, ldq, ldv, mxiter,
     &           n, nev, np
      Double precision   
     &           tol
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    ipntr(13)
      Complex*16 
     &           bounds(nev+np), h(ldh,nev+np), q(ldq,nev+np), 
     &           resid(n), ritz(nev+np),  v(ldv,nev+np), 
     &           workd(3*n), workl( (nev+np)*(nev+np+3) )
       Double precision   
     &           rwork(nev+np)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Complex*16 
     &           one, zero
      Double precision 
     &           rzero
      parameter (one = (1.0D+0, 0.0D+0) , zero = (0.0D+0, 0.0D+0) ,
     &           rzero = 0.0D+0 )
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      logical    cnorm , getv0, initv , update, ushift
      integer    ierr  , iter , kplusp, msglvl, nconv, 
     &           nevbef, nev0 , np0   , nptemp, i    ,
     &           j    
      Complex*16 
     &           cmpnorm
      Double precision 
     &           rnorm , eps23, rtemp
      character  wprime*2
c
      save       cnorm,  getv0, initv , update, ushift, 
     &           rnorm,  iter , kplusp, msglvl, nconv ,
     &           nevbef, nev0 , np0   , eps23
c
c
c     %-----------------------%
c     | Local array arguments |
c     %-----------------------%
c
      integer    kp(3)
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   zcopy , zgetv0 , znaitr , zneigh , zngets , znapps ,
     &           zsortc , zswap , zmout , zvout , ivout, second
c
c     %--------------------%
c     | External functions |
c     %--------------------%
c
      Complex*16 
     &           zdotc 
      Double precision   
     &           dznrm2 , dlamch , dlapy2 
      external   zdotc , dznrm2 , dlamch , dlapy2 
c
c     %---------------------%
c     | Intrinsic Functions |
c     %---------------------%
c
      intrinsic  dimag , dble , min, max
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      if (ido .eq. 0) then
c 
         call second (t0)
c 
         msglvl = mcaup2
c 
         nev0   = nev
         np0    = np
c
c        %-------------------------------------%
c        | kplusp is the bound on the largest  |
c        |        Lanczos factorization built. |
c        | nconv is the current number of      |
c        |        "converged" eigenvalues.     |
c        | iter is the counter on the current  |
c        |      iteration step.                |
c        %-------------------------------------%
c
         kplusp = nev + np
         nconv  = 0
         iter   = 0
c 
c        %---------------------------------%
c        | Get machine dependent constant. |
c        %---------------------------------%
c
         eps23 = dlamch ('Epsilon-Machine')
         eps23 = eps23**(2.0D+0  / 3.0D+0 )
c
c        %---------------------------------------%
c        | Set flags for computing the first NEV |
c        | steps of the Arnoldi factorization.   |
c        %---------------------------------------%
c
         getv0    = .true.
         update   = .false.
         ushift   = .false.
         cnorm    = .false.
c
         if (info .ne. 0) then
c
c           %--------------------------------------------%
c           | User provides the initial residual vector. |
c           %--------------------------------------------%
c
            initv = .true.
            info  = 0
         else
            initv = .false.
         end if
      end if
c 
c     %---------------------------------------------%
c     | Get a possibly random starting vector and   |
c     | force it into the range of the operator OP. |
c     %---------------------------------------------%
c
   10 continue
c
      if (getv0) then
         call zgetv0  (ido, bmat, 1, initv, n, 1, v, ldv, resid, rnorm,
     &                ipntr, workd, info)
c
         if (ido .ne. 99) go to 9000
c
         if (rnorm .eq. rzero) then
c
c           %-----------------------------------------%
c           | The initial vector is zero. Error exit. | 
c           %-----------------------------------------%
c
            info = -9
            go to 1100
         end if
         getv0 = .false.
         ido  = 0
      end if
c 
c     %-----------------------------------%
c     | Back from reverse communication : |
c     | continue with update step         |
c     %-----------------------------------%
c
      if (update) go to 20
c
c     %-------------------------------------------%
c     | Back from computing user specified shifts |
c     %-------------------------------------------%
c
      if (ushift) go to 50
c
c     %-------------------------------------%
c     | Back from computing residual norm   |
c     | at the end of the current iteration |
c     %-------------------------------------%
c
      if (cnorm)  go to 100
c 
c     %----------------------------------------------------------%
c     | Compute the first NEV steps of the Arnoldi factorization |
c     %----------------------------------------------------------%
c
      call znaitr  (ido, bmat, n, 0, nev, mode, resid, rnorm, v, ldv, 
     &             h, ldh, ipntr, workd, info)
c
      if (ido .ne. 99) go to 9000
c
      if (info .gt. 0) then
         np   = info
         mxiter = iter
         info = -9999
         go to 1200
      end if
c 
c     %--------------------------------------------------------------%
c     |                                                              |
c     |           M A I N  ARNOLDI  I T E R A T I O N  L O O P       |
c     |           Each iteration implicitly restarts the Arnoldi     |
c     |           factorization in place.                            |
c     |                                                              |
c     %--------------------------------------------------------------%
c 
 1000 continue
c
         iter = iter + 1
c
         if (msglvl .gt. 0) then
            call ivout (logfil, 1, iter, ndigit, 
     &           '_naup2: **** Start of major iteration number ****')
         end if
c 
c        %-----------------------------------------------------------%
c        | Compute NP additional steps of the Arnoldi factorization. |
c        | Adjust NP since NEV might have been updated by last call  |
c        | to the shift application routine znapps .                  |
c        %-----------------------------------------------------------%
c
         np  = kplusp - nev
c
         if (msglvl .gt. 1) then
            call ivout (logfil, 1, nev, ndigit, 
     &     '_naup2: The length of the current Arnoldi factorization')
            call ivout (logfil, 1, np, ndigit, 
     &           '_naup2: Extend the Arnoldi factorization by')
         end if
c
c        %-----------------------------------------------------------%
c        | Compute NP additional steps of the Arnoldi factorization. |
c        %-----------------------------------------------------------%
c
         ido = 0
   20    continue
         update = .true.
c
         call znaitr (ido, bmat, n, nev, np,    mode,  resid, rnorm,
     &               v  , ldv , h, ldh, ipntr, workd, info)
c
         if (ido .ne. 99) go to 9000
c
         if (info .gt. 0) then
            np = info
            mxiter = iter
            info = -9999
            go to 1200
         end if
         update = .false.
c
         if (msglvl .gt. 1) then
            call dvout  (logfil, 1, rnorm, ndigit, 
     &           '_naup2: Corresponding B-norm of the residual')
         end if
c 
c        %--------------------------------------------------------%
c        | Compute the eigenvalues and corresponding error bounds |
c        | of the current upper Hessenberg matrix.                |
c        %--------------------------------------------------------%
c
         call zneigh  (rnorm, kplusp, h, ldh, ritz, bounds,
     &                q, ldq, workl, rwork,  ierr)
c
         if (ierr .ne. 0) then
            info = -8
            go to 1200
         end if
c
c        %---------------------------------------------------%
c        | Select the wanted Ritz values and their bounds    |
c        | to be used in the convergence test.               |
c        | The wanted part of the spectrum and corresponding |
c        | error bounds are in the last NEV loc. of RITZ,    |
c        | and BOUNDS respectively.                          | 
c        %---------------------------------------------------%
c
         nev = nev0
         np = np0
c
c        %--------------------------------------------------%
c        | Make a copy of Ritz values and the corresponding |
c        | Ritz estimates obtained from zneigh .             |
c        %--------------------------------------------------%
c
         call zcopy (kplusp,ritz,1,workl(kplusp**2+1),1)
         call zcopy (kplusp,bounds,1,workl(kplusp**2+kplusp+1),1)
c
c        %---------------------------------------------------%
c        | Select the wanted Ritz values and their bounds    |
c        | to be used in the convergence test.               |
c        | The wanted part of the spectrum and corresponding |
c        | bounds are in the last NEV loc. of RITZ           |
c        | BOUNDS respectively.                              |
c        %---------------------------------------------------%
c
         call zngets  (ishift, which, nev, np, ritz, bounds)
c 
c        %------------------------------------------------------------%
c        | Convergence test: currently we use the following criteria. |
c        | The relative accuracy of a Ritz value is considered        |
c        | acceptable if:                                             |
c        |                                                            |
c        | error_bounds(i) .le. tol*max(eps23, magnitude_of_ritz(i)). |
c        |                                                            |
c        %------------------------------------------------------------%
c
         nconv  = 0
c
         do 25 i = 1, nev
            rtemp = max( eps23, dlapy2 ( dble (ritz(np+i)),
     &                                  dimag (ritz(np+i)) ) ) 
            if ( dlapy2 (dble (bounds(np+i)),dimag (bounds(np+i))) 
     &                 .le. tol*rtemp ) then
               nconv = nconv + 1
            end if
   25    continue
c 
         if (msglvl .gt. 2) then
            kp(1) = nev
            kp(2) = np
            kp(3) = nconv
            call ivout (logfil, 3, kp, ndigit, 
     &                  '_naup2: NEV, NP, NCONV are')
            call zvout  (logfil, kplusp, ritz, ndigit,
     &           '_naup2: The eigenvalues of H')
            call zvout  (logfil, kplusp, bounds, ndigit, 
     &          '_naup2: Ritz estimates of the current NCV Ritz values')
         end if
c
c        %---------------------------------------------------------%
c        | Count the number of unwanted Ritz values that have zero |
c        | Ritz estimates. If any Ritz estimates are equal to zero |
c        | then a leading block of H of order equal to at least    |
c        | the number of Ritz values with zero Ritz estimates has  |
c        | split off. None of these Ritz values may be removed by  |
c        | shifting. Decrease NP the number of shifts to apply. If |
c        | no shifts may be applied, then prepare to exit          |
c        %---------------------------------------------------------%
c
         nptemp = np
         do 30 j=1, nptemp
            if (bounds(j) .eq. zero) then
               np = np - 1
               nev = nev + 1
            end if
 30      continue
c     
         if ( (nconv .ge. nev0) .or. 
     &        (iter .gt. mxiter) .or.
     &        (np .eq. 0) ) then
c
            if (msglvl .gt. 4) then
               call zvout (logfil, kplusp, workl(kplusp**2+1), ndigit,
     &             '_naup2: Eigenvalues computed by _neigh:')
               call zvout (logfil, kplusp, workl(kplusp**2+kplusp+1),
     &                     ndigit,
     &             '_naup2: Ritz estimates computed by _neigh:')
            end if
c     
c           %------------------------------------------------%
c           | Prepare to exit. Put the converged Ritz values |
c           | and corresponding bounds in RITZ(1:NCONV) and  |
c           | BOUNDS(1:NCONV) respectively. Then sort. Be    |
c           | careful when NCONV > NP                        |
c           %------------------------------------------------%
c
c           %------------------------------------------%
c           |  Use h( 3,1 ) as storage to communicate  |
c           |  rnorm to zneupd  if needed               |
c           %------------------------------------------%

            h(3,1) = dcmplx (rnorm,rzero)
c
c           %----------------------------------------------%
c           | Sort Ritz values so that converged Ritz      |
c           | values appear within the first NEV locations |
c           | of ritz and bounds, and the most desired one |
c           | appears at the front.                        |
c           %----------------------------------------------%
c
            if (which .eq. 'LM') wprime = 'SM'
            if (which .eq. 'SM') wprime = 'LM'
            if (which .eq. 'LR') wprime = 'SR'
            if (which .eq. 'SR') wprime = 'LR'
            if (which .eq. 'LI') wprime = 'SI'
            if (which .eq. 'SI') wprime = 'LI'
c
            call zsortc (wprime, .true., kplusp, ritz, bounds)
c
c           %--------------------------------------------------%
c           | Scale the Ritz estimate of each Ritz value       |
c           | by 1 / max(eps23, magnitude of the Ritz value).  |
c           %--------------------------------------------------%
c
            do 35 j = 1, nev0 
                rtemp = max( eps23, dlapy2 ( dble (ritz(j)),
     &                                       dimag (ritz(j)) ) )
                bounds(j) = bounds(j)/rtemp
 35         continue
c
c           %---------------------------------------------------%
c           | Sort the Ritz values according to the scaled Ritz |
c           | estimates.  This will push all the converged ones |
c           | towards the front of ritz, bounds (in the case    |
c           | when NCONV < NEV.)                                |
c           %---------------------------------------------------%
c
            wprime = 'LM'
            call zsortc (wprime, .true., nev0, bounds, ritz)
c
c           %----------------------------------------------%
c           | Scale the Ritz estimate back to its original |
c           | value.                                       |
c           %----------------------------------------------%
c
            do 40 j = 1, nev0
                rtemp = max( eps23, dlapy2 ( dble (ritz(j)),
     &                                       dimag (ritz(j)) ) )
                bounds(j) = bounds(j)*rtemp
 40         continue
c
c           %-----------------------------------------------%
c           | Sort the converged Ritz values again so that  |
c           | the "threshold" value appears at the front of |
c           | ritz and bound.                               |
c           %-----------------------------------------------%
c
            call zsortc (which, .true., nconv, ritz, bounds)
c
            if (msglvl .gt. 1) then
               call zvout  (logfil, kplusp, ritz, ndigit,
     &            '_naup2: Sorted eigenvalues')
               call zvout  (logfil, kplusp, bounds, ndigit,
     &            '_naup2: Sorted ritz estimates.')
            end if
c
c           %------------------------------------%
c           | Max iterations have been exceeded. | 
c           %------------------------------------%
c
            if (iter .gt. mxiter .and. nconv .lt. nev0) info = 1
c
c           %---------------------%
c           | No shifts to apply. | 
c           %---------------------%
c
            if (np .eq. 0 .and. nconv .lt. nev0)  info = 2
c
            np = nconv
            go to 1100
c
         else if ( (nconv .lt. nev0) .and. (ishift .eq. 1) ) then
c     
c           %-------------------------------------------------%
c           | Do not have all the requested eigenvalues yet.  |
c           | To prevent possible stagnation, adjust the size |
c           | of NEV.                                         |
c           %-------------------------------------------------%
c
            nevbef = nev
            nev = nev + min(nconv, np/2)
            if (nev .eq. 1 .and. kplusp .ge. 6) then
               nev = kplusp / 2
            else if (nev .eq. 1 .and. kplusp .gt. 3) then
               nev = 2
            end if
            np = kplusp - nev
c     
c           %---------------------------------------%
c           | If the size of NEV was just increased |
c           | resort the eigenvalues.               |
c           %---------------------------------------%
c     
            if (nevbef .lt. nev) 
     &         call zngets  (ishift, which, nev, np, ritz, bounds)
c
         end if              
c     
         if (msglvl .gt. 0) then
            call ivout (logfil, 1, nconv, ndigit, 
     &           '_naup2: no. of "converged" Ritz values at this iter.')
            if (msglvl .gt. 1) then
               kp(1) = nev
               kp(2) = np
               call ivout (logfil, 2, kp, ndigit, 
     &              '_naup2: NEV and NP are')
               call zvout  (logfil, nev, ritz(np+1), ndigit,
     &              '_naup2: "wanted" Ritz values ')
               call zvout  (logfil, nev, bounds(np+1), ndigit,
     &              '_naup2: Ritz estimates of the "wanted" values ')
            end if
         end if
c
         if (ishift .eq. 0) then
c
c           %-------------------------------------------------------%
c           | User specified shifts: pop back out to get the shifts |
c           | and return them in the first 2*NP locations of WORKL. |
c           %-------------------------------------------------------%
c
            ushift = .true.
            ido = 3
            go to 9000
         end if
   50    continue
         ushift = .false.
c
         if ( ishift .ne. 1 ) then
c 
c            %----------------------------------%
c            | Move the NP shifts from WORKL to |
c            | RITZ, to free up WORKL           |
c            | for non-exact shift case.        |
c            %----------------------------------%
c
             call zcopy  (np, workl, 1, ritz, 1)
         end if
c
         if (msglvl .gt. 2) then 
            call ivout (logfil, 1, np, ndigit, 
     &                  '_naup2: The number of shifts to apply ')
            call zvout  (logfil, np, ritz, ndigit,
     &                  '_naup2: values of the shifts')
            if ( ishift .eq. 1 ) 
     &          call zvout  (logfil, np, bounds, ndigit,
     &                  '_naup2: Ritz estimates of the shifts')
         end if
c
c        %---------------------------------------------------------%
c        | Apply the NP implicit shifts by QR bulge chasing.       |
c        | Each shift is applied to the whole upper Hessenberg     |
c        | matrix H.                                               |
c        | The first 2*N locations of WORKD are used as workspace. |
c        %---------------------------------------------------------%
c
         call znapps  (n, nev, np, ritz, v, ldv, 
     &                h, ldh, resid, q, ldq, workl, workd)
c
c        %---------------------------------------------%
c        | Compute the B-norm of the updated residual. |
c        | Keep B*RESID in WORKD(1:N) to be used in    |
c        | the first step of the next call to znaitr .  |
c        %---------------------------------------------%
c
         cnorm = .true.
         call second (t2)
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            call zcopy  (n, resid, 1, workd(n+1), 1)
            ipntr(1) = n + 1
            ipntr(2) = 1
            ido = 2
c 
c           %----------------------------------%
c           | Exit in order to compute B*RESID |
c           %----------------------------------%
c 
            go to 9000
         else if (bmat .eq. 'I') then
            call zcopy  (n, resid, 1, workd, 1)
         end if
c 
  100    continue
c 
c        %----------------------------------%
c        | Back from reverse communication; |
c        | WORKD(1:N) := B*RESID            |
c        %----------------------------------%
c
         if (bmat .eq. 'G') then
            call second (t3)
            tmvbx = tmvbx + (t3 - t2)
         end if
c 
         if (bmat .eq. 'G') then         
            cmpnorm = zdotc  (n, resid, 1, workd, 1)
            rnorm = sqrt(dlapy2 (dble (cmpnorm),dimag (cmpnorm)))
         else if (bmat .eq. 'I') then
            rnorm = dznrm2 (n, resid, 1)
         end if
         cnorm = .false.
c
         if (msglvl .gt. 2) then
            call dvout  (logfil, 1, rnorm, ndigit, 
     &      '_naup2: B-norm of residual for compressed factorization')
            call zmout  (logfil, nev, nev, h, ldh, ndigit,
     &        '_naup2: Compressed upper Hessenberg matrix H')
         end if
c 
      go to 1000
c
c     %---------------------------------------------------------------%
c     |                                                               |
c     |  E N D     O F     M A I N     I T E R A T I O N     L O O P  |
c     |                                                               |
c     %---------------------------------------------------------------%
c
 1100 continue
c
      mxiter = iter
      nev = nconv
c     
 1200 continue
      ido = 99
c
c     %------------%
c     | Error Exit |
c     %------------%
c
      call second (t1)
      tcaup2 = t1 - t0
c     
 9000 continue
c
c     %---------------%
c     | End of znaup2  |
c     %---------------%
c
      return
      end


c
c\SCCS Information: @(#)
c FILE: statn.F   SID: 2.2   DATE OF SID: 4/20/96   RELEASE: 2
c
c     %---------------------------------------------%
c     | Initialize statistic and timing information |
c     | for complex nonsymmetric Arnoldi code.      |
c     %---------------------------------------------%

      subroutine zstatn
c
c     %--------------------------------%
c     | See stat.doc for documentation |
c     %--------------------------------%
c
      include   'stat.h'
 
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%

      nopx   = 0
      nbx    = 0
      nrorth = 0
      nitref = 0
      nrstrt = 0
 
      tcaupd = 0.0D+0
      tcaup2 = 0.0D+0
      tcaitr = 0.0D+0
      tceigh = 0.0D+0
      tcgets = 0.0D+0
      tcapps = 0.0D+0
      tcconv = 0.0D+0
      titref = 0.0D+0
      tgetv0 = 0.0D+0
      trvec  = 0.0D+0
 
c     %----------------------------------------------------%
c     | User time including reverse communication overhead |
c     %----------------------------------------------------%
      tmvopx = 0.0D+0
      tmvbx  = 0.0D+0
 
      return
c
c     %---------------%
c     | End of zstatn |
c     %---------------%
c
      end



c\BeginDoc
c
c\Name: zngets
c
c\Description: 
c  Given the eigenvalues of the upper Hessenberg matrix H,
c  computes the NP shifts AMU that are zeros of the polynomial of 
c  degree NP which filters out components of the unwanted eigenvectors
c  corresponding to the AMU's based on some given criteria.
c
c  NOTE: call this even in the case of user specified shifts in order
c  to sort the eigenvalues, and error bounds of H for later use.
c
c\Usage:
c  call zngets
c      ( ISHIFT, WHICH, KEV, NP, RITZ, BOUNDS )
c
c\Arguments
c  ISHIFT  Integer.  (INPUT)
c          Method for selecting the implicit shifts at each iteration.
c          ISHIFT = 0: user specified shifts
c          ISHIFT = 1: exact shift with respect to the matrix H.
c
c  WHICH   Character*2.  (INPUT)
c          Shift selection criteria.
c          'LM' -> want the KEV eigenvalues of largest magnitude.
c          'SM' -> want the KEV eigenvalues of smallest magnitude.
c          'LR' -> want the KEV eigenvalues of largest REAL part.
c          'SR' -> want the KEV eigenvalues of smallest REAL part.
c          'LI' -> want the KEV eigenvalues of largest imaginary part.
c          'SI' -> want the KEV eigenvalues of smallest imaginary part.
c
c  KEV     Integer.  (INPUT)
c          The number of desired eigenvalues.
c
c  NP      Integer.  (INPUT)
c          The number of shifts to compute.
c
c  RITZ    Complex*16 array of length KEV+NP.  (INPUT/OUTPUT)
c          On INPUT, RITZ contains the the eigenvalues of H.
c          On OUTPUT, RITZ are sorted so that the unwanted
c          eigenvalues are in the first NP locations and the wanted
c          portion is in the last KEV locations.  When exact shifts are 
c          selected, the unwanted part corresponds to the shifts to 
c          be applied. Also, if ISHIFT .eq. 1, the unwanted eigenvalues
c          are further sorted so that the ones with largest Ritz values
c          are first.
c
c  BOUNDS  Complex*16 array of length KEV+NP.  (INPUT/OUTPUT)
c          Error bounds corresponding to the ordering in RITZ.
c
c  
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  Complex*16
c
c\Routines called:
c     zsortc  ARPACK sorting routine.
c     ivout   ARPACK utility routine that prints integers.
c     second  ARPACK utility routine for timing.
c     zvout   ARPACK utility routine that prints vectors.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics 
c     Rice University           
c     Houston, Texas 
c
c\SCCS Information: @(#)
c FILE: ngets.F   SID: 2.2   DATE OF SID: 4/20/96   RELEASE: 2
c
c\Remarks
c     1. This routine does not keep complex conjugate pairs of
c        eigenvalues together.
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine zngets ( ishift, which, kev, np, ritz, bounds)
c
c     %----------------------------------------------------%
c     | Include files for debugging and timing information |
c     %----------------------------------------------------%
c
      include   'debug.h'
      include   'stat.h'
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character*2 which
      integer    ishift, kev, np
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Complex*16
     &           bounds(kev+np), ritz(kev+np)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Complex*16
     &           one, zero
      parameter (one = (1.0D+0, 0.0D+0), zero = (0.0D+0, 0.0D+0))
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer    msglvl
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   zvout,  zsortc, second
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
c     %-------------------------------%
c     | Initialize timing statistics  |
c     | & message level for debugging |
c     %-------------------------------%
c 
      call second (t0)
      msglvl = mcgets
c 
      call zsortc (which, .true., kev+np, ritz, bounds)
c     
      if ( ishift .eq. 1 ) then
c     
c        %-------------------------------------------------------%
c        | Sort the unwanted Ritz values used as shifts so that  |
c        | the ones with largest Ritz estimates are first        |
c        | This will tend to minimize the effects of the         |
c        | forward instability of the iteration when the shifts  |
c        | are applied in subroutine znapps.                     |
c        | Be careful and use 'SM' since we want to sort BOUNDS! |
c        %-------------------------------------------------------%
c     
         call zsortc ( 'SM', .true., np, bounds, ritz )
c
      end if
c     
      call second (t1)
      tcgets = tcgets + (t1 - t0)
c
      if (msglvl .gt. 0) then
         call ivout (logfil, 1, kev, ndigit, '_ngets: KEV is')
         call ivout (logfil, 1, np, ndigit, '_ngets: NP is')
         call zvout (logfil, kev+np, ritz, ndigit,
     &        '_ngets: Eigenvalues of current H matrix ')
         call zvout (logfil, kev+np, bounds, ndigit, 
     &      '_ngets: Ritz estimates of the current KEV+NP Ritz values')
      end if
c     
      return
c     
c     %---------------%
c     | End of zngets |
c     %---------------%
c     
      end



c\BeginDoc
c
c\Name: zgetv0
c
c\Description: 
c  Generate a random initial residual vector for the Arnoldi process.
c  Force the residual vector to be in the range of the operator OP.  
c
c\Usage:
c  call zgetv0
c     ( IDO, BMAT, ITRY, INITV, N, J, V, LDV, RESID, RNORM, 
c       IPNTR, WORKD, IERR )
c
c\Arguments
c  IDO     Integer.  (INPUT/OUTPUT)
c          Reverse communication flag.  IDO must be zero on the first
c          call to zgetv0.
c          -------------------------------------------------------------
c          IDO =  0: first call to the reverse communication interface
c          IDO = -1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c                    This is for the initialization phase to force the
c                    starting vector into the range of OP.
c          IDO =  2: compute  Y = B * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c          IDO = 99: done
c          -------------------------------------------------------------
c
c  BMAT    Character*1.  (INPUT)
c          BMAT specifies the type of the matrix B in the (generalized)
c          eigenvalue problem A*x = lambda*B*x.
c          B = 'I' -> standard eigenvalue problem A*x = lambda*x
c          B = 'G' -> generalized eigenvalue problem A*x = lambda*B*x
c
c  ITRY    Integer.  (INPUT)
c          ITRY counts the number of times that zgetv0 is called.  
c          It should be set to 1 on the initial call to zgetv0.
c
c  INITV   Logical variable.  (INPUT)
c          .TRUE.  => the initial residual vector is given in RESID.
c          .FALSE. => generate a random initial residual vector.
c
c  N       Integer.  (INPUT)
c          Dimension of the problem.
c
c  J       Integer.  (INPUT)
c          Index of the residual vector to be generated, with respect to
c          the Arnoldi process.  J > 1 in case of a "restart".
c
c  V       Complex*16 N by J array.  (INPUT)
c          The first J-1 columns of V contain the current Arnoldi basis
c          if this is a "restart".
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling 
c          program.
c
c  RESID   Complex*16 array of length N.  (INPUT/OUTPUT)
c          Initial residual vector to be generated.  If RESID is 
c          provided, force RESID into the range of the operator OP.
c
c  RNORM   Double precision scalar.  (OUTPUT)
c          B-norm of the generated residual.
c
c  IPNTR   Integer array of length 3.  (OUTPUT)
c
c  WORKD   Complex*16 work array of length 2*N.  (REVERSE COMMUNICATION).
c          On exit, WORK(1:N) = B*RESID to be used in SSAITR.
c
c  IERR    Integer.  (OUTPUT)
c          =  0: Normal exit.
c          = -1: Cannot generate a nontrivial restarted residual vector
c                in the range of the operator OP.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  Complex*16
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c
c\Routines called:
c     second  ARPACK utility routine for timing.
c     zvout   ARPACK utility routine that prints vectors.
c     zlarnv  LAPACK routine for generating a random vector. 
c     zgemv   Level 2 BLAS routine for matrix vector multiplication.
c     zcopy   Level 1 BLAS that copies one vector to another.
c     zdotc   Level 1 BLAS that computes the scalar product of two vectors.
c     dznrm2  Level 1 BLAS that computes the norm of a vector. 
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics 
c     Rice University           
c     Houston, Texas            
c
c\SCCS Information: @(#)
c FILE: getv0.F   SID: 2.3   DATE OF SID: 08/27/96   RELEASE: 2
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine zgetv0 
     &   ( ido, bmat, itry, initv, n, j, v, ldv, resid, rnorm, 
     &     ipntr, workd, ierr )
c 
c     %----------------------------------------------------%
c     | Include files for debugging and timing information |
c     %----------------------------------------------------%
c
      include   'debug.h'
      include   'stat.h'
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character  bmat*1
      logical    initv
      integer    ido, ierr, itry, j, ldv, n
      Double precision
     &           rnorm
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    ipntr(3)
      Complex*16
     &           resid(n), v(ldv,j), workd(2*n)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Complex*16
     &           one, zero
      Double precision
     &           rzero
      parameter  (one = (1.0D+0, 0.0D+0), zero = (0.0D+0, 0.0D+0),
     &            rzero = 0.0D+0)
c
c     %------------------------%
c     | Local Scalars & Arrays |
c     %------------------------%
c
      logical    first, inits, orth
      integer    idist, iseed(4), iter, msglvl, jj
      Double precision
     &           rnorm0
      Complex*16
     &           cnorm
      save       first, iseed, inits, iter, msglvl, orth, rnorm0
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   zcopy, zgemv, zlarnv, zvout, second
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision 
     &           dznrm2, dlapy2
      Complex*16
     &           zdotc
      external   zdotc, dznrm2, dlapy2
c
c     %-----------------%
c     | Data Statements |
c     %-----------------%
c
      data       inits /.true./
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
c
c     %-----------------------------------%
c     | Initialize the seed of the LAPACK |
c     | random number generator           |
c     %-----------------------------------%
c
      if (inits) then
          iseed(1) = 1
          iseed(2) = 3
          iseed(3) = 5
          iseed(4) = 7
          inits = .false.
      end if
c
      if (ido .eq.  0) then
c 
c        %-------------------------------%
c        | Initialize timing statistics  |
c        | & message level for debugging |
c        %-------------------------------%
c
         call second (t0)
         msglvl = mgetv0
c 
         ierr   = 0
         iter   = 0
         first  = .FALSE.
         orth   = .FALSE.
c
c        %-----------------------------------------------------%
c        | Possibly generate a random starting vector in RESID |
c        | Use a LAPACK random number generator used by the    |
c        | matrix generation routines.                         |
c        |    idist = 1: uniform (0,1)  distribution;          |
c        |    idist = 2: uniform (-1,1) distribution;          |
c        |    idist = 3: normal  (0,1)  distribution;          |
c        %-----------------------------------------------------%
c
         if (.not.initv) then
            idist = 2
            call zlarnv (idist, iseed, n, resid)
         end if
c 
c        %----------------------------------------------------------%
c        | Force the starting vector into the range of OP to handle |
c        | the generalized problem when B is possibly (singular).   |
c        %----------------------------------------------------------%
c
         call second (t2)
         if (bmat .eq. 'G') then
            nopx = nopx + 1
            ipntr(1) = 1
            ipntr(2) = n + 1
            call zcopy (n, resid, 1, workd, 1)
            ido = -1
            go to 9000
         end if
      end if
c 
c     %----------------------------------------%
c     | Back from computing B*(initial-vector) |
c     %----------------------------------------%
c
      if (first) go to 20
c
c     %-----------------------------------------------%
c     | Back from computing B*(orthogonalized-vector) |
c     %-----------------------------------------------%
c
      if (orth)  go to 40
c 
      call second (t3)
      tmvopx = tmvopx + (t3 - t2)
c 
c     %------------------------------------------------------%
c     | Starting vector is now in the range of OP; r = OP*r; |
c     | Compute B-norm of starting vector.                   |
c     %------------------------------------------------------%
c
      call second (t2)
      first = .TRUE.
      if (bmat .eq. 'G') then
         nbx = nbx + 1
         call zcopy (n, workd(n+1), 1, resid, 1)
         ipntr(1) = n + 1
         ipntr(2) = 1
         ido = 2
         go to 9000
      else if (bmat .eq. 'I') then
         call zcopy (n, resid, 1, workd, 1)
      end if
c 
   20 continue
c
      if (bmat .eq. 'G') then
         call second (t3)
         tmvbx = tmvbx + (t3 - t2)
      end if
c 
      first = .FALSE.
      if (bmat .eq. 'G') then
          cnorm  = zdotc (n, resid, 1, workd, 1)
          rnorm0 = sqrt(dlapy2(dble(cnorm),dimag(cnorm)))
      else if (bmat .eq. 'I') then
           rnorm0 = dznrm2(n, resid, 1)
      end if
      rnorm  = rnorm0
c
c     %---------------------------------------------%
c     | Exit if this is the very first Arnoldi step |
c     %---------------------------------------------%
c
      if (j .eq. 1) go to 50
c 
c     %----------------------------------------------------------------
c     | Otherwise need to B-orthogonalize the starting vector against |
c     | the current Arnoldi basis using Gram-Schmidt with iter. ref.  |
c     | This is the case where an invariant subspace is encountered   |
c     | in the middle of the Arnoldi factorization.                   |
c     |                                                               |
c     |       s = V^{T}*B*r;   r = r - V*s;                           |
c     |                                                               |
c     | Stopping criteria used for iter. ref. is discussed in         |
c     | Parlett's book, page 107 and in Gragg & Reichel TOMS paper.   |
c     %---------------------------------------------------------------%
c
      orth = .TRUE.
   30 continue
c
      call zgemv ('C', n, j-1, one, v, ldv, workd, 1, 
     &            zero, workd(n+1), 1)
      call zgemv ('N', n, j-1, -one, v, ldv, workd(n+1), 1, 
     &            one, resid, 1)
c 
c     %----------------------------------------------------------%
c     | Compute the B-norm of the orthogonalized starting vector |
c     %----------------------------------------------------------%
c
      call second (t2)
      if (bmat .eq. 'G') then
         nbx = nbx + 1
         call zcopy (n, resid, 1, workd(n+1), 1)
         ipntr(1) = n + 1
         ipntr(2) = 1
         ido = 2
         go to 9000
      else if (bmat .eq. 'I') then
         call zcopy (n, resid, 1, workd, 1)
      end if
c 
   40 continue
c
      if (bmat .eq. 'G') then
         call second (t3)
         tmvbx = tmvbx + (t3 - t2)
      end if
c 
      if (bmat .eq. 'G') then
         cnorm = zdotc (n, resid, 1, workd, 1)
         rnorm = sqrt(dlapy2(dble(cnorm),dimag(cnorm)))
      else if (bmat .eq. 'I') then
         rnorm = dznrm2(n, resid, 1)
      end if
c
c     %--------------------------------------%
c     | Check for further orthogonalization. |
c     %--------------------------------------%
c
      if (msglvl .gt. 2) then
          call dvout (logfil, 1, rnorm0, ndigit, 
     &                '_getv0: re-orthonalization ; rnorm0 is')
          call dvout (logfil, 1, rnorm, ndigit, 
     &                '_getv0: re-orthonalization ; rnorm is')
      end if
c
      if (rnorm .gt. 0.717*rnorm0) go to 50
c 
      iter = iter + 1
      if (iter .le. 1) then
c
c        %-----------------------------------%
c        | Perform iterative refinement step |
c        %-----------------------------------%
c
         rnorm0 = rnorm
         go to 30
      else
c
c        %------------------------------------%
c        | Iterative refinement step "failed" |
c        %------------------------------------%
c
         do 45 jj = 1, n
            resid(jj) = zero
   45    continue
         rnorm = rzero
         ierr = -1
      end if
c 
   50 continue
c
      if (msglvl .gt. 0) then
         call dvout (logfil, 1, rnorm, ndigit,
     &        '_getv0: B-norm of initial / restarted starting vector')
      end if
      if (msglvl .gt. 2) then
         call zvout (logfil, n, resid, ndigit,
     &        '_getv0: initial / restarted starting vector')
      end if
      ido = 99
c 
      call second (t1)
      tgetv0 = tgetv0 + (t1 - t0)
c 
 9000 continue
      return
c
c     %---------------%
c     | End of zgetv0 |
c     %---------------%
c
      end



c\BeginDoc
c
c\Name: znaitr
c
c\Description: 
c  Reverse communication interface for applying NP additional steps to 
c  a K step nonsymmetric Arnoldi factorization.
c
c  Input:  OP*V_{k}  -  V_{k}*H = r_{k}*e_{k}^T
c
c          with (V_{k}^T)*B*V_{k} = I, (V_{k}^T)*B*r_{k} = 0.
c
c  Output: OP*V_{k+p}  -  V_{k+p}*H = r_{k+p}*e_{k+p}^T
c
c          with (V_{k+p}^T)*B*V_{k+p} = I, (V_{k+p}^T)*B*r_{k+p} = 0.
c
c  where OP and B are as in znaupd.  The B-norm of r_{k+p} is also
c  computed and returned.
c
c\Usage:
c  call znaitr
c     ( IDO, BMAT, N, K, NP, NB, RESID, RNORM, V, LDV, H, LDH, 
c       IPNTR, WORKD, INFO )
c
c\Arguments
c  IDO     Integer.  (INPUT/OUTPUT)
c          Reverse communication flag.
c          -------------------------------------------------------------
c          IDO =  0: first call to the reverse communication interface
c          IDO = -1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORK for X,
c                    IPNTR(2) is the pointer into WORK for Y.
c                    This is for the restart phase to force the new
c                    starting vector into the range of OP.
c          IDO =  1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORK for X,
c                    IPNTR(2) is the pointer into WORK for Y,
c                    IPNTR(3) is the pointer into WORK for B * X.
c          IDO =  2: compute  Y = B * X  where
c                    IPNTR(1) is the pointer into WORK for X,
c                    IPNTR(2) is the pointer into WORK for Y.
c          IDO = 99: done
c          -------------------------------------------------------------
c          When the routine is used in the "shift-and-invert" mode, the
c          vector B * Q is already available and do not need to be
c          recomputed in forming OP * Q.
c
c  BMAT    Character*1.  (INPUT)
c          BMAT specifies the type of the matrix B that defines the
c          semi-inner product for the operator OP.  See znaupd.
c          B = 'I' -> standard eigenvalue problem A*x = lambda*x
c          B = 'G' -> generalized eigenvalue problem A*x = lambda*M**x
c
c  N       Integer.  (INPUT)
c          Dimension of the eigenproblem.
c
c  K       Integer.  (INPUT)
c          Current size of V and H.
c
c  NP      Integer.  (INPUT)
c          Number of additional Arnoldi steps to take.
c
c  NB      Integer.  (INPUT)
c          Blocksize to be used in the recurrence.          
c          Only work for NB = 1 right now.  The goal is to have a 
c          program that implement both the block and non-block method.
c
c  RESID   Complex*16 array of length N.  (INPUT/OUTPUT)
c          On INPUT:  RESID contains the residual vector r_{k}.
c          On OUTPUT: RESID contains the residual vector r_{k+p}.
c
c  RNORM   Double precision scalar.  (INPUT/OUTPUT)
c          B-norm of the starting residual on input.
c          B-norm of the updated residual r_{k+p} on output.
c
c  V       Complex*16 N by K+NP array.  (INPUT/OUTPUT)
c          On INPUT:  V contains the Arnoldi vectors in the first K 
c          columns.
c          On OUTPUT: V contains the new NP Arnoldi vectors in the next
c          NP columns.  The first K columns are unchanged.
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling 
c          program.
c
c  H       Complex*16 (K+NP) by (K+NP) array.  (INPUT/OUTPUT)
c          H is used to store the generated upper Hessenberg matrix.
c
c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling 
c          program.
c
c  IPNTR   Integer array of length 3.  (OUTPUT)
c          Pointer to mark the starting locations in the WORK for 
c          vectors used by the Arnoldi iteration.
c          -------------------------------------------------------------
c          IPNTR(1): pointer to the current operand vector X.
c          IPNTR(2): pointer to the current result vector Y.
c          IPNTR(3): pointer to the vector B * X when used in the 
c                    shift-and-invert mode.  X is the current operand.
c          -------------------------------------------------------------
c          
c  WORKD   Complex*16 work array of length 3*N.  (REVERSE COMMUNICATION)
c          Distributed array to be used in the basic Arnoldi iteration
c          for reverse communication.  The calling program should not 
c          use WORKD as temporary workspace during the iteration !!!!!!
c          On input, WORKD(1:N) = B*RESID and is used to save some 
c          computation at the first step.
c
c  INFO    Integer.  (OUTPUT)
c          = 0: Normal exit.
c          > 0: Size of the spanning invariant subspace of OP found.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  Complex*16
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
c     Restarted Arnoldi Iteration", Rice University Technical Report
c     TR95-13, Department of Computational and Applied Mathematics.
c
c\Routines called:
c     zgetv0  ARPACK routine to generate the initial vector.
c     ivout   ARPACK utility routine that prints integers.
c     second  ARPACK utility routine for timing.
c     zmout   ARPACK utility routine that prints matrices
c     zvout   ARPACK utility routine that prints vectors.
c     zlanhs  LAPACK routine that computes various norms of a matrix.
c     zlascl  LAPACK routine for careful scaling of a matrix.
c     dlabad  LAPACK routine for defining the underflow and overflow
c             limits.
c     dlamch  LAPACK routine that determines machine constants.
c     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     zgemv   Level 2 BLAS routine for matrix vector multiplication.
c     zaxpy   Level 1 BLAS that computes a vector triad.
c     zcopy   Level 1 BLAS that copies one vector to another .
c     zdotc   Level 1 BLAS that computes the scalar product of two vectors. 
c     zscal   Level 1 BLAS that scales a vector.
c     zdscal  Level 1 BLAS that scales a complex vector by a real number. 
c     dznrm2  Level 1 BLAS that computes the norm of a vector.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas 
c     Applied Mathematics 
c     Rice University           
c     Houston, Texas 
c 
c\SCCS Information: @(#)
c FILE: naitr.F   SID: 2.3   DATE OF SID: 8/27/96   RELEASE: 2
c
c\Remarks
c  The algorithm implemented is:
c  
c  restart = .false.
c  Given V_{k} = [v_{1}, ..., v_{k}], r_{k}; 
c  r_{k} contains the initial residual vector even for k = 0;
c  Also assume that rnorm = || B*r_{k} || and B*r_{k} are already 
c  computed by the calling program.
c
c  betaj = rnorm ; p_{k+1} = B*r_{k} ;
c  For  j = k+1, ..., k+np  Do
c     1) if ( betaj < tol ) stop or restart depending on j.
c        ( At present tol is zero )
c        if ( restart ) generate a new starting vector.
c     2) v_{j} = r(j-1)/betaj;  V_{j} = [V_{j-1}, v_{j}];  
c        p_{j} = p_{j}/betaj
c     3) r_{j} = OP*v_{j} where OP is defined as in znaupd
c        For shift-invert mode p_{j} = B*v_{j} is already available.
c        wnorm = || OP*v_{j} ||
c     4) Compute the j-th step residual vector.
c        w_{j} =  V_{j}^T * B * OP * v_{j}
c        r_{j} =  OP*v_{j} - V_{j} * w_{j}
c        H(:,j) = w_{j};
c        H(j,j-1) = rnorm
c        rnorm = || r_(j) ||
c        If (rnorm > 0.717*wnorm) accept step and go back to 1)
c     5) Re-orthogonalization step:
c        s = V_{j}'*B*r_{j}
c        r_{j} = r_{j} - V_{j}*s;  rnorm1 = || r_{j} ||
c        alphaj = alphaj + s_{j};   
c     6) Iterative refinement step:
c        If (rnorm1 > 0.717*rnorm) then
c           rnorm = rnorm1
c           accept step and go back to 1)
c        Else
c           rnorm = rnorm1
c           If this is the first time in step 6), go to 5)
c           Else r_{j} lies in the span of V_{j} numerically.
c              Set r_{j} = 0 and rnorm = 0; go to 1)
c        EndIf 
c  End Do
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine znaitr
     &   (ido, bmat, n, k, np, nb, resid, rnorm, v, ldv, h, ldh, 
     &    ipntr, workd, info)
c
c     %----------------------------------------------------%
c     | Include files for debugging and timing information |
c     %----------------------------------------------------%
c
      include   'debug.h'
      include   'stat.h'
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character  bmat*1
      integer    ido, info, k, ldh, ldv, n, nb, np
      Double precision
     &           rnorm
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    ipntr(3)
      Complex*16
     &           h(ldh,k+np), resid(n), v(ldv,k+np), workd(3*n)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Complex*16
     &           one, zero
      Double precision
     &           rone, rzero
      parameter (one = (1.0D+0, 0.0D+0), zero = (0.0D+0, 0.0D+0), 
     &           rone = 1.0D+0, rzero = 0.0D+0)
c
c     %--------------%
c     | Local Arrays |
c     %--------------%
c
      Double precision
     &           rtemp(2)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      logical    first, orth1, orth2, rstart, step3, step4
      integer    ierr, i, infol, ipj, irj, ivj, iter, itry, j, msglvl,
     &           jj
      Double precision            
     &           ovfl, smlnum, tst1, ulp, unfl, betaj,
     &           temp1, rnorm1, wnorm
      Complex*16
     &           cnorm
c
      save       first, orth1, orth2, rstart, step3, step4,
     &           ierr, ipj, irj, ivj, iter, itry, j, msglvl, ovfl,
     &           betaj, rnorm1, smlnum, ulp, unfl, wnorm
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   zaxpy, zcopy, zscal, zdscal, zgemv, zgetv0, 
     &           dlabad, zvout, zmout, ivout, second
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Complex*16
     &           zdotc 
      Double precision            
     &           dlamch,  dznrm2, zlanhs, dlapy2
      external   zdotc, dznrm2, zlanhs, dlamch, dlapy2
c
c     %---------------------%
c     | Intrinsic Functions |
c     %---------------------%
c
      intrinsic  dimag, dble, max, sqrt 
c
c     %-----------------%
c     | Data statements |
c     %-----------------%
c
      data       first / .true. /
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      if (first) then
c
c        %-----------------------------------------%
c        | Set machine-dependent constants for the |
c        | the splitting and deflation criterion.  |
c        | If norm(H) <= sqrt(OVFL),               |
c        | overflow should not occur.              |
c        | REFERENCE: LAPACK subroutine zlahqr     |
c        %-----------------------------------------%
c
         unfl = dlamch( 'safe minimum' )
         ovfl = dble(one / unfl)
         call dlabad( unfl, ovfl )
         ulp = dlamch( 'precision' )
         smlnum = unfl*( n / ulp )
         first = .false.
      end if
c
      if (ido .eq. 0) then
c 
c        %-------------------------------%
c        | Initialize timing statistics  |
c        | & message level for debugging |
c        %-------------------------------%
c
         call second (t0)
         msglvl = mcaitr
c 
c        %------------------------------%
c        | Initial call to this routine |
c        %------------------------------%
c
         info   = 0
         step3  = .false.
         step4  = .false.
         rstart = .false.
         orth1  = .false.
         orth2  = .false.
         j      = k + 1
         ipj    = 1
         irj    = ipj   + n
         ivj    = irj   + n
      end if
c 
c     %-------------------------------------------------%
c     | When in reverse communication mode one of:      |
c     | STEP3, STEP4, ORTH1, ORTH2, RSTART              |
c     | will be .true. when ....                        |
c     | STEP3: return from computing OP*v_{j}.          |
c     | STEP4: return from computing B-norm of OP*v_{j} |
c     | ORTH1: return from computing B-norm of r_{j+1}  |
c     | ORTH2: return from computing B-norm of          |
c     |        correction to the residual vector.       |
c     | RSTART: return from OP computations needed by   |
c     |         zgetv0.                                 |
c     %-------------------------------------------------%
c
      if (step3)  go to 50
      if (step4)  go to 60
      if (orth1)  go to 70
      if (orth2)  go to 90
      if (rstart) go to 30
c
c     %-----------------------------%
c     | Else this is the first step |
c     %-----------------------------%
c
c     %--------------------------------------------------------------%
c     |                                                              |
c     |        A R N O L D I     I T E R A T I O N     L O O P       |
c     |                                                              |
c     | Note:  B*r_{j-1} is already in WORKD(1:N)=WORKD(IPJ:IPJ+N-1) |
c     %--------------------------------------------------------------%
 
 1000 continue
c
         if (msglvl .gt. 1) then
            call ivout (logfil, 1, j, ndigit, 
     &                  '_naitr: generating Arnoldi vector number')
            call dvout (logfil, 1, rnorm, ndigit, 
     &                  '_naitr: B-norm of the current residual is')
         end if
c 
c        %---------------------------------------------------%
c        | STEP 1: Check if the B norm of j-th residual      |
c        | vector is zero. Equivalent to determine whether   |
c        | an exact j-step Arnoldi factorization is present. |
c        %---------------------------------------------------%
c
         betaj = rnorm
         if (rnorm .gt. rzero) go to 40
c
c           %---------------------------------------------------%
c           | Invariant subspace found, generate a new starting |
c           | vector which is orthogonal to the current Arnoldi |
c           | basis and continue the iteration.                 |
c           %---------------------------------------------------%
c
            if (msglvl .gt. 0) then
               call ivout (logfil, 1, j, ndigit,
     &                     '_naitr: ****** RESTART AT STEP ******')
            end if
c 
c           %---------------------------------------------%
c           | ITRY is the loop variable that controls the |
c           | maximum amount of times that a restart is   |
c           | attempted. NRSTRT is used by stat.h         |
c           %---------------------------------------------%
c 
            betaj  = rzero
            nrstrt = nrstrt + 1
            itry   = 1
   20       continue
            rstart = .true.
            ido    = 0
   30       continue
c
c           %--------------------------------------%
c           | If in reverse communication mode and |
c           | RSTART = .true. flow returns here.   |
c           %--------------------------------------%
c
            call zgetv0 (ido, bmat, itry, .false., n, j, v, ldv, 
     &                   resid, rnorm, ipntr, workd, ierr)
            if (ido .ne. 99) go to 9000
            if (ierr .lt. 0) then
               itry = itry + 1
               if (itry .le. 3) go to 20
c
c              %------------------------------------------------%
c              | Give up after several restart attempts.        |
c              | Set INFO to the size of the invariant subspace |
c              | which spans OP and exit.                       |
c              %------------------------------------------------%
c
               info = j - 1
               call second (t1)
               tcaitr = tcaitr + (t1 - t0)
               ido = 99
               go to 9000
            end if
c 
   40    continue
c
c        %---------------------------------------------------------%
c        | STEP 2:  v_{j} = r_{j-1}/rnorm and p_{j} = p_{j}/rnorm  |
c        | Note that p_{j} = B*r_{j-1}. In order to avoid overflow |
c        | when reciprocating a small RNORM, test against lower    |
c        | machine bound.                                          |
c        %---------------------------------------------------------%
c
         call zcopy (n, resid, 1, v(1,j), 1)
         if ( rnorm .ge. unfl) then
             temp1 = rone / rnorm
             call zdscal (n, temp1, v(1,j), 1)
             call zdscal (n, temp1, workd(ipj), 1)
         else
c
c            %-----------------------------------------%
c            | To scale both v_{j} and p_{j} carefully |
c            | use LAPACK routine zlascl               |
c            %-----------------------------------------%
c
             call zlascl ('General', i, i, rnorm, rone,
     &                    n, 1, v(1,j), n, infol)
             call zlascl ('General', i, i, rnorm, rone,  
     &                    n, 1, workd(ipj), n, infol)
         end if
c
c        %------------------------------------------------------%
c        | STEP 3:  r_{j} = OP*v_{j}; Note that p_{j} = B*v_{j} |
c        | Note that this is not quite yet r_{j}. See STEP 4    |
c        %------------------------------------------------------%
c
         step3 = .true.
         nopx  = nopx + 1
         call second (t2)
         call zcopy (n, v(1,j), 1, workd(ivj), 1)
         ipntr(1) = ivj
         ipntr(2) = irj
         ipntr(3) = ipj
         ido = 1
c 
c        %-----------------------------------%
c        | Exit in order to compute OP*v_{j} |
c        %-----------------------------------%
c 
         go to 9000 
   50    continue
c 
c        %----------------------------------%
c        | Back from reverse communication; |
c        | WORKD(IRJ:IRJ+N-1) := OP*v_{j}   |
c        | if step3 = .true.                |
c        %----------------------------------%
c
         call second (t3)
         tmvopx = tmvopx + (t3 - t2)
 
         step3 = .false.
c
c        %------------------------------------------%
c        | Put another copy of OP*v_{j} into RESID. |
c        %------------------------------------------%
c
         call zcopy (n, workd(irj), 1, resid, 1)
c 
c        %---------------------------------------%
c        | STEP 4:  Finish extending the Arnoldi |
c        |          factorization to length j.   |
c        %---------------------------------------%
c
         call second (t2)
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            step4 = .true.
            ipntr(1) = irj
            ipntr(2) = ipj
            ido = 2
c 
c           %-------------------------------------%
c           | Exit in order to compute B*OP*v_{j} |
c           %-------------------------------------%
c 
            go to 9000
         else if (bmat .eq. 'I') then
            call zcopy (n, resid, 1, workd(ipj), 1)
         end if
   60    continue
c 
c        %----------------------------------%
c        | Back from reverse communication; |
c        | WORKD(IPJ:IPJ+N-1) := B*OP*v_{j} |
c        | if step4 = .true.                |
c        %----------------------------------%
c
         if (bmat .eq. 'G') then
            call second (t3)
            tmvbx = tmvbx + (t3 - t2)
         end if
c 
         step4 = .false.
c
c        %-------------------------------------%
c        | The following is needed for STEP 5. |
c        | Compute the B-norm of OP*v_{j}.     |
c        %-------------------------------------%
c
         if (bmat .eq. 'G') then  
             cnorm = zdotc (n, resid, 1, workd(ipj), 1)
             wnorm = sqrt( dlapy2(dble(cnorm),dimag(cnorm)) )
         else if (bmat .eq. 'I') then
             wnorm = dznrm2(n, resid, 1)
         end if
c
c        %-----------------------------------------%
c        | Compute the j-th residual corresponding |
c        | to the j step factorization.            |
c        | Use Classical Gram Schmidt and compute: |
c        | w_{j} <-  V_{j}^T * B * OP * v_{j}      |
c        | r_{j} <-  OP*v_{j} - V_{j} * w_{j}      |
c        %-----------------------------------------%
c
c
c        %------------------------------------------%
c        | Compute the j Fourier coefficients w_{j} |
c        | WORKD(IPJ:IPJ+N-1) contains B*OP*v_{j}.  |
c        %------------------------------------------%
c 
         call zgemv ('C', n, j, one, v, ldv, workd(ipj), 1,
     &               zero, h(1,j), 1)
c
c        %--------------------------------------%
c        | Orthogonalize r_{j} against V_{j}.   |
c        | RESID contains OP*v_{j}. See STEP 3. | 
c        %--------------------------------------%
c
         call zgemv ('N', n, j, -one, v, ldv, h(1,j), 1,
     &               one, resid, 1)
c
         if (j .gt. 1) h(j,j-1) = dcmplx(betaj, rzero)
c
         call second (t4)
c 
         orth1 = .true.
c 
         call second (t2)
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            call zcopy (n, resid, 1, workd(irj), 1)
            ipntr(1) = irj
            ipntr(2) = ipj
            ido = 2
c 
c           %----------------------------------%
c           | Exit in order to compute B*r_{j} |
c           %----------------------------------%
c 
            go to 9000
         else if (bmat .eq. 'I') then
            call zcopy (n, resid, 1, workd(ipj), 1)
         end if 
   70    continue
c 
c        %---------------------------------------------------%
c        | Back from reverse communication if ORTH1 = .true. |
c        | WORKD(IPJ:IPJ+N-1) := B*r_{j}.                    |
c        %---------------------------------------------------%
c
         if (bmat .eq. 'G') then
            call second (t3)
            tmvbx = tmvbx + (t3 - t2)
         end if
c 
         orth1 = .false.
c
c        %------------------------------%
c        | Compute the B-norm of r_{j}. |
c        %------------------------------%
c
         if (bmat .eq. 'G') then         
            cnorm = zdotc (n, resid, 1, workd(ipj), 1)
            rnorm = sqrt( dlapy2(dble(cnorm),dimag(cnorm)) )
         else if (bmat .eq. 'I') then
            rnorm = dznrm2(n, resid, 1)
         end if
c 
c        %-----------------------------------------------------------%
c        | STEP 5: Re-orthogonalization / Iterative refinement phase |
c        | Maximum NITER_ITREF tries.                                |
c        |                                                           |
c        |          s      = V_{j}^T * B * r_{j}                     |
c        |          r_{j}  = r_{j} - V_{j}*s                         |
c        |          alphaj = alphaj + s_{j}                          |
c        |                                                           |
c        | The stopping criteria used for iterative refinement is    |
c        | discussed in Parlett's book SEP, page 107 and in Gragg &  |
c        | Reichel ACM TOMS paper; Algorithm 686, Dec. 1990.         |
c        | Determine if we need to correct the residual. The goal is |
c        | to enforce ||v(:,1:j)^T * r_{j}|| .le. eps * || r_{j} ||  |
c        | The following test determines whether the sine of the     |
c        | angle between  OP*x and the computed residual is less     |
c        | than or equal to 0.717.                                   |
c        %-----------------------------------------------------------%
c
         if ( rnorm .gt. 0.717*wnorm ) go to 100
c
         iter  = 0
         nrorth = nrorth + 1
c 
c        %---------------------------------------------------%
c        | Enter the Iterative refinement phase. If further  |
c        | refinement is necessary, loop back here. The loop |
c        | variable is ITER. Perform a step of Classical     |
c        | Gram-Schmidt using all the Arnoldi vectors V_{j}  |
c        %---------------------------------------------------%
c 
   80    continue
c
         if (msglvl .gt. 2) then
            rtemp(1) = wnorm
            rtemp(2) = rnorm
            call dvout (logfil, 2, rtemp, ndigit, 
     &      '_naitr: re-orthogonalization; wnorm and rnorm are')
            call zvout (logfil, j, h(1,j), ndigit,
     &                  '_naitr: j-th column of H')
         end if
c
c        %----------------------------------------------------%
c        | Compute V_{j}^T * B * r_{j}.                       |
c        | WORKD(IRJ:IRJ+J-1) = v(:,1:J)'*WORKD(IPJ:IPJ+N-1). |
c        %----------------------------------------------------%
c
         call zgemv ('C', n, j, one, v, ldv, workd(ipj), 1, 
     &               zero, workd(irj), 1)
c
c        %---------------------------------------------%
c        | Compute the correction to the residual:     |
c        | r_{j} = r_{j} - V_{j} * WORKD(IRJ:IRJ+J-1). |
c        | The correction to H is v(:,1:J)*H(1:J,1:J)  |
c        | + v(:,1:J)*WORKD(IRJ:IRJ+J-1)*e'_j.         |
c        %---------------------------------------------%
c
         call zgemv ('N', n, j, -one, v, ldv, workd(irj), 1, 
     &               one, resid, 1)
         call zaxpy (j, one, workd(irj), 1, h(1,j), 1)
c 
         orth2 = .true.
         call second (t2)
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            call zcopy (n, resid, 1, workd(irj), 1)
            ipntr(1) = irj
            ipntr(2) = ipj
            ido = 2
c 
c           %-----------------------------------%
c           | Exit in order to compute B*r_{j}. |
c           | r_{j} is the corrected residual.  |
c           %-----------------------------------%
c 
            go to 9000
         else if (bmat .eq. 'I') then
            call zcopy (n, resid, 1, workd(ipj), 1)
         end if 
   90    continue
c
c        %---------------------------------------------------%
c        | Back from reverse communication if ORTH2 = .true. |
c        %---------------------------------------------------%
c
         if (bmat .eq. 'G') then
            call second (t3)
            tmvbx = tmvbx + (t3 - t2)
         end if 
c
c        %-----------------------------------------------------%
c        | Compute the B-norm of the corrected residual r_{j}. |
c        %-----------------------------------------------------%
c 
         if (bmat .eq. 'G') then         
             cnorm  = zdotc (n, resid, 1, workd(ipj), 1)
             rnorm1 = sqrt( dlapy2(dble(cnorm),dimag(cnorm)) )
         else if (bmat .eq. 'I') then
             rnorm1 = dznrm2(n, resid, 1)
         end if
c 
         if (msglvl .gt. 0 .and. iter .gt. 0 ) then
            call ivout (logfil, 1, j, ndigit,
     &           '_naitr: Iterative refinement for Arnoldi residual')
            if (msglvl .gt. 2) then
                rtemp(1) = rnorm
                rtemp(2) = rnorm1
                call dvout (logfil, 2, rtemp, ndigit,
     &           '_naitr: iterative refinement ; rnorm and rnorm1 are')
            end if
         end if
c
c        %-----------------------------------------%
c        | Determine if we need to perform another |
c        | step of re-orthogonalization.           |
c        %-----------------------------------------%
c
         if ( rnorm1 .gt. 0.717*rnorm ) then
c
c           %---------------------------------------%
c           | No need for further refinement.       |
c           | The cosine of the angle between the   |
c           | corrected residual vector and the old |
c           | residual vector is greater than 0.717 |
c           | In other words the corrected residual |
c           | and the old residual vector share an  |
c           | angle of less than arcCOS(0.717)      |
c           %---------------------------------------%
c
            rnorm = rnorm1
c 
         else
c
c           %-------------------------------------------%
c           | Another step of iterative refinement step |
c           | is required. NITREF is used by stat.h     |
c           %-------------------------------------------%
c
            nitref = nitref + 1
            rnorm  = rnorm1
            iter   = iter + 1
            if (iter .le. 1) go to 80
c
c           %-------------------------------------------------%
c           | Otherwise RESID is numerically in the span of V |
c           %-------------------------------------------------%
c
            do 95 jj = 1, n
               resid(jj) = zero
  95        continue 
            rnorm = rzero
         end if
c 
c        %----------------------------------------------%
c        | Branch here directly if iterative refinement |
c        | wasn't necessary or after at most NITER_REF  |
c        | steps of iterative refinement.               |
c        %----------------------------------------------%
c 
  100    continue
c 
         rstart = .false.
         orth2  = .false.
c 
         call second (t5)
         titref = titref + (t5 - t4)
c 
c        %------------------------------------%
c        | STEP 6: Update  j = j+1;  Continue |
c        %------------------------------------%
c
         j = j + 1
         if (j .gt. k+np) then
            call second (t1)
            tcaitr = tcaitr + (t1 - t0)
            ido = 99
            do 110 i = max(1,k), k+np-1
c     
c              %--------------------------------------------%
c              | Check for splitting and deflation.         |
c              | Use a standard test as in the QR algorithm |
c              | REFERENCE: LAPACK subroutine zlahqr        |
c              %--------------------------------------------%
c     
               tst1 = dlapy2(dble(h(i,i)),dimag(h(i,i)))
     &              + dlapy2(dble(h(i+1,i+1)), dimag(h(i+1,i+1)))
               if( tst1.eq.dble(zero) )
     &              tst1 = zlanhs( '1', k+np, h, ldh, workd(n+1) )
               if( dlapy2(dble(h(i+1,i)),dimag(h(i+1,i))) .le. 
     &                    max( ulp*tst1, smlnum ) ) 
     &             h(i+1,i) = zero
 110        continue
c     
            if (msglvl .gt. 2) then
               call zmout (logfil, k+np, k+np, h, ldh, ndigit, 
     &          '_naitr: Final upper Hessenberg matrix H of order K+NP')
            end if
c     
            go to 9000
         end if
c
c        %--------------------------------------------------------%
c        | Loop back to extend the factorization by another step. |
c        %--------------------------------------------------------%
c
      go to 1000
c 
c     %---------------------------------------------------------------%
c     |                                                               |
c     |  E N D     O F     M A I N     I T E R A T I O N     L O O P  |
c     |                                                               |
c     %---------------------------------------------------------------%
c
 9000 continue
      return
c
c     %---------------%
c     | End of znaitr |
c     %---------------%
c
      end


c\BeginDoc
c
c\Name: zneigh
c
c\Description:
c  Compute the eigenvalues of the current upper Hessenberg matrix
c  and the corresponding Ritz estimates given the current residual norm.
c
c\Usage:
c  call zneigh
c     ( RNORM, N, H, LDH, RITZ, BOUNDS, Q, LDQ, WORKL, RWORK, IERR )
c
c\Arguments
c  RNORM   Double precision scalar.  (INPUT)
c          Residual norm corresponding to the current upper Hessenberg 
c          matrix H.
c
c  N       Integer.  (INPUT)
c          Size of the matrix H.
c
c  H       Complex*16 N by N array.  (INPUT)
c          H contains the current upper Hessenberg matrix.
c
c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling
c          program.
c
c  RITZ    Complex*16 array of length N.  (OUTPUT)
c          On output, RITZ(1:N) contains the eigenvalues of H.
c
c  BOUNDS  Complex*16 array of length N.  (OUTPUT)
c          On output, BOUNDS contains the Ritz estimates associated with
c          the eigenvalues held in RITZ.  This is equal to RNORM 
c          times the last components of the eigenvectors corresponding 
c          to the eigenvalues in RITZ.
c
c  Q       Complex*16 N by N array.  (WORKSPACE)
c          Workspace needed to store the eigenvectors of H.
c
c  LDQ     Integer.  (INPUT)
c          Leading dimension of Q exactly as declared in the calling
c          program.
c
c  WORKL   Complex*16 work array of length N**2 + 3*N.  (WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.  This is needed to keep the full Schur form
c          of H and also in the calculation of the eigenvectors of H.
c
c  RWORK   Double precision  work array of length N (WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end. 
c
c  IERR    Integer.  (OUTPUT)
c          Error exit flag from zlahqr or ztrevc.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  Complex*16
c
c\Routines called:
c     ivout   ARPACK utility routine that prints integers.
c     second  ARPACK utility routine for timing.
c     zmout   ARPACK utility routine that prints matrices
c     zvout   ARPACK utility routine that prints vectors.
c     dvout   ARPACK utility routine that prints vectors.
c     zlacpy  LAPACK matrix copy routine.
c     zlahqr  LAPACK routine to compute the Schur form of an
c             upper Hessenberg matrix.
c     zlaset  LAPACK matrix initialization routine.
c     ztrevc  LAPACK routine to compute the eigenvectors of a matrix
c             in upper triangular form
c     zcopy   Level 1 BLAS that copies one vector to another. 
c     zdscal  Level 1 BLAS that scales a complex vector by a real number.
c     dznrm2  Level 1 BLAS that computes the norm of a vector.
c     
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics 
c     Rice University           
c     Houston, Texas 
c
c\SCCS Information: @(#)
c FILE: neigh.F   SID: 2.2   DATE OF SID: 4/20/96   RELEASE: 2
c
c\Remarks
c     None
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine zneigh (rnorm, n, h, ldh, ritz, bounds, 
     &                   q, ldq, workl, rwork, ierr)
c
c     %----------------------------------------------------%
c     | Include files for debugging and timing information |
c     %----------------------------------------------------%
c
      include   'debug.h'
      include   'stat.h'
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      integer    ierr, n, ldh, ldq
      Double precision     
     &           rnorm
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Complex*16     
     &           bounds(n), h(ldh,n), q(ldq,n), ritz(n),
     &           workl(n*(n+3)) 
      Double precision 
     &           rwork(n)
c 
c     %------------%
c     | Parameters |
c     %------------%
c
      Complex*16     
     &           one, zero
      Double precision
     &           rone
      parameter  (one = (1.0D+0, 0.0D+0), zero = (0.0D+0, 0.0D+0),
     &           rone = 1.0D+0)
c 
c     %------------------------%
c     | Local Scalars & Arrays |
c     %------------------------%
c
      logical    select(1)
      integer    j,  msglvl
      Complex*16     
     &           vl(1)
      Double precision
     &           temp
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   zlacpy, zlahqr, ztrevc, zcopy, 
     &           zdscal, zmout, zvout, second
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision 
     &           dznrm2
      external   dznrm2
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
c     %-------------------------------%
c     | Initialize timing statistics  |
c     | & message level for debugging |
c     %-------------------------------%
c
      call second (t0)
      msglvl = mceigh
c 
      if (msglvl .gt. 2) then
          call zmout (logfil, n, n, h, ldh, ndigit, 
     &         '_neigh: Entering upper Hessenberg matrix H ')
      end if
c 
c     %----------------------------------------------------------%
c     | 1. Compute the eigenvalues, the last components of the   |
c     |    corresponding Schur vectors and the full Schur form T |
c     |    of the current upper Hessenberg matrix H.             |
c     |    zlahqr returns the full Schur form of H               | 
c     |    in WORKL(1:N**2), and the Schur vectors in q.         |
c     %----------------------------------------------------------%
c
      call zlacpy ('All', n, n, h, ldh, workl, n)
      call zlaset ('All', n, n, zero, one, q, ldq)
      call zlahqr (.true., .true., n, 1, n, workl, ldh, ritz,
     &             1, n, q, ldq, ierr)
      if (ierr .ne. 0) go to 9000
c
      call zcopy (n, q(n-1,1), ldq, bounds, 1)
      if (msglvl .gt. 1) then
         call zvout (logfil, n, bounds, ndigit,
     &              '_neigh: last row of the Schur matrix for H')
      end if
c
c     %----------------------------------------------------------%
c     | 2. Compute the eigenvectors of the full Schur form T and |
c     |    apply the Schur vectors to get the corresponding      |
c     |    eigenvectors.                                         |
c     %----------------------------------------------------------%
c
      call ztrevc ('Right', 'Back', select, n, workl, n, vl, n, q, 
     &             ldq, n, n, workl(n*n+1), rwork, ierr)
c
      if (ierr .ne. 0) go to 9000
c
c     %------------------------------------------------%
c     | Scale the returning eigenvectors so that their |
c     | Euclidean norms are all one. LAPACK subroutine |
c     | ztrevc returns each eigenvector normalized so  |
c     | that the element of largest magnitude has      |
c     | magnitude 1; here the magnitude of a complex   |
c     | number (x,y) is taken to be |x| + |y|.         |
c     %------------------------------------------------%
c
      do 10 j=1, n
            temp = dznrm2( n, q(1,j), 1 )
            call zdscal ( n, rone / temp, q(1,j), 1 )
   10 continue
c
      if (msglvl .gt. 1) then
         call zcopy(n, q(n,1), ldq, workl, 1)
         call zvout (logfil, n, workl, ndigit,
     &              '_neigh: Last row of the eigenvector matrix for H')
      end if
c
c     %----------------------------%
c     | Compute the Ritz estimates |
c     %----------------------------%
c
      call zcopy(n, q(n,1), n, bounds, 1)
      call zdscal(n, rnorm, bounds, 1)
c
      if (msglvl .gt. 2) then
         call zvout (logfil, n, ritz, ndigit,
     &              '_neigh: The eigenvalues of H')
         call zvout (logfil, n, bounds, ndigit,
     &              '_neigh: Ritz estimates for the eigenvalues of H')
      end if
c
      call second(t1)
      tceigh = tceigh + (t1 - t0)
c
 9000 continue
      return
c
c     %---------------%
c     | End of zneigh |
c     %---------------%
c
      end



c\BeginDoc
c
c\Name: zsortc
c
c\Description:
c  Sorts the Complex*16 array in X into the order 
c  specified by WHICH and optionally applies the permutation to the
c  Double precision  array Y. 
c
c\Usage:
c  call zsortc
c     ( WHICH, APPLY, N, X, Y )
c
c\Arguments
c  WHICH   Character*2.  (Input)
c          'LM' -> sort X into increasing order of magnitude.
c          'SM' -> sort X into decreasing order of magnitude.
c          'LR' -> sort X with real(X) in increasing algebraic order 
c          'SR' -> sort X with real(X) in decreasing algebraic order
c          'LI' -> sort X with imag(X) in increasing algebraic order
c          'SI' -> sort X with imag(X) in decreasing algebraic order
c
c  APPLY   Logical.  (Input)
c          APPLY = .TRUE.  -> apply the sorted order to array Y.
c          APPLY = .FALSE. -> do not apply the sorted order to array Y.
c
c  N       Integer.  (INPUT)
c          Size of the arrays.
c
c  X       Complex*16 array of length N.  (INPUT/OUTPUT)
c          This is the array to be sorted.
c
c  Y       Complex*16 array of length N.  (INPUT/OUTPUT)
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Routines called:
c     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics 
c     Rice University           
c     Houston, Texas 
c
c     Adapted from the sort routine in LANSO.
c
c\SCCS Information: @(#)
c FILE: sortc.F   SID: 2.2   DATE OF SID: 4/20/96   RELEASE: 2
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine zsortc (which, apply, n, x, y)
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character*2 which
      logical    apply
      integer    n
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Complex*16     
     &           x(0:n-1), y(0:n-1)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer    i, igap, j
      Complex*16     
     &           temp
      Double precision 
     &           temp1, temp2
c
c     %--------------------%
c     | External functions |
c     %--------------------%
c
      Double precision
     &           dlapy2
c
c     %--------------------%
c     | Intrinsic Functions |
c     %--------------------%
       Intrinsic
     &           dble, dimag
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      igap = n / 2
c 
      if (which .eq. 'LM') then
c
c        %--------------------------------------------%
c        | Sort X into increasing order of magnitude. |
c        %--------------------------------------------%
c
   10    continue
         if (igap .eq. 0) go to 9000
c
         do 30 i = igap, n-1
            j = i-igap
   20       continue
c
            if (j.lt.0) go to 30
c
            temp1 = dlapy2(dble(x(j)),dimag(x(j)))
            temp2 = dlapy2(dble(x(j+igap)),dimag(x(j+igap)))
c
            if (temp1.gt.temp2) then
                temp = x(j)
                x(j) = x(j+igap)
                x(j+igap) = temp
c
                if (apply) then
                    temp = y(j)
                    y(j) = y(j+igap)
                    y(j+igap) = temp
                end if
            else
                go to 30
            end if
            j = j-igap
            go to 20
   30    continue
         igap = igap / 2
         go to 10
c
      else if (which .eq. 'SM') then
c
c        %--------------------------------------------%
c        | Sort X into decreasing order of magnitude. |
c        %--------------------------------------------%
c
   40    continue
         if (igap .eq. 0) go to 9000
c
         do 60 i = igap, n-1
            j = i-igap
   50       continue
c
            if (j .lt. 0) go to 60
c
            temp1 = dlapy2(dble(x(j)),dimag(x(j)))
            temp2 = dlapy2(dble(x(j+igap)),dimag(x(j+igap)))
c
            if (temp1.lt.temp2) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
c 
               if (apply) then
                  temp = y(j)
                  y(j) = y(j+igap)
                  y(j+igap) = temp
               end if
            else
               go to 60
            endif
            j = j-igap
            go to 50
   60    continue
         igap = igap / 2
         go to 40
c 
      else if (which .eq. 'LR') then
c
c        %------------------------------------------------%
c        | Sort XREAL into increasing order of algebraic. |
c        %------------------------------------------------%
c
   70    continue
         if (igap .eq. 0) go to 9000
c
         do 90 i = igap, n-1
            j = i-igap
   80       continue
c
            if (j.lt.0) go to 90
c
            if (dble(x(j)).gt.dble(x(j+igap))) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
c 
               if (apply) then
                  temp = y(j)
                  y(j) = y(j+igap)
                  y(j+igap) = temp
               end if
            else
               go to 90
            endif
            j = j-igap
            go to 80
   90    continue
         igap = igap / 2
         go to 70
c 
      else if (which .eq. 'SR') then
c
c        %------------------------------------------------%
c        | Sort XREAL into decreasing order of algebraic. |
c        %------------------------------------------------%
c
  100    continue
         if (igap .eq. 0) go to 9000
         do 120 i = igap, n-1
            j = i-igap
  110       continue
c
            if (j.lt.0) go to 120
c
            if (dble(x(j)).lt.dble(x(j+igap))) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
c 
               if (apply) then
                  temp = y(j)
                  y(j) = y(j+igap)
                  y(j+igap) = temp
               end if
            else
               go to 120
            endif
            j = j-igap
            go to 110
  120    continue
         igap = igap / 2
         go to 100
c 
      else if (which .eq. 'LI') then
c
c        %--------------------------------------------%
c        | Sort XIMAG into increasing algebraic order |
c        %--------------------------------------------%
c
  130    continue
         if (igap .eq. 0) go to 9000
         do 150 i = igap, n-1
            j = i-igap
  140       continue
c
            if (j.lt.0) go to 150
c
            if (dimag(x(j)).gt.dimag(x(j+igap))) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
c
               if (apply) then
                  temp = y(j)
                  y(j) = y(j+igap)
                  y(j+igap) = temp
               end if
            else
               go to 150
            endif
            j = j-igap
            go to 140
  150    continue
         igap = igap / 2
         go to 130
c 
      else if (which .eq. 'SI') then
c
c        %---------------------------------------------%
c        | Sort XIMAG into decreasing algebraic order  |
c        %---------------------------------------------%
c
  160    continue
         if (igap .eq. 0) go to 9000
         do 180 i = igap, n-1
            j = i-igap
  170       continue
c
            if (j.lt.0) go to 180
c
            if (dimag(x(j)).lt.dimag(x(j+igap))) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
c 
               if (apply) then
                  temp = y(j)
                  y(j) = y(j+igap)
                  y(j+igap) = temp
               end if
            else
               go to 180
            endif
            j = j-igap
            go to 170
  180    continue
         igap = igap / 2
         go to 160
      end if
c 
 9000 continue
      return
c
c     %---------------%
c     | End of zsortc |
c     %---------------%
c
      end



c\BeginDoc
c
c\Name: znapps
c
c\Description:
c  Given the Arnoldi factorization
c
c     A*V_{k} - V_{k}*H_{k} = r_{k+p}*e_{k+p}^T,
c
c  apply NP implicit shifts resulting in
c
c     A*(V_{k}*Q) - (V_{k}*Q)*(Q^T* H_{k}*Q) = r_{k+p}*e_{k+p}^T * Q
c
c  where Q is an orthogonal matrix which is the product of rotations
c  and reflections resulting from the NP bulge change sweeps.
c  The updated Arnoldi factorization becomes:
c
c     A*VNEW_{k} - VNEW_{k}*HNEW_{k} = rnew_{k}*e_{k}^T.
c
c\Usage:
c  call znapps
c     ( N, KEV, NP, SHIFT, V, LDV, H, LDH, RESID, Q, LDQ, 
c       WORKL, WORKD )
c
c\Arguments
c  N       Integer.  (INPUT)
c          Problem size, i.e. size of matrix A.
c
c  KEV     Integer.  (INPUT/OUTPUT)
c          KEV+NP is the size of the input matrix H.
c          KEV is the size of the updated matrix HNEW. 
c
c  NP      Integer.  (INPUT)
c          Number of implicit shifts to be applied.
c
c  SHIFT   Complex*16 array of length NP.  (INPUT)
c          The shifts to be applied.
c
c  V       Complex*16 N by (KEV+NP) array.  (INPUT/OUTPUT)
c          On INPUT, V contains the current KEV+NP Arnoldi vectors.
c          On OUTPUT, V contains the updated KEV Arnoldi vectors
c          in the first KEV columns of V.
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling
c          program.
c
c  H       Complex*16 (KEV+NP) by (KEV+NP) array.  (INPUT/OUTPUT)
c          On INPUT, H contains the current KEV+NP by KEV+NP upper 
c          Hessenberg matrix of the Arnoldi factorization.
c          On OUTPUT, H contains the updated KEV by KEV upper Hessenberg
c          matrix in the KEV leading submatrix.
c
c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling
c          program.
c
c  RESID   Complex*16 array of length N.  (INPUT/OUTPUT)
c          On INPUT, RESID contains the the residual vector r_{k+p}.
c          On OUTPUT, RESID is the update residual vector rnew_{k} 
c          in the first KEV locations.
c
c  Q       Complex*16 KEV+NP by KEV+NP work array.  (WORKSPACE)
c          Work array used to accumulate the rotations and reflections
c          during the bulge chase sweep.
c
c  LDQ     Integer.  (INPUT)
c          Leading dimension of Q exactly as declared in the calling
c          program.
c
c  WORKL   Complex*16 work array of length (KEV+NP).  (WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.
c
c  WORKD   Complex*16 work array of length 2*N.  (WORKSPACE)
c          Distributed array used in the application of the accumulated
c          orthogonal matrix Q.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  Complex*16
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c
c\Routines called:
c     ivout   ARPACK utility routine that prints integers.
c     second  ARPACK utility routine for timing.
c     zmout   ARPACK utility routine that prints matrices
c     zvout   ARPACK utility routine that prints vectors.
c     zlacpy  LAPACK matrix copy routine.
c     zlanhs  LAPACK routine that computes various norms of a matrix.
c     zlartg  LAPACK Givens rotation construction routine.
c     zlaset  LAPACK matrix initialization routine.
c     dlabad  LAPACK routine for defining the underflow and overflow
c             limits.
c     dlamch  LAPACK routine that determines machine constants.
c     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     zgemv   Level 2 BLAS routine for matrix vector multiplication.
c     zaxpy   Level 1 BLAS that computes a vector triad.
c     zcopy   Level 1 BLAS that copies one vector to another.
c     zscal   Level 1 BLAS that scales a vector.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics 
c     Rice University           
c     Houston, Texas 
c
c\SCCS Information: @(#)
c FILE: napps.F   SID: 2.3   DATE OF SID: 3/28/97   RELEASE: 2
c
c\Remarks
c  1. In this version, each shift is applied to all the sublocks of
c     the Hessenberg matrix H and not just to the submatrix that it
c     comes from. Deflation as in LAPACK routine zlahqr (QR algorithm
c     for upper Hessenberg matrices ) is used.
c     Upon output, the subdiagonals of H are enforced to be non-negative
c     real numbers.
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine znapps
     &   ( n, kev, np, shift, v, ldv, h, ldh, resid, q, ldq, 
     &     workl, workd )
c
c     %----------------------------------------------------%
c     | Include files for debugging and timing information |
c     %----------------------------------------------------%
c
      include   'debug.h'
      include   'stat.h'
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      integer    kev, ldh, ldq, ldv, n, np
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Complex*16
     &           h(ldh,kev+np), resid(n), shift(np), 
     &           v(ldv,kev+np), q(ldq,kev+np), workd(2*n), workl(kev+np)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Complex*16
     &           one, zero
      Double precision
     &           rzero
      parameter (one = (1.0D+0, 0.0D+0), zero = (0.0D+0, 0.0D+0),
     &           rzero = 0.0D+0)
c
c     %------------------------%
c     | Local Scalars & Arrays |
c     %------------------------%
c
      integer    i, iend, istart, j, jj, kplusp, msglvl
      logical    first
      Complex*16
     &           cdum, f, g, h11, h21, r, s, sigma, t
      Double precision             
     &           c,  ovfl, smlnum, ulp, unfl, tst1
      save       first, ovfl, smlnum, ulp, unfl 
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   zaxpy, zcopy, zgemv, zscal, zlacpy, zlartg, 
     &           zvout, zlaset, dlabad, zmout, second, ivout
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision                 
     &           zlanhs, dlamch, dlapy2
      external   zlanhs, dlamch, dlapy2
c
c     %----------------------%
c     | Intrinsics Functions |
c     %----------------------%
c
      intrinsic  abs, dimag, conjg, dcmplx, max, min, dble
c
c     %---------------------%
c     | Statement Functions |
c     %---------------------%
c
      Double precision     
     &           zabs1
      zabs1( cdum ) = abs( dble( cdum ) ) + abs( dimag( cdum ) )
c
c     %----------------%
c     | Data statments |
c     %----------------%
c
      data       first / .true. /
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      if (first) then
c
c        %-----------------------------------------------%
c        | Set machine-dependent constants for the       |
c        | stopping criterion. If norm(H) <= sqrt(OVFL), |
c        | overflow should not occur.                    |
c        | REFERENCE: LAPACK subroutine zlahqr           |
c        %-----------------------------------------------%
c
         unfl = dlamch( 'safe minimum' )
         ovfl = dble(one / unfl)
         call dlabad( unfl, ovfl )
         ulp = dlamch( 'precision' )
         smlnum = unfl*( n / ulp )
         first = .false.
      end if
c
c     %-------------------------------%
c     | Initialize timing statistics  |
c     | & message level for debugging |
c     %-------------------------------%
c
      call second (t0)
      msglvl = mcapps
c 
      kplusp = kev + np 
c 
c     %--------------------------------------------%
c     | Initialize Q to the identity to accumulate |
c     | the rotations and reflections              |
c     %--------------------------------------------%
c
      call zlaset ('All', kplusp, kplusp, zero, one, q, ldq)
c
c     %----------------------------------------------%
c     | Quick return if there are no shifts to apply |
c     %----------------------------------------------%
c
      if (np .eq. 0) go to 9000
c
c     %----------------------------------------------%
c     | Chase the bulge with the application of each |
c     | implicit shift. Each shift is applied to the |
c     | whole matrix including each block.           |
c     %----------------------------------------------%
c
      do 110 jj = 1, np
         sigma = shift(jj)
c
         if (msglvl .gt. 2 ) then
            call ivout (logfil, 1, jj, ndigit, 
     &               '_napps: shift number.')
            call zvout (logfil, 1, sigma, ndigit, 
     &               '_napps: Value of the shift ')
         end if
c
         istart = 1
   20    continue
c
         do 30 i = istart, kplusp-1
c
c           %----------------------------------------%
c           | Check for splitting and deflation. Use |
c           | a standard test as in the QR algorithm |
c           | REFERENCE: LAPACK subroutine zlahqr    |
c           %----------------------------------------%
c
            tst1 = zabs1( h( i, i ) ) + zabs1( h( i+1, i+1 ) )
            if( tst1.eq.rzero )
     &         tst1 = zlanhs( '1', kplusp-jj+1, h, ldh, workl )
            if ( abs(dble(h(i+1,i))) 
     &           .le. max(ulp*tst1, smlnum) )  then
               if (msglvl .gt. 0) then
                  call ivout (logfil, 1, i, ndigit, 
     &                 '_napps: matrix splitting at row/column no.')
                  call ivout (logfil, 1, jj, ndigit, 
     &                 '_napps: matrix splitting with shift number.')
                  call zvout (logfil, 1, h(i+1,i), ndigit, 
     &                 '_napps: off diagonal element.')
               end if
               iend = i
               h(i+1,i) = zero
               go to 40
            end if
   30    continue
         iend = kplusp
   40    continue
c
         if (msglvl .gt. 2) then
             call ivout (logfil, 1, istart, ndigit, 
     &                   '_napps: Start of current block ')
             call ivout (logfil, 1, iend, ndigit, 
     &                   '_napps: End of current block ')
         end if
c
c        %------------------------------------------------%
c        | No reason to apply a shift to block of order 1 |
c        | or if the current block starts after the point |
c        | of compression since we'll discard this stuff  |
c        %------------------------------------------------%
c
         if ( istart .eq. iend .or. istart .gt. kev) go to 100
c
         h11 = h(istart,istart)
         h21 = h(istart+1,istart)
         f = h11 - sigma
         g = h21
c 
         do 80 i = istart, iend-1
c
c           %------------------------------------------------------%
c           | Construct the plane rotation G to zero out the bulge |
c           %------------------------------------------------------%
c
            call zlartg (f, g, c, s, r)
            if (i .gt. istart) then
               h(i,i-1) = r
               h(i+1,i-1) = zero
            end if
c
c           %---------------------------------------------%
c           | Apply rotation to the left of H;  H <- G'*H |
c           %---------------------------------------------%
c
            do 50 j = i, kplusp
               t        =  c*h(i,j) + s*h(i+1,j)
               h(i+1,j) = -conjg(s)*h(i,j) + c*h(i+1,j)
               h(i,j)   = t   
   50       continue
c
c           %---------------------------------------------%
c           | Apply rotation to the right of H;  H <- H*G |
c           %---------------------------------------------%
c
            do 60 j = 1, min(i+2,iend)
               t        =  c*h(j,i) + conjg(s)*h(j,i+1)
               h(j,i+1) = -s*h(j,i) + c*h(j,i+1)
               h(j,i)   = t   
   60       continue
c
c           %-----------------------------------------------------%
c           | Accumulate the rotation in the matrix Q;  Q <- Q*G' |
c           %-----------------------------------------------------%
c
            do 70 j = 1, min(i+jj, kplusp)
               t        =   c*q(j,i) + conjg(s)*q(j,i+1)
               q(j,i+1) = - s*q(j,i) + c*q(j,i+1)
               q(j,i)   = t   
   70       continue
c
c           %---------------------------%
c           | Prepare for next rotation |
c           %---------------------------%
c
            if (i .lt. iend-1) then
               f = h(i+1,i)
               g = h(i+2,i)
            end if
   80    continue
c
c        %-------------------------------%
c        | Finished applying the shift.  |
c        %-------------------------------%
c 
  100    continue
c
c        %---------------------------------------------------------%
c        | Apply the same shift to the next block if there is any. |
c        %---------------------------------------------------------%
c
         istart = iend + 1
         if (iend .lt. kplusp) go to 20
c
c        %---------------------------------------------%
c        | Loop back to the top to get the next shift. |
c        %---------------------------------------------%
c
  110 continue
c
c     %---------------------------------------------------%
c     | Perform a similarity transformation that makes    |
c     | sure that the compressed H will have non-negative |
c     | real subdiagonal elements.                        |
c     %---------------------------------------------------%
c
      do 120 j=1,kev
         if ( dble( h(j+1,j) ) .lt. rzero .or.
     &        dimag( h(j+1,j) ) .ne. rzero ) then
            t = h(j+1,j) / dlapy2(dble(h(j+1,j)),dimag(h(j+1,j)))
            call zscal( kplusp-j+1, conjg(t), h(j+1,j), ldh )
            call zscal( min(j+2, kplusp), t, h(1,j+1), 1 )
            call zscal( min(j+np+1,kplusp), t, q(1,j+1), 1 )
            h(j+1,j) = dcmplx( dble( h(j+1,j) ), rzero )
         end if
  120 continue
c
      do 130 i = 1, kev
c
c        %--------------------------------------------%
c        | Final check for splitting and deflation.   |
c        | Use a standard test as in the QR algorithm |
c        | REFERENCE: LAPACK subroutine zlahqr.       |
c        | Note: Since the subdiagonals of the        |
c        | compressed H are nonnegative real numbers, |
c        | we take advantage of this.                 |
c        %--------------------------------------------%
c
         tst1 = zabs1( h( i, i ) ) + zabs1( h( i+1, i+1 ) )
         if( tst1 .eq. rzero )
     &       tst1 = zlanhs( '1', kev, h, ldh, workl )
         if( dble( h( i+1,i ) ) .le. max( ulp*tst1, smlnum ) ) 
     &       h(i+1,i) = zero
 130  continue
c
c     %-------------------------------------------------%
c     | Compute the (kev+1)-st column of (V*Q) and      |
c     | temporarily store the result in WORKD(N+1:2*N). |
c     | This is needed in the residual update since we  |
c     | cannot GUARANTEE that the corresponding entry   |
c     | of H would be zero as in exact arithmetic.      |
c     %-------------------------------------------------%
c
      if ( dble( h(kev+1,kev) ) .gt. rzero )
     &   call zgemv ('N', n, kplusp, one, v, ldv, q(1,kev+1), 1, zero, 
     &                workd(n+1), 1)
c 
c     %----------------------------------------------------------%
c     | Compute column 1 to kev of (V*Q) in backward order       |
c     | taking advantage of the upper Hessenberg structure of Q. |
c     %----------------------------------------------------------%
c
      do 140 i = 1, kev
         call zgemv ('N', n, kplusp-i+1, one, v, ldv,
     &               q(1,kev-i+1), 1, zero, workd, 1)
         call zcopy (n, workd, 1, v(1,kplusp-i+1), 1)
  140 continue
c
c     %-------------------------------------------------%
c     |  Move v(:,kplusp-kev+1:kplusp) into v(:,1:kev). |
c     %-------------------------------------------------%
c
      call zlacpy ('A', n, kev, v(1,kplusp-kev+1), ldv, v, ldv)
c 
c     %--------------------------------------------------------------%
c     | Copy the (kev+1)-st column of (V*Q) in the appropriate place |
c     %--------------------------------------------------------------%
c
      if ( dble( h(kev+1,kev) ) .gt. rzero )
     &   call zcopy (n, workd(n+1), 1, v(1,kev+1), 1)
c 
c     %-------------------------------------%
c     | Update the residual vector:         |
c     |    r <- sigmak*r + betak*v(:,kev+1) |
c     | where                               |
c     |    sigmak = (e_{kev+p}'*Q)*e_{kev}  |
c     |    betak = e_{kev+1}'*H*e_{kev}     |
c     %-------------------------------------%
c
      call zscal (n, q(kplusp,kev), resid, 1)
      if ( dble( h(kev+1,kev) ) .gt. rzero )
     &   call zaxpy (n, h(kev+1,kev), v(1,kev+1), 1, resid, 1)
c
      if (msglvl .gt. 1) then
         call zvout (logfil, 1, q(kplusp,kev), ndigit,
     &        '_napps: sigmak = (e_{kev+p}^T*Q)*e_{kev}')
         call zvout (logfil, 1, h(kev+1,kev), ndigit,
     &        '_napps: betak = e_{kev+1}^T*H*e_{kev}')
         call ivout (logfil, 1, kev, ndigit, 
     &               '_napps: Order of the final Hessenberg matrix ')
         if (msglvl .gt. 2) then
            call zmout (logfil, kev, kev, h, ldh, ndigit,
     &      '_napps: updated Hessenberg matrix H for next iteration')
         end if
c
      end if
c
 9000 continue
      call second (t1)
      tcapps = tcapps + (t1 - t0)
c 
      return
c
c     %---------------%
c     | End of znapps |
c     %---------------%
c
      end
