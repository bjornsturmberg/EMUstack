C     Calculate absorption (A = 1 - t - r) and write Transmission,
C          Refelction and Absorption for each Lambda to file
C                 
      subroutine A_and_W_Lambda_sub(TLambda, RLambda, neq_PW,  
     * numberprop_S, numberprop_S_b, lambda, freq, pol,
     * Zeroth_Order_inv, debug, d_in_nm, 
     * incident, what4incident, out4incident, Checks, h)
C
      implicit none
C
      integer*8 neq_PW, debug, d_in_nm
      complex*16 TLambda(2*neq_PW,2*neq_PW) 
      complex*16 RLambda(2*neq_PW,2*neq_PW)
      complex*16 Lambda_t, Lambda_r, Lambda_a, tot_0
      integer*8 incident, what4incident, inc
      integer*8 i, j, out4incident
      integer*8 numberprop_S, numberprop_S_b
      integer*8 pol, Zeroth_Order_inv
      double precision lambda, freq, lambda_nm, freq_nm, h
      integer*8 Checks
C
CCCCCCCCCCC Start Program CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      if (incident .eq. 0) then
        inc = Zeroth_Order_inv
      else
        inc = incident
      endif
C
      if (debug .eq. 1) then
        write(*,*) "pol", pol,"incident = ", inc
      endif
C
C  Calculate the desired flavour of total transmission and
C   reflection as specified in the parameter/output file
      tot_0 = 1.0d0
      Lambda_t = 0.0d0
      Lambda_r = 0.0d0
      lambda_nm = lambda*d_in_nm
      freq_nm = freq/d_in_nm
C
C CHECK IF TLambda, RLambda output ScatMat and input A_and_W_Lambda_sub diff!?!?
C
      if (Checks .eq. 2) then
      open (unit=341, file='Matrices/TLambda-AW-2.txt',
     *         status='unknown')
      open (unit=342, file='Matrices/RLambda-AW-2.txt',
     *         status='unknown')
      do j=1,2*neq_PW
        do i=1,2*neq_PW
          write(341,102) i, j, TLambda(i,j), abs(TLambda(i,j))**2
          write(342,102) i, j, RLambda(i,j), abs(RLambda(i,j))**2
        enddo
      enddo       
      close(341)      
      close(342)
      endif

102   format(I4,I4,3(G25.17))
C
C
      if (what4incident .eq. 1) then
        if (pol .eq. 1) then
          do i=1,numberprop_S
            do j=1,numberprop_S
C  For TE polarisation:
            Lambda_t = Lambda_t + ABS(TLambda(i,j))**2
     *          + ABS(TLambda(neq_PW+i,j))**2
            Lambda_r = Lambda_r + ABS(RLambda(i,j))**2
     *           + ABS(RLambda(neq_PW+i,j))**2
            enddo
          enddo
        endif
C
C  For TM polarisation:
        if (pol .eq. 2) then
      do i=1,numberprop_S
        do j=1,numberprop_S
            Lambda_t = Lambda_t + ABS(TLambda(i,neq_PW+j))**2
     *           + ABS(TLambda(neq_PW+i,neq_PW+j))**2
            Lambda_r = Lambda_r + ABS(RLambda(i,neq_PW+j))**2
     *           + ABS(RLambda(neq_PW+i,neq_PW+j))**2
        enddo
      enddo
        endif
      Lambda_a = tot_0*numberprop_S
     *                - Lambda_t - Lambda_r
      write(643,101) lambda_nm, freq_nm, REAL(Lambda_t), h*d_in_nm
      write(644,101) lambda_nm, freq_nm, REAL(Lambda_r), h*d_in_nm
      write(645,101) lambda_nm, freq_nm, REAL(Lambda_a), h*d_in_nm
C
C
C
      elseif (what4incident .eq. 2) then
        if (pol .eq. 1) then
C  For TE polarisation:
      do i=1,numberprop_S_b
          Lambda_t = Lambda_t + ABS(TLambda(i,inc))**2
     *         + ABS(TLambda(neq_PW+i,inc))**2
      enddo
      do i=1,numberprop_S
          Lambda_r = Lambda_r + ABS(RLambda(i,inc))**2
     *         + ABS(RLambda(neq_PW+i,inc))**2
      enddo

        endif
C  For TM polarisation:
        if (pol .eq. 2) then
C
        if (Checks .eq. 2) then
        open (unit=341, file='Normed/AW_Lambda_t-2.txt',
     *           status='unknown')
        open (unit=342, file='Normed/AW_Lambda_r-2.txt',
     *           status='unknown')
        endif
C
      do i=1,numberprop_S_b
          Lambda_t = Lambda_t + ABS(TLambda(i,neq_PW+inc))**2
     *         + ABS(TLambda(neq_PW+i,neq_PW+inc))**2
        if (Checks .eq. 2) then
            write(341,*) i, neq_PW, inc, incident, 
     *     TLambda(i,neq_PW+inc), TLambda(neq_PW+i,neq_PW+inc),
     *     Lambda_t
        endif
      enddo
       do i=1,numberprop_S
          Lambda_r = Lambda_r + ABS(RLambda(i,neq_PW+inc))**2
     *         + ABS(RLambda(neq_PW+i,neq_PW+inc))**2
         if (Checks .eq. 2) then
            write(342,*) i, neq_PW, inc, incident, 
     *     RLambda(i,neq_PW+inc), RLambda(neq_PW+i,neq_PW+inc),
     *     Lambda_r
         endif
       enddo
        endif

        Lambda_a = tot_0 - Lambda_t - Lambda_r
C
      if (Checks .eq. 2) then
      write(341,*) tot_0
      write(341,*) pol
      write(341,*) Lambda_a
      write(342,*) tot_0
      write(342,*) pol
      write(342,*) Lambda_a
      close(341)      
      close(342)
      endif
C
      write(643,101) lambda_nm, freq_nm, REAL(Lambda_t), h*d_in_nm
      write(644,101) lambda_nm, freq_nm, REAL(Lambda_r), h*d_in_nm
      write(645,101) lambda_nm, freq_nm, REAL(Lambda_a), h*d_in_nm    
C
C
C
      else
C  For TE polarisation:
        if (pol .eq. 1) then
          Lambda_t = Lambda_t + ABS(TLambda(out4incident,inc))**2
     *         + ABS(TLambda(neq_PW+out4incident,inc))**2
          Lambda_r = Lambda_r + ABS(RLambda(out4incident,inc))**2
     *         + ABS(RLambda(neq_PW+out4incident,inc))**2
        endif
C  For TM polarisation:
        if (pol .eq. 2) then
          Lambda_t = Lambda_t 
     *      + ABS(TLambda(out4incident,neq_PW+inc))**2
     *         + ABS(TLambda(neq_PW+out4incident,neq_PW+inc))**2
          Lambda_r = Lambda_r 
     *      + ABS(RLambda(out4incident,neq_PW+inc))**2
     *         + ABS(RLambda(neq_PW+out4incident,neq_PW+inc))**2
        endif
        Lambda_a = tot_0 - Lambda_t - Lambda_r
      write(643,101) lambda_nm, freq_nm, REAL(Lambda_t), h*d_in_nm
      write(644,101) lambda_nm, freq_nm, REAL(Lambda_r), h*d_in_nm
      write(645,101) lambda_nm, freq_nm, REAL(Lambda_a), h*d_in_nm
      endif
C
101    format(2(f18.10),(g25.17),(f18.10))
C
      return
      end 
