C  Calculate the Projection integral of the conjugate Plane Waves & Bloch Modes
C
      subroutine J_overlap_sub( lat_vecs, lambda, freq, n_eff_0,
     *  eps_eff_sub, neq_PW, bloch_vec, X_mat, numberprop_S,
     *  index_pw_inv, PrintAll, ordre_ls, k_0)
C
      implicit none 
C  input output parameters
      double precision lat_vecs(2,2), n_eff_0, eps_eff_sub
      double precision lambda, freq
      integer*8 neq_PW, PrintAll, ordre_ls
      double precision bloch_vec(2), k_0, d
      integer*8 numberprop_S
      integer*8 index_pw_inv(neq_PW)
C  local parameters - purely internal
      integer*8 i, j, s, s2
      integer*8 px, py, PW_max
      parameter (PW_max = 200)
      complex*16 X_mat(2*neq_PW,2*neq_PW)  
      complex*16 beta_z_pw(PW_max), test
      complex*16 z_tmp, X_tmp
      double precision vec_kx, vec_ky, k_1
      double precision bloch1, bloch2, pi, alpha, beta, norm
C
CCCCCCCCCCC Start Program CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      d = lat_vecs(1,1)  
      pi = 3.141592653589793d0
      bloch1 = bloch_vec(1)
      bloch2 = bloch_vec(2)
      vec_kx = 2.0d0*pi/d
      vec_ky = 2.0d0*pi/d
      k_1 = 2.0d0*pi*n_eff_0*freq
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCC

C
C  Search for number of propagating plane wave orders

      s = 1
      do px = -ordre_ls, ordre_ls
        do py = -ordre_ls, ordre_ls
          if (px**2 + py**2 .le. ordre_ls**2) then
            alpha = bloch1 + vec_kx*px  ! Bloch vector along x
            beta  = bloch2 + vec_ky*py  ! Bloch vector along y
            z_tmp = k_1**2 - alpha**2 - beta**2
            beta_z_pw(s) = sqrt(z_tmp)
            s = s + 1
          endif
        enddo
      enddo
C  
      numberprop_S = 0
C        write(6,*) "J_o_sub k1", k_1
      do i=1,neq_PW
        test = (beta_z_pw(i)**2)
C        write(6,*) "J_o_sub beta_z", beta_z_pw(i)
        if (REAL(test) .gt. 0.0d-5) then
          numberprop_S = numberprop_S + 1
        endif
      enddo
C
C     set up empty X_mat
      do j=1,2*neq_PW        
        do i=1,2*neq_PW
          X_mat(i,j) = 0.0D0
        enddo 
      enddo 

C     Set up X_mat
      s = 1
      do px = -ordre_ls, ordre_ls
        do py = -ordre_ls, ordre_ls
          if (px**2 + py**2 .le. ordre_ls**2) then
            alpha = bloch1 + vec_kx*px  ! Bloch vector along x
            beta  = bloch2 + vec_ky*py  ! Bloch vector along y
            s2 = index_pw_inv(s)
            X_tmp = (k_1**2 - alpha**2 - beta**2)
            X_tmp = SQRT(X_tmp)
            X_mat(s2,s2) = X_tmp/k_0
            X_mat(neq_PW+s2,neq_PW+s2) = eps_eff_sub * k_0/X_tmp
            s = s + 1
          endif
        enddo
      enddo

C
CCCCCCCCCCCCCCCCC	  save results   CCCCCCCCCCCCCCCC
C
      if (PrintAll .eq. 1) then
C
      open (unit=34, file="Matrices/X_mat_b.txt", status='unknown')
      write(34,131) lambda, freq    
      do s = 1, 2*neq_PW
        write(34,13) X_mat(s,s), abs(X_mat(s,s))     
      enddo
      close(34)
C
131   format(2(f12.4))
13    format(2(g25.17),g18.10)
      endif
C
C
      return
      end 
