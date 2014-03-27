c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      x=[x0,y0] = coordinatate in the reference tetrahedron
c      x_g = Corresponding coordinatate in the actual tetrahedron
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine jacobian_p1_2d (x, xel, nnodes, 
     *  x_g, det_jacobian, mat_B_0, mat_T)
c
      implicit none
      integer*8 nnodes
      double precision x(2), xel(2,nnodes)
      double precision mat_B(2,2), mat_T(2,2), mat_B_0(2,2)
      double precision det_jacobian, x_g(2)
      integer*8 i, j
      double precision phi, phi_x
c
      double precision ZERO, ONE
      parameter ( ZERO = 0.0D0, ONE = 1.0D0)
c     32-but integers for BLAS and LAPACK
      integer*4 INFO_32, LDB_32, LDT_32
      integer*4 IPIV_32(2), NRHS_32, N_32
c
c
c
c    Compute the Affine mappings from the current Tetrahedron to the
c     reference unit Tetrahedron N. Integration will be performed on that
c     Tetrahedron by Gaussian quadrature.
c
c    X_g = B X + (x_0, y_0, z_0)^t
c
      do i=1,2
        do j=1,2
          mat_B(j,i) = xel(j,i+1) - xel(j,1)
        enddo
      enddo

      do i=1,2
        do j=1,2
          mat_B_0(j,i) = mat_B(j,i)
        enddo
      enddo
c
ccccccccccccccccccccccccccccc
c        mat_T = mat_B
c
      N_32 = 2 ! The order of the matrix mat_B
      NRHS_32 = N_32 ! The number of right hand sides
      LDB_32 = N_32 ! The leading dimension of the array mat_B
      LDT_32 = N_32 ! The leading dimension of the array mat_T
c
c     X_g = B X + (x_0, y_0, z_0)^t
c     Initialisation for DGEMV X_g = xel(1:2,1)
      do i=1,2
        X_g(i) = xel(i,1)
      enddo
c     DGEMV  performs one of the matrix-vector operations
      call DGEMV('No transpose', N_32, N_32, ONE, mat_B, 
     *  LDB_32, X, 1, ONE, X_g, 1)
c
c     Initialisation for DGESV: mat_T = identity
      do i=1,2
        do j=1,2
          mat_T(i,j) = 0.0d0
        enddo
          mat_T(i,i) = 1.0d0
      enddo
      call DGESV( N_32, NRHS_32, mat_B, LDB_32, IPIV_32, 
     *            mat_T, LDT_32, INFO_32 )
c
      if(INFO_32 .ne. 0) then
        write(*,*) 
        write(*,*) "jacobian_p1_2d: TRANSF_MAT: ATTENTION, INFO_32 = ",
     *   INFO_32
        stop
      endif

c     The value determinant can be obtained from the factorization P*L*U
      det_jacobian = 1
      do i=1,2
        if( (IPIV_32(i)-i) .eq. 0) then
          det_jacobian = det_jacobian*mat_B(i,i)
        else
          det_jacobian = -det_jacobian*mat_B(i,i)
        endif
      enddo

       if(abs(det_jacobian) .lt. 1.0d-10) then
         write(*,*)
         write(*,*) "   ???"
         write(*,*) "jacobian_p1_2d: det = 0 : det = ", det_jacobian
         do i=1,nnodes
           write(*,*) "i, x, y, x, = ", i, (xel(j,i),j=1,2)
         enddo
         write(*,*) "jacobian_p1_2d: Aborting..."
         stop
       endif

c
      return
      end
