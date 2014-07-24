c
c***********************************************************************
c
      subroutine periodic_N_E_F (n_ddl, 
     *      type_N_E_F, x_N_E_F, ip_period_N_E_F, 
     *      n_period_N_E_F, lat_vecs)
c
c***********************************************************************
c
c   type_N_E_F(1,i) = 0  => interiour point
c   type_N_E_F(1,i) != 0 => boundary point
c   type_N_E_F(2,i)  = 0, 1, 2 => diemnsion of the domain
c
c
c***********************************************************************
c
      implicit none
      integer*8 n_ddl
      integer*8 type_N_E_F(2,n_ddl), ip_period_N_E_F(n_ddl)
      integer*8 n_period_N_E_F(n_ddl)
      double precision lat_vecs(2,2)
      double precision x_N_E_F(2,n_ddl)
      integer*8 i, j, i1, j1, i_not_periodic
      integer*8 i_boundary1, i_dim1, i_boundary2, i_dim2
      double precision tmp1, tmp2, tol
      double precision delta_v(2),  vec(2)
      integer*8 ix, iy, test_lattice
      integer*8 k, debug
c
      debug = 0
      tol = 1.0d-6
      do i=1,n_ddl
        ip_period_N_E_F(i) = 0
        n_period_N_E_F(i) = 0
      enddo
      do i=1,n_ddl-1
        i_boundary1 = type_N_E_F(1,i)
        i_dim1 = type_N_E_F(2,i)
        if(i_boundary1 .ne. 0) then
          do j=i+1,n_ddl
            i_boundary2 = type_N_E_F(1,j)
            i_dim2 = type_N_E_F(2,j)
            if(i_boundary2 .ne. 0  .and. i_dim2 .eq. i_dim1 
     *          .and. ip_period_N_E_F(j) .eq. 0) then
              delta_v(1) = x_N_E_F(1,j) - x_N_E_F(1,i)
              delta_v(2) = x_N_E_F(2,j) - x_N_E_F(2,i)
              test_lattice = 0
              do ix=-1,1
                do iy=-1,1
                  do k=1,2
                    vec(k) = ix*lat_vecs(k,1)
     *                + iy*lat_vecs(k,2)
                  enddo
                  tmp1 = 0.0d0
                  do k=1,2
                    tmp2 = delta_v(k) - vec(k)
                    tmp1 = tmp1 + abs(tmp2)
                  enddo
                  if (tmp1 .lt. tol) then
                    test_lattice = 1
                    goto 10
                  endif
                enddo
              enddo
10            continue
c
              if(test_lattice .eq. 1) then
                n_period_N_E_F(i) = n_period_N_E_F(i) + 1
                i1 = ip_period_N_E_F(i)
                j1 = ip_period_N_E_F(j)
                if(i1 .eq. 0) then
                  ip_period_N_E_F(i) = i
                  ip_period_N_E_F(j) = i
                elseif(i1 .ne. 0 .and. j1 .eq. 0) then
                  ip_period_N_E_F(j) = i1
                else
                  write(*,*)
                  write(*,*) "  ???"
                  write(*,*) "periodic_N_E_F: for i, j = ", i, j
                  write(*,*) "periodic_N_E_F: ip_period_N_E_F : ",i1,j1
                  write(*,*) "periodic_N_E_F: Aborting..."
                  stop
c                  ip_period_N_E_F(i) = min(i,i1)
c                  ip_period_N_E_F(j) = min(i,i1)
                endif
              endif
            endif
          enddo
        endif
      enddo
c
      i_not_periodic = 0
      do i=1,n_ddl
        i_boundary1 = type_N_E_F(1,i)
        i_dim1 = type_N_E_F(2,i)
        if(i_boundary1 .ne. 0 .and. ip_period_N_E_F(i) .eq. 0) then
          if(i_not_periodic .eq. 0) then
            open (unit=11,file="not_period.txt",status='unknown')
          endif
        i_not_periodic = i_not_periodic + 1
        write(11,*) i, i_boundary1, i_dim1, ip_period_N_E_F(i), 
     *       n_period_N_E_F(i)
        endif
      enddo

      if(i_not_periodic .gt. 0) then
        close ( unit = 11)
        write(*,*)
        write(*,*) "  ???"
        write(*,*) "periodic_N_E_F: the FEM mesh is not periodic"
        write(*,*) "periodic_N_E_F: see the file not_period.txt"
        write(*,*) "periodic_N_E_F: Aborting..."
        write(*,*)
        stop
      endif
c
      if (debug .eq. 1) then
      open (unit = 10, file="ip_period_N_E_F.txt",status='unknown')
        do i=1,n_ddl
          i_boundary1 = type_N_E_F(1,i)
          i_dim1 = type_N_E_F(2,i)
          write(10,*) i, i_boundary1, i_dim1, ip_period_N_E_F(i), 
     *       n_period_N_E_F(i), "     ", (x_N_E_F(j,i),j=1,2)
        enddo
      close ( unit = 10)
      endif
c
cccccccccccccccccccccccccccccccccccccccccccc
c
      return
      end
