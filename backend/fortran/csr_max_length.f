
c
c***********************************************************************
c
      subroutine csr_max_length (nel, n_ddl, neq, nnodes, 
     *  table_N_E_F, ineq, lb, nonz)
c
      implicit none
      integer*8 nel, neq, n_ddl, nnodes, nonz
      integer*8 table_N_E_F (14,nel)
      integer*8 ineq(3,n_ddl), lb(neq+1)

c     Local variables
      integer*8 nddl_0
      parameter (nddl_0 = 14)

      integer*8 i, k, iel, ind_ip, ip
      integer*8 k_copy1, k_copy2
c
      if ( nnodes .ne. 6 ) then
        write(*,*) "csr_max_length: problem nnodes = ", nnodes
        write(*,*) "csr_max_length: nnodes should be equal to 6 !"
        write(*,*) "csr_max_length: Aborting..."
        stop
      endif
c
      do i=1,neq+1
         lb(i) = 0
      enddo
c
c   Determination of the bandwidths
c
      do 20 iel=1,nel
        do i=1,nddl_0
          ip = table_N_E_F(i,iel)
          do k=1,3
            ind_ip = ineq(k,ip)
            if (ind_ip .ne. 0) lb(ind_ip) = lb(ind_ip)+1
          enddo
        enddo
20    continue
c
      nonz = 0
      do i=1,neq
        nonz = nonz + 3*nddl_0 + 3*(nddl_0-1)*(lb(i)-1)
      enddo
c
c      print*
c      do i=1,neq
c        print*, "csr_max_length: i, lb(i) = ", i, lb(i)
c      enddo
c      print*, "csr_max_length: neq, n_ddl nonz, = ", neq, n_ddl, nonz
c      print*
c

c     Compressed Row Storage (CRS): determine the row pointer
c
      k_copy1 = lb(1)
      lb(1) = 1
      do i=2,neq+1
        k_copy2 = lb(i)
        lb(i) = lb(i-1) + 3*nddl_0 + 3*(nddl_0-1)*(k_copy1-1)
        k_copy1 = k_copy2
      enddo

      nonz = lb(neq+1) - 1
c
c
      return
      end
