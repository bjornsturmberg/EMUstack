      SUBROUTINE z_indexx(n,arr,indx)

      INTEGER*8 n,indx(n),M,NSTACK
      complex*16 arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER*8 i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      double precision a
      integer :: allocate_status=0
      double precision, dimension(:), allocatable :: arr_0
c
C
      allocate(arr_0(n), STAT=allocate_status)
      if (allocate_status /= 0) then
        write(*,*) "The allocation is unsuccessful"
        write(*,*) "allocate_status = ", allocate_status
        write(*,*) "z_indexx: Not enough memory for arr_0"
        write(*,*) "nval = ", n
        write(*,*) "Aborting..."
        stop
      endif
C
      do 11 j=1,n
        indx(j) = j
11    continue
c
      do j=1,n
        arr_0(j) = -1*dble(arr(j)**2)
      enddo

c
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a = arr_0(indxt)
          do 12 i=j-1,l,-1
            if( arr_0(indx(i)) .le. a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=l-1
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if( arr_0(indx(l)) .gt. arr_0(indx(ir))) then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if( arr_0(indx(l+1)) .gt. arr_0(indx(ir))) then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if( arr_0(indx(l)) .gt. arr_0(indx(l+1))) then
          itemp=indx(l)
          indx(l)=indx(l+1)
          indx(l+1)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l+1)
        a = arr_0(indxt)
3       continue
          i=i+1
        if( arr_0(indx(i)) .lt. a)goto 3
4       continue
          j=j-1
        if( arr_0(indx(j)) .gt. a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l+1)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
C        if(jstack.gt.NSTACK)pause 'NSTACK too small in indexx'
        if(jstack.gt.NSTACK) then
          write(*,*) 'NSTACK too small in indexx'
          write(*,*) "z_indexx: Aborting..."
          stop
        endif
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
C
      deallocate(arr_0)
      END
