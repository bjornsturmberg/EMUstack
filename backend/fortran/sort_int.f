c
c *********************************************************
c
c     Sort an array in increasing order
c
c    Sorts an array arr(1:n) into ascending numerical order using the Quicksort algorithm. n
c    is input; arr is replaced on output by its sorted rearrangement.
c    Parameters: M is the size of subarrays sorted by straight insertion and NSTACK is the required
c    auxiliary storage.
c
c    Reference 1: Numerical Recipes in Fortran 77; The Art of Scientific Computing
c                 Second Edition
c                 See "Chapter 8. Sorting", "Section 8.2 Quicksort", page 324
c
c    Reference 2: Sedgewick, R. 1978, Communications of the ACM, vol. 21, pp. 847–857. [1]
c
c *********************************************************
c
c
      SUBROUTINE sort_int (n, arr, indx, istack)

c      istack is a working array
      implicit none
      INTEGER*8 n, indx(n), M, NSTACK
      INTEGER*8 arr(n), istack(n)
      PARAMETER (M=7)
      INTEGER*8 i,indxt,ir,itemp,j,jstack,k,l
      double precision a


      NSTACK = n
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,l,-1
            if(arr(indx(i)).le.a)goto 2
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
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(l+1)))then
          itemp=indx(l)
          indx(l)=indx(l+1)
          indx(l+1)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l+1)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l+1)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK) then
          write(*,*) 'NSTACK too small in sort_int'
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
      END
