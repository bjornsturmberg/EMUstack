
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine type_arete(i, i1, i2, ne_d1, nu_d1, typ_el_d1)
c
      implicit none
      integer i, i1, i2, ne_d1
      integer nu_d1(3,ne_d1), typ_el_d1(ne_d1)
      integer j, k, k1, k2
c
      i = 0
      do j=1,ne_d1
        k = typ_el_d1(j)
        k1 = nu_d1(1,j)
        k2 = nu_d1(2,j)
        if(k1 .eq. i1 .and. k2 .eq. i2) then
          i = k
          return
        elseif(k1 .eq. i2 .and. k2 .eq. i1) then
          i = k
          return
        endif
      enddo
c      print*, '?? type_arete: Resultat negatif'
c      print*, 'i1, i2 = ', i1, i2
c      stop
c
      return
      end
c


c
c###############################################
c


c
      subroutine renumerote(npt, ne, idfn, nu, x, y, ui)
c
      implicit none
      integer ne, npt, ui
      integer nu(6,ne), idfn(npt)
      double precision x(npt), y(npt)
c
      integer i, j, k, i1, ip(2,3)
      integer long_adj
c      xadj(npt+1),
      integer npt_max, long_adj_max
      parameter(npt_max=250000,long_adj_max=2500000 )
      integer visite(npt_max), lb(npt_max)
      integer xadj(npt_max+1), adjncy(long_adj_max)
      integer nut(6)
      character file_ui*100
      common/imp_file/file_ui
c
c ip(1,i) = i+1 MOD 3
c ip(2,i) = i+2 MOD 3
c
      ip(1,1) = 2
      ip(1,2) = 3
      ip(1,3) = 1
c
      ip(2,1) = 3
      ip(2,2) = 1
      ip(2,3) = 2
c
      if(npt_max .lt. npt) then
        open (unit=ui,file=file_ui)
        write(*,*) 'Renumerote (in conv_gmsh.f): ',
     *    ' npt_max < npt ', npt_max, npt
        close(ui)
        stop
      endif
c
      do i = 1,npt
        visite(i) = 0
      enddo
c
      do i = 1,npt
        lb(i)=0
      enddo
c
      do j = 1,ne
c
        do i=1,6
          nut(i) = nu(i,j)
        enddo
        do i=1,6
          k = nut(i)
          lb(k) = lb(k) + 6 - 1
        enddo
        do i=1,3
          i1 = nut(i+3)
          if(visite(i1) .eq. 0) then
            visite(i1) = 1
          else  ! les points sur larete i ont deja ete compte
            k = nut(i)
            lb(k) = lb(k) - 2
            k = nut(ip(1,i))
            lb(k) = lb(k) - 2
            k = nut(i+3)
            lb(k) = lb(k) - 2
          endif
        enddo
      enddo
c
      xadj(1) = 1
      do i = 1,npt
        xadj(i+1) = xadj(i) + lb(i)
      enddo
      long_adj = xadj(npt+1)
c      print*, 'renumerote: long_adj = ', long_adj
c
      if(long_adj_max .lt. long_adj) then
        open (unit=ui,file=file_ui)
        write(*,*) 'Renumerote (in conv_gmsh.f): ',
     *      'long_adj_max < long_adj ',
     *      long_adj_max, long_adj
        close(ui)
        stop
      endif
c
c      print*, 'Appel de mat_ad'
      call mat_adj(npt, ne, long_adj, nu,
     *      xadj, adjncy, x, y, idfn, ui)
c
      return
      end

c
c###############################################
c


Cc***********************************************************************
Cc
Cc
C      subroutine test_border(npt, type_nod0, x, y, i_err,
C     *                       file1_mesh, ui)
Cc
C      implicit none
C      integer npt, i_err
C      integer type_nod0(npt)
C      double precision x(npt), y(npt)
Cc
C      integer type_nod, n_typ2
C      double precision d_period, tmp, tmp2
Cc
C      integer  n_border1, n_border2, n_typ1, ui
C      integer i, j, i1, i2, j1, j2, max_period
C      parameter(max_period=2000)
C      double precision x2(2), x1(2), y_min, y_max
C      double precision delta_x, tol
C      integer period1(max_period), period2(max_period)
C      integer period3(2,max_period)
C      character*(*) file1_mesh
C      character file_ui*100
C      common/imp_file/file_ui
Cc
Cc      character fichier_mail*20
Cc
Cc     debut de la programmation de la sous-routine maillage
Cc     -----------------------------------------------------
Cc
C      i_err = 0
C      tol = 1.0d-6
C      y_min = y(1)
C      y_max = y(1)
C      do i=1,npt
C        if(y(i) .lt. y_min) y_min = y(i)
C        if(y(i) .gt. y_max) y_max = y(i)
C      enddo
C      i1 = 0
C      i2 = 0 
C      do i=1,npt
C        if(abs(y(i)-y_max) .lt. tol) then
C          if(i1 .eq. 0) then
C            i1 = 1
C            x1(1) = x(i)
C            x1(2) = x(i)
C          endif
C          if(x(i) .lt. x1(1)) x1(1) = x(i)
C          if(x(i) .gt. x1(2)) x1(2) = x(i)
C        endif
C        if(abs(y(i)-y_min) .lt. tol) then
C          if(i2 .eq. 0) then
C            i2 = 1
C            x2(1) = x(i)
C            x2(2) = x(i)
C          endif
C          if(x(i) .lt. x2(1)) x2(1) = x(i)
C          if(x(i) .gt. x2(2)) x2(2) = x(i)
C        endif
C      enddo
C      delta_x = x1(1) - x2(1)
C      d_period = x1(2)-x1(1)
Cc
C      n_border1 = 0
C      n_border2 = 0
C      n_typ1 = 0
C      n_typ2 = 0
Cc
C      do i=1,npt
C        type_nod = type_nod0(i)
C        if(type_nod .eq. 1) then
C          n_typ1 = n_typ1 + 1
C        elseif(type_nod .eq. 2) then
C          n_typ2 = n_typ2 + 1
C        elseif(type_nod .eq. 3) then
C          n_border1 = n_border1 + 1
C         period1(n_border1) = i
C        elseif(type_nod .eq. 4) then
C          n_border2 = n_border2 + 1
C         period2(n_border2) = i
C        endif
C       if(n_border2 .ge. max_period) then
C        open (unit=ui,file=file_ui)
C           write(*,*) 'get_border (in conv_gmsh.f): attention'
C           write(*,*) 'n_border2 >= max_period = ', 
C     *        n_border2, max_period
C        close(ui)
C         i_err = 1
Cc         stop
C       endif
C      enddo
Cc
C      if(n_border1 .ne. n_border2) then
C        open (unit=ui,file=file_ui)
C          write(*,*) 'get_border (in conv_gmsh.f): ',
C     *    ' n_border1 .ne. n_border2'
C        close(ui)
C          i_err = 2
Cc        stop
C      endif
C      do i=1,n_border1
C        i1 = period1(i)
C        period3(1,i) = i1
C        period3(2,i) = 0
C        do j=1,n_border1
C          j1=period2(j)
C          tmp=dsqrt((x(j1)-x(i1))**2+(y(j1)-y(i1))**2)
C          tmp2 = abs(y(j1)-y(i1))
C        if(abs(tmp-d_period) .lt. 5.0d-6 .and. tmp2 .lt. 5.0d-6) then
C            if(period3(2,i) .ne. 0) then
C           open (unit=ui,file=file_ui)
C              write(*,*) 'get_border  (in conv_gmsh.f): ',
C     *        'i,j = ', i, j, i1, j1
C              write(*,*) 'probleme avec period3(2,i) = ', period3(2,i)
C              write(*,*) x(i1), y(i1), d_period
C              write(*,*) x(j1), y(j1), tmp, (tmp-d_period)
C              j2 = period3(2,i)
C              tmp=dsqrt((x(j2)-x(i1))**2+(y(j2)-y(i1))**2)
C              write(*,*) x(j2), y(j2), tmp, (tmp-d_period)
C              i_err = 3
C              close(ui)
Cc              stop
C            endif
C            period3(2,i) = j1
C          endif
C        enddo
C      enddo
C      do i=1,n_border1
C        if(period3(2,i) .eq. 0) then
C           open (unit=ui,file=file_ui)
C              write(*,*) 'get_border (in conv_gmsh.f): ',
C     *        ' period3(2,i) = 0', i
C           close(ui)
C          i_err = 4
Cc          stop
C        endif
C      enddo
Cc
C      return
C      end
Cc
C
c
c###############################################
c

c     ------------------------------------------------------------------

      SUBROUTINE  DEGREE (ROOT,XADJ,ADJNCY,MASK,DEG,CCSIZE,LS)

c     ------------------------------------------------------------------

      INTEGER ADJNCY(*), DEG(*), LS(*), MASK(*)
      INTEGER XADJ(*), CCSIZE, I, IDEG, J, JSTOP, JSTRT
      INTEGER LBEGIN, LVLEND, LVSIZE, NBR, NODE, ROOT

c     ------------------------------------------------------------------

      LS(1) = ROOT
      XADJ(ROOT) = -XADJ(ROOT)
      LVLEND = 0
      CCSIZE = 1
100   LBEGIN = LVLEND + 1
      LVLEND = CCSIZE
      DO 400 I = LBEGIN, LVLEND
         NODE = LS(I)
         JSTRT = -XADJ(NODE)
         JSTOP = IABS(XADJ(NODE + 1)) - 1
         IDEG = 0
         IF (JSTOP.LT.JSTRT) GO TO 300
         DO 200 J = JSTRT, JSTOP
            NBR = ADJNCY(J)
            IF (MASK(NBR).EQ.0) GO TO 200
            IDEG = IDEG + 1
            IF (XADJ(NBR).LT.0) GO TO 200
            XADJ(NBR) = -XADJ(NBR)
            CCSIZE = CCSIZE + 1
            LS(CCSIZE) = NBR
200      CONTINUE
300      DEG(NODE) = IDEG
400   CONTINUE
      LVSIZE = CCSIZE - LVLEND
      IF (LVSIZE.GT.0) GO TO 100
      DO 500 I = 1, CCSIZE
         NODE = LS(I)
         XADJ(NODE) = -XADJ(NODE)
500   CONTINUE


      END

c
c###############################################
c

c     ------------------------------------------------------------------

      SUBROUTINE FNROOT (ROOT,XADJ,ADJNCY,MASK,NLVL,XLS,LS)

c     ------------------------------------------------------------------

      INTEGER ADJNCY(*), LS(*), MASK(*), XLS(*)
      INTEGER XADJ(*), CCSIZE, J, JSTRT, K, KSTOP, KSTRT
      INTEGER MINDEG, NABOR, NDEG, NLVL, NODE, NUNLVL, ROOT

c     ------------------------------------------------------------------

      CALL ROOTLS (ROOT,XADJ,ADJNCY,MASK,NLVL,XLS,LS)
      CCSIZE = XLS(NLVL+1) - 1
      IF (NLVL.EQ.1 .OR. NLVL.EQ.CCSIZE) RETURN
100   JSTRT = XLS(NLVL)
      MINDEG = CCSIZE
      ROOT = LS(JSTRT)
      IF (CCSIZE.EQ.JSTRT) GO TO 400
      DO 300 J = JSTRT, CCSIZE
         NODE = LS(J)
         NDEG = 0
         KSTRT = XADJ(NODE)
         KSTOP = XADJ(NODE+1) - 1
         DO 200 K = KSTRT, KSTOP
            NABOR = ADJNCY(K)
            IF (MASK(NABOR).GT.0) NDEG = NDEG + 1
200      CONTINUE
         IF (NDEG.GE.MINDEG) GO TO 300
         ROOT = NODE
         MINDEG = NDEG
300   CONTINUE
400   CALL ROOTLS (ROOT,XADJ,ADJNCY,MASK,NUNLVL,XLS,LS)
      IF (NUNLVL.LE.NLVL) RETURN
      NLVL = NUNLVL
      IF (NLVL.LT.CCSIZE) GO TO 100


      END

c
c###############################################
c

c     ------------------------------------------------------------------

      SUBROUTINE  GENRCM (NEQNS,XADJ,ADJNCY,PERM,MASK,XLS) 

c     ------------------------------------------------------------------

      INTEGER ADJNCY(*), MASK(*), PERM(*), XLS(*)
      INTEGER XADJ(*), CCSIZE, I, NEQNS, NLVL, NUM, ROOT                                                
c     ------------------------------------------------------------------


      DO 100 I = 1, NEQNS
         MASK(I) = 1
100   CONTINUE
      NUM = 1
      DO 200 I = 1, NEQNS
         IF (MASK(I).EQ.0) GO TO 200
         ROOT = I
         CALL FNROOT (ROOT,XADJ,ADJNCY,MASK,NLVL,XLS,PERM(NUM))
         CALL RCM (ROOT,XADJ,ADJNCY,MASK,PERM(NUM),CCSIZE,XLS)
         NUM = NUM + CCSIZE
         IF (NUM.GT.NEQNS) RETURN
200   CONTINUE


      END
c
c###############################################
c


c     ------------------------------------------------------------------
c
      subroutine mailp2(tcp2,maxvis,ne,npt,ns,liste,nb,numero)

c     ------------------------------------------------------------------

      integer flag, i, ip(3), is1, is2, j, jel, jj, ne, npt, ns, temp
      integer liste(maxvis,ns), maxvis, nb(ns), numero(maxvis,ns)
      integer tcp2(6,ne)
c
c     ------------------------------------------------------------------
c
c      print*, 'MAILP2: 0: npt = ', npt
      ip(1) = 2
      ip(2) = 3
      ip(3) = 1
c
      do i = 1, ns
         nb(i) = 0
      enddo
c
      npt = ns
      do jel = 1, ne
         do i = 1, 3
c
            is1 = tcp2(i,jel)
            is2 = tcp2(ip(i),jel)
            if (is1.gt.is2) then
               temp = is1
               is1  = is2
               is2  = temp
            endif
c
            flag = 0
c
            if (nb(is1).eq.0) go to 1
c
            do j = 1, nb(is1)
               if (liste(j,is1).eq.is2) then
                  flag = 1
                  jj   = j
                  goto 1
               endif
            enddo
c
1           continue
            if (flag.eq.0) then
c
c              l'arete n'est pas dans la liste
c              -------------------------------
c
               npt = npt+1
               tcp2(i+3,jel) = npt
               nb(is1) = nb(is1)+1
               liste(nb(is1),is1) = is2
               numero(nb(is1),is1) = npt
            else
c
c              l'arete est deja dans la liste
c              ------------------------------
c
               tcp2(i+3,jel) = numero(jj,is1)
            endif
c
         enddo
c
      enddo
c
      end

c
c###############################################
c


c
      subroutine mat_adj(npt, ne, long_adj, nu,
     *      xadj, adjncy, x, y, idfn, ui)
c
      implicit none
      integer ne, npt, long_adj, ui
      integer nu(6,ne), idfn(npt)
      integer xadj(npt+1), adjncy(long_adj)
      double precision x(npt), y(npt)
c
      integer i, j, k, i1, j2, j3, k1, ind1, ind2
      integer npt_max, ip(2,3), m, m1
      parameter(npt_max=250000)
      integer lb2(npt_max)
      integer nut(6)
      integer mask(npt_max), perm(npt_max)
      integer xls(npt_max), invperm(npt_max)
      integer idfn_r(npt_max)
      double precision x_r(npt_max), y_r(npt_max)
      character file_ui*100
      common/imp_file/file_ui
c
c ip(1,i) = i+1 MOD 3
c ip(2,i) = i+2 MOD 3
c
      ip(1,1) = 2
      ip(1,2) = 3
      ip(1,3) = 1
c
      ip(2,1) = 3
      ip(2,2) = 1
      ip(2,3) = 2
c
      if(npt_max .lt. npt) then
        open (unit=ui,file=file_ui)
          write(*,*) 'mat_adj (in conv_gmsh.f): ',
     *    ' npt_max < npt', npt_max, npt
        close(ui)
        stop
      endif
c
      do i = 1,npt
        lb2(i)=0
      enddo
c
c
      do j = 1,ne
        do i=1,6
          nut(i) = nu(i,j)
        enddo
        do i=1,6
          i1=nut(i)
          ind1=xadj(i1)
          do k=1,6
            k1=nut(k)
            ind2=0
            if(k .eq. i) then
              ind2=1
              goto 5
            endif
            do m=1,lb2(i1)
              m1=adjncy(ind1+m-1)
              if(m1 .eq. k1) then
                ind2=1
                goto 5
              endif
            enddo
5           continue
            if(ind2 .eq. 0) then
              lb2(i1) = lb2(i1)+1
              adjncy(ind1+lb2(i1)-1)=k1
            endif
          enddo
        enddo
      enddo
c
c     execution de l'algorithme RCM
c     -----------------------------

      call genrcm(npt,xadj,adjncy,perm,mask,xls)
c
      do i = 1,long_adj-1
        if(adjncy(i) .eq. 0) then
        open (unit=ui,file=file_ui)
          write(*,*) 'MAT_ADJ (in conv_gmsh.f): ',
     *    ' Attention, adjncy(i)=0'
          write(*,*) 'i, adjncy(i) = ', i, adjncy(i)
        close(ui)
          stop
        endif
      enddo
c
c     creation de la nouvelle table tcp2
c     ----------------------------------
c
      if(.true.) then   ! ##########################
      do i = 1, npt
        invperm(perm(i)) = i
      enddo
c
      do i = 1, ne
         do j = 1, 6
           nu(j,i) = invperm(nu(j,i))
         enddo
      enddo
c
      do i = 1, npt
        idfn_r(i) = idfn(perm(i))
        x_r(i) = x(perm(i))
        y_r(i) = y(perm(i))
      enddo
c
c     $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
      do i = 1, npt
        idfn(i) = idfn_r(i)
        x(i) = x_r(i)
        y(i) = y_r(i)
      enddo
c
c     $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
      endif   ! ##########################
c
      return
      end

c
c###############################################
c


c     ------------------------------------------------------------------
c
      subroutine prepmailp2(tcp1,maxvis,ne,ns,visite)

c     ------------------------------------------------------------------

      integer i, ii, jel, maxvis, ne, ns, tcp1(6,ne), visite(ns)

c     ------------------------------------------------------------------


      do i = 1, ns
         visite(i) = 0
      enddo

      do jel = 1, ne
         do i = 1, 3
            ii = tcp1(i,jel)
            visite(ii) = visite(ii)+1
         enddo
      enddo

      maxvis = 0
      do i = 1, ns
         maxvis = max(maxvis,visite(i))
      enddo


      end

c
c###############################################
c

c     ------------------------------------------------------------------

      SUBROUTINE  RCM (ROOT,XADJ,ADJNCY,MASK,PERM,CCSIZE,DEG)

c     ------------------------------------------------------------------

      INTEGER ADJNCY(*),DEG(*),MASK(*),PERM(*)
      INTEGER XADJ(*),CCSIZE,FNBR,I,J,JSTOP,JSTRT,K,L,LBEGIN,LNBR
      INTEGER LPERM,LVLEND,NBR,NODE,ROOT

c     ------------------------------------------------------------------

      CALL DEGREE (ROOT,XADJ,ADJNCY,MASK,DEG,CCSIZE,PERM)
      MASK(ROOT) = 0
      IF (CCSIZE.LE.1) RETURN
      LVLEND = 0
      LNBR = 1
100   LBEGIN = LVLEND + 1
      LVLEND = LNBR
      DO 600 I = LBEGIN, LVLEND
         NODE = PERM(I)
         JSTRT = XADJ(NODE)
         JSTOP = XADJ(NODE+1) - 1
         FNBR = LNBR + 1
         DO 200 J = JSTRT, JSTOP
            NBR = ADJNCY(J)
            IF (MASK(NBR).EQ.0) GO TO 200
            LNBR = LNBR + 1
            MASK(NBR) = 0
            PERM(LNBR) = NBR
200      CONTINUE
         IF (FNBR.GE.LNBR) GO TO 600
         K = FNBR
300      L = K
         K = K + 1
         NBR = PERM(K)
400      IF (L.LT.FNBR) GO TO 500
         LPERM = PERM(L)
         IF (DEG(LPERM).LE.DEG(NBR)) GO TO 500
         PERM(L+1) = LPERM
         L = L - 1
         GO TO 400
500      PERM(L+1) = NBR
         IF (K.LT.LNBR) GO TO 300
600   CONTINUE
      IF (LNBR.GT.LVLEND) GO TO 100
      K = CCSIZE/2
      L = CCSIZE
      DO 700 I = 1, K
         LPERM = PERM(L)
         PERM(L) = PERM(I)
         PERM(I) = LPERM
         L = L - 1
700   CONTINUE


      END

c
c###############################################
c

c     ------------------------------------------------------------------

      SUBROUTINE  ROOTLS (ROOT,XADJ,ADJNCY,MASK,NLVL,XLS,LS)

c     ------------------------------------------------------------------

      INTEGER ADJNCY(*), LS(*), MASK(*), XLS(*)
      INTEGER XADJ(*), I, J, JSTOP, JSTRT, LBEGIN
      INTEGER CCSIZE, LVLEND, LVSIZE, NBR, NLVL, NODE, ROOT                                               
c     ------------------------------------------------------------------

      MASK(ROOT) = 0
      LS(1) = ROOT
      NLVL = 0
      LVLEND = 0
      CCSIZE = 1
200   LBEGIN = LVLEND + 1
      LVLEND = CCSIZE
      NLVL = NLVL + 1
      XLS(NLVL) = LBEGIN
      DO 400 I = LBEGIN, LVLEND
         NODE = LS(I)
         JSTRT = XADJ(NODE)
         JSTOP = XADJ(NODE + 1) - 1
         IF (JSTOP.LT.JSTRT) GO TO 400
         DO 300 J = JSTRT, JSTOP
            NBR = ADJNCY(J)
            IF (MASK(NBR).EQ. 0) GO TO 300
            CCSIZE = CCSIZE + 1
            LS(CCSIZE) = NBR
            MASK(NBR) = 0
300      CONTINUE
400   CONTINUE
      LVSIZE = CCSIZE - LVLEND
      IF (LVSIZE.GT.0) GO TO 200
      XLS(NLVL+1) = LVLEND + 1
      DO 500 I = 1, CCSIZE
         NODE = LS(I)
         MASK(NODE) = 1
500   CONTINUE


      END

c
c###############################################
c

c
c###############################################
c
      subroutine symmetry(npt, ne, 
     *      max_ne, max_npt, idfn, nu, typ_el, 
     *      x, y, i_sym)

c*******************************************************
c
c     symmetry: symmetrize an FEM mesh
c
c*******************************************************
c
      implicit none
      integer i_sym, max_ne, max_npt

      integer ne, npt
      integer nu(6,max_ne), typ_el(max_ne)
      integer idfn(max_npt)
      double precision x(max_npt), y(max_npt)

c     Local data

      integer max_ne_0, max_npt_0
      parameter(max_npt_0=250000, max_ne_0=120000)
      integer ne_0, npt_0, idfn_0(max_npt_0)
      integer nu_0(6,max_ne_0), typ_el_0(max_ne_0)
      double precision x_0(max_npt_0),  y_0(max_npt_0)
      integer tab_ne(max_ne_0), tab_npt(max_npt_0,3)
c
c
      integer i, j, k
c
      integer debug, ui
      character file_ui*100
c
      common/imp/ui, debug
      common/imp_file/file_ui
c
ccccccccccccccccccccccccc
c
        npt_0 = npt
        ne_0 = ne
        do i=1,npt_0
          x_0(i) = x(i)
          y_0(i) = y(i)
          idfn_0(i) = idfn(i)
        enddo
c
        do i=1,ne_0
          do j=1,6
            nu_0(j,i) = nu(j,i)
          enddo
          typ_el_0(i) = typ_el(i)
        enddo
c
ccccccccccc
      if(i_sym .eq. 1) then
        call y_symmetry(npt, ne, ne_0, npt_0, 
     *      max_ne, max_npt, idfn, nu, typ_el, 
     *      idfn_0, nu_0, typ_el_0, x, y, x_0, y_0,
     *      tab_ne, tab_npt)
cccccc
      elseif(i_sym .eq. 2) then
        call x_symmetry(npt, ne, ne_0, npt_0, 
     *      max_ne, max_npt, idfn, nu, typ_el, 
     *      idfn_0, nu_0, typ_el_0, x, y, x_0, y_0,
     *      tab_ne, tab_npt)
cccccc
      elseif(i_sym .eq. 3) then
        call y_symmetry(npt, ne, ne_0, npt_0, 
     *      max_ne, max_npt, idfn, nu, typ_el, 
     *      idfn_0, nu_0, typ_el_0, x, y, x_0, y_0,
     *      tab_ne, tab_npt)

        npt_0 = npt
        ne_0 = ne
        do i=1,npt_0
          x_0(i) = x(i)
          y_0(i) = y(i)
          idfn_0(i) = idfn(i)
        enddo
c
        do i=1,ne_0
          do j=1,6
            nu_0(j,i) = nu(j,i)
          enddo
          typ_el_0(i) = typ_el(i)
        enddo

        call x_symmetry(npt, ne, ne_0, npt_0, 
     *      max_ne, max_npt, idfn, nu, typ_el, 
     *      idfn_0, nu_0, typ_el_0, x, y, x_0, y_0,
     *      tab_ne, tab_npt)

cccccc
      else
        return
      endif

ccccccccccc
c
c
      return
      end
c
c###############################################
c


c
c###############################################
c
      subroutine y_symmetry(npt, ne, ne_0, npt_0, 
     *      max_ne, max_npt, idfn, nu, typ_el, 
     *      idfn_0, nu_0, typ_el_0, x, y, x_0, y_0,
     *      tab_ne, tab_npt)
c
      implicit none
      integer max_ne, max_npt

      integer ne_0, npt_0, idfn_0(max_npt)
      integer nu_0(6,max_ne), typ_el_0(max_ne)

      integer ne, npt, idfn(max_npt)
      integer nu(6,max_ne), typ_el(max_ne)

      integer tab_ne(max_ne), tab_npt(max_npt,3)
      double precision x(max_npt), y(max_npt)
      double precision x_0(max_npt),  y_0(max_npt)
c
c     Local variables
      integer i, i1, i2, i_a, i_b, j, j1, j2, k
      integer ne_1, npt_1, npt_2
      integer nut_0(6), nut_a(6), nut_b(6)
      double precision tol, y_min, y_max, y_mid
      double precision x_a, y_a, x_b, y_b

      integer ui, debug
      character file_ui*100
      common/imp/ui, debug
      common/imp_file/file_ui
c
      y_min = y_0(1)
      y_max = y_0(1)
      do i=1,npt_0
        if(y_0(i) .lt. y_min) y_min = y_0(i)
        if(y_0(i) .gt. y_max) y_max = y_0(i)
      enddo

c     Selecting points in upper half of the domain
      tol = 1.0d-7
      y_mid = (y_min+y_max)/2.0d0
      i1 = 0
      i2 = 0
      do i=1,npt_0
        if(y_0(i) .ge. y_mid-tol) then
          i1 = i1 + 1
          tab_npt(i,1) = i1
          if(Abs(y_0(i)-y_mid) .le. tol) then
c           No duplication for a point on the symmetry line
            i2 = i2 + 1
            tab_npt(i,2) = 1
          else
            tab_npt(i,2) = 2
          endif
        else
          tab_npt(i,1) = 0
          tab_npt(i,2) = 0
        endif
      enddo
      npt_1 = i1
      npt_2 = i2
      npt = 2*i1 - i2
c
c     Selecting triangles in upper half of the domain
      i1 = 0
      do i=1,ne_0
        i2 = 0
        do j=1,3
          j1 = nu_0(j,i)
          j2 = tab_npt(j1,1)
          if(j2 .gt. 0) i2 = i2 + 1
        enddo
        if(i2 .eq. 3) then
          i1 = i1 + 1
          tab_ne(i) = i1
        else
          tab_ne(i) = 0
        endif
      enddo
      ne_1 = i1
      ne = 2*i1
c
c     Generating the symmetrized FEM mesh
c
      i_b = npt_1
      do i=npt_0,1,-1
        tab_npt(i,3) = 0
        x_a = x_0(i)
        y_a = y_0(i)

        x_b = x_a
        y_b = 2.0d0*y_mid - y_a

        i_a = tab_npt(i,1)
c        i_b = npt - i_a + 1
        i1 = tab_npt(i,2)

        if(i_a .gt. 0) then
          x(i_a) = x_a
          y(i_a) = y_a
          idfn(i_a) = idfn_0(i)
          if(i1 .eq. 2) then
            i_b = i_b + 1
            tab_npt(i,3) = i_b
            x(i_b) = x_b
            y(i_b) = y_b
            idfn(i_b) = idfn_0(i)
            if(idfn(i_b) .eq. 1) idfn(i_b) = 2
          endif
        endif
      enddo
c
ccccccccccccccccccccccc
c
      i_b = ne_1
      do i=ne_0,1,-1
        i_a = tab_ne(i)
        if(i_a .gt. 0) then
          i_b = i_b + 1
          do j=1,6
            nut_0(j) = nu_0(j,i)
          enddo

          do j=1,6
            j1 = tab_npt(nut_0(j),1)
            nut_a(j) = j1
            j2 = tab_npt(nut_0(j),2)

            if(j2 .eq. 1) then
              nut_b(j) = j1
            elseif(j2 .eq. 2) then
              nut_b(j) = tab_npt(nut_0(j),3)
c              nut_b(j) = npt - j1 + 1
            else
              open (unit=ui,file=file_ui)
              write(*,*) 'SYMMETRY: tab_npt(i,2) = ', j2
              write(*,*) 'i, tab_npt(i,1) = ', nut_0(j), j1
              stop
            endif
          enddo
          do j=1,6
            nu(j,i_a) = nut_a(j)
          enddo
          typ_el(i_a) = typ_el_0(i)

c          i_b = ne - i_a + 1
          if(i_a .gt. ne/2) stop
          do j=1,3
            nu(j,i_b) = nut_b(3-j+1)
            nu(j+3,i_b) = nut_b(6-j+1)
          enddo
c
c       Symmetry reverses the orientation
c       so we must reverse the numbering to get the positive orientation

          nu(1,i_b) = nut_b(1)
          nu(2,i_b) = nut_b(3)
          nu(3,i_b) = nut_b(2)
          nu(4,i_b) = nut_b(6)
          nu(5,i_b) = nut_b(5)
          nu(6,i_b) = nut_b(4)

          typ_el(i_b) = typ_el_0(i)
        endif
      enddo
c
      return
      end
c
c###############################################
c

c
c###############################################
c
      subroutine x_symmetry(npt, ne, ne_0, npt_0, 
     *      max_ne, max_npt, idfn, nu, typ_el, 
     *      idfn_0, nu_0, typ_el_0, x, y, x_0, y_0,
     *      tab_ne, tab_npt)
c
      implicit none
      integer max_ne, max_npt

      integer ne_0, npt_0, idfn_0(max_npt)
      integer nu_0(6,max_ne), typ_el_0(max_ne)

      integer ne, npt, idfn(max_npt)
      integer nu(6,max_ne), typ_el(max_ne)

      integer tab_ne(max_ne), tab_npt(max_npt,3)
      double precision x(max_npt), y(max_npt)
      double precision x_0(max_npt),  y_0(max_npt)
c
c     Local variables
      integer i, i1, i2, i_a, i_b, j, j1, j2, k
      integer ne_1, npt_1, npt_2
      integer nut_0(6), nut_a(6), nut_b(6)
      double precision tol, x_min, x_max, x_mid
      double precision x_a, y_a, x_b, y_b

      integer ui, debug
      character file_ui*100
      common/imp/ui, debug
      common/imp_file/file_ui
c
      x_min = x_0(1)
      x_max = x_0(1)
      do i=1,npt_0
        if(x_0(i) .lt. x_min) x_min = x_0(i)
        if(x_0(i) .gt. x_max) x_max = x_0(i)
      enddo

c     Selecting points in upper half of the domain
      tol = 1.0d-7
      x_mid = (x_min+x_max)/2.0d0
      i1 = 0
      i2 = 0
      do i=1,npt_0
        if(x_0(i) .le. x_mid+tol) then
          i1 = i1 + 1
          tab_npt(i,1) = i1
          if(Abs(x_0(i)-x_mid) .le. tol) then
c           No duplication for a point on the symmetry line
            i2 = i2 + 1
            tab_npt(i,2) = 1
          else
            tab_npt(i,2) = 2
          endif
        else
          tab_npt(i,1) = 0
          tab_npt(i,2) = 0
        endif
      enddo
      npt_1 = i1
      npt_2 = i2
      npt = 2*i1 - i2
c
c     Selecting triangles in left half of the domain
      i1 = 0
      do i=1,ne_0
        i2 = 0
        do j=1,3
          j1 = nu_0(j,i)
          j2 = tab_npt(j1,1)
          if(j2 .gt. 0) i2 = i2 + 1
        enddo
        if(i2 .eq. 3) then
          i1 = i1 + 1
          tab_ne(i) = i1
        else
          tab_ne(i) = 0
        endif
      enddo
      ne_1 = i1
      ne = 2*i1
c
c     Generating the symmetrized FEM mesh
c
      i_b = npt_1
      do i=npt_0,1,-1
        tab_npt(i,3) = 0
        x_a = x_0(i)
        y_a = y_0(i)

        x_b = 2.0d0*x_mid - x_a
        y_b = y_a

        i_a = tab_npt(i,1)
        i1 = tab_npt(i,2)

        if(i_a .gt. 0) then
          x(i_a) = x_a
          y(i_a) = y_a
          idfn(i_a) = idfn_0(i)
          if(i1 .eq. 2) then
            i_b = i_b + 1
            tab_npt(i,3) = i_b
            x(i_b) = x_b
            y(i_b) = y_b
            idfn(i_b) = idfn_0(i)
            if(idfn(i_b) .eq. 3) idfn(i_b) = 4
          endif
        endif
      enddo
c
ccccccccccccccccccccccc
c
      i_b = ne_1
      do i=ne_0,1,-1
        i_a = tab_ne(i)
        if(i_a .gt. 0) then
          i_b = i_b + 1
          do j=1,6
            nut_0(j) = nu_0(j,i)
          enddo

          do j=1,6
            j1 = tab_npt(nut_0(j),1)
            nut_a(j) = j1
            j2 = tab_npt(nut_0(j),2)

            if(j2 .eq. 1) then
              nut_b(j) = j1
            elseif(j2 .eq. 2) then
              nut_b(j) = tab_npt(nut_0(j),3)
c              nut_b(j) = npt - j1 + 1
            else
              open (unit=ui,file=file_ui)
              write(*,*) 'SYMMETRY_X: tab_npt(i,2) = ', j2
              write(*,*) 'i, tab_npt(i,1) = ', nut_0(j), j1
              stop
            endif
          enddo
          do j=1,6
            nu(j,i_a) = nut_a(j)
          enddo
          typ_el(i_a) = typ_el_0(i)

c          i_b = ne - i_a + 1
          if(i_a .gt. ne/2) stop
          do j=1,3
            nu(j,i_b) = nut_b(3-j+1)
            nu(j+3,i_b) = nut_b(6-j+1)
          enddo
c
c       Symmetry reverses the orientation
c       so we must reverse the numbering to get the positive orientation

          nu(1,i_b) = nut_b(1)
          nu(2,i_b) = nut_b(3)
          nu(3,i_b) = nut_b(2)
          nu(4,i_b) = nut_b(6)
          nu(5,i_b) = nut_b(5)
          nu(6,i_b) = nut_b(4)

          typ_el(i_b) = typ_el_0(i)
        endif
      enddo
c

      return
      end
c
c###############################################
c
