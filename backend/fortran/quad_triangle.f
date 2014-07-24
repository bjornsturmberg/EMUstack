c
c***********************************************************************
c
      subroutine quad_triangle (nquad, nquad_max, wq, xq, yq)
c
c
c***********************************************************************
c
c                      Integration numerique
c                      ---------------------
c       Evalue les integrales elementaires des composantes convectives
c       sur chaque triangle. on utilise ici la methode de hammer a
c       seize points de gauss qui integre exactement des polynomes du
c       huitieme degre.
c
c       Reference
c       J. N. Lyness and D. Jespersen
c       "Moderate Degree Symmetric Quadrature Rules for the Triangle"
c       J. Inst. Math. Appl., 1975, 15(1), pp. 19-32
c
c       "J. Inst. Math. Appl." is now Continued as "IMA J. Appl. Math."
c       J. Inst. Math. Appl. = Journal of the Institute of Mathematics and its Applications
c       IMA J. Appl. Math.   = IMA Journal of Applied Mathematics
c
c***********************************************************************
c
      implicit none
      integer*8 nquad, nquad_max
      double precision wq(nquad_max)
      double precision xq(nquad_max), yq(nquad_max)

c     Local variables
      integer*8 i
      double precision poidbar, coorbar
      double precision poid1, coor1grp1, coor2grp1
      double precision poid2, coor1grp2, coor2grp2
      double precision poid3, coor1grp3, coor2grp3
      double precision poid4, coor1grp4, coor2grp4, coor3grp4
c
c     Check the number of quadrature points
c
      nquad = 16
      if (nquad .gt. nquad_max) then
        write(*,*)
        write(*,*) "   ???"
        write(*,*) "quad_triangle: nquad > nquad_max : ",
     *  nquad, nquad_max
        write(*,*) "quad_triangle: Aborting..."
        stop
      endif
c

c
c_____________________________________________________________________
c
      poidbar = 1.443156076777862d-1 / 2.d0
      poid1   = 2.852749028018549d-1 / 6.d0
      poid2   = 9.737549286959440d-2 / 6.d0
      poid3   = 3.096521116041552d-1 / 6.d0
      poid4   = 1.633818850466092d-1 / 1.2d1
c
c_____________________________________________________________________
      coorbar = 1.d0 / 3.d0
c
      coor1grp1 = 4.592925882927229d-1
      coor2grp1 = 8.141482341455413d-2
c
      coor1grp2 = 5.054722831703103d-2
      coor2grp2 = 8.989055433659379d-1
c
      coor1grp3 = 1.705693077517601d-1
      coor2grp3 = 6.588613844964797d-1
c
      coor1grp4 = 7.284923929554041d-1
      coor2grp4 = 2.63112829634638689d-1
      coor3grp4 = 8.394777409957211d-3
c_____________________________________________________________________
      i = 1
        xq(i) = coorbar
        yq(i) = coorbar
        wq(i) = poidbar
c_____________________________________________________________________
      i = 2
        xq(i) = coor1grp1
        yq(i) = coor1grp1
        wq(i) = poid1
      i = 3
        xq(i) = coor1grp1
        yq(i) = coor2grp1
        wq(i) = poid1
      i = 4
        xq(i) = coor2grp1
        yq(i) = coor1grp1
        wq(i) = poid1
c_____________________________________________________________________
      i = 5
        xq(i) = coor1grp2
        yq(i) = coor1grp2
        wq(i) = poid2
      i = 6
        xq(i) = coor1grp2
        yq(i) = coor2grp2
        wq(i) = poid2
      i = 7
        xq(i) = coor2grp2
        yq(i) = coor1grp2
        wq(i) = poid2
c_____________________________________________________________________
      i = 8
        xq(i) = coor1grp3
        yq(i) = coor1grp3
        wq(i) = poid3
      i = 9
        xq(i) = coor1grp3
        yq(i) = coor2grp3
        wq(i) = poid3
      i = 10
        xq(i) = coor2grp3
        yq(i) = coor1grp3
        wq(i) = poid3
c_____________________________________________________________________
      i = 11
        xq(i) = coor1grp4
        yq(i) = coor2grp4
        wq(i) = poid4
      i = 12
        xq(i) = coor2grp4
        yq(i) = coor1grp4
        wq(i) = poid4
      i = 13
        xq(i) = coor2grp4
        yq(i) = coor3grp4
        wq(i) = poid4
      i = 14
        xq(i) = coor3grp4
        yq(i) = coor2grp4
        wq(i) = poid4
      i = 15
        xq(i) = coor3grp4
        yq(i) = coor1grp4
        wq(i) = poid4
      i = 16
        xq(i) = coor1grp4
        yq(i) = coor3grp4
        wq(i) = poid4
c_____________________________________________________________________
c
      return
      end
