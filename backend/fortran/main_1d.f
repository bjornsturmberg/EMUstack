      program main_1d

C************************************************************************
C
C  Program:
C    Finite Element Method for 1D array of planear waveguides
C
C************************************************************************
C

      implicit none
C  Local parameters:
      integer*8 nval, npt, nel, nb_typ_el
      double precision lambda

c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

C

      call calc_modes_1d(lambda, nval, npt, nel, nb_typ_el)

C
      stop
      end 

