c
C  Calculate  Int(unit cell) ||E||^2 and Int(over cylinder) epsilon ||E||^2
C  mode_pol(1) contains Int(unit cell) epsilon |E_x|^2 / Int(unit cell) epsilon ||E||^2
C  mode_pol(2) contains Int(unit cell) epsilon |E_y|^2 / Int(unit cell) epsilon ||E||^2
C  mode_pol(3) contains Int(unit cell) epsilon |E_z|^2 / Int(unit cell) epsilon ||E||^2
C  mode_pol(4) contains Int(over cylinder) epsilon ||E||^2 / Int(unit cell) epsilon ||E||^2
c  Note that we use E_z = i * \hat{E}_z * (i*beta) (because of the change of variable)
c  It is assumed that the type number triangle in the cylinders is : 
c                           typ_e=n_core(1) or typ_e=n_core(2)
c
      subroutine mode_energy (nval, nel, npt, n_ddl, nnodes, 
     *     n_core, table_nod, type_el, nb_typ_el, eps_eff, 
     *     x, sol, beta1, mode_pol)

      implicit none
      integer*8 nval, nel, npt, n_ddl, nnodes
      integer*8 nb_typ_el, n_core(2), type_el(nel)
      integer*8 table_nod(nnodes,nel)
      complex*16 x(2,npt), sol(3,nnodes+7,nval,nel)
      complex*16 eps_eff(nb_typ_el), mode_pol(4,nval)
      complex*16 beta1(nval)
c
C  local parameters - purely internal
c
C  variables for quadrature interpolation
      integer*8 nquad, nquad_max, iq
      parameter (nquad_max = 25)
      double precision wq(nquad_max)
      double precision xq(nquad_max), yq(nquad_max)
      double precision xx(2), xx_g(2), ww, det
      double precision mat_B(2,2), mat_T(2,2)

      integer*8 nnodes_0
      parameter (nnodes_0 = 6)
      integer*8 nod_el_p(nnodes_0)
      double precision xel(2,nnodes_0)
      complex*16 vec_phi(3)

      double precision phi2_list(6)

      integer*8 i, j, iel, ival, ival2, typ_e
      integer*8 inode, global, trans
      integer*8 debug, ui

      complex*16 ii, z_tmp, coeff_1
C
CCCCCCCCCCC Start Program CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      ui = 6
      debug = 1
c
      if ( nnodes .ne. 6 ) then
        write(ui,*) "overlap_J: problem nnodes = ", nnodes
        write(ui,*) "overlap_J: nnodes should be equal to 6!"
        write(ui,*) "overlap_J: Aborting..."
        stop
      endif 
c
      do ival=1,nval
        do j=1,4
          mode_pol(j,ival) = 0.0d0
        enddo
      enddo
c
      call quad_triangle(nquad,nquad_max,wq,xq,yq)
C
CCCCCCCCCCCCCCCCC	loop over all elements	CCCCCCCCCCCCCCCC
C
      do iel=1,nel
        typ_e = type_el(iel)
        do inode=1,nnodes
          global = table_nod(inode,iel)
          nod_el_p(inode) = global
          xel(1,inode) = x(1,global)
          xel(2,inode) = x(2,global)
        enddo
C
        do iq=1,nquad
          xx(1) = xq(iq)
          xx(2) = yq(iq)
          ww = wq(iq)
          call phi2_2d_mat_J(xx, phi2_list)
          call jacobian_p1_2d(xx, xel, nnodes, 
     *               xx_g, det, mat_B, mat_T)
c
          coeff_1 = ww * ABS(det) * eps_eff(typ_e)
          do ival=1,nval
            do trans=1,3
              vec_phi(trans) = 0.0d0
            enddo
            do inode=1,nnodes
              do trans=1,3		! transverse field components
                vec_phi(trans) = vec_phi(trans) + 
     *            sol(trans,inode,ival,iel) * phi2_list(inode)
              enddo
            enddo
            vec_phi(3) = vec_phi(3) * beta1(ival)  ! E_z = i * \hat{E}_z * (i*beta) (because of the change of variable)
            do trans=1,3
              z_tmp = coeff_1 * abs(vec_phi(trans))**2
              mode_pol(trans,ival) = mode_pol(trans,ival) + z_tmp
            enddo
            if (typ_e .eq. n_core(1) .or. typ_e .eq. n_core(2)) then
              do trans=1,3
                z_tmp = coeff_1 * abs(vec_phi(trans))**2
                mode_pol(4,ival) = mode_pol(4,ival) + z_tmp
              enddo
            endif
          enddo
        enddo
      enddo
c
c       Total energy and normalization
      do ival=1,nval
        z_tmp = mode_pol(1,ival) + mode_pol(2,ival) 
     *        + mode_pol(3,ival)
        if (abs(z_tmp) .lt. 1.0d-10) then
          write(*,*) "mode_energy: the total energy ",
     *       "is too small : ", z_tmp
          write(*,*) "mode_energy: ival = ", ival
          write(*,*) "mode_energy: zero eigenvector; aborting..."
          stop
        endif
        do j=1,4
          mode_pol(j,ival) = mode_pol(j,ival) / z_tmp
        enddo
      enddo
c
      return
      end
