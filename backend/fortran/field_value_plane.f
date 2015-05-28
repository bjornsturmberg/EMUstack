c  Write an ASCII file that can be read into GMSH.
c  P1-element is used to represent the  3D vector field.
c  P2-element is used to represent each component of the 3D vector field.
c
      subroutine field_value_plane (nval,
     *     nel, npt, nnodes, nb_typ_el, table_nod, type_el,
     *     eps_eff, x, beta, sol, vec_coef, h, hz,m_E)
c
      implicit none
      integer*8 nval, nel, npt, nnodes
      integer*8 nb_typ_el
      double precision h, hz
      integer*8 table_nod(nnodes,nel), type_el(nel)
      double precision x(2,npt)
      complex*16 eps_eff(nb_typ_el)
      complex*16 sol(3,nnodes+7,nval,nel)
      complex*16 vec_coef(2*nval)
      complex*16 beta(nval)
      complex*16 m_E(nel,nnodes,3)

      !Local variables

      integer*8 nnodes_0
      parameter (nnodes_0 = 6)

      character*20 dir_name
      integer alloc_stat
      double precision xel(3,nnodes_0), xel_p1(3,3)
      complex*16 sol_el(3,nnodes_0), sol_max(4)
      double precision sol_el_abs2(nnodes_0)
      double precision ls_index(nnodes_0), r_index, zz

      complex*16 P_down, P_up, coef_down, coef_up, coef_t, coef_z

      integer*8 i, j, i1, iel, ival, namelen, typ_e
      integer*8 plot_imag, plot_real, plot_abs
      integer*8 debug, ui
      complex*16 z_tmp1
      complex*16 ii
      character*100 tchar
      character*1 tE_H
      integer*8 namelength


Cf2py intent(in) nval, nel, npt, nnodes nb_typ_el, table_nod,
Cf2py intent(in) type_el, eps_eff, x, beta, sol, vec_coef, h, hz
Cf2py intent(out) m_E

Cf2py depend(table_nod) nnodes, nel
Cf2py depend(type_el) nel
Cf2py depend(x) npt
Cf2py depend(eps_eff) nb_typ_el
Cf2py depend(sol) nnodes, nval, nel
Cf2py depend(vec_coef) nval
Cf2py depend(beta) nval

      ! imaginary unit
      ii = cmplx(0.0d0, 1.0d0)

      ! echoing the inputs
    !   WRITE(*,*) "nval",nval
    !   WRITE(*,*) "nel",nel
    !   WRITE(*,*) "npt",npt
    !   WRITE(*,*) "nnodes",nnodes
    !   WRITE(*,*) "nb_typ_el",nb_typ_el
    !   WRITE(*,*) "table_nod(1,1)",table_nod(1,1)
    !   WRITE(*,*) "type_el(1)",type_el(1)
    !   WRITE(*,*) "eps_eff(1)",eps_eff(1)
    !   WRITE(*,*) "x(1,1)",x(1,1)
    !   WRITE(*,*) "beta(1)",beta(1)
    !   WRITE(*,*) "sol(1,1,1,1)",sol(1,1,1,1)
    !   WRITE(*,*) "vec_coef(1)",vec_coef(1)
    !   WRITE(*,*) "h",h
    !   WRITE(*,*) "hz",hz
    !   WRITE(*,*)

      ! checking if mesh nodes are ok
      ui = 6
      debug = 0
      if ( nnodes .ne. 6 ) then
          write(ui,*) "gmsh_plot_field: problem nnodes = ", nnodes
          write(ui,*) "gmsh_plot_field: nnodes should be equal to 6 !"
          write(ui,*) "gmsh_plot_field: Aborting..."
          stop
      endif

      ! Is the z values ok? If not aborting
      if (hz .lt. 0 .or. hz .gt. h ) then
          write(ui,*)
          write(ui,*) "gmsh_plot_field: invalid value for hz"
          write(ui,*) "hz should be in the interval [0,h]"
          write(ui,*) "hz, h = ", hz, h
          write(ui,*) "gmsh_plot_field: Aborting..."
          stop
      endif

      ! Initializing the array contaning the field
      do iel=1,nel
          do i=1,nnodes
              do j=1,3
                  m_E(iel,i,j) = 0.0d0
              enddo
          enddo
      enddo

      ! Stacking all the block modes to get the whole field
      do ival=1,nval

          ! Preparing the propagation in z
          P_down = EXP(ii*beta(ival)*hz)   !  Introduce Propagation in -z
          P_up = EXP(ii*beta(ival)*(h-hz)) !  Introduce Propagation in +z
          coef_down = vec_coef(ival) * P_down
          coef_up = vec_coef(ival+nval) * P_up

          coef_t = coef_up + coef_down
          coef_z = (coef_up - coef_down)/beta(ival) ! Taking into accout the change of variable for Ez

          ! Looping over elements, and nodes to compute the total field
          do iel=1,nel
              do i=1,nnodes
                  do j=1,2
                      z_tmp1 = sol(j,i,ival,iel) * coef_t
                      m_E(iel,i,j) = m_E(iel,i,j) + z_tmp1
                  enddo
                  j=3
                  z_tmp1 = sol(j,i,ival,iel) * coef_z
                  m_E(iel,i,j) = m_E(iel,i,j) + z_tmp1
              enddo
          enddo
      enddo
      return
      end
