C  
C 
      subroutine cmp_ref_ind_lambda (hole_wire, Lambda_Res, eps_eff,
     *    n_eff, Complex_refract1, Complex_refract2,
     *    Complex_refract3, Complex_refract4,
     *    i_lambda, substrate)
c
      implicit none
C
      integer*8 hole_wire, Lambda_Res, i_lambda, max_typ_el
      integer*8 substrate
      parameter (max_typ_el=10)
      complex*16 eps_eff(max_typ_el), n_eff(max_typ_el)
      complex*16 Complex_refract1(818), Complex_refract2(1635)
      complex*16 Complex_refract3(818), Complex_refract4(1635)
C
CCCCCCCCCCCC Update for each wavelength CCCCCCCCCCCCCCCCC
C
C substrate
      if (substrate .eq. 1) then
        if (Lambda_Res .eq. 1) then
          eps_eff(2) = Complex_refract3(i_lambda)**2
          n_eff(2) = Complex_refract3(i_lambda)
        elseif (Lambda_Res .eq. 2) then
          eps_eff(2) = Complex_refract4(i_lambda)**2
          n_eff(2) = Complex_refract4(i_lambda)
        endif
      endif
C
C  air holes in background
      if (hole_wire .eq. 0) then
        if (Lambda_Res .eq. 1) then
          eps_eff(3) = Complex_refract1(i_lambda)**2
          n_eff(3) = Complex_refract1(i_lambda)
        elseif (Lambda_Res .eq. 2) then
          eps_eff(3) = Complex_refract2(i_lambda)**2
          n_eff(3) = Complex_refract2(i_lambda)
        endif
C
C  nanowires in air background
      elseif (hole_wire .eq. 1) then
        if (Lambda_Res .eq. 1) then
          eps_eff(4) = Complex_refract1(i_lambda)**2
          n_eff(4) = Complex_refract1(i_lambda)
        elseif (Lambda_Res .eq. 2) then
          eps_eff(4) = Complex_refract2(i_lambda)**2
          n_eff(4) = Complex_refract2(i_lambda)
        endif
C
C  thin film (uniform slab)
      elseif (hole_wire .eq. 2) then
        if (Lambda_Res .eq. 1) then
          eps_eff(3) = Complex_refract1(i_lambda)**2
          n_eff(3) = Complex_refract1(i_lambda)
          eps_eff(4) = Complex_refract1(i_lambda)**2
          n_eff(4) = Complex_refract1(i_lambda)
        elseif (Lambda_Res .eq. 2) then
          eps_eff(3) = Complex_refract2(i_lambda)**2
          n_eff(3) = Complex_refract2(i_lambda)
          eps_eff(4) = Complex_refract2(i_lambda)**2
          n_eff(4) = Complex_refract2(i_lambda)
        endif
C
C  constant refractive index
      elseif (hole_wire .eq. 3) then
        write(*,*) "cmp_ref_ind_lambda: ",
     * "using constant refractive index!"
      else
        write(*,*) "cmp_ref_ind_lambda: ",
     * "need to chose either nanowires, holes"
        write(*,*) "cmp_ref_ind_lambda: or n constant"
        stop
      endif
C      
C      write(*,*) 'Re^2', Complex_refract1(i_lambda)**2
C      
      return
      end 
