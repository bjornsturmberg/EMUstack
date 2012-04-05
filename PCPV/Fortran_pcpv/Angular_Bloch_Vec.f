C  Calculate the bloch vector components for a given wavelength
C     and incident angles theta (polar) and phi (azimuthal)
C
      subroutine Angular_Bloch_Vec(theta, phi, bloch_vec0, k_0)
C
      implicit none
C
      double precision theta, phi, bloch_vec0(2), k_0
      double precision eps_eff_0, pi, common_factor
C  constants - free space refractive index and pi
      pi = 3.141592653589793d0
C  calculate bloch vectors - for derivation see Bjorns Theory Notes
C      common_factor = sin(theta*pi/180.0d0)*2.0d0*pi*eps_eff_0/lambda
      common_factor = sin(theta*pi/180.0d0)*k_0
      bloch_vec0(1) = common_factor*cos(phi*pi/180.0d0)
      bloch_vec0(2) = common_factor*sin(phi*pi/180.0d0)
C  make sure not both are equal zero to avoid FEM problem
      if (abs(bloch_vec0(1)) + abs(bloch_vec0(2))
     *           .lt. 0.000000001d0) then
        bloch_vec0(1) = 0.000000001d0
      endif
C
      return
      end 