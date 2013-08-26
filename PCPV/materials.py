"""
    materials.py is a subroutine of PCPV that defines Material objects,
    these represent dispersive lossy refractive indices and possess 
    methods to interpolate n from tabulated data.

    Copyright (C) 2013  Bjorn Sturmberg

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import numpy as np
from scipy.interpolate import UnivariateSpline, interp1d

data_location = '../PCPV/Data/'

class Material(object):
    """ Represents a material with a refractive index n.

        If the material is dispersive, the refractive index at a given
        wavelength is calculated by linear interpolation from the 
        initially given data `n`.

        INPUTS:

        - `n`: Either a scalar refractive index, or an array
                of values `(wavelength, n)`, or 
                `(wavelength, real(n), imag(n))`.
    """
    def __init__(self, n):
        if () == np.shape(n):
            # n is a scalar, the medium is non-dispersive.
            self._n = lambda x: n
            self.data_wls = None
            self.data_ns = n
        else:
            self.data_wls = n[:,0]
            if len(n[0]) == 2:
                # n is an array of wavelengths and (possibly-complex)
                # refractive indices.
                self.data_ns = n[:,1]
            elif len(n[0]) == 3:
                self.data_ns = n[:,1] + 1j * n[:,2]
            else:
                raise ValueError
            # Do cubic interpolation if we get the chance
            # if len(self.data_wls) > 3:
            #     self._n = interp1d(self.data_wls, self.data_ns, 'cubic')
            # else:
            self._n = interp1d(self.data_wls, self.data_ns)

    def n(self, wl_nm):
        """ Return n for the specified wavelength."""
        return self._n(wl_nm)

    def n_drude(self, wavelength, plasma_wl, gamma_wl):
        # http://www.wave-scattering.com/drudefit.html
        c = 2.98e8*1e9
        numer = c*wavelength**2
        denom = c*plasma_wl**2 - 1j*2*np.pi*gamma_wl*wavelength*plasma_wl**2
        drude_n = 1 - numer/denom
        return np.sqrt(drude_n)

    def __getstate__(self):
        """ Can't pickle self._n, so remove it from what is pickled."""
        d = self.__dict__.copy()
        d.pop('_n')
        return d

    def __setstate__(self, d):
        """ Recreate self._n when unpickling."""
        self.__dict__ = d
        if None == self.data_wls:
            self._n = lambda x: self.data_ns
        else:
            self._n = interp1d(self.data_wls, self.data_ns)


Air    = Material(np.loadtxt('%sAir.txt'% data_location))
# Transparent oxides
TiO2   = Material(np.loadtxt('%sTiO2.txt'% data_location))
ITO    = Material(np.loadtxt('%sITO.txt'% data_location))           #Filmetrics.com
# Semiconductors
Si_c   = Material(np.loadtxt('%sSi_c.txt'% data_location))
Si_a   = Material(np.loadtxt('%sSi_a.txt'% data_location))
SiO2_a = Material(np.loadtxt('%sSiO2_a.txt'% data_location))
CuO    = Material(np.loadtxt('%sCuO.txt'% data_location))
CdTe   = Material(np.loadtxt('%sCdTe.txt'% data_location))
FeS2   = Material(np.loadtxt('%sFeS2.txt'% data_location))
Zn3P2  = Material(np.loadtxt('%sZn3P2.txt'% data_location))
Sb2S3  = Material(np.loadtxt('%sSb2S3.txt'% data_location))
AlGaAs = Material(np.loadtxt('%sAlGaAs.txt'% data_location))
GaAs   = Material(np.loadtxt('%sGaAs.txt'% data_location))
Si3N4  = Material(np.loadtxt('%sSi3N4.txt'% data_location))
InP    = Material(np.loadtxt('%sInP.txt'% data_location))
InAs   = Material(np.loadtxt('%sInAs.txt'% data_location))          #Filmetrics.com
GaP    = Material(np.loadtxt('%sGaP.txt'% data_location))           #Filmetrics.com
# Metals
Au     = Material(np.loadtxt('%sAu.txt'% data_location))
Ag     = Material(np.loadtxt('%sAg.txt'% data_location))


# import matplotlib
# matplotlib.use('pdf')
# import matplotlib.pyplot as plt
# def n_plot(spectra_name, wavelengths, data):
#     fig = plt.figure(num=None, figsize=(9, 12), dpi=80, facecolor='w', edgecolor='k')
#     ax1 = fig.add_subplot(2,1,1)
#     ax1.plot(wavelengths, np.real(data))
#     ax1.set_xlabel('Wavelength (nm)')
#     ax1.set_ylabel('n')    
#     ax1 = fig.add_subplot(2,1,2)
#     ax1.plot(wavelengths, np.imag(data))
#     ax1.set_xlabel('Wavelength (nm)')
#     ax1.set_ylabel('k')
#     # plt.axis([wavelengths[0], wavelengths[-1], 0, 1])
#     plt.suptitle(spectra_name)
#     plt.savefig(spectra_name)

