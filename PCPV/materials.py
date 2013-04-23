# doc

import numpy as np
from scipy.interpolate import UnivariateSpline
# import copy



data_location = '../PCPV/Data/'

class Material(object):
    """docstring for """
    def __init__(self, n_data):

    # if n_data.isdigit() == True:
    #     print 'hi babes'
    #     # n_inc = dict(zip(wavelengths, inclusion_a))
    # else:
        self.data_wls   = n_data[:,0]
        self.data_re_ns = n_data[:,1]
        self.data_im_ns = n_data[:,2]
        self.stored_ns  = {}

    def n_spline(self, wavelengths):
        if wavelengths.min() < self.data_wls.min():
            raise ValueError, "Input wavelength: %(in)f smaller than \
             data min of %(d)f" % {
         'in' : wavelengths.min(), 'd' : self.data_wls.min()
            }
        if wavelengths.max() > self.data_wls.max():
            raise ValueError, "Input wavelength: %(in)f larger than \
             data max of %(d)f" % {
         'in' : wavelengths.max(), 'd' : self.data_wls.max()
            }
        # for i in range(len(self.data_wls)-1):
        #     if self.data_wls[i] == self.data_wls[i+1]:
        #         print self.data_wls[i]
        spline_fit  = UnivariateSpline(self.data_wls, self.data_re_ns, s = 1)
        interp_re_n = spline_fit(wavelengths)
        spline_fit  = UnivariateSpline(self.data_wls, self.data_im_ns, s = 1)
        interp_im_n = spline_fit(wavelengths)
        self.interp_data = interp_re_n + 1j*interp_im_n
        # for single wavelength simo turn self.interp_data into array for zipping
        if isinstance(interp_re_n, float): 
            self.interp_data = np.array([interp_re_n + 1j*interp_im_n])
        # print repr(self.interp_data)
        self.stored_ns   = dict(zip(wavelengths, self.interp_data))
        # print self.stored_ns

    def n_interp(self, wavelengths):
        if wavelengths.min() < self.data_wls.min():
            raise ValueError, "Input wavelength: %(in)f smaller than \
             data min of %(d)f" % {
         'in' : wavelengths.min(), 'd' : self.data_wls.min()
            }
        if wavelengths.max() > self.data_wls.max():
            raise ValueError, "Input wavelength: %(in)f larger than \
             data max of %(d)f" % {
         'in' : wavelengths.max(), 'd' : self.data_wls.max()
            }
        interp_re_n = np.interp(wavelengths, 
            self.data_wls, self.data_re_ns)
        interp_im_n = np.interp(wavelengths, 
            self.data_wls, self.data_im_ns)
        self.interp_data = interp_re_n + 1j*interp_im_n
        self.stored_ns   = dict(zip(wavelengths, self.interp_data))
        # print self.stored_ns


    def n_drude(self, wavelengths, plasma_wl, gamma_wl):
        # http://www.wave-scattering.com/drudefit.html
        c = 2.98e8*1e9
        numer = c*wavelengths**2
        denom = c*plasma_wl**2 - 1j*2*np.pi*gamma_wl*wavelengths*plasma_wl**2
        drude_n = 1 - numer/denom
        self.interp_data = np.sqrt(drude_n)
        self.stored_ns   = dict(zip(wavelengths, self.interp_data))


    def n(self, wl):
        """ Return n at the specified wavelength."""
        # if self.stored_ns.has_key(wl):
        return self.stored_ns[wl]
        # else:
        #     except KeyError:
        #     raise  NotImplementedError
        #     "Cannot interpolate for individual wavelengths"

    def max_n(self):
        """ Return maximum of Re(n) over the interpolation range."""
        return self.interp_data.real.max()

Air    = Material(np.loadtxt('%sAir.txt'% data_location))
Si_c   = Material(np.loadtxt('%sSi_c.txt'% data_location))
Si_a   = Material(np.loadtxt('%sSi_a.txt'% data_location))
SiO2_a = Material(np.loadtxt('%sSiO2_a.txt'% data_location))
CuO    = Material(np.loadtxt('%sCuO.txt'% data_location))
CdTe   = Material(np.loadtxt('%sCdTe.txt'% data_location))
FeS2   = Material(np.loadtxt('%sFeS2.txt'% data_location))
Zn3P2  = Material(np.loadtxt('%sZn3P2.txt'% data_location))
Au     = Material(np.loadtxt('%sAu.txt'% data_location))
Ag     = Material(np.loadtxt('%sAg.txt'% data_location))
Sb2S3  = Material(np.loadtxt('%sSb2S3.txt'% data_location))
AlGaAs = Material(np.loadtxt('%sAlGaAs.txt'% data_location))
GaAs   = Material(np.loadtxt('%sGaAs.txt'% data_location))
Si3N4  = Material(np.loadtxt('%sSi3N4.txt'% data_location))
TiO2   = Material(np.loadtxt('%sTiO2.txt'% data_location))
InP    = Material(np.loadtxt('%sInP.txt'% data_location))


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

def interp_needed(wavelengths, material=Air):
    if material == Air:
        # Air.n_spline(wavelengths)
        Air.n_interp(wavelengths)
        # n_plot('Air',wavelengths, Air.interp_data)
    if material == Si_c:
        # Si_c.n_spline(wavelengths)
        Si_c.n_interp(wavelengths)
        # n_plot('Si_c',wavelengths, Si_c.interp_data)
    if material == Si_a:
        Si_a.n_spline(wavelengths)
        # Si_a.n_interp(wavelengths)
        # n_plot('Si_a',wavelengths, Si_a.interp_data)
    if material == SiO2_a:
        # SiO2_a.n_spline(wavelengths)
        SiO2_a.n_interp(wavelengths)
        # n_plot('SiO2_a',wavelengths, SiO2_a.interp_data)
    if material == CuO:
        CuO.n_spline(wavelengths)#.n_interp(wavelengths)
        # n_plot('CuO',wavelengths, CuO.interp_data)
    if material == CdTe:
        CdTe.n_spline(wavelengths)#.n_interp(wavelengths)
        # n_plot('CdTe',wavelengths, CdTe.interp_data)
    if material == FeS2:
        FeS2.n_spline(wavelengths)#.n_interp(wavelengths)
        # n_plot('FeS2',wavelengths, FeS2.interp_data)
    if material == Zn3P2:
        Zn3P2.n_spline(wavelengths)
        # Zn3P2.n_interp(wavelengths)
        # n_plot('Zn3P2',wavelengths, Zn3P2.interp_data)
    if material == Au:
        # Au.n_spline(wavelengths)
        Au.n_interp(wavelengths)
        # plasma_wl_in_nm = 137.36e1*np.ones(len(wavelengths))
        # gamma_for_wl_nm = 0#215*1e-7
        # Au.n_drude(wavelengths, plasma_wl_in_nm, gamma_for_wl_nm)
        # n_plot('Au',wavelengths, Au.interp_data)
    if material == Ag:
        # Ag.n_spline(wavelengths)
        Ag.n_interp(wavelengths)
        # n_plot('Ag',wavelengths, Ag.interp_data)
    if material == Sb2S3:
        Sb2S3.n_spline(wavelengths)#.n_interp(wavelengths)
        # n_plot('Sb2S3',wavelengths, Sb2S3.interp_data)
    if material == AlGaAs:
        # AlGaAs.n_spline(wavelengths)
        AlGaAs.n_interp(wavelengths)
        # n_plot('AlGaAs',wavelengths, AlGaAs.interp_data)
    if material == GaAs:
        # GaAs.n_spline(wavelengths)
        GaAs.n_interp(wavelengths)
        # n_plot('GaAs',wavelengths, GaAs.interp_data)
    if material == Si3N4:
        # Si3N4.n_spline(wavelengths)
        Si3N4.n_interp(wavelengths)
        # n_plot('Si3N4',wavelengths, Si3N4.interp_data)
    if material == TiO2:
        # TiO2.n_spline(wavelengths)
        TiO2.n_interp(wavelengths)
        # n_plot('TiO2',wavelengths, TiO2.interp_data)
    if material == InP:
        # InP.n_spline(wavelengths)
        InP.n_interp(wavelengths)
        # n_plot('InP',wavelengths, InP.interp_data)
