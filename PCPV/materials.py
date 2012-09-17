# doc

import numpy as np
import copy

data_location = '../PCPV/Data/'

class Material(object):
    """docstring for """
    def __init__(self, n_data):
         self.data_wls   = n_data[:,0]
         self.data_re_ns = n_data[:,1]
         self.data_im_ns = n_data[:,2]
         self.stored_ns  = {}

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
         # spline_fit = interp1d(self.data_wls, self.data_re_ns, kind = 2)
         # interp_re_n = spline_fit(wavelengths)
         interp_re_n = np.interp(wavelengths, 
             self.data_wls, self.data_re_ns)
         interp_im_n = np.interp(wavelengths, 
             self.data_wls, self.data_im_ns)
         self.interp_data = interp_re_n + 1j*interp_im_n
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
Sb2S3  = Material(np.loadtxt('%sSb2S3.txt'% data_location))

test   = Material(np.loadtxt('%sTuniz.txt'% data_location))#copy.deepcopy(Air)

# def interp_all(wavelengths):
#     Air.n_interp(wavelengths)
#     Si_c.n_interp(wavelengths)
#     Si_a.n_interp(wavelengths)
#     SiO2_a.n_interp(wavelengths)
#     CuO.n_interp(wavelengths)
#     CdTe.n_interp(wavelengths)
#     FeS2.n_interp(wavelengths)
#     Zn3P2.n_interp(wavelengths)

def interp_needed(wavelengths, inclusion_a, inclusion_b, background, superstrate, substrate):
    if inclusion_a == Air or inclusion_b == Air or background == Air or superstrate == Air or substrate == Air:
        Air.n_interp(wavelengths)
        # n_plot('Air',wavelengths, Air.interp_data)
    if inclusion_a == Si_c or inclusion_b == Si_c or background == Si_c or superstrate == Si_c or substrate == Si_c:
        Si_c.n_interp(wavelengths)
        # n_plot('Si_c',wavelengths, Si_c.interp_data)
    if inclusion_a == Si_a or inclusion_b == Si_a or background == Si_a or superstrate == Si_a or substrate == Si_a:
        Si_a.n_interp(wavelengths)
        # n_plot('Si_a',wavelengths, Si_a.interp_data)
    if inclusion_a == SiO2_a or inclusion_b == SiO2_a or background == SiO2_a or superstrate == SiO2_a or substrate == SiO2_a:
        SiO2_a.n_interp(wavelengths)
        # n_plot('SiO2_a',wavelengths, SiO2_a.interp_data)
    if inclusion_a == CuO or inclusion_b == CuO or background == CuO or superstrate == CuO or substrate == CuO:
        CuO.n_interp(wavelengths)
        # n_plot('CuO',wavelengths, CuO.interp_data)
    if inclusion_a == CdTe or inclusion_b == CdTe or background == CdTe or superstrate == CdTe or substrate == CdTe:
        CdTe.n_interp(wavelengths)
        # n_plot('CdTe',wavelengths, CdTe.interp_data)
    if inclusion_a == FeS2 or inclusion_b == FeS2 or background == FeS2 or superstrate == FeS2 or substrate == FeS2:
        FeS2.n_interp(wavelengths)
        # n_plot('FeS2',wavelengths, FeS2.interp_data)
    if inclusion_a == Zn3P2 or inclusion_b == Zn3P2 or background == Zn3P2 or superstrate == Zn3P2 or substrate == Zn3P2:
        Zn3P2.n_interp(wavelengths)
        # n_plot('Zn3P2',wavelengths, Zn3P2.interp_data)
    if inclusion_a == Au or inclusion_b == Au or background == Au or superstrate == Au or substrate == Au:
        Au.n_interp(wavelengths)
        # n_plot('Au',wavelengths, Au.interp_data)
    if inclusion_a == Sb2S3 or inclusion_b == Sb2S3 or background == Sb2S3 or superstrate == Sb2S3 or substrate == Sb2S3:
        Sb2S3.n_interp(wavelengths)
        # n_plot('Sb2S3',wavelengths, Sb2S3.interp_data)
    if inclusion_a == test or inclusion_b == test or background == test or superstrate == test or substrate == test:
        test.n_interp(wavelengths)

# from matplotlib.mlab import griddata
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