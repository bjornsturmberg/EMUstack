# doc

import numpy as np

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

air    = Material(np.loadtxt('%sair.txt'% data_location))
Si_c   = Material(np.loadtxt('%sSi_c.txt'% data_location))
Si_a   = Material(np.loadtxt('%sSi_a.txt'% data_location))
SiO2_a = Material(np.loadtxt('%sSiO2_a.txt'% data_location))
CuO    = Material(np.loadtxt('%sCuO.txt'% data_location))
CdTe   = Material(np.loadtxt('%sCdTe.txt'% data_location))
FeS2   = Material(np.loadtxt('%sFeS2.txt'% data_location))
Zn2P3  = Material(np.loadtxt('%sZn2P3.txt'% data_location))

# def interp_all(wavelengths):
#     air.n_interp(wavelengths)
#     Si_c.n_interp(wavelengths)
#     Si_a.n_interp(wavelengths)
#     SiO2_a.n_interp(wavelengths)
#     CuO.n_interp(wavelengths)
#     # CdTe.n_interp(wavelengths)
#     FeS2.n_interp(wavelengths)
#     # Zn2P3.n_interp(wavelengths)

def interp_needed(wavelengths, inclusion, background, superstrate, substrate):
    if inclusion == air or background == air or superstrate == air or substrate == air:
        air.n_interp(wavelengths)
    if inclusion == Si_c or background == Si_c or superstrate == Si_c or substrate == Si_c:
        Si_c.n_interp(wavelengths)
    if inclusion == Si_a or background == Si_a or superstrate == Si_a or substrate == Si_a:
      Si_a.n_interp(wavelengths)
    if inclusion == SiO2_a or background == SiO2_a or superstrate == SiO2_a or substrate == SiO2_a:
       SiO2_a.n_interp(wavelengths)
    if inclusion == CuO or background == CuO or superstrate == CuO or substrate == CuO:
       CuO.n_interp(wavelengths)
    if inclusion == CdTe or background == CdTe or superstrate == CdTe or substrate == CdTe:
       CdTe.n_interp(wavelengths)
    if inclusion == FeS2 or background == FeS2 or superstrate == FeS2 or substrate == FeS2:
       FeS2.n_interp(wavelengths)
    if inclusion == Zn2P3 or background == Zn2P3 or superstrate == Zn2P3 or substrate == Zn2P3:
       Zn2P3.n_interp(wavelengths)