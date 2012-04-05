# doc

import numpy as np
# import data

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


# air  = Material(data.air_data)
air  = Material(np.loadtxt('%sair_r_index_1.txt'% data_location))
cSi  = Material(np.loadtxt('%scSi_r_index_1.txt'% data_location))
aSi  = Material(np.loadtxt('%saSi_r_index_1.txt'% data_location))
SiO2 = Material(np.loadtxt('%sSiO2_r_index_1.txt'% data_location))

def interp_all(wavelengths):
    air.n_interp(wavelengths)
    cSi.n_interp(wavelengths)
    aSi.n_interp(wavelengths)
    SiO2.n_interp(wavelengths)
