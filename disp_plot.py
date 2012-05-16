#

import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt

def omega_plot():#spectra_name, spec_list, solar_cell, light, max_num_BMs, max_order_PWs, Efficiency):
	fig = plt.figure(num=None, figsize=(12, 9), dpi=80, facecolor='w', edgecolor='k')
	ax1 = fig.add_subplot(1,1,1)
	# ax2 = fig.add_subplot(2,1,2)
	# data  = np.loadtxt('%s.txt' % 'omega2')
	data        = np.genfromtxt('%s.txt' % 'omega', usecols=(0, 2))
	wavelengths = np.genfromtxt('%s.txt' % 'omega', usecols=(2))
	num_BMs     = np.genfromtxt('%s.txt' % 'omega', usecols=(1))

	# h  = 6.626068e-34;
	# c  = 299792458;
	# eV = 1.60217646e-19;

	for i in range(len(num_BMs)):
		prop = []
		re = np.genfromtxt('%s.txt' % 'omega', usecols=(5+2*i), invalid_raise=False)
		im = np.genfromtxt('%s.txt' % 'omega', usecols=(5+2*i+1), invalid_raise=False)
		for j in range(len(re)):
			if re[j] > im[j]:
				prop.append(re[j])
		# print np.size(re)
		# print max(re)
		# print re[0]
		# print re[-1]

		trim_wls = wavelengths[0:len(prop)]
		ax1.plot(trim_wls, prop, 'ro')
		ax1.set_xlabel('Wavelength (nm)')
		ax1.set_ylabel(r'Re(k$_z$)')
		# e = (h*c/(trim_wls*1e-9))/eV
		# ax2.plot(e, prop, 'ro')
		# ax2.set_xlabel('E (eV)')
		# ax2.set_ylabel(r'Re(k$_z$)')
	plt.xlim((wavelengths[0], wavelengths[-1]))
	plt.savefig('Disp_Diagram')


if __name__ == "__main__":
    import sys
    omega_plot()