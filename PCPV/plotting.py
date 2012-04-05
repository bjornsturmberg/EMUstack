#

import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt




def irradiance(Absorptance, Irradiance, Weighted):
	a_data         = np.loadtxt('%s.txt' % Absorptance)
	wavelengths    = a_data[:,0]
	a_spec         = a_data[:,2]
	i_data         = np.loadtxt('%s.txt' % Irradiance)
	i_spec         = np.interp(wavelengths, i_data[:,0], i_data[:,2]) 
	weighted_abs   = a_data[:,2]*i_spec/i_spec.max()
	weighted_array = zip(wavelengths, a_data[:,1], weighted_abs)
	np.savetxt('%s.txt' % Weighted, weighted_array, fmt = '%18.12f')

	#  Total solar irradiance - integral of I(lambda) from 310nm-4000nm
	#  intergral done in Mathematica (OtherCode/Silicon_ASTM/ASTMG173.nb)
	tot_irradiance = 900.084
	bandgap_wl     = wavelengths[-1] #have as property of material.
	expression     = i_spec*a_spec*wavelengths
	integral_tmp   = np.trapz(expression, x=wavelengths)
	Efficiency     = integral_tmp/(bandgap_wl*tot_irradiance)
	np.savetxt('Efficiency.txt', [Efficiency], fmt = '%12.8f')
	return Efficiency


def tra_plot(spec_list, solar_cell, light, max_num_BMs, max_order_PWs):
	fig = plt.figure(num=None, figsize=(8, 12), dpi=80, facecolor='w', edgecolor='k')
	for i in range(len(spec_list)):
		ax1 = fig.add_subplot(4,1,i+1, adjustable='box', aspect=400)
		spec_name = spec_list.pop(0)
		s_data  = np.loadtxt('%s.txt' % spec_name)
		wavelengths = s_data[:,0]
		spectrum = s_data[:,2]
		ax1.plot(wavelengths, spectrum)
		ax1.set_xlabel('Wavelength (nm)')
		ax1.set_ylabel(spec_name)
		plt.axis([wavelengths[0], wavelengths[-1], 0, 1])
		# ax2 = ax1.twiny()
		# frequencies = s_data[:,1]
		# ax2.plot(frequencies, spectrum)
		# ax2.set_ylabel('Frequency')
		# this is an inset axes over the main axes
		# create some data to use for the plot
	tmp1 = 'a1 = %(radius)d, a2 = %(rad)d, d = %(period)d, ff = %(ff)4.2f, h = %(height)d, '% {
	'radius' 	    : solar_cell.radius1,
	'rad' 	        : solar_cell.radius2,
	'period' 	    : solar_cell.period,
	'ff' 	  		: solar_cell.ff, 
	'height' 	    : solar_cell.height, }
	tmp2 = r'$\theta$ = '
	tmp3 = '%(theta)6.2f, '% {
	'theta' 	    : light.theta, }
	tmp4 = r'$\phi$ = '
	tmp5 = '%(phi)6.2f, '% {
	'phi'	 	    : light.phi, }
	tmp6 = '\nmax_BMs = %(max_num_BMs)d, max_PW_order = %(max_order_PWs)d'% {
	'max_num_BMs'	: max_num_BMs,
	'max_order_PWs'	: max_order_PWs, }

	imp_facts = tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6
	plt.suptitle(imp_facts)
	plt.savefig('tra_spectra')