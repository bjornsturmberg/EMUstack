#

import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt




def irradiance(Irradiance, Absorb, W_Absorb, Trans, W_Trans, Reflec, W_Reflec, radius1, radius2,
		period, ff):
	a_data         = np.loadtxt('%s.txt' % Absorb)
	wavelengths    = a_data[:,0]
	a_spec         = a_data[:,2]
	i_data         = np.loadtxt('%s.txt' % Irradiance)
	i_spec         = np.interp(wavelengths, i_data[:,0], i_data[:,2])

	#  Total solar irradiance - integral of I(lambda) from 310nm-4000nm
	#  intergral done in Mathematica (OtherCode/Silicon_ASTM/ASTMG173.nb)
	tot_irradiance = 900.084
	bandgap_wl     = wavelengths[-1] #have as property of material.
	expression     = i_spec*a_spec*wavelengths
	integral_tmp   = np.trapz(expression, x=wavelengths)
	Efficiency     = integral_tmp/(bandgap_wl*tot_irradiance)
	np.savetxt('Efficiency.txt', [Efficiency, radius1, radius2,
		period, ff], fmt = '%12.8f')

	# weighted absorptance
	weighting(i_spec, wavelengths, Absorb, W_Absorb)
	# weighted transmittance
	weighting(i_spec, wavelengths, Trans, W_Trans)
	# weighted reflectance
	weighting(i_spec, wavelengths, Reflec, W_Reflec)
	return Efficiency

def weighting(weight_by, wavelengths, spectrum, w_spectum):
	spec_data      = np.loadtxt('%s.txt' % spectrum)
	spec           = spec_data[:,2]
	weighted       = spec*weight_by/weight_by.max()
	weigthed_array = zip(wavelengths, spec_data[:,1], weighted)
	np.savetxt('%s.txt' % w_spectum, weigthed_array, fmt = '%18.12f')


def tra_plot(spectra_name, spec_list, solar_cell, light, max_num_BMs, max_order_PWs, Efficiency):
	fig = plt.figure(num=None, figsize=(8, 12), dpi=80, facecolor='w', edgecolor='k')
	for i in range(len(spec_list)):
		ax1 = fig.add_subplot(3,1,i+1, adjustable='box', aspect=400)
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
	tmp7 = '\n' r'$\eta$ = %(Efficiency)6.2f'% {
	'Efficiency'	: Efficiency*100, }

	imp_facts = tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6 + tmp7
	plt.suptitle(imp_facts)
	plt.savefig(spectra_name)
