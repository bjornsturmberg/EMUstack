#

import numpy as np
from matplotlib.mlab import griddata
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import os

Irrad_spec_file = '../PCPV/Data/ASTM_1_5_spectrum'

def average_spec(spec_name,av_spec_name,num_wl,num_h):
	data      = np.loadtxt('%s.txt' % spec_name)
	av_wl     = []
	av_spec   = []
	av_h      = []
	for i in np.linspace(0,num_wl-1,num_wl):
		av_wl   = np.append(av_wl,data[i*num_h,0])
		av_tmp  = np.mean(data[i*num_h:(i+1)*num_h,1])
		av_spec = np.append(av_spec,av_tmp)
		av_tmp  = np.mean(data[i*num_h:(i+1)*num_h,2])
		av_h    = np.append(av_h,av_tmp)
	# Save averages to file
	av_array = zip(av_wl, av_spec, av_h)
	np.savetxt('%s.txt' % av_spec_name, av_array, fmt = '%18.12f')


def efficiency_h(spec_name,h_spec_name,wavelengths,num_wl,num_h, Animate):
	data   = np.loadtxt('%s.txt' % spec_name)
	spec   = data[:,1] 		#NEEED CHANGE
	i_data = np.loadtxt('%s.txt' % Irrad_spec_file)
	i_spec = np.interp(wavelengths, i_data[:,0], i_data[:,2])
	h_Efficiency = []
	h_data       = []
	for i in np.linspace(0,num_h-1,num_h):
		h_wl = []
		h = data[i,2]		#NEEED CHANGE
		for j in np.linspace(0,num_wl-1,num_wl):
			tmp  = spec[i + j*(num_h)]
			h_wl = np.append(h_wl,tmp)
		eta_calc = efficiency(wavelengths,i_spec,h_wl)
		h_Efficiency = np.append(h_Efficiency,eta_calc)
		h_data       = np.append(h_data,data[i,2])
		if Animate == True:
			spectra_h(wavelengths,h_wl, h, eta_calc, i)
	h_array = zip(h_Efficiency, h_data)
	np.savetxt('%s.txt' % h_spec_name, h_array, fmt = '%18.12f')
	if Animate == True:
		ps = os.system("convert -delay 3 +dither -layers Optimize -colors 16 Animated/*.pdf Animated/Abs_Spectra.gif")
		p  = os.system("gifsicle -O2 Animated/Abs_Spectra.gif -o Animated/Abs_Spectra-opt.gif")

def spectra_h(wavelengths,h_wl, h, eta_calc,i):
	fig = plt.figure(num=None, figsize=(9, 4), dpi=80, facecolor='w', edgecolor='k')
	ax1 = fig.add_subplot(1,1,1)#, adjustable='box', aspect=400)
	ax1.plot(wavelengths, h_wl)
	ax1.set_xlabel('Wavelength (nm)')
	ax1.set_ylabel('Absorptance')
	plt.axis([wavelengths[0], wavelengths[-1], 0, 1])
	tmp1 = 'h = %(h)d (nm), '% {
	'h' 	    : h, }
	tmp7 = r'$\eta$ = %(eta_calc)6.2f'% {
	'eta_calc'	: eta_calc*100, }
	tmp8 = ' %'

	imp_facts = tmp1 + tmp7 + tmp8
	plt.suptitle(imp_facts)
	if i < 10:
		plt.savefig('Animated/000%d' % i)
	elif i < 100:
		plt.savefig('Animated/00%d' % i)
	elif i < 1000:
		plt.savefig('Animated/0%d' % i)
	else:
		plt.savefig('Animated/%d' % i)


def irradiance(Absorb, W_Absorb, Trans, W_Trans, Reflec, W_Reflec, radius1, radius2,
		period, ff):
	data         = np.loadtxt('%s.txt' % Absorb)
	wavelengths  = data[:,0]
	spec         = data[:,1]
	h_av         = data[:,2] 
	i_data       = np.loadtxt('%s.txt' % Irrad_spec_file)
	i_spec       = np.interp(wavelengths, i_data[:,0], i_data[:,2])

	# call efficiency function 
	Efficiency = efficiency(wavelengths,i_spec,spec)
	np.savetxt('Efficiency.txt', [Efficiency, radius1, radius2,
		period, ff], fmt = '%12.8f')

	# weighted absorptance
	weighting(i_spec, wavelengths, spec, h_av, W_Absorb)
	# weighted transmittance
	data         = np.loadtxt('%s.txt' % Trans)
	spec         = data[:,1] 
	weighting(i_spec, wavelengths, spec, h_av, W_Trans)
	# weighted reflectance
	data         = np.loadtxt('%s.txt' % Reflec)
	spec         = data[:,1] 
	weighting(i_spec, wavelengths, spec, h_av, W_Reflec)
	return Efficiency

def efficiency(wavelengths,i_spec,spec):
	#  Total solar irradiance - integral of I(lambda) from 310nm-4000nm
	#  intergral done in Mathematica (OtherCode/Silicon_ASTM/ASTMG173.nb)
	tot_irradiance = 900.084
	bandgap_wl     = wavelengths[-1] #have as property of material.
	expression     = i_spec*spec*wavelengths
	integral_tmp   = np.trapz(expression, x=wavelengths)
	Efficiency     = integral_tmp/(bandgap_wl*tot_irradiance)
	return Efficiency

def weighting(weight_by, wavelengths, spectrum, h_av, w_spectum):
	# spec_data      = np.loadtxt('%s.txt' % spectrum)
	# spec           = spec_data[:,1]
	weighted       = spectrum*weight_by/weight_by.max()
	weigthed_array = zip(wavelengths, weighted, h_av)
	np.savetxt('%s.txt' % w_spectum, weigthed_array, fmt = '%18.12f')


def tra_plot(spectra_name, spec_list, solar_cell, light, max_num_BMs, max_order_PWs, Efficiency):
	fig = plt.figure(num=None, figsize=(9, 12), dpi=80, facecolor='w', edgecolor='k')
	for i in range(len(spec_list)):
		ax1 = fig.add_subplot(3,1,i+1)#, adjustable='box', aspect=400)
		spec_name = spec_list.pop(0)
		s_data  = np.loadtxt('%s.txt' % spec_name)
		wavelengths = s_data[:,0]
		spectrum = s_data[:,1]
		ax1.plot(wavelengths, spectrum)
		ax1.set_xlabel('Wavelength (nm)')
		ax1.set_ylabel(spec_name)
		plt.axis([wavelengths[0], wavelengths[-1], 0, 1])
	tmp1 = 'a1 = %(radius)d, a2 = %(rad)d, d = %(period)d, ff = %(ff)4.2f, h_1 = %(h_one)d, h_2 = %(h_two)d, num_h = %(num_h)d, \n'% {
	'radius' 	    : solar_cell.radius1,
	'rad' 	        : solar_cell.radius2,
	'period' 	    : solar_cell.period,
	'ff' 	  		: solar_cell.ff, 
	'h_one' 	    : solar_cell.height_1,
	'h_two' 	    : solar_cell.height_2,
	'num_h' 	    : solar_cell.num_h, }
	tmp2 = r'$\theta$ = '
	tmp3 = '%(theta)6.2f, '% {
	'theta' 	    : light.theta, }
	tmp4 = r'$\phi$ = '
	tmp5 = '%(phi)6.2f, '% {
	'phi'	 	    : light.phi, }
	tmp6 = 'max_BMs = %(max_num_BMs)d, max_PW_order = %(max_order_PWs)d'% {
	'max_num_BMs'	: max_num_BMs,
	'max_order_PWs'	: max_order_PWs, }
	tmp7 = '\n' r'$\eta$ = %(Efficiency)6.2f'% {
	'Efficiency'	: Efficiency*100, }
	tmp8 = ' %'

	imp_facts = tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6 + tmp7 + tmp8
	plt.suptitle(imp_facts)
	plt.savefig(spectra_name)


def overlay_plot(spectra_name, spec_list, solar_cell, light, max_num_BMs, max_order_PWs, Efficiency):
	fig = plt.figure(num=None, figsize=(9, 12), dpi=80, facecolor='w', edgecolor='k')
	ax1 = fig.add_subplot(3,1,2)
	line_list = ['r-.','b--','k-']
	for i in range(len(spec_list)):
		spec_name   = spec_list.pop(0)
		line_name   = line_list.pop(0)
		s_data      = np.loadtxt('%s.txt' % spec_name)
		wavelengths = s_data[:,0]
		spectrum    = s_data[:,1]
		ax1.plot(wavelengths, spectrum, line_name)
		ax1.set_xlabel('Wavelength (nm)')
		ax1.set_ylabel(spec_name)
		plt.axis([wavelengths[0], wavelengths[-1], 0, 1])
	plt.legend( ('a = 70 nm', 'a = 123 nm','supercell') )
	tmp1 = 'a1 = %(radius)d, a2 = %(rad)d, d = %(period)d, ff = %(ff)4.2f, h_1 = %(h_one)d, h_2 = %(h_two)d, num_h = %(num_h)d, '% {
	'radius' 	    : solar_cell.radius1,
	'rad' 	        : solar_cell.radius2,
	'period' 	    : solar_cell.period,
	'ff' 	  		: solar_cell.ff, 
	'h_one' 	    : solar_cell.height_1,
	'h_two' 	    : solar_cell.height_2,
	'num_h' 	    : solar_cell.num_h, }
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
	tmp8 = ' %'

	imp_facts = tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6 + tmp7 + tmp8
	plt.suptitle(imp_facts)
	plt.savefig(spectra_name)


def height_plot(name_out, name_in, solar_cell, light, max_num_BMs, max_order_PWs, Efficiency_h,
	Efficiency,num_h):
	fig   = plt.figure(num=None, figsize=(12, 9), dpi=80, facecolor='w', edgecolor='k')
	ax1   = fig.add_subplot(111)

	data    = np.loadtxt('%s.txt' % name_in)
	wl_data = data[:,0]
	s_data  = data[:,1]
	h_data  = data[:,2]

	wl_1    = wl_data[0]
	wl_2    = wl_data[-1]
	d_wl    = 1
	h_1     = h_data[0]
	h_2     = h_data[-1]
	h_int   = (h_2-h_1)/(num_h-1)

	# Interpolate the data
	int_wl,int_h = np.mgrid[wl_1:wl_2+d_wl:d_wl, h_1:h_2+h_int:h_int]
	int_Eff    = griddata(wl_data, h_data, s_data, int_wl, int_h)
	int_Eff_t  = int_Eff.T

	# CS = plt.contourf(int_wl,int_h,int_Eff,100,cmap=plt.cm.jet)
	CS = plt.imshow(int_Eff_t, cmap=plt.cm.hot, norm=None, aspect='auto', interpolation=None,
	  vmin=0, vmax=s_data.max(), origin='lower', extent=(wl_1, wl_2, h_1, h_2))
	ax1.set_xlabel('Wavelength (nm)')
	ax1.set_ylabel('Height (nm)')
	tick_array = []
	for i in range(0,10,1):
		tick = float(i)/10
		if tick < s_data.max():
			tick_array.append(tick)
	tick_array.append(round(s_data.max(),4))
	cbar = plt.colorbar(extend='neither',ticks=tick_array)
	cbar.set_ticklabels(tick_array)
	cbar.ax.set_ylabel('Absorptance')

	tmp1 = 'a1 = %(radius)d, a2 = %(rad)d, d = %(period)d, ff = %(ff)4.2f, h_1 = %(h_one)d, h_2 = %(h_two)d, num_h = %(num_h)d, '% {
	'radius' 	    : solar_cell.radius1,
	'rad' 	        : solar_cell.radius2,
	'period' 	    : solar_cell.period,
	'ff' 	  		: solar_cell.ff, 
	'h_one' 	    : solar_cell.height_1,
	'h_two' 	    : solar_cell.height_2,
	'num_h' 	    : solar_cell.num_h, }
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
	tmp8 = ' %'
	imp_facts = tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6 + tmp7 + tmp8
	plt.suptitle(imp_facts)
	plt.savefig(name_out)

# Plot vs Energy
	fig  = plt.figure(num=None, figsize=(12, 9), dpi=80, facecolor='w', edgecolor='k')
	ax1  = fig.add_subplot(111)
	h  = 6.626068e-34;
	c  = 299792458;
	eV = 1.60217646e-19;
	e_data  = (h*c/(wl_data*1e-9))/eV
	e_1    = e_data[0]
	e_2    = e_data[-1]
	d_e    = (e_2-e_1)/(500-1)

	# Interpolate the data
	int_e,int_h = np.mgrid[e_1:e_2+d_e:d_e, h_1:h_2+h_int:h_int]
	int_Eff_e   = griddata(e_data, h_data, s_data, int_e, int_h)
	int_Eff_e_t = np.fliplr(int_Eff_e.T)

	CS = plt.imshow(int_Eff_e_t, cmap=plt.cm.hot, norm=None, aspect='auto', interpolation=None,
	  alpha=None, vmin=0, vmax=s_data.max(), origin='lower', extent=(e_2, e_1, h_1, h_2))#, shape=None, filternorm=1, filterrad=4.0, imlim=None, resample=None, url=None, hold=None
	ax1.set_xlabel('Energy (eV)')
	ax1.set_ylabel('Height (nm)')
	cbar = plt.colorbar(extend='neither',ticks=tick_array)
	cbar.set_ticklabels(tick_array)
	cbar.ax.set_ylabel('Absorptance')
	plt.suptitle(imp_facts)
	plt.savefig('%s_energy' % name_out)

# Plot efficiency vs height
	fig2= plt.figure(num=None, figsize=(12, 9), dpi=80, facecolor='w', edgecolor='k')
	ax2 = fig2.add_subplot(153)
	data = np.loadtxt('%s.txt' % Efficiency_h)
	es = data[:,0]*100
	hs = data[:,1]
	ax2.plot(es, hs,'r')
	ax2.set_xlabel(r'$\eta$' ' (%)')
	ax2.set_ylabel('Height (nm)')
	plt.axis([es[-1], 0, hs[0], hs[-1]])
	plt.xticks([round(es.max(),1), 0])
	plt.fill_betweenx(hs,es,color='blue', alpha=0.5)
	plt.suptitle(imp_facts)
	plt.savefig('%s_eta_bar' % name_out)


def omega_plot(solar_cell, light, max_num_BMs, max_order_PWs, Efficiency):
	fig = plt.figure(num=None, figsize=(12, 9), dpi=80, facecolor='w', edgecolor='k')
	ax1 = fig.add_subplot(1,1,1)
	wavelengths = np.genfromtxt('%s.txt' % 'omega', usecols=(2))
	num_BMs     = np.genfromtxt('%s.txt' % 'omega', usecols=(1))

	for i in range(len(num_BMs)):
		prop = []
		re = np.genfromtxt('%s.txt' % 'omega', usecols=(5+2*i), invalid_raise=False)
		im = np.genfromtxt('%s.txt' % 'omega', usecols=(5+2*i+1), invalid_raise=False)
		for j in range(len(re)):
			if re[j] > im[j]:
				prop.append(re[j])

		trim_wls = wavelengths[0:len(prop)]
		ax1.plot(trim_wls, prop, 'ro')
		ax1.set_xlabel('Wavelength (nm)')
		ax1.set_ylabel(r'Re(k$_z$)')
	plt.xlim((wavelengths[0], wavelengths[-1]))

	tmp1 = 'a1 = %(radius)d, a2 = %(rad)d, d = %(period)d, ff = %(ff)4.2f, h_1 = %(h_one)d, h_2 = %(h_two)d, num_h = %(num_h)d, '% {
	'radius' 	    : solar_cell.radius1,
	'rad' 	   		: solar_cell.radius2,
	'period' 	    : solar_cell.period,
	'ff' 	  		: solar_cell.ff, 
	'h_one' 	    : solar_cell.height_1,
	'h_two' 	    : solar_cell.height_2,
	'num_h' 	    : solar_cell.num_h, }
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
	tmp8 = ' %'
	imp_facts = tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6 + tmp7 + tmp8
	plt.suptitle(imp_facts)
	plt.savefig('Disp_Diagram')