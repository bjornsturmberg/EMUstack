"""
    plotting.py is a subroutine of PCPV that contains numerous plotting
    routines. These were developed during simulations for photovoltaics,
    hence the efficiency calculations.

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

import objects
import mode_calcs

import numpy as np
import subprocess
from matplotlib.mlab import griddata
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt





#######################################################################################
def clear_previous(file_type):
    """ Delete all files of specified type 'file_type', eg '.txt' """

    files_rm = 'rm *%s'% file_type
    subprocess.call(files_rm, shell = True)
#######################################################################################



#######################################################################################
def gen_params_string(param_layer, light, max_num_BMs=0):
    """ Generate the string of simulation info that is to be printed at the top of plots. """

    # Plot t,r,a for each layer & total, then save each to text files
    if isinstance(param_layer,objects.NanoStruct):
        params_2_print = 'ff = %5.3f, '% param_layer.ff
        params_2_print += 'd = %(period)d, a1 = %(diameter)d, '% {
        'period'        : param_layer.period, 'diameter' : param_layer.diameter1,}
        if param_layer.diameter2 != 0: params_2_print += 'a2 = %(rad)d '% {'rad' : param_layer.diameter2,}
        if param_layer.geometry == 'NW_array':
            if param_layer.diameter3  != 0: params_2_print += 'a3 = %(rad)d '%    {'rad' : param_layer.diameter3,}
            if param_layer.diameter4  != 0: params_2_print += 'a4 = %(rad)d '%    {'rad' : param_layer.diameter4,}
            if param_layer.diameter5  != 0: params_2_print += '\na5 = %(rad)d '%  {'rad' : param_layer.diameter5,}
            if param_layer.diameter6  != 0: params_2_print += 'a6 = %(rad)d '%    {'rad' : param_layer.diameter6,}
            if param_layer.diameter7  != 0: params_2_print += 'a7 = %(rad)d '%    {'rad' : param_layer.diameter7,}
            if param_layer.diameter8  != 0: params_2_print += 'a8 = %(rad)d '%    {'rad' : param_layer.diameter8,}
            if param_layer.diameter9  != 0: params_2_print += 'a9 = %(rad)d \n'%  {'rad' : param_layer.diameter9,}
            if param_layer.diameter10 != 0: params_2_print += 'a10 = %(rad)d '%   {'rad' : param_layer.diameter10,}
            if param_layer.diameter11 != 0: params_2_print += 'a11 = %(rad)d '%   {'rad' : param_layer.diameter11,}
            if param_layer.diameter12 != 0: params_2_print += 'a12 = %(rad)d '%   {'rad' : param_layer.diameter12,}
            if param_layer.diameter13 != 0: params_2_print += 'a13 = %(rad)d '%   {'rad' : param_layer.diameter13,}
            if param_layer.diameter14 != 0: params_2_print += 'a14 = %(rad)d '%   {'rad' : param_layer.diameter14,}
            if param_layer.diameter15 != 0: params_2_print += 'a15 = %(rad)d '%   {'rad' : param_layer.diameter15,}
            if param_layer.diameter16 != 0: params_2_print += 'a16 = %(rad)d \n'% {'rad' : param_layer.diameter16,}
            if param_layer.square == True: params_2_print += '\nSquare NWs '
            if param_layer.ellipticity == True: params_2_print += '\nEllipticity = %(rad)5.3f '% {'rad' : param_layer.ellipticity,}
        elif param_layer.geometry == '1D_grating':
            params_2_print += ''
        params_2_print += 'max_BMs = %(max_num_BMs)d, max_PW_order = %(max_order_PWs)d, '% {
        'max_num_BMs'   : max_num_BMs,'max_order_PWs' : light.max_order_PWs, }
    else:
        params_2_print = 'max_PW_order = %(max_order_PWs)d, '% {'max_order_PWs' : light.max_order_PWs, }


    # k_pll = light.k_pll * param_layer.period
    # params_2_print += r'$k_\parallel d$ = '
    # tmp3 = '(%(kx)1.4f, %(ky)1.4f), '% {'kx' : k_pll[0], 'ky' : k_pll[1]}
    # params_2_print += r'$\theta$ = %(theta)6.2f, $\phi$ = %(phi)6.2f, '% {
    # 'theta' : light.theta,'phi' : light.phi, } 

    return params_2_print
#######################################################################################



def zeros_int_str(zero_int):
    """ Convert integer into string with '0' in place of ' ' """
    string = '%4.0f' % zero_int
    fmt_string = string.replace(' ','0')
    return fmt_string

def tick_function(energies):
    """ Convert energy in eV into wavelengths in nm """
    h  = 6.626068e-34;
    c  = 299792458;
    eV = 1.60217646e-19;
    wls = h*c/(energies*eV)*1e9
    return wls

#######################################################################################
def t_r_a_plots(stack_wl_list, wavelengths, params_2_print, active_layer_nu=0, stack_label=1, add_name=''):
    """ Plot t,r,a for each layer & total, then save each to text files. """

    height_list = stack_wl_list[0].heights_nm()[::-1]
    params_2_print += '\n'r'$h_t,...,h_b$ = '
    params_2_print += ''.join('%4d, ' % num for num in height_list)

    stack_label = zeros_int_str(stack_label)

    # wavelengths = np.array([s.layers[0].light.wl_nm for s in stack_wl_list]) #look at first layer to find wls.
    a_list = []
    t_list = []
    r_list = []
    for stack in stack_wl_list:
        a_list.extend(stack.a_list)
        t_list.extend(stack.t_list)
        r_list.extend(stack.r_list)

    layers_steps = len(stack_wl_list[0].layers) - 1
    active_abs = []
    a_tot      = []
    t_tot      = []
    r_tot      = []
    for i in range(len(wavelengths)): 
        active_abs.append(float(a_list[active_layer_nu + i*layers_steps]))
        a_tot.append(float(a_list[layers_steps-1+(i*layers_steps)]))
        t_tot.append(float(t_list[layers_steps-1+(i*layers_steps)]))
        r_tot.append(float(r_list[i]))
    Efficiency, Irradiance = ult_efficiency(active_abs, wavelengths, params_2_print, stack_label,add_name)
    params_2_print += r'$\eta$ = %(Efficiency)6.2f'% {'Efficiency' : Efficiency*100, }
    params_2_print += ' %'

    total_h = sum(stack_wl_list[0].heights_nm()) #look at first wl result to find h.
    layers_plot('Lay_Absorb', a_list, wavelengths, total_h, params_2_print, stack_label, add_name)
    layers_plot('Lay_Trans',  t_list, wavelengths, total_h, params_2_print, stack_label, add_name)
    layers_plot('Lay_Reflec', r_list, wavelengths, total_h, params_2_print, stack_label, add_name)
    # Also plot total t,r,a on a single plot
    plot_name = 'Total_Spectra'
    total_tra_plot(plot_name, a_tot, t_tot, r_tot, wavelengths, params_2_print, stack_label, add_name)
    # Also plot totals weighted by solar irradiance
    bandgap_wl   = wavelengths[-1] 
    weighting = Irradiance/Irradiance.max()*(wavelengths/bandgap_wl)
    a_weighted = a_tot*weighting
    t_weighted = t_tot*weighting
    r_weighted = r_tot*weighting
    plot_name = 'Weighted_Total_Spectra'
    total_tra_plot(plot_name, a_weighted, t_weighted, r_weighted, wavelengths, params_2_print, stack_label, add_name)
    # extinction_plot(t_tot, wavelengths, params_2_print, stack_label, add_name)
    return Efficiency

def ult_efficiency(active_abs, wavelengths, params_2_print, stack_label,add_name):
    """ Calculate the photovoltaic ultimate efficiency achieved in the specified active layer. 
        For definition see Sturmberg et al., Optics Express, Vol. 19, Issue S5, pp. A1067-A1081 (2011)
        http://dx.doi.org/10.1364/OE.19.0A1067
    """
    # TODO make E_g a property of material, not just longest wl included.

    Irrad_spec_file = '../PCPV/Data/ASTM_1_5_spectrum'
    #  Total solar irradiance - integral of I(lambda) from 310nm-4000nm
    #  intergral done in Mathematica (OtherCode/Silicon_ASTM/ASTMG173.nb)
    tot_irradiance = 900.084
    i_data       = np.loadtxt('%s.txt' % Irrad_spec_file)
    i_spec       = np.interp(wavelengths, i_data[:,0], i_data[:,2])
    bandgap_wl   = wavelengths[-1] #have as property of material.
    expression   = i_spec*active_abs*wavelengths
    integral_tmp = np.trapz(expression, x=wavelengths)
    Efficiency   = integral_tmp/(bandgap_wl*tot_irradiance)   
    nums_2_print = params_2_print.split()
    if len(nums_2_print) >= 8:
        eta_string   = '%8.6f \n'% Efficiency + nums_2_print[5].replace(',','\n') + \
          nums_2_print[8].replace(',','\n') #save params in easy to read in fmt
    else:
        eta_string   = '%8.6f \n'% Efficiency
    np.savetxt('Efficiency_stack%(bon)s_%(add)s.txt'% {'bon' : stack_label,'add' : add_name}, np.array([eta_string]), fmt = '%s')
    return Efficiency, i_spec


def layers_plot(spectra_name, spec_list, wavelengths, total_h, params_2_print, stack_label, add_name):
    """ The plotting function for one type of spectrum across all layers. """

    fig = plt.figure(num=None, figsize=(9, 12), dpi=80, facecolor='w', edgecolor='k')
    nu_layers = len(spec_list)/len(wavelengths)
    h_array = np.ones(len(wavelengths))*total_h
    for i in range(nu_layers):
        layer_spec = []
        for wl in range(len(wavelengths)):
            layer_spec = np.append(layer_spec,spec_list[wl*nu_layers + i])
        ax1 = fig.add_subplot(nu_layers,1,i+1)
        ax1.plot(wavelengths, layer_spec)
        ax2 = ax1.twiny()
        new_tick_values = np.linspace(10,0.5,20)
        new_tick_locations = tick_function(new_tick_values)
        new_tick_labels = ["%.1f" % z for z in new_tick_values]
        ax2.set_xticks(new_tick_locations)
        ax2.set_xticklabels(new_tick_labels)
        ax2.set_xlim((wavelengths[0], wavelengths[-1]))  
        if i == 0:
            ax2.set_xlabel('Energy (eV)')
        elif i == nu_layers-1:
            ax1.set_xlabel('Wavelength (nm)')
            ax1.set_ylabel('Total')
        if i != 0:
            ax2.set_xticklabels( () )
        if i != nu_layers-1:
            ax1.set_xticklabels( () )
        if spectra_name == 'Lay_Absorb':
            if i == 0: ax1.set_ylabel('Top Layer')
            if i == nu_layers-2: ax1.set_ylabel('Bottom Layer')
            suptitle_w_params = 'Absorptance in each layer'+add_name+'\n'+params_2_print
            plt.suptitle(suptitle_w_params)
            lay_spec_name = 'Lay_Absorb'
            if i == nu_layers-1: 
                ax1.set_ylabel('Total')
                lay_spec_name = 'Absorptance'
        elif spectra_name == 'Lay_Trans':
            if i == 0: ax1.set_ylabel('Top Layer')
            if i == nu_layers-2: ax1.set_ylabel('Bottom Layer')
            suptitle_w_params = 'Transmittance in each layer'+add_name+'\n'+params_2_print
            plt.suptitle(suptitle_w_params)
            lay_spec_name = 'Lay_Trans'
            if i == nu_layers-1: 
                ax1.set_ylabel('Total')
                lay_spec_name = 'Transmittance'
        elif spectra_name == 'Lay_Reflec':
            if i == 0: ax1.set_ylabel('Top Layer')
            if i == nu_layers-3: ax1.set_ylabel('Bottom Layer')
            if i == nu_layers-2: ax1.set_ylabel('Substrate')
            suptitle_w_params = 'Reflectance in each layer'+add_name+'\n'+params_2_print
            plt.suptitle(suptitle_w_params)
            lay_spec_name = 'Lay_Reflec'
            if i == nu_layers-1: 
                ax1.set_ylabel('Total')
                lay_spec_name = 'Reflectance'
        av_array = zip(wavelengths, layer_spec, h_array)
        ax1.set_xlim((wavelengths[0], wavelengths[-1]))
        plt.ylim((0, 1))

        if i != nu_layers-1: 
            np.savetxt('%(s)s_%(i)i_stack%(bon)s_%(add)s.txt'% {'s' : lay_spec_name, 'i' : i, 
                'bon' : stack_label,'add' : add_name}, av_array, fmt = '%18.11f')
        else:
            np.savetxt('%(s)s_stack%(bon)s_%(add)s.txt'% {'s' : lay_spec_name, 
                'bon' : stack_label,'add' : add_name}, av_array, fmt = '%18.11f')

        plt.savefig('%(s)s_stack%(bon)s_%(add)s'% {'s' : spectra_name, 'bon' : stack_label,'add' : add_name})

def total_tra_plot(plot_name, a_spec, t_spec, r_spec, wavelengths, params_2_print, stack_label, add_name):
    """ The plotting function for total t,r,a spectra on one plot. """

    fig = plt.figure(num=None, figsize=(9, 12), dpi=80, facecolor='w', edgecolor='k')
    ax1 = fig.add_subplot(3,1,1)
    ax1.plot(wavelengths, a_spec)
    ax1.set_xlabel('Wavelength (nm)')
    ax1.set_ylabel('Absorptance')
    ax1.set_xlim((wavelengths[0], wavelengths[-1]))
    ax2 = ax1.twiny()
    new_tick_values = np.linspace(10,0.5,20)
    new_tick_locations = tick_function(new_tick_values)
    new_tick_labels = ["%.1f" % z for z in new_tick_values]
    ax2.set_xticks(new_tick_locations)
    ax2.set_xticklabels(new_tick_labels)
    ax2.set_xlim((wavelengths[0], wavelengths[-1])) 
    ax2.set_xlabel('Energy (eV)')
    plt.ylim((0, 1))
    ax1 = fig.add_subplot(3,1,2)
    ax1.plot(wavelengths, t_spec)
    ax1.set_xlabel('Wavelength (nm)')
    ax1.set_ylabel('Transmittance')
    ax1.set_xlim((wavelengths[0], wavelengths[-1]))
    ax2 = ax1.twiny()
    ax2.set_xticklabels( () )
    ax2.set_xticks(new_tick_locations)
    ax2.set_xlim((wavelengths[0], wavelengths[-1])) 
    plt.ylim((0, 1))
    ax1 = fig.add_subplot(3,1,3)
    ax1.plot(wavelengths, r_spec)
    ax1.set_xlabel('Wavelength (nm)')
    ax1.set_ylabel('Reflectance')
    ax1.set_xlim((wavelengths[0], wavelengths[-1]))
    ax2 = ax1.twiny()
    ax2.set_xticklabels( () )
    ax2.set_xticks(new_tick_locations)
    ax2.set_xlim((wavelengths[0], wavelengths[-1])) 
    plt.ylim((0, 1))
    plt.suptitle(plot_name+add_name+'\n'+params_2_print)
    plt.savefig('%(s)s_stack%(bon)s_%(add)s'% {'s' : plot_name, 'bon' : stack_label,'add' : add_name})



def t_r_a_plots_subs(stack_wl_list, wavelengths, params_2_print, period, sub_n, active_layer_nu=0, stack_label=1, add_name=''):
    """ Plot t,r,a for each layer & total, then save each to text files. """

    height_list = stack_wl_list[0].heights_nm()[::-1]
    params_2_print += '\n'r'$h_t,...,h_b$ = '
    params_2_print += ''.join('%4d, ' % num for num in height_list)

    stack_label = zeros_int_str(stack_label)

    # wavelengths = np.array([s.layers[0].light.wl_nm for s in stack_wl_list]) #look at first layer to find wls.
    a_list = []
    t_list = []
    r_list = []
    for stack in stack_wl_list:
        a_list.extend(stack.a_list)
        t_list.extend(stack.t_list)
        r_list.extend(stack.r_list)

    layers_steps = len(stack_wl_list[0].layers) - 1
    active_abs = []
    a_tot      = []
    t_tot      = []
    r_tot      = []
    for i in range(len(wavelengths)): 
        active_abs.append(float(a_list[active_layer_nu + i*layers_steps]))
        a_tot.append(float(a_list[layers_steps-1+(i*layers_steps)]))
        t_tot.append(float(t_list[layers_steps-1+(i*layers_steps)]))
        r_tot.append(float(r_list[i]))
    Efficiency, Irradiance = ult_efficiency(active_abs, wavelengths, params_2_print, stack_label)
    params_2_print += r'$\eta$ = %(Efficiency)6.2f'% {'Efficiency' : Efficiency*100, }
    params_2_print += ' %'

    total_h = sum(stack_wl_list[0].heights_nm()) #look at first wl result to find h.
    layers_plot('Lay_Absorb', a_list, wavelengths, total_h, params_2_print, stack_label, add_name)
    layers_plot('Lay_Trans',  t_list, wavelengths, total_h, params_2_print, stack_label, add_name)
    layers_plot('Lay_Reflec', r_list, wavelengths, total_h, params_2_print, stack_label, add_name)
    # Also plot total t,r,a on a single plot
    plot_name = 'Total_Spectra'
    total_tra_plot_subs(plot_name, a_tot, t_tot, r_tot, wavelengths, params_2_print,
      stack_label, add_name, period, sub_n)
    # Also plot totals weighted by solar irradiance
    weighting = Irradiance/Irradiance.max()
    a_weighted = a_tot*weighting
    t_weighted = t_tot*weighting
    r_weighted = r_tot*weighting
    plot_name = 'Weighted_Total_Spectra'
    total_tra_plot_subs(plot_name, a_weighted, t_weighted, r_weighted, wavelengths, params_2_print,
      stack_label, add_name, period, sub_n)
    return Efficiency
def total_tra_plot_subs(plot_name, a_spec, t_spec, r_spec, wavelengths, params_2_print, stack_label, add_name, period, sub_n):
    """ The plotting function for total t,r,a spectra on one plot. """

    t_WA_01 = sub_n * period
    t_WA_11 = sub_n * period / np.sqrt(2)
    t_WA_20 = sub_n * period / 2
    t_WA_21 = sub_n * period / np.sqrt(5)
    t_WA_22 = sub_n * period / np.sqrt(8)
    t_WA_30 = sub_n * period / 3

    r_WA_01 = period
    r_WA_11 = period / np.sqrt(2)
    r_WA_20 = period / 2
    r_WA_21 = period / np.sqrt(5)
    r_WA_22 = period / np.sqrt(8)

    sub_line01 = [t_WA_01, t_WA_01]
    sub_line11 = [t_WA_11, t_WA_11]
    sub_line20 = [t_WA_20, t_WA_20]
    sub_line21 = [t_WA_21, t_WA_21]
    sub_line22 = [t_WA_22, t_WA_22]
    sub_line30 = [t_WA_30, t_WA_30]
    sup_line01 = [r_WA_01, r_WA_01]
    sup_line11 = [r_WA_11, r_WA_11]
    sup_line20 = [r_WA_20, r_WA_20]
    sup_line21 = [r_WA_21, r_WA_21]
    sup_line22 = [r_WA_22, r_WA_22]
    v_line = [0, 1]

    fig = plt.figure(num=None, figsize=(9, 12), dpi=80, facecolor='w', edgecolor='k')
    ax1 = fig.add_subplot(3,1,1)
    ax1.plot(wavelengths, a_spec)
    ax1.plot(sub_line01, v_line, 'k')
    ax1.plot(sub_line11, v_line, 'k')
    ax1.plot(sub_line20, v_line, 'k')
    ax1.plot(sub_line21, v_line, 'k')
    ax1.plot(sub_line22, v_line, 'k')
    ax1.plot(sub_line30, v_line, 'k')
    ax1.plot(sup_line01, v_line, 'r')
    ax1.plot(sup_line11, v_line, 'r')
    ax1.plot(sup_line20, v_line, 'r')
    ax1.plot(sup_line21, v_line, 'r')
    ax1.plot(sup_line22, v_line, 'r')
    ax1.set_xlabel('Wavelength (nm)')
    ax1.set_ylabel('Absorptance')
    ax1.set_xlim((wavelengths[0], wavelengths[-1]))
    ax2 = ax1.twiny()
    new_tick_values = np.linspace(10,0.5,20)
    new_tick_locations = tick_function(new_tick_values)
    new_tick_labels = ["%.1f" % z for z in new_tick_values]
    ax2.set_xticks(new_tick_locations)
    ax2.set_xticklabels(new_tick_labels)
    ax2.set_xlim((wavelengths[0], wavelengths[-1])) 
    ax2.set_xlabel('Energy (eV)')
    plt.ylim((0, 1))
    ax1 = fig.add_subplot(3,1,2)
    ax1.plot(wavelengths, t_spec)
    ax1.plot(sub_line01, v_line, 'k')
    ax1.plot(sub_line11, v_line, 'k')
    ax1.plot(sub_line20, v_line, 'k')
    ax1.plot(sub_line21, v_line, 'k')
    ax1.plot(sub_line22, v_line, 'k')
    ax1.plot(sub_line30, v_line, 'k')
    ax1.plot(sup_line01, v_line, 'r')
    ax1.plot(sup_line11, v_line, 'r')
    ax1.plot(sup_line20, v_line, 'r')
    ax1.plot(sup_line21, v_line, 'r')
    ax1.plot(sup_line22, v_line, 'r')
    ax1.set_xlabel('Wavelength (nm)')
    ax1.set_ylabel('Transmittance')
    ax1.set_xlim((wavelengths[0], wavelengths[-1]))
    ax2 = ax1.twiny()
    ax2.set_xticklabels( () )
    ax2.set_xticks(new_tick_locations)
    ax2.set_xlim((wavelengths[0], wavelengths[-1])) 
    plt.ylim((0, 1))
    ax1 = fig.add_subplot(3,1,3)
    ax1.plot(wavelengths, r_spec)
    ax1.plot(sub_line01, v_line, 'k')
    ax1.plot(sub_line11, v_line, 'k')
    ax1.plot(sub_line20, v_line, 'k')
    ax1.plot(sub_line21, v_line, 'k')
    ax1.plot(sub_line22, v_line, 'k')
    ax1.plot(sub_line30, v_line, 'k')
    ax1.plot(sup_line01, v_line, 'r')
    ax1.plot(sup_line11, v_line, 'r')
    ax1.plot(sup_line20, v_line, 'r')
    ax1.plot(sup_line21, v_line, 'r')
    ax1.plot(sup_line22, v_line, 'r')
    ax1.set_xlabel('Wavelength (nm)')
    ax1.set_ylabel('Reflectance')
    ax1.set_xlim((wavelengths[0], wavelengths[-1]))
    ax2 = ax1.twiny()
    ax2.set_xticklabels( () )
    ax2.set_xticks(new_tick_locations)
    ax2.set_xlim((wavelengths[0], wavelengths[-1])) 
    plt.ylim((0, 1))
    plt.suptitle(plot_name+add_name+'\n'+params_2_print)
    plt.savefig('%(s)s_stack%(bon)s__%(add)s'% {'s' : plot_name, 'bon' : stack_label,'add' : add_name})



def t_r_a_NO_plots(stack_wl_list, wavelengths, params_2_print, active_layer_nu=0, stack_label=1, add_name=''):
    """ Plot t,r,a for each layer & total, then save each to text files. """

    height_list = stack_wl_list[0].heights_nm()[::-1]
    params_2_print += '\n'r'$h_t,...,h_b$ = '
    params_2_print += ''.join('%4d, ' % num for num in height_list)

    stack_label = zeros_int_str(stack_label)

    # wavelengths = np.array([s.layers[0].light.wl_nm for s in stack_wl_list]) #look at first layer to find wls.
    a_list = []
    t_list = []
    r_list = []
    for stack in stack_wl_list:
        a_list.extend(stack.a_list)
        t_list.extend(stack.t_list)
        r_list.extend(stack.r_list)

    layers_steps = len(stack_wl_list[0].layers) - 1
    active_abs = []
    a_tot      = []
    t_tot      = []
    r_tot      = []
    for i in range(len(wavelengths)): 
        active_abs.append(float(a_list[active_layer_nu + i*layers_steps]))
        a_tot.append(float(a_list[layers_steps-1+(i*layers_steps)]))
        t_tot.append(float(t_list[layers_steps-1+(i*layers_steps)]))
        r_tot.append(float(r_list[i]))
    Efficiency, Irradiance = ult_efficiency(active_abs, wavelengths, params_2_print, stack_label)
    return Efficiency





def extinction_plot(t_spec, wavelengths, params_2_print, stack_label, add_name):
    """ The plotting function for total t,r,a spectra on one plot. """

    fig = plt.figure(num=None, figsize=(9, 12), dpi=80, facecolor='w', edgecolor='k')
    ax1 = fig.add_subplot(3,1,2)
    extinciton = np.log10(1.0/np.array(t_spec))
    ax1.plot(wavelengths, extinciton)
    plt.ylim((0, 4.5))
    ax1.set_xlabel('Wavelength (nm)')
    ax1.set_ylabel('Extinciton')
    plot_name = 'extinciton'
    plt.suptitle(plot_name+add_name+'\n'+params_2_print)
    plt.savefig('%(s)s_stack%(bon)s_%(add)s'% {'s' : plot_name, 'bon' : stack_label,'add' : add_name})




def EOT_plot(stack_wl_list, wavelengths, params_2_print, add_name=''):
    """ Plot T_{00} as in Martin-Moreno PRL 86 2001. """

    height_list = stack_wl_list[0].heights_nm()[::-1]
    params_2_print += r'$h_t,...,h_b$ = '
    params_2_print += ''.join('%4d, ' % num for num in height_list)

    T_00 = []
    for i in range(len(wavelengths)): 
        t00  = stack_wl_list[i].T_net[0,0]
        t00c = t00.conjugate()
        T_00.append(np.real(t00*t00c))

    fig = plt.figure()#num=None, figsize=(9, 12), dpi=80, facecolor='w', edgecolor='k')
    ax1 = fig.add_subplot(2,1,1)
    ax1.plot(wavelengths, T_00)
    # ax1.set_xlabel('Wavelength (nm)')
    ax1.set_xlabel(r'$\lambda/a$')
    ax1.set_ylabel(r'T$_{00}$')
    # plt.ylim((0, 0.3))


    R_00 = []
    for i in range(len(wavelengths)): 
        r00  = stack_wl_list[i].R_net[0,0]
        r00c = r00.conjugate()
        R_00.append(np.real(r00*r00c))

    ax1 = fig.add_subplot(2,1,2)
    ax1.plot(wavelengths, R_00)
    # ax1.set_xlabel('Wavelength (nm)')
    ax1.set_xlabel(r'$\lambda/a$')
    ax1.set_ylabel(r'R$_{00}$')
    plt.ylim((0.2, 1.0))
    plot_name = 'EOT'
    plt.suptitle(params_2_print)
    plt.savefig('%(s)s_%(add)s'% {'s' : plot_name, 'add' : add_name})






#######################################################################################



#######################################################################################
def omega_plot(stack_wl_list, wavelengths, params_2_print, stack_label=1):
    """ Plots the dispersion diagram of each layer in one plot. """

    stack_label = zeros_int_str(stack_label)
    period = np.array(stack_wl_list[0].layers[0].structure.period)
    normed_omegas = 1/wavelengths*period

    num_layers = len(stack_wl_list[0].layers)
    fig1 = plt.figure(num=None, figsize=(9, 12), dpi=80, facecolor='w', edgecolor='k')
    fig2 = plt.figure(num=None, figsize=(9, 12), dpi=80, facecolor='w', edgecolor='k')
    fig3 = plt.figure(num=None, figsize=(9, 12), dpi=80, facecolor='w', edgecolor='k')
    fig4 = plt.figure(num=None, figsize=(9, 12), dpi=80, facecolor='w', edgecolor='k')
    for l in range(num_layers):
        ax1 = fig1.add_subplot(num_layers,1,num_layers-l)
        ax2 = fig2.add_subplot(num_layers,1,num_layers-l)
        ax3 = fig3.add_subplot(num_layers,1,num_layers-l)
        ax4 = fig4.add_subplot(num_layers,1,num_layers-l)
        for i in range(len(wavelengths)):
            k_zs = stack_wl_list[i].layers[l].k_z
            real_k_zs = []
            imag_k_zs = []
            for k in k_zs:
                if np.real(k) > np.imag(k): #alternatively np.imag(k)< small
                    real_k_zs.append(np.real(k))
                    imag_k_zs.append(np.imag(k))
            wl = np.ones(len(real_k_zs))*wavelengths[i]
            ax1.plot(wl,real_k_zs,'bo')
            wl = np.ones(len(imag_k_zs))*wavelengths[i]
            ax2.plot(wl,imag_k_zs,'ro')
            om = np.ones(len(real_k_zs))*normed_omegas[i]
            ax3.plot(real_k_zs,om,'bo')
            om = np.ones(len(imag_k_zs))*normed_omegas[i]
            ax4.plot(imag_k_zs, om,'ro')
        ax1.set_ylabel(r'Real $k_zd$'), ax2.set_ylabel(r'Imaginary $k_zd$')
        ax3.set_ylabel(r'Frequency ($\omega$d/2$\pi$c)'), ax4.set_ylabel(r'Frequency ($\omega$d/2$\pi$c)')
        if l == 0: 
            ax1.set_ylabel('Bottom Layer'), ax2.set_ylabel('Bottom Layer')
            ax3.set_ylabel('Bottom Layer'), ax4.set_ylabel('Bottom Layer')
            ax3.set_xlabel(r'Real $k_zd$'), ax4.set_xlabel(r'Imaginary $k_zd$')
            ax1.set_xlabel('Wavelength (nm)'), ax2.set_xlabel('Wavelength (nm)')
        else:
            ax1.set_xticklabels( () )
            ax2.set_xticklabels( () )
            ax3.set_xticklabels( () )
            ax4.set_xticklabels( () )
        if l == num_layers-1:
            ax1.set_ylabel('Top Layer'), ax2.set_ylabel('Top Layer')
            ax3.set_ylabel('Top Layer'), ax4.set_ylabel('Top Layer')
        ax2.set_xlim((wavelengths[0], wavelengths[-1]))
        ax1.set_xlim((wavelengths[0], wavelengths[-1]))
        # Plot the (dispersive) light line in homogeneous layers.
        if isinstance(stack_wl_list[0].layers[l], mode_calcs.Anallo):
            ns = [stack_wl_list[i].layers[l].n() for i in range(len(wavelengths))]
            ax3.plot(np.real(ns)*normed_omegas, normed_omegas,'k',linewidth=2)
    fig1.suptitle(r'Real $k_z$'+params_2_print+'\n')
    fig2.suptitle(r'Imaginary $k_z$'+params_2_print+'\n')
    fig3.suptitle(r'Real $k_z$'+params_2_print+'\n')
    fig4.suptitle(r'Imaginary $k_z$'+params_2_print+'\n')
    fig1.savefig('Disp_Diagram_Re_stack%(bon)s'% {'bon' : stack_label})
    fig2.savefig('Disp_Diagram_Im_stack%(bon)s'% {'bon' : stack_label})
    fig3.savefig('Disp_Diagram_w(k)_Re_stack%(bon)s'% {'bon' : stack_label})
    fig4.savefig('Disp_Diagram_w(k)_Im_stack%(bon)s'% {'bon' : stack_label})
    # Uncomment if you wish to save the dispersion data of a simulation to file.
    # np.savetxt('Disp_Data_stack%(bon)i.txt'% {'bon' : stack_label}, av_array, fmt = '%18.11f')

def E_conc_plot(stack_wl_list, which_layer, which_modes, wavelengths, params_2_print, stack_label=1):
    """ Plots the energy concentration (epsilon |E_{cyl}| / epsilon |E_{cell}|) of given layer. """

    stack_label = zeros_int_str(stack_label)
    if isinstance(stack_wl_list[0].layers[which_layer], mode_calcs.Simmo):
        num_layers = len(stack_wl_list[0].layers)
        fig1 = plt.figure(num=None, figsize=(9, 4), dpi=80, facecolor='w', edgecolor='k')
        ax1 = fig1.add_subplot(1,1,1)
        # Uncomment if you wish to have a continuous curve plotted.
        # ax1 = fig1.add_subplot(1,1,1)       
        # ax2 = fig1.add_subplot(2,1,2)
        # E_conc = []
        for i in range(len(wavelengths)):
            for mode in which_modes:
                E_conc_tmp = np.real(stack_wl_list[i].layers[which_layer].mode_pol[3,mode])
                # E_conc.append(E_conc_tmp)
                ax1.plot(wavelengths[i],E_conc_tmp,'bo')
        plt.xlim((wavelengths[0], wavelengths[-1]))
        ax1.set_xlabel('Wavelength (nm)')
        ax1.set_ylabel(r'$\epsilon |E_{cyl}| / \epsilon |E_{cell}|$')
        # ax2.plot(wavelengths,E_conc,'k')
        # ax2.set_xlabel('Wavelength (nm)')
        # ax2.set_ylabel(r'$E_{cyl} / E_{cell}$')
        fig1.suptitle('Energy Concentration = '+r'$E_{cyl} / E_{cell}$'+'\n'+params_2_print)
        fig1.savefig('Energy_Concentration_stack%(bon)s'% {'bon' : stack_label})
    else:
        print "\n ERROR: plotting.E_conc_plot; \n Can only calculate energy concentration in NanoStruct layers."
        print repr(stack_wl_list[0].layers[which_layer])
#######################################################################################






#######################################################################################
def vis_scat_mats(scat_mat,wl=0,extra_title=''):
    """ Plot given scattering matrix as grayscale images. """
    image = np.real(abs(scat_mat))
    plt.matshow(image,cmap=plt.cm.gray)
    print np.shape(scat_mat)
    scat_mat_dim_x = np.shape(scat_mat)[0]
    scat_mat_dim_y = np.shape(scat_mat)[1]
    half_dim_x = scat_mat_dim_x/2-0.5
    half_dim_y = scat_mat_dim_y/2-0.5
    plt.plot([-0.5, scat_mat_dim_y-0.5],[half_dim_x,half_dim_x],'r', linewidth=1)
    plt.plot([half_dim_y,half_dim_y],[-0.5, scat_mat_dim_x-0.5],'r', linewidth=1)
    plt.axis([-0.5, scat_mat_dim_y-0.5, scat_mat_dim_x-0.5,  -0.5])
    plt.xticks([half_dim_y/2, scat_mat_dim_y-half_dim_y/2],['TE', 'TM'])
    plt.yticks([half_dim_x/2, scat_mat_dim_x-half_dim_x/2],['TE', 'TM'])
    cbar = plt.colorbar(extend='neither')
    plt.xlabel('Incoming Orders')
    plt.ylabel('Outgoing Orders')
    plt.suptitle('%s Scattering Matrix' % extra_title)
    plt.savefig('Scat_mat-%(title)s_wl%(wl)i' % {'title' : extra_title, 'wl' : wl})
#######################################################################################






#######################################################################################
def k_plot(t_func_k,nu_PW_pols):
    fig = plt.figure(num=None, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
    ax1 = fig.add_subplot(1,1,1)
    ks = range(nu_PW_pols)
    np.savetxt('%(s)s.txt'% {'s' : 't_func_k',}, t_func_k, fmt = '%18.12f')
    ax1.plot(ks,t_func_k)
    ax1.set_xlabel('k vector')
    ax1.set_ylabel(r'$|E|_{trans}$')
    plt.savefig('t_func_k')

def t_func_k_plot(stack_list):
    for stack in stack_list:
        # print stack.T_net[0,:]
        nu_PW_pols = 2*stack.layers[0].structure.num_pw_per_pol
        trans_k = np.abs(stack.trans_vector)
        trans_k_array = np.array(trans_k).reshape(-1,)
        tot_trans_k_array = trans_k_array#**2
        # print repr(trans_k)
        # print repr(tot_trans_k_array)
        k_plot(tot_trans_k_array,nu_PW_pols)
#######################################################################################









# np.genfromtxt can deal with incomplete data!
#     num_BMs     = np.genfromtxt(file_name, usecols=(1))




# def efficiency_h(spec_name,h_spec_name,wavelengths,num_wl,num_h, Animate):
#     data   = np.loadtxt('%s.txt' % spec_name)
#     spec   = data[:,1]
#     i_data = np.loadtxt('%s.txt' % Irrad_spec_file)
#     i_spec = np.interp(wavelengths, i_data[:,0], i_data[:,2])
#     h_Efficiency = []
#     h_data       = []
#     for i in np.linspace(0,num_h-1,num_h):
#         h_wl = []
#         h = data[i,2]
#         for j in np.linspace(0,num_wl-1,num_wl):
#             tmp  = spec[i + j*(num_h)]
#             h_wl = np.append(h_wl,tmp)
#         eta_calc = efficiency(wavelengths,i_spec,h_wl)
#         h_Efficiency = np.append(h_Efficiency,eta_calc)
#         h_data       = np.append(h_data,data[i,2])
#         if Animate == True:
#             spectra_h(wavelengths,h_wl, h, eta_calc, i)
#     h_array = zip(h_Efficiency, h_data)
#     np.savetxt('%s.txt' % h_spec_name, h_array, fmt = '%18.12f')
#     if Animate == True:
#         ps = os.system("convert -delay 3 +dither -layers Optimize -colors 16 Animated/*.pdf Animated/Abs_Spectra.gif")
#         p  = os.system("gifsicle -O2 Animated/Abs_Spectra.gif -o Animated/Abs_Spectra-opt.gif")


# def height_plot(name_out, name_in, layer, light, max_num_BMs, max_order_PWs, Efficiency_h,
#     Efficiency,num_h):
#     fig   = plt.figure(num=None, figsize=(12, 9), dpi=80, facecolor='w', edgecolor='k')
#     ax1   = fig.add_subplot(111)

#     data    = np.loadtxt('%s.txt' % name_in)
#     wl_data = data[:,0]
#     s_data  = data[:,1]
#     h_data  = data[:,2]

#     wl_1    = wl_data[0]
#     wl_2    = wl_data[-1]
#     d_wl    = 1
#     h_1     = h_data[0]
#     h_2     = h_data[-1]
#     h_int   = (h_2-h_1)/(num_h-1)

#     # Interpolate the data
#     int_wl,int_h = np.mgrid[wl_1:wl_2+d_wl:d_wl, h_1:h_2+h_int:h_int]
#     int_Eff    = griddata(wl_data, h_data, s_data, int_wl, int_h)
#     int_Eff_t  = int_Eff.T

#     # CS = plt.contourf(int_wl,int_h,int_Eff,100,cmap=plt.cm.jet)
#     CS = plt.imshow(int_Eff_t, cmap=plt.cm.hot, norm=None, aspect='auto', interpolation=None,
#       vmin=0, vmax=s_data.max(), origin='lower', extent=(wl_1, wl_2, h_1, h_2))
#     ax1.set_xlabel('Wavelength (nm)')
#     ax1.set_ylabel('Height (nm)')
#     tick_array = []
#     for i in range(0,10,1):
#         tick = float(i)/10
#         if tick < s_data.max():
#             tick_array.append(tick)
#     tick_array.append(round(s_data.max(),4))
#     cbar = plt.colorbar(extend='neither',ticks=tick_array)
#     cbar.set_ticklabels(tick_array)
#     cbar.ax.set_ylabel('Absorptance')



# # Plot efficiency vs height
#     fig2= plt.figure(num=None, figsize=(12, 9), dpi=80, facecolor='w', edgecolor='k')
#     ax2 = fig2.add_subplot(153)
#     data = np.loadtxt('%s.txt' % Efficiency_h)
#     es = data[:,0]*100
#     hs = data[:,1]
#     ax2.plot(es, hs,'r')
#     ax2.set_xlabel(r'$\eta$' ' (%)')
#     ax2.set_ylabel('Height (nm)')
#     plt.axis([es[-1], 0, hs[0], hs[-1]])
#     plt.xticks([round(es.max(),1), 0])
#     plt.fill_betweenx(hs,es,color='blue', alpha=0.5)
#     plt.suptitle(imp_facts)
#     plt.savefig('%s_eta_bar' % name_out)









# ================    ===============================
# character           description
# ================    ===============================
# ``'-'``             solid line style
# ``'--'``            dashed line style
# ``'-.'``            dash-dot line style
# ``':'``             dotted line style
# ``'.'``             point marker
# ``','``             pixel marker
# ``'o'``             circle marker
# ``'v'``             triangle_down marker
# ``'^'``             triangle_up marker
# ``'<'``             triangle_left marker
# ``'>'``             triangle_right marker
# ``'1'``             tri_down marker
# ``'2'``             tri_up marker
# ``'3'``             tri_left marker
# ``'4'``             tri_right marker
# ``'s'``             square marker
# ``'p'``             pentagon marker
# ``'*'``             star marker
# ``'h'``             hexagon1 marker
# ``'H'``             hexagon2 marker
# ``'+'``             plus marker
# ``'x'``             x marker
# ``'D'``             diamond marker
# ``'d'``             thin_diamond marker
# ``'|'``             vline marker
# ``'_'``             hline marker
# ================    ===============================