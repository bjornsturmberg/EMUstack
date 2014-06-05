"""
    plotting.py is a subroutine of EMUstack that contains numerous plotting
    routines. These were developed during simulations for photovoltaics,
    hence the efficiency calculations.

    Copyright (C) 2013  Bjorn Sturmberg, Kokou Dossou, Felix Lawrence

    EMUstack is free software: you can redistribute it and/or modify
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
from scipy import sqrt
import subprocess
from matplotlib.mlab import griddata
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import os

# font = {'family' : 'normal',
#         'weight' : 'bold',
#         'size'   : 18}

font = {'size'   : 18}
matplotlib.rc('font', **font)
linesstrength = 2.5
title_font = 10


####Natural constants##################################################################
ASTM15_tot_I   = 900.084            # Integral of ASTM 1.5 solar irradiance in W/m**2
Plancks_h      = 6.62606957*1e-34   # Planck's constant
speed_c        = 299792458          # Speed of light in vacuum
charge_e       = 1.602176565*1e-19  # Charge of an electron
#######################################################################################


####Short utility functions############################################################
def clear_previous(file_type):
    """ Delete all files of specified type 'file_type', eg '.txt' """

    files_rm = 'rm *%s'% file_type
    subprocess.call(files_rm, shell = True)

def zeros_int_str(zero_int):
    """ Convert integer into string with '0' in place of ' ' """
    # if zero_int == 0:
    #     fmt_string = '0000'
    # else:
    #     string = '%4.0f' % zero_int
    #     fmt_string = string.replace(' ','0')
    string = '%4.0f' % zero_int
    fmt_string = string.replace(' ','0')
    return fmt_string

def tick_function(energies):
    """ Convert energy in eV into wavelengths in nm """
    wls = Plancks_h*speed_c/(energies*charge_e)*1e9
    return wls

def max_n(stacks_list):
    ns = []
    for s in stacks_list:
        for l in s.layers:
            if isinstance(l,objects.Anallo):
                ns.append(l.n())
            if isinstance(l,objects.Simmo):
                wl = l.light.wl_nm
                ns.append(l.structure.background.n(wl))
                ns.append(l.structure.inclusion_a.n(wl))
                ns.append(l.structure.inclusion_b.n(wl))
    return np.real(np.max(ns))

def gen_params_string(stack, layer=1):
    """ Generate the string of simulation info that is to be printed at the top of plots. """
    param_layer = stack[0].layers[layer].structure
    # Plot t,r,a for each layer & total, then save each to text files
    if isinstance(param_layer,objects.NanoStruct):
        params_2_print = 'ff = %5.3f, '% param_layer.ff
        params_2_print += 'd = %(period)d, a1 = %(diameter)d, '% {
        'period'        : param_layer.period, 'diameter' : param_layer.diameter1,}
        if param_layer.diameter2 != 0: params_2_print += 'a2 = %(rad)d '% {'rad' : param_layer.diameter2,}
        if param_layer.geometry == '2D_array':
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
        elif param_layer.geometry == '1D_array':
            params_2_print += ''
        params_2_print += '%(BMs)dBMs, PW_radius = %(PWs)d, '% {
        'BMs' : stack[0].layers[layer].num_BM,'PWs' : stack[0].layers[layer].max_order_PWs, }
    else:
        params_2_print = 'PW_radius = %(PWs)d, '% {'PWs' : stack[0].layers[layer].max_order_PWs, }

    # light = stack[0].layers[layer].light
    # k_pll = light.k_pll * param_layer.period
    # params_2_print += r'$k_\parallel d$ = '
    # tmp3 = '(%(kx)1.4f, %(ky)1.4f), '% {'kx' : k_pll[0], 'ky' : k_pll[1]}
    # params_2_print += r'$\theta$ = %(theta)6.2f, $\phi$ = %(phi)6.2f, '% {
    # 'theta' : light.theta,'phi' : light.phi, } 

    return params_2_print
#######################################################################################


####Standard plotting of spectra#######################################################
def t_r_a_plots(stacks_list, xvalues=None, params_layer=1, active_layer_nu=0,\
    stack_label=1, ult_eta=False, J_sc=False, weight_spec=False, extinct=False,\
    add_height=0, add_name='', force_txt_save=False):
    """ Plot t,r,a for each layer & total, then save each to text files. """

    height_list = stacks_list[0].heights_nm()[::-1]
    params_2_print = gen_params_string(stacks_list, params_layer)
    params_2_print += '\n'r'$h_t,...,h_b$ = '
    params_2_print += ''.join('%4d, ' % num for num in height_list)

    if stacks_list[0].layers[0].light.wl_nm != stacks_list[-1].layers[0].light.wl_nm:
        xvalues = [s.layers[0].light.wl_nm for s in stacks_list]
        xlabel = r'$\lambda$ (nm)'
    elif set(stacks_list[0].layers[0].light.k_pll) != set(stacks_list[-1].layers[0].light.k_pll):
        xvalues = [np.sqrt(s.layers[0].light.k_pll[0]**2 + s.layers[0].light.k_pll[1]**2) for s in stacks_list]
        xlabel = r'$|k_\parallel|$'
    elif xvalues==None: 
        "t_r_a_plots cannot guess what to plot on x-axis, specify with xvalues input."
        return


    if add_height!=0: add_name += zeros_int_str(add_height)
    stack_label = zeros_int_str(stack_label)
    
    a_list = []
    t_list = []
    r_list = []
    for stack in stacks_list:
        a_list.extend(stack.a_list)
        t_list.extend(stack.t_list)
        r_list.extend(stack.r_list)

    layers_steps = len(stacks_list[0].layers) - 1
    active_abs = []
    a_tot      = []
    t_tot      = []
    r_tot      = []
    for i in range(len(xvalues)): 
        active_abs.append(float(a_list[active_layer_nu + i*layers_steps]))
        a_tot.append(float(a_list[layers_steps-1+(i*layers_steps)]))
        t_tot.append(float(t_list[layers_steps-1+(i*layers_steps)]))
        r_tot.append(float(r_list[i]))
    
    if ult_eta == True:
        Efficiency = ult_efficiency(active_abs, xvalues, params_2_print, stack_label, add_name)
        params_2_print += r'$\eta$ = %(Efficiency)6.2f'% {'Efficiency' : Efficiency*100, }
        params_2_print += ' %'

    if J_sc == True:
        J = J_short_circuit(active_abs, xvalues, params_2_print, stack_label, add_name)
        params_2_print += r'$J_{sc}$ = %(J)6.2f'% {'J' : J, }
        params_2_print += r' mA/cm$^2$'

    total_h = sum(stacks_list[0].heights_nm()) # look at first wl result to find h.
    # Plot t,r,a for each layer
    layers_plot('Lay_Absorb', a_list, xvalues, xlabel,total_h,params_2_print,stack_label,add_name,force_txt_save)
    layers_plot('Lay_Trans',  t_list, xvalues, xlabel,total_h,params_2_print,stack_label,add_name,force_txt_save)
    layers_plot('Lay_Reflec', r_list, xvalues, xlabel,total_h,params_2_print,stack_label,add_name,force_txt_save)
    # Also plot total t,r,a on a single plot
    plot_name = 'Total_Spectra'
    total_tra_plot(plot_name, a_tot, t_tot, r_tot, xvalues, xlabel, params_2_print, stack_label, add_name)

    if weight_spec == True:
        # Also plot totals weighted by solar irradiance
        Irrad_spec_file = '../backend/data/ASTMG173'
        i_data       = np.loadtxt('%s.txt' % Irrad_spec_file)
        i_spec       = np.interp(xvalues, i_data[:,0], i_data[:,3])
        bandgap_wl   = xvalues[-1]
        weighting = i_spec/i_spec.max()*(xvalues/bandgap_wl)
        a_weighted = a_tot*weighting
        t_weighted = t_tot*weighting
        r_weighted = r_tot*weighting
        plot_name = 'Weighted_Total_Spectra'
        total_tra_plot(plot_name, a_weighted, t_weighted, r_weighted, xvalues, 
            xlabel, params_2_print, stack_label, add_name)

    if  extinct == True:
        extinction_plot(t_tot, xvalues, params_2_print, stack_label, add_name)
    return

def layers_plot(spectra_name, spec_list, xvalues, xlabel, total_h, params_2_print,\
    stack_label, add_name, force_txt_save):
    """ The plotting function for one type of spectrum across all layers. """

    fig = plt.figure(num=None, figsize=(9, 12), dpi=80, facecolor='w', edgecolor='k')
    nu_layers = len(spec_list)/len(xvalues)
    h_array = np.ones(len(xvalues))*total_h
    for i in range(nu_layers):
        layer_spec = []
        for wl in range(len(xvalues)):
            layer_spec = np.append(layer_spec,spec_list[wl*nu_layers + i])
        ax1 = fig.add_subplot(nu_layers,1,i+1)
        ax1.plot(xvalues, layer_spec, linewidth=linesstrength)
        ax2 = ax1.twiny()
        new_tick_values = np.linspace(10,0.5,20)
        new_tick_locations = tick_function(new_tick_values)
        new_tick_labels = ["%.1f" % z for z in new_tick_values]
        ax2.set_xticks(new_tick_locations)
        ax2.set_xticklabels(new_tick_labels)
        ax2.set_xlim((xvalues[0], xvalues[-1]))  
        if i == 0:
            if xlabel == r'$\lambda$ (nm)': ax2.set_xlabel('Energy (eV)')
        elif i == nu_layers-1:
            ax1.set_xlabel(xlabel)
            ax1.set_ylabel('Total')
        if i != 0:
            ax2.set_xticklabels( () )
        if i != nu_layers-1:
            ax1.set_xticklabels( () )
        if spectra_name == 'Lay_Absorb':
            if i == 0: ax1.set_ylabel('Top Layer')
            if i == nu_layers-2: ax1.set_ylabel('Bottom Layer')
            suptitle_w_params = 'Absorptance in each layer'+add_name+'\n'+params_2_print
            plt.suptitle(suptitle_w_params,fontsize=title_font)
            lay_spec_name = 'Lay_Absorb'
            if i == nu_layers-1: 
                ax1.set_ylabel('Total')
                lay_spec_name = 'Absorptance'
        elif spectra_name == 'Lay_Trans':
            if i == 0: ax1.set_ylabel('Top Layer')
            if i == nu_layers-2: ax1.set_ylabel('Bottom Layer')
            suptitle_w_params = 'Transmittance in each layer'+add_name+'\n'+params_2_print
            plt.suptitle(suptitle_w_params,fontsize=title_font)
            lay_spec_name = 'Lay_Trans'
            if i == nu_layers-1: 
                ax1.set_ylabel('Total')
                lay_spec_name = 'Transmittance'
        elif spectra_name == 'Lay_Reflec':
            if i == 0: ax1.set_ylabel('Top Layer')
            if i == nu_layers-3: ax1.set_ylabel('Bottom Layer')
            if i == nu_layers-2: ax1.set_ylabel('Substrate')
            suptitle_w_params = 'Reflectance in each layer'+add_name+'\n'+params_2_print
            plt.suptitle(suptitle_w_params,fontsize=title_font)
            lay_spec_name = 'Lay_Reflec'
            if i == nu_layers-1: 
                ax1.set_ylabel('Total')
                lay_spec_name = 'Reflectance'
        av_array = zip(xvalues, layer_spec, h_array)
        ax1.set_xlim((xvalues[0], xvalues[-1]))
        plt.ylim((0, 1))

        if force_txt_save == True:
            if i != nu_layers-1: 
                np.savetxt('%(s)s_%(i)i_stack%(bon)s%(add)s.txt'% {'s' : lay_spec_name, 'i' : i, 
                    'bon' : stack_label,'add' : add_name}, av_array, fmt = '%18.11f')
            else:
                np.savetxt('%(s)s_stack%(bon)s%(add)s.txt'% {'s' : lay_spec_name, 
                    'bon' : stack_label,'add' : add_name}, av_array, fmt = '%18.11f')

        plt.savefig('%(s)s_stack%(bon)s_%(add)s'% {'s': spectra_name, 'bon': stack_label,\
            'add' : add_name})

def total_tra_plot(plot_name, a_spec, t_spec, r_spec, xvalues, xlabel, params_2_print,\
    stack_label, add_name):
    """ The plotting function for total t,r,a spectra on one plot. """

    fig = plt.figure(num=None, figsize=(9, 12), dpi=80, facecolor='w', edgecolor='k')
    
    ax1 = fig.add_subplot(3,1,1)
    ax1.plot(xvalues, a_spec, linewidth=linesstrength)
    ax1.set_ylabel('Absorptance')
    ax1.set_xlim((xvalues[0], xvalues[-1]))
    ax2 = ax1.twiny()
    new_tick_values = np.linspace(10,0.5,20)
    new_tick_locations = tick_function(new_tick_values)
    new_tick_labels = ["%.1f" % z for z in new_tick_values]
    ax2.set_xticks(new_tick_locations)
    ax2.set_xticklabels(new_tick_labels)
    ax2.set_xlim((xvalues[0], xvalues[-1])) 
    if xlabel == r'$\lambda$ (nm)': ax2.set_xlabel('Energy (eV)')
    plt.ylim((0, 1))

    ax1 = fig.add_subplot(3,1,2)
    ax1.plot(xvalues, t_spec, linewidth=linesstrength)
    ax1.set_ylabel('Transmittance')
    ax1.set_xlim((xvalues[0], xvalues[-1]))
    ax2 = ax1.twiny()
    ax2.set_xticklabels( () )
    ax2.set_xticks(new_tick_locations)
    ax2.set_xlim((xvalues[0], xvalues[-1])) 
    plt.ylim((0, 1))

    ax1 = fig.add_subplot(3,1,3)
    ax1.plot(xvalues, r_spec, linewidth=linesstrength)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel('Reflectance')
    ax1.set_xlim((xvalues[0], xvalues[-1]))
    ax2 = ax1.twiny()
    ax2.set_xticklabels( () )
    ax2.set_xticks(new_tick_locations)
    ax2.set_xlim((xvalues[0], xvalues[-1])) 
    plt.ylim((0, 1))

    plt.suptitle(params_2_print,fontsize=title_font)
    # plt.suptitle(plot_name+add_name+'\n'+params_2_print)
    plt.savefig('%(s)s_stack%(bon)s_%(add)s'% {'s' : plot_name, 'bon' : stack_label,'add' : add_name})
#######################################################################################


####Plot spectra indicating Wood anomalies in substrate################################
def t_r_a_plots_subs(stacks_list, wavelengths, period, sub_n, params_layer=1, \
    active_layer_nu=0, stack_label=1, ult_eta=False, J_sc=False, weight_spec=False,\
    extinct=False, add_height=0, add_name=''):
    """ Plot t,r,a indicating Wood anomalies in substrate for each layer & total,
        then save each to text files. """

    height_list = stacks_list[0].heights_nm()[::-1]
    params_2_print = gen_params_string(stacks_list, params_layer)
    params_2_print += '\n'r'$h_t,...,h_b$ = '
    params_2_print += ''.join('%4d, ' % num for num in height_list)

    stack_label = zeros_int_str(stack_label)

    a_list = []
    t_list = []
    r_list = []
    for stack in stacks_list:
        a_list.extend(stack.a_list)
        t_list.extend(stack.t_list)
        r_list.extend(stack.r_list)

    layers_steps = len(stacks_list[0].layers) - 1
    active_abs = []
    a_tot      = []
    t_tot      = []
    r_tot      = []
    for i in range(len(wavelengths)): 
        active_abs.append(float(a_list[active_layer_nu + i*layers_steps]))
        a_tot.append(float(a_list[layers_steps-1+(i*layers_steps)]))
        t_tot.append(float(t_list[layers_steps-1+(i*layers_steps)]))
        r_tot.append(float(r_list[i]))

    if ult_eta == True:
        Efficiency = ult_efficiency(active_abs, wavelengths, params_2_print, stack_label,add_name)
        params_2_print += r'$\eta$ = %(Efficiency)6.2f'% {'Efficiency' : Efficiency*100, }
        params_2_print += ' %'

    if J_sc == True:
        J = J_short_circuit(active_abs, wavelengths, params_2_print, stack_label, add_name)
        params_2_print += r'$J_{sc}$ = %(J)6.2f'% {'J' : J, }
        params_2_print += r' mA/cm$^2$'

    # Plot total t,r,a on a single plot indicating Wood anomalies
    plot_name = 'Total_Spectra_subs'
    total_tra_plot_subs(plot_name, a_tot, t_tot, r_tot, wavelengths, params_2_print,
      stack_label, add_name, period, sub_n)

    if weight_spec == True:
        # Also plot totals weighted by solar irradiance
        Irrad_spec_file = '../backend/data/ASTMG173'
        i_data       = np.loadtxt('%s.txt' % Irrad_spec_file)
        i_spec       = np.interp(wavelengths, i_data[:,0], i_data[:,3])
        bandgap_wl   = wavelengths[-1]
        weighting = i_spec/i_spec.max()*(wavelengths/bandgap_wl)
        a_weighted = a_tot*weighting
        t_weighted = t_tot*weighting
        r_weighted = r_tot*weighting
        plot_name = 'Weighted_Total_Spectra_subs'
        total_tra_plot_subs(plot_name, a_weighted, t_weighted, r_weighted, wavelengths, 
            params_2_print, stack_label, add_name, period, sub_n)

    if  extinct == True:
        extinction_plot(t_tot, wavelengths, params_2_print, stack_label, add_name)
    return

def total_tra_plot_subs(plot_name, a_spec, t_spec, r_spec, wavelengths, \
    params_2_print, stack_label, add_name, period, sub_n):
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
    ax1.plot(sub_line01, v_line, 'k', linewidth=linesstrength)
    ax1.plot(sub_line11, v_line, 'k', linewidth=linesstrength)
    ax1.plot(sub_line20, v_line, 'k', linewidth=linesstrength)
    ax1.plot(sub_line21, v_line, 'k', linewidth=linesstrength)
    ax1.plot(sub_line22, v_line, 'k', linewidth=linesstrength)
    ax1.plot(sub_line30, v_line, 'k', linewidth=linesstrength)
    ax1.plot(sup_line01, v_line, 'r', linewidth=linesstrength)
    ax1.plot(sup_line11, v_line, 'r', linewidth=linesstrength)
    ax1.plot(sup_line20, v_line, 'r', linewidth=linesstrength)
    ax1.plot(sup_line21, v_line, 'r', linewidth=linesstrength)
    ax1.plot(sup_line22, v_line, 'r', linewidth=linesstrength)
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
    ax1.plot(sub_line01, v_line, 'k', linewidth=linesstrength)
    ax1.plot(sub_line11, v_line, 'k', linewidth=linesstrength)
    ax1.plot(sub_line20, v_line, 'k', linewidth=linesstrength)
    ax1.plot(sub_line21, v_line, 'k', linewidth=linesstrength)
    ax1.plot(sub_line22, v_line, 'k', linewidth=linesstrength)
    ax1.plot(sub_line30, v_line, 'k', linewidth=linesstrength)
    ax1.plot(sup_line01, v_line, 'r', linewidth=linesstrength)
    ax1.plot(sup_line11, v_line, 'r', linewidth=linesstrength)
    ax1.plot(sup_line20, v_line, 'r', linewidth=linesstrength)
    ax1.plot(sup_line21, v_line, 'r', linewidth=linesstrength)
    ax1.plot(sup_line22, v_line, 'r', linewidth=linesstrength)
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
    ax1.plot(sub_line01, v_line, 'k', linewidth=linesstrength)
    ax1.plot(sub_line11, v_line, 'k', linewidth=linesstrength)
    ax1.plot(sub_line20, v_line, 'k', linewidth=linesstrength)
    ax1.plot(sub_line21, v_line, 'k', linewidth=linesstrength)
    ax1.plot(sub_line22, v_line, 'k', linewidth=linesstrength)
    ax1.plot(sub_line30, v_line, 'k', linewidth=linesstrength)
    ax1.plot(sup_line01, v_line, 'r', linewidth=linesstrength)
    ax1.plot(sup_line11, v_line, 'r', linewidth=linesstrength)
    ax1.plot(sup_line20, v_line, 'r', linewidth=linesstrength)
    ax1.plot(sup_line21, v_line, 'r', linewidth=linesstrength)
    ax1.plot(sup_line22, v_line, 'r', linewidth=linesstrength)
    ax1.set_xlabel('Wavelength (nm)')
    ax1.set_ylabel('Reflectance')
    ax1.set_xlim((wavelengths[0], wavelengths[-1]))
    ax2 = ax1.twiny()
    ax2.set_xticklabels( () )
    ax2.set_xticks(new_tick_locations)
    ax2.set_xlim((wavelengths[0], wavelengths[-1])) 
    plt.ylim((0, 1))
    plt.suptitle(params_2_print,fontsize=title_font)
    # plt.suptitle(plot_name+add_name+'\n'+params_2_print)
    plt.savefig('%(s)s_stack%(bon)s__%(add)s'% {'s' : plot_name, 'bon' : stack_label,'add' : add_name})
#######################################################################################


####Save J_sc & ult efficiency w/o spectra#############################################
def J_sc_eta_NO_plots(stacks_list, wavelengths, params_layer=1, active_layer_nu=0,\
    stack_label=1, add_name=''):
    """ Calculate J_sc & efficiency but do not save or plot spectra. """

    height_list = stacks_list[0].heights_nm()[::-1]
    params_2_print = gen_params_string(stacks_list, params_layer)
    params_2_print += '\n'r'$h_t,...,h_b$ = '
    params_2_print += ''.join('%4d, ' % num for num in height_list)

    stack_label = zeros_int_str(stack_label)

    a_list = []
    for stack in stacks_list:
        a_list.extend(stack.a_list)

    layers_steps = len(stacks_list[0].layers) - 1
    active_abs = []
    for i in range(len(wavelengths)): 
        active_abs.append(float(a_list[active_layer_nu + i*layers_steps]))

    Efficiency = ult_efficiency(active_abs, wavelengths, params_2_print, stack_label, add_name)

    J = J_short_circuit(active_abs, wavelengths, params_2_print, stack_label, add_name)
    return
#######################################################################################


####Saving spectra to files############################################################
def layers_print(spectra_name, spec_list, wavelengths, total_h, stack_label=1,\
    add_name=''):
    """ Save spectra to text files. """

    nu_layers = len(spec_list)/len(wavelengths)
    h_array = np.ones(len(wavelengths))*total_h
    for i in range(nu_layers):
        layer_spec = []
        for wl in range(len(wavelengths)):
            layer_spec = np.append(layer_spec,spec_list[wl*nu_layers + i])
        if spectra_name == 'Lay_Absorb':
            lay_spec_name = 'Lay_Absorb'
            if i == nu_layers-1: 
                lay_spec_name = 'Absorptance'
        elif spectra_name == 'Lay_Trans':
            lay_spec_name = 'Lay_Trans'
            if i == nu_layers-1: 
                lay_spec_name = 'Transmittance'
        elif spectra_name == 'Lay_Reflec':
            lay_spec_name = 'Lay_Reflec'
            if i == nu_layers-1: 
                lay_spec_name = 'Reflectance'
        av_array = zip(wavelengths, layer_spec, h_array)

        if i != nu_layers-1: 
            np.savetxt('%(s)s_%(i)i_stack%(bon)s%(add)s.txt'% {'s' : lay_spec_name, 'i' : i, 
                'bon' : stack_label,'add' : add_name}, av_array, fmt = '%18.11f')
        else:
            np.savetxt('%(s)s_stack%(bon)s%(add)s.txt'% {'s' : lay_spec_name, 
                'bon' : stack_label,'add' : add_name}, av_array, fmt = '%18.11f')

def t_r_a_write_files(stacks_list, wavelengths, stack_label=1, add_name=''):
    """ Save t,r,a for each layer & total in text files. """
    stack_label = zeros_int_str(stack_label)

    a_list = []
    t_list = []
    r_list = []
    for stack in stacks_list:
        a_list.extend(stack.a_list)
        t_list.extend(stack.t_list)
        r_list.extend(stack.r_list)

    layers_steps = len(stacks_list[0].layers) - 1
    a_tot      = []
    t_tot      = []
    r_tot      = []
    for i in range(len(wavelengths)): 
        a_tot.append(float(a_list[layers_steps-1+(i*layers_steps)]))
        t_tot.append(float(t_list[layers_steps-1+(i*layers_steps)]))
        r_tot.append(float(r_list[i]))

    total_h = sum(stacks_list[0].heights_nm()) #look at first wl result to find h.
    layers_print('Lay_Absorb', a_list, wavelengths, total_h, stack_label)
    layers_print('Lay_Trans',  t_list, wavelengths, total_h, stack_label)
    layers_print('Lay_Reflec', r_list, wavelengths, total_h, stack_label)
#######################################################################################


####Plot spectra on other scales#######################################################
def extinction_plot(t_spec, wavelengths, params_2_print, stack_label, add_name):
    """ Plot extinction ratio in transmission extinct = log_10(1/t). """

    fig = plt.figure(num=None, figsize=(9, 12), dpi=80, facecolor='w', edgecolor='k')
    ax1 = fig.add_subplot(1,1,1)
    extinciton = np.log10(1.0/np.array(t_spec))
    ax1.plot(wavelengths, extinciton, linewidth=linesstrength)
    ax1.set_xlabel('Wavelength (nm)')
    ax1.set_ylabel('Extinction')
    plot_name = 'Extinction'
    plt.suptitle(plot_name+add_name+'\n'+params_2_print,fontsize=title_font)
    plt.savefig('%(s)s_stack%(bon)s_%(add)s'% {'s' : plot_name, 'bon' : stack_label,'add' : add_name})

def EOT_plot(stacks_list, wavelengths, params_layer=1, num_pw_per_pol=0, add_name=''):
    """ Plot T_{00} as in Martin-Moreno PRL 86 2001. 
        To plot {9,0} component of TM polarisation set num_pw_per_pol = num_pw_per_pol.
    """

    height_list = stacks_list[0].heights_nm()[::-1]
    params_2_print = gen_params_string(stacks_list, params_layer)
    params_2_print += r'$h_t,...,h_b$ = '
    params_2_print += ''.join('%4d, ' % num for num in height_list)

    T_00 = []
    for i in range(len(wavelengths)): 
        t00  = stacks_list[i].T_net[num_pw_per_pol,num_pw_per_pol]
        t00c = t00.conjugate()
        T_00.append(np.real(t00*t00c))

    fig = plt.figure()#num=None, figsize=(9, 12), dpi=80, facecolor='w', edgecolor='k')
    ax1 = fig.add_subplot(2,1,1)
    ax1.plot(wavelengths, T_00, linewidth=linesstrength)
    # ax1.set_xlabel('Wavelength (nm)')
    ax1.set_xlabel(r'$\lambda/a$')
    ax1.set_ylabel(r'T$_{00}$')
    # plt.ylim((0, 0.3))

    R_00 = []
    for i in range(len(wavelengths)): 
        r00  = stacks_list[i].R_net[num_pw_per_pol,num_pw_per_pol]
        r00c = r00.conjugate()
        R_00.append(np.real(r00*r00c))

    ax1 = fig.add_subplot(2,1,2)
    ax1.plot(wavelengths, R_00, linewidth=linesstrength)
    # ax1.set_xlabel('Wavelength (nm)')
    ax1.set_xlabel(r'$\lambda/a$')
    ax1.set_ylabel(r'R$_{00}$')
    # plt.ylim((0.2, 1.0))
    plot_name = 'EOT'
    plt.suptitle(params_2_print,fontsize=title_font)
    plt.savefig('%(s)s_%(add)s'% {'s' : plot_name, 'add' : add_name})
#######################################################################################


####Calculate short circuit current and ultimate efficiency############################
def J_short_circuit(active_abs, wavelengths, params_2_print, stack_label, add_name):
    """ Calculate the short circuit current J_sc under ASTM 1.5 illumination.
        Assuming every absorbed photon produces a pair of charge carriers.
    """

    Irrad_spec_file = '../backend/data/ASTMG173'
    i_data       = np.loadtxt('%s.txt' % Irrad_spec_file)
    i_spec       = np.interp(wavelengths, i_data[:,0], i_data[:,3])
    expression   = i_spec*active_abs*wavelengths
    integral_tmp = np.trapz(expression, x=wavelengths)
    J = (charge_e/(Plancks_h*speed_c)) * integral_tmp *1e-10 # in mA/cm^2  
    nums_2_print = params_2_print.split()
    if len(nums_2_print) >= 8:
        eta_string   = '%8.6f \n'% J + nums_2_print[5].replace(',','\n') + \
          nums_2_print[8].replace(',','\n') #save params in easy to read in fmt
    else:
        eta_string   = '%8.6f \n'% J
    np.savetxt('J_sc_stack%(bon)s%(add)s.txt'% {'bon' : stack_label,'add' : add_name}, np.array([eta_string]), fmt = '%s')
    return J

def ult_efficiency(active_abs, wavelengths, params_2_print, stack_label,add_name):
    """ Calculate the photovoltaic ultimate efficiency achieved in the specified active layer. 
        For definition see Sturmberg et al., Optics Express, Vol. 19, Issue S5, pp. A1067-A1081 (2011)
        http://dx.doi.org/10.1364/OE.19.0A1067
    """
    # TODO make E_g a property of material, not just longest wl included.

    Irrad_spec_file = '../backend/data/ASTMG173'
    i_data       = np.loadtxt('%s.txt' % Irrad_spec_file)
    i_spec       = np.interp(wavelengths, i_data[:,0], i_data[:,3])
    bandgap_wl   = wavelengths[-1] #have as property of material.
    expression   = i_spec*active_abs*wavelengths
    integral_tmp = np.trapz(expression, x=wavelengths)
    Efficiency   = integral_tmp/(bandgap_wl*ASTM15_tot_I)   
    nums_2_print = params_2_print.split()
    if len(nums_2_print) >= 8:
        eta_string   = '%8.6f \n'% Efficiency + nums_2_print[5].replace(',','\n') + \
          nums_2_print[8].replace(',','\n') #save params in easy to read in fmt
    else:
        eta_string   = '%8.6f \n'% Efficiency
    np.savetxt('Efficiency_stack%(bon)s%(add)s.txt'% {'bon' : stack_label,'add' : add_name}, np.array([eta_string]), fmt = '%s')
    return Efficiency
#######################################################################################


####Plot dispersion diagrams & field concentrations as function of wavelength##########
def omega_plot(stacks_list, wavelengths, params_layer=1, stack_label=1):
    """ Plots the dispersion diagram of each layer in one plot. """

    params_2_print = gen_params_string(stacks_list, params_layer)
    stack_label = zeros_int_str(stack_label)
    period = np.array(stacks_list[0].layers[0].structure.period)
    normed_omegas = 1/wavelengths*period

    num_layers = len(stacks_list[0].layers)
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
            k_zs = stacks_list[i].layers[l].k_z
            real_k_zs = []
            imag_k_zs = []
            for k in k_zs:
                if np.real(k) > np.imag(k): #alternatively np.imag(k)< small
                    real_k_zs.append(np.real(k))
                    imag_k_zs.append(np.imag(k))
            wl = np.ones(len(real_k_zs))*wavelengths[i]
            ax1.plot(wl,real_k_zs,'bo', linewidth=linesstrength)
            wl = np.ones(len(imag_k_zs))*wavelengths[i]
            ax2.plot(wl,imag_k_zs,'ro', linewidth=linesstrength)
            om = np.ones(len(real_k_zs))*normed_omegas[i]
            ax3.plot(real_k_zs,om,'bo', linewidth=linesstrength)
            om = np.ones(len(imag_k_zs))*normed_omegas[i]
            ax4.plot(imag_k_zs, om,'ro', linewidth=linesstrength)
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
        if isinstance(stacks_list[0].layers[l], mode_calcs.Anallo):
            ns = [stacks_list[i].layers[l].n() for i in range(len(wavelengths))]
            ax3.plot(np.real(ns)*normed_omegas, normed_omegas,'k', linewidth=linesstrength)
    fig1.suptitle(r'Real $k_z$'+params_2_print+'\n',fontsize=title_font)
    fig2.suptitle(r'Imaginary $k_z$'+params_2_print+'\n',fontsize=title_font)
    fig3.suptitle(r'Real $k_z$'+params_2_print+'\n',fontsize=title_font)
    fig4.suptitle(r'Imaginary $k_z$'+params_2_print+'\n',fontsize=title_font)
    fig1.savefig('Disp_Diagram_Re_stack%(bon)s'% {'bon' : stack_label}, bbox_inches='tight')
    fig2.savefig('Disp_Diagram_Im_stack%(bon)s'% {'bon' : stack_label}, bbox_inches='tight')
    fig3.savefig('Disp_Diagram_w(k)_Re_stack%(bon)s'% {'bon' : stack_label}, bbox_inches='tight')
    fig4.savefig('Disp_Diagram_w(k)_Im_stack%(bon)s'% {'bon' : stack_label}, bbox_inches='tight')
    # Uncomment if you wish to save the dispersion data of a simulation to file.
    # np.savetxt('Disp_Data_stack%(bon)i.txt'% {'bon' : stack_label}, av_array, fmt = '%18.11f')

def E_conc_plot(stacks_list, which_layer, which_modes, wavelengths, params_layer=1,\
    stack_label=1):
    """ Plots the energy concentration (epsilon |E_{cyl}| / epsilon |E_{cell}|) of given layer. """

    params_2_print = gen_params_string(stacks_list, params_layer)
    stack_label = zeros_int_str(stack_label)
    if isinstance(stacks_list[0].layers[which_layer], mode_calcs.Simmo):
        num_layers = len(stacks_list[0].layers)
        fig1 = plt.figure(num=None, figsize=(9, 4), dpi=80, facecolor='w', edgecolor='k')
        ax1 = fig1.add_subplot(1,1,1)
        # Uncomment if you wish to have a continuous curve plotted.
        # ax1 = fig1.add_subplot(1,1,1)       
        # ax2 = fig1.add_subplot(2,1,2)
        # E_conc = []
        for i in range(len(wavelengths)):
            for mode in which_modes:
                E_conc_tmp = np.real(stacks_list[i].layers[which_layer].mode_pol[3,mode])
                # E_conc.append(E_conc_tmp)
                ax1.plot(wavelengths[i],E_conc_tmp,'bo', linewidth=linesstrength)
        plt.xlim((wavelengths[0], wavelengths[-1]))
        ax1.set_xlabel('Wavelength (nm)')
        ax1.set_ylabel(r'$\epsilon |E_{cyl}| / \epsilon |E_{cell}|$')
        # ax2.plot(wavelengths,E_conc,'k', linewidth=linesstrength)
        # ax2.set_xlabel('Wavelength (nm)')
        # ax2.set_ylabel(r'$E_{cyl} / E_{cell}$')
        fig1.suptitle('Energy Concentration = '+r'$E_{cyl} / E_{cell}$'+'\n'+params_2_print,fontsize=title_font)
        fig1.savefig('Energy_Concentration_stack%(bon)s'% {'bon' : stack_label}, bbox_inches='tight')
    else:
        print "\n ERROR: plotting.E_conc_plot; \n Can only calculate energy concentration in NanoStruct layers."
        print repr(stacks_list[0].layers[which_layer])
#######################################################################################


####Visualise scattering matrices######################################################
def vis_scat_mats(scat_mat,nu_prop_PWs=0,wl=None,extra_title=None):
    """ Plot given scattering matrix as greyscale images. """
    fig = plt.figure(num=None, figsize=(14, 6), dpi=80, facecolor='w', edgecolor='k')
    if wl != None and extra_title == None: title = '-wl%(wl)i' % {'wl' : wl}
    elif wl == None and extra_title != None: title = '-%(ti)s' % {'ti' : extra_title}
    elif wl != None and extra_title != None: title = '-wl%(wl)i-%(ti)s' % {'wl' : wl, 'ti' : extra_title}
    else: title = ''

    for i in [1,2]:
        ax1 = fig.add_subplot(1,2,i)
        if i==0: image = abs(np.abs(scat_mat))
        if i==1: image = abs(np.real(scat_mat))
        mat = ax1.matshow(image,cmap=plt.cm.gray)
        scat_mat_dim_x = np.shape(scat_mat)[0]
        scat_mat_dim_y = np.shape(scat_mat)[1]
        half_dim_x = scat_mat_dim_x/2-0.5
        half_dim_y = scat_mat_dim_y/2-0.5
        ax1.plot([-0.5, scat_mat_dim_y-0.5],[half_dim_x,half_dim_x],'w', linewidth=2)
        ax1.plot([half_dim_y,half_dim_y],[-0.5, scat_mat_dim_x-0.5],'w', linewidth=2)
        ax1.axis([-0.5, scat_mat_dim_y-0.5, scat_mat_dim_x-0.5,  -0.5])
        ax1.set_xticks([half_dim_y/2, scat_mat_dim_y-half_dim_y/2],['TE', 'TM'])
        ax1.set_yticks([half_dim_x/2, scat_mat_dim_x-half_dim_x/2],['TE', 'TM'])
        # proping = []
        # half_k = k_array[0:len(k_array)/2]
        # for i in range(len(half_k)):
        #     if np.real(half_k[i])>0: proping.append(i)
        # print max(proping)
        ax1.plot([-0.5, scat_mat_dim_y-0.5],[nu_prop_PWs-0.5,nu_prop_PWs-0.5],'r', linewidth=1)
        ax1.plot([nu_prop_PWs-0.5,nu_prop_PWs-0.5],[-0.5, scat_mat_dim_x-0.5],'r', linewidth=1)
        ax1.plot([-0.5, scat_mat_dim_y-0.5],[half_dim_x+nu_prop_PWs,half_dim_x+nu_prop_PWs],'r', linewidth=1)
        ax1.plot([half_dim_y+nu_prop_PWs,half_dim_y+nu_prop_PWs],[-0.5, scat_mat_dim_x-0.5],'r', linewidth=1)
        cbar = fig.colorbar(mat,extend='neither')
        # if i==1: cbar.ax.set_ylabel('|Re(matrix)|')
        # if i==2: cbar.ax.set_ylabel('|Imag(matrix)|')
        if i==1: ax1.set_title('|Re(matrix)|')
        if i==2: ax1.set_title('|Imag(matrix)|')
        ax1.set_xticklabels('')
        ax1.set_yticklabels('')
        ax1.set_xlabel('Incoming Orders')
        ax1.set_ylabel('Outgoing Orders')
        

    plt.suptitle('Scattering Matrices' + title)
    plt.savefig('Scat_mat' + title)
#######################################################################################


####Plot PW amplitudes function k-vector###############################################
def t_func_k_plot_1D(stacks_list, lay_interest=0, pol='TE'):
    """ Plot PW amplitudes as a function of their in-plane k-vector. """
    fig = plt.figure(num=None, dpi=80, facecolor='w', edgecolor='k')
    ax1 = fig.add_subplot(1,1,1)

    # Create arrays of grating order indexes (-p, ..., p)
    max_ords = stacks_list[0].layers[-1].max_order_PWs
    pxs = np.arange(-max_ords, max_ords + 1)

    # vec_coef sorted from top of stack, everything else is sorted from bottom
    vec_index = len(stacks_list[0].layers) - lay_interest - 1

    store_alphas = []
    store_k_trans = []
    for stack in stacks_list:
        k0 = stack.layers[0].k()
        n_PW_p_pols = stack.layers[0].structure.num_pw_per_pol
        # Calculate k_x that correspond to k_y = beta0 = 0 (in normalized units)
        alpha0, beta0 = stack.layers[0].k_pll_norm()
        alphas = alpha0 + pxs * 2 * np.pi
        on_axis_kzs = sqrt(k0**2 - alphas**2 - beta0**2)
        full_k_space = stack.layers[0].k_z
        # consider singular polarization
        one_pol_k_space = full_k_space[0:n_PW_p_pols]

        axis_indices = []
        for a in on_axis_kzs:
            ix = np.in1d(one_pol_k_space.ravel(), a).reshape(one_pol_k_space.shape)
            axis_indices = np.append(axis_indices, np.ravel(np.array(np.where(ix))))
        axis_indices = axis_indices.astype(int)

        # Outgoing TE polarisation 
        if pol=='TE': trans_k = np.abs(stack.vec_coef_down[vec_index][0:n_PW_p_pols]).reshape(-1,)
        # Outgoing TM polarisation
        if pol=='TM': trans_k = np.abs(stack.vec_coef_down[vec_index][n_PW_p_pols-1:-1]).reshape(-1,) 
        trans_k_array = np.array(trans_k).reshape(-1,)

        select_trans = trans_k_array[axis_indices]
        store_alphas = np.append(store_alphas,alphas)
        store_k_trans = np.append(store_k_trans,select_trans)

    sort_indices = np.argsort(store_alphas)
    plot_alphas = store_alphas[sort_indices]
    plot_k_trans = store_k_trans[sort_indices]

    ax1.plot(plot_alphas,plot_k_trans, linewidth=linesstrength)

    min_k_label = np.max(np.abs(plot_alphas))
    k0 = abs(k0)
    min_k_l_k0  = np.max(np.abs(plot_alphas))/k0

    n_H = max_n(stacks_list)
    new_tick_values = [-min_k_label, -n_H*k0, -k0, 0, k0, n_H*k0, min_k_label]
    new_tick_labels = [r"$-%ik_0$"%min_k_l_k0,r'$-n_Hk_0$',r'$-k_0$',r'0',
        r'$k_0$',r'$n_Hk_0$',r"$%ik_0$"%min_k_l_k0]
    ax1.set_xticks(new_tick_values)
    ax1.set_xticklabels(new_tick_labels)
    ax1.set_xlim(-min_k_label,min_k_label)

    ax1.set_xlabel(r'$k_\parallel$')
    ax1.set_ylabel(r'$|E|$')
    plt.savefig('k_vector_excitation-lay_%s' % lay_interest, bbox_inches='tight')
#######################################################################################


####Plot amplitudes of PW orders#######################################################
def amps_of_orders(stacks_list, xvalues=None, chosen_PW_order=None,\
    lay_interest=0, add_height=None, add_title=None):
    """ Plot the amplitudes of plane wave orders in given layer. """
    fig = plt.figure(num=None, dpi=80, facecolor='w', edgecolor='k')
    ax1 = fig.add_subplot(1,1,1)
    # vec_coef sorted from top of stack, everything else is sorted from bottom
    vec_index = len(stacks_list[-1].layers) - lay_interest - 1
    if chosen_PW_order == None:
        # Create arrays of grating order indexes (-p, ..., p)
        max_ords = stacks_list[0].layers[-1].max_order_PWs
        chosen_PW_order = np.arange(-max_ords, max_ords + 1)

    if stacks_list[0].layers[0].light.wl_nm != stacks_list[-1].layers[0].light.wl_nm:
        xvalues = [s.layers[0].light.wl_nm for s in stacks_list]
        xlabel = r'$\lambda$ (nm)'
    elif set(stacks_list[0].layers[0].light.k_pll) != set(stacks_list[-1].layers[0].light.k_pll):
        xvalues = [np.sqrt(s.layers[0].light.k_pll[0]**2 + s.layers[0].light.k_pll[1]**2) for s in stacks_list]
        xlabel = r'$|k_\parallel|$'
    elif xvalues==None: 
        "amps_of_orders cannot guess what to plot on x-axis, specify with xvalues input."
        return

    for pxs in chosen_PW_order:
        store_trans = []
        for stack in stacks_list:
            assert isinstance(stack.layers[lay_interest],objects.Anallo), \
            "amps_of_orders only works in ThinFilm layers, change lay_interest input."
            k0 = stack.layers[0].k()
            n_PW_p_pols = stack.layers[0].structure.num_pw_per_pol
            # Calculate k_x that correspond to k_y = beta0 = 0 (in normalized units)
            alpha0, beta0 = stack.layers[0].k_pll_norm()
            alphas = alpha0 + pxs * 2 * np.pi
            on_axis_kzs = sqrt(k0**2 - alphas**2 - beta0**2)
            full_k_space = stack.layers[0].k_z
            # consider only transmission into singular polarization
            one_pol_k_space = full_k_space[0:n_PW_p_pols]

            ix = np.in1d(one_pol_k_space.ravel(), on_axis_kzs).reshape(one_pol_k_space.shape)
            axis_indices = np.ravel(np.array(np.where(ix))).astype(int)

            trans = np.abs(stack.vec_coef_down[vec_index][axis_indices]).reshape(-1,) # Outgoing TE polarisation
            trans += np.abs(stack.vec_coef_down[vec_index][n_PW_p_pols+axis_indices]).reshape(-1,) # Outgoing TM polarisation
            store_trans = np.append(store_trans,trans)

        ax1.plot(xvalues,store_trans, label="m = %s" %str(pxs))#, linewidth=linesstrength)

    handles, labels = ax1.get_legend_handles_labels()
    lgd = ax1.legend(handles, labels, loc='center left', bbox_to_anchor=(1.0,0.5))
    ax1.set_ylabel(r'$|E|_{trans}$')
    ax1.set_xlabel(xlabel)
    if add_title != None: plt.suptitle('%s' % add_title)
    if add_height!= None: add_title += zeros_int_str(add_height)
    if add_title != None: plt.savefig('PW_orders-lay_%s' % lay_interest + add_title, \
        bbox_extra_artists=(lgd,), bbox_inches='tight')
    else: plt.savefig('PW_orders-lay_%s' % lay_interest, bbox_extra_artists=(lgd,), \
        bbox_inches='tight')

def evanescent_merit(stacks_list, xvalues=None, chosen_PW_order=None,\
    lay_interest=0, add_height=None, add_title=None):
    """ Create a figure of merit for the 'evanescent-ness' of excited fields. """
    fig = plt.figure(num=None, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
    ax1 = fig.add_subplot(1,1,1)
    # vec_coef sorted from top of stack, everything else is sorted from bottom
    vec_index = len(stacks_list[-1].layers) - lay_interest - 1
    if chosen_PW_order == None:
        # Create arrays of grating order indexes (-p, ..., p)
        max_ords = stacks_list[0].layers[-1].max_order_PWs
        chosen_PW_order = np.arange(-max_ords, max_ords + 1)
    n_H = max_n(stacks_list)

    if stacks_list[0].layers[0].light.wl_nm != stacks_list[-1].layers[0].light.wl_nm:
        xvalues = [s.layers[0].light.wl_nm for s in stacks_list]
        xlabel = r'$\lambda$ (nm)'
    elif set(stacks_list[0].layers[0].light.k_pll) != set(stacks_list[-1].layers[0].light.k_pll):
        xvalues = [np.sqrt(s.layers[0].light.k_pll[0]**2 + s.layers[0].light.k_pll[1]**2) for s in stacks_list]
        xlabel = r'$|k_\parallel|$'
    elif xvalues==None: 
        "evanescent_merit cannot guess what to plot on x-axis, specify with xvalues input."
        return

    store_m_p  = []
    store_m_ne = []
    store_m_fe = []
    store_x_p  = []
    store_x_ne = []
    store_x_fe = []
    s = 0
    for stack in stacks_list:
        assert isinstance(stack.layers[lay_interest],objects.Anallo), \
        "evanescent_merit only works in ThinFilm layers, change lay_interest input."
        merit_prop    = 0.0
        merit_near_ev = 0.0
        merit_far_ev  = 0.0
        for pxs in chosen_PW_order:
            k0 = stack.layers[lay_interest].k()
            n_PW_p_pols = stack.layers[lay_interest].structure.num_pw_per_pol
            # Calculate k_x that correspond to k_y = beta0 = 0 (in normalized units)
            alpha0, beta0 = stack.layers[lay_interest].k_pll_norm()
            alphas = alpha0 + pxs * 2 * np.pi
            this_k_pll2 = alphas**2 + beta0**2
            on_axis_kzs = sqrt(k0**2 - alphas**2 - beta0**2)
            full_k_space = stack.layers[lay_interest].k_z
            # consider only transmission into singular polarization
            one_pol_k_space = full_k_space[0:n_PW_p_pols]

            ix = np.in1d(one_pol_k_space.ravel(), on_axis_kzs).reshape(one_pol_k_space.shape)
            axis_indices = np.ravel(np.array(np.where(ix))).astype(int)

            trans = np.abs(stack.vec_coef_down[vec_index][axis_indices]).reshape(-1,) # Outgoing TE polarisation
            trans += np.abs(stack.vec_coef_down[vec_index][n_PW_p_pols+axis_indices]).reshape(-1,) # Outgoing TM polarisation

            if this_k_pll2 < k0**2: merit_prop += trans
            if k0**2 < this_k_pll2 < (n_H*k0)**2: merit_near_ev += trans
            if this_k_pll2 > (n_H*k0)**2: merit_far_ev += trans

        if merit_prop != 0.0:
            store_m_p = np.append(store_m_p,merit_prop)
            store_x_p = np.append(store_x_p,xvalues[s])
        if merit_near_ev != 0.0:
            store_m_ne = np.append(store_m_ne,merit_near_ev)
            store_x_ne = np.append(store_x_ne,xvalues[s])
        if merit_far_ev != 0.0:
            store_m_fe = np.append(store_m_fe,merit_far_ev)
            store_x_fe = np.append(store_x_fe,xvalues[s])
        s+=1

    if len(store_m_p)  != 0: ax1.plot(store_x_p,store_m_p, label="prop")
    if len(store_m_ne) != 0: ax1.plot(store_x_ne,store_m_ne, label="near ev")
    if len(store_m_fe) != 0: ax1.plot(store_x_fe,store_m_fe, label="far ev")

    handles, labels = ax1.get_legend_handles_labels()
    lgd = ax1.legend(handles, labels, ncol=3, loc='upper center', bbox_to_anchor=(0.5,1.2))
    ax1.set_ylabel('Ev FoM')
    ax1.set_xlabel(xlabel)
    if add_title != None: plt.suptitle('%s' % add_title)
    if add_height!= None: add_title += zeros_int_str(add_height)
    if add_title != None: plt.savefig('evanescent_merit-lay_%s' % lay_interest + add_title, \
        bbox_extra_artists=(lgd,), bbox_inches='tight')
    else: plt.savefig('evanescent_merit-lay_%s' % lay_interest, bbox_extra_artists=(lgd,),\
     bbox_inches='tight')


    fig = plt.figure(num=None, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
    ax1 = fig.add_subplot(1,1,1)

    store_m_p  = []
    store_m_ne = []
    store_m_fe = []
    for stack in stacks_list:
        merit_prop    = 0.0
        merit_near_ev = 0.0
        for pxs in chosen_PW_order:
            k0 = stack.layers[lay_interest].k()
            n_PW_p_pols = stack.layers[lay_interest].structure.num_pw_per_pol
            # Calculate k_x that correspond to k_y = beta0 = 0 (in normalized units)
            alpha0, beta0 = stack.layers[lay_interest].k_pll_norm()
            alphas = alpha0 + pxs * 2 * np.pi
            this_k_pll2 = alphas**2 + beta0**2
            on_axis_kzs = sqrt(k0**2 - alphas**2 - beta0**2)
            full_k_space = stack.layers[lay_interest].k_z
            # consider only transmission into singular polarization
            one_pol_k_space = full_k_space[0:n_PW_p_pols]

            ix = np.in1d(one_pol_k_space.ravel(), on_axis_kzs).reshape(one_pol_k_space.shape)
            axis_indices = np.ravel(np.array(np.where(ix))).astype(int)

            trans = np.abs(stack.vec_coef_down[vec_index][axis_indices]).reshape(-1,) # Outgoing TE polarisation
            trans += np.abs(stack.vec_coef_down[vec_index][n_PW_p_pols+axis_indices]).reshape(-1,) # Outgoing TM polarisation

            merit_prop += np.abs(pxs) * trans
            merit_near_ev += trans

        store_m_p  = np.append(store_m_p,merit_prop)
        store_m_ne = np.append(store_m_ne,merit_prop/merit_near_ev)

    s = 0
    for stack in stacks_list:
        merit_far_ev  = 0.0
        for pxs in chosen_PW_order:
            k0 = stack.layers[lay_interest].k()
            n_PW_p_pols = stack.layers[lay_interest].structure.num_pw_per_pol
            # Calculate k_x that correspond to k_y = beta0 = 0 (in normalized units)
            alpha0, beta0 = stack.layers[lay_interest].k_pll_norm()
            alphas = alpha0 + pxs * 2 * np.pi
            this_k_pll2 = alphas**2 + beta0**2
            on_axis_kzs = sqrt(k0**2 - alphas**2 - beta0**2)
            full_k_space = stack.layers[lay_interest].k_z
            # consider only transmission into singular polarization
            one_pol_k_space = full_k_space[0:n_PW_p_pols]

            ix = np.in1d(one_pol_k_space.ravel(), on_axis_kzs).reshape(one_pol_k_space.shape)
            axis_indices = np.ravel(np.array(np.where(ix))).astype(int)

            trans = np.abs(stack.vec_coef_down[vec_index][axis_indices]).reshape(-1,) # Outgoing TE polarisation
            trans += np.abs(stack.vec_coef_down[vec_index][n_PW_p_pols+axis_indices]).reshape(-1,) # Outgoing TM polarisation

            merit_far_ev += np.abs(np.abs(pxs) - np.abs(store_m_ne[s]))**2 * trans

        store_m_fe = np.append(store_m_fe, merit_far_ev)
        s +=1

    ax1.plot(xvalues,store_m_p, label=r'$\Sigma |p| |a_p|$')
    ax1.plot(xvalues,store_m_ne, label=r'$\Sigma|p| |a_p| / |a_p|$')
    ax1.plot(xvalues,store_m_fe, label=r'$\Sigma||p| - |\overline{p}||^2 |a_p|$')

    handles, labels = ax1.get_legend_handles_labels()
    lgd = ax1.legend(handles, labels, ncol=3, loc='upper center', bbox_to_anchor=(0.5,1.2))
    ax1.set_ylabel('Ev FoM')
    if np.max(xvalues) < 90: ax1.set_xlabel(r'$\theta$') # hack way to tell what plotting against
    else: ax1.set_xlabel(r'$\lambda$ (nm)')
    if add_title != None: plt.suptitle('%s' % add_title)
    if add_height!= None: add_title += zeros_int_str(add_height)
    if add_title != None: plt.savefig('evanescent_merit-2-lay_%s' % lay_interest + add_title, \
        bbox_extra_artists=(lgd,), bbox_inches='tight')
    else: plt.savefig('evanescent_merit-2-lay_%s' % lay_interest, bbox_extra_artists=(lgd,),\
     bbox_inches='tight')
#######################################################################################


####Field plotting routines############################################################
def fields_2d(pstack, wl, Struc_lay = 1, TF_lay=0):
    """
    """
    from fortran import EMUstack

    dir_name = "2D_Fields"
    if os.path.exists(dir_name):
        subprocess.call('rm %s -r' %dir_name, shell = True)
    os.mkdir(dir_name)

    #
    # plot fields inside nanostructure (Bloch Mode Basis)
    #
    meat = pstack.layers[Struc_lay]
    gmsh_file_pos = meat.structure.mesh_file[0:-5]
    eps_eff = meat.n_effs**2
    h_normed = float(meat.structure.height_nm)/float(meat.structure.period)
    wl_normed = pstack.layers[Struc_lay].wl_norm()

    nnodes=6
    q_average = 0 # at a discontinuity, use average value if q_average = 1
    if meat.E_H_field == 1:
        EH_name = "E_"
    else:
        EH_name = "H_"
    extra_name = EH_name + "Lay" + zeros_int_str(Struc_lay)

    if Struc_lay < TF_lay: select_h = 0.0   # top interface
    else: select_h = h_normed               # bottom interface

    # vec_coef sorted from top of stack, everything else is sorted from bottom
    vec_index = len(pstack.layers) - Struc_lay - 1
    vec_coef = np.concatenate((pstack.vec_coef_down[vec_index],pstack.vec_coef_up[vec_index]))

    EMUstack.gmsh_plot_field (meat.num_BM, 
        meat.n_msh_el, meat.n_msh_pts, nnodes, meat.nb_typ_el, meat.table_nod, meat.type_el,
        eps_eff, meat.x_arr, meat.k_z, meat.sol1, vec_coef, h_normed, select_h, 
        gmsh_file_pos, q_average, meat.structure.plot_real, 
        meat.structure.plot_imag, meat.structure.plot_abs, extra_name)


    #
    # plot fields in Thin Film (Plane Wave Basis)
    #
    extra_name = EH_name + "Lay" + zeros_int_str(TF_lay)
    select_h = 0.0
    lat_vec = [[1.0, 0.0], [0.0, 1.0]]
    bun = pstack.layers[TF_lay] # superstrate # substrate
    n_eff_0 = bun.n()
    neq_PW = bun.structure.num_pw_per_pol
    bloch_vec = bun.k_pll_norm()
    ordre_ls = bun.max_order_PWs
    index_pw = bun.sort_order
    index_pw_inv = np.zeros(shape=(np.shape(index_pw)),dtype='int')
    for s in range(len(index_pw)):
        s2 = index_pw[s]
        index_pw_inv[s2] = s
    index_pw_inv += 1 # convert to fortran indices (starting from 1)

    # vec_coef sorted from top of stack, everything else is sorted from bottom
    vec_index = len(pstack.layers) - TF_lay - 1
    vec_coef_down = pstack.vec_coef_down[vec_index]
    if TF_lay == 0: vec_coef_up = np.zeros(shape=(np.shape(vec_coef_down)),dtype='complex128')
    else: vec_coef_up = pstack.vec_coef_up[vec_index]

    EMUstack.gmsh_plot_pw(meat.n_msh_el, meat.n_msh_pts, nnodes, neq_PW,
        bloch_vec, meat.table_nod, meat.x_arr, lat_vec, wl_normed,
        n_eff_0, vec_coef_down, vec_coef_up, 
        index_pw_inv, ordre_ls, select_h,
        gmsh_file_pos, q_average, meat.structure.plot_real, meat.structure.plot_imag,
        meat.structure.plot_abs, extra_name)

    # # Semi-inf case
    # vec_coef_up = meat.R12[:,0] # TE polarisation
    # vec_coef_up = meat.R12[:,neq_PW] # TM polarisation

    # vec_coef_down = np.zeros(shape=(np.shape(vec_coef_up)),dtype='complex128')
    # vec_coef_down[neq_PW] = 1.0

def E_substrate2(stacks_list, wavelengths, period, r_a, pw, bm):
    alpha_unsrt = np.array(stacks_list[0].layers[0].alphas)
    beta_unsrt = np.array(stacks_list[0].layers[0].betas)
    gamma_unsrt = np.array(stacks_list[0].layers[0].k_z_unsrt)
    k_array = np.column_stack((alpha_unsrt,beta_unsrt,gamma_unsrt))
    k_sort = k_array [np.argsort(-1*np.real(k_array [:,2])+np.imag(k_array [:,2]))]
    alpha = np.array(k_sort[:,0])
    beta = np.array(k_sort[:,1])
    gamma = np.array(k_sort[:,2])
    n = stacks_list[0].layers[0].n()
    PWordtot = stacks_list[0].layers[0].structure.num_pw_per_pol
    Tnet = stacks_list[0].T_net*stacks_list[0].layers[-1].specular_incidence(pol='TE')
    Tnet_TE = np.array(Tnet[0:PWordtot]).flatten()
    Tnet_TM = np.array(Tnet[PWordtot::]).flatten()
    nu_calc_pts = 50
    xrange = np.linspace(0,1.0,nu_calc_pts)
    yrange = np.linspace(0,1.0,nu_calc_pts)
    zrange = np.array([0])
    norm = np.sqrt(alpha**2+beta**2)
    k = np.sqrt(alpha**2+beta**2+gamma**2)
    chi_TE = np.sqrt((n*gamma)/k)
    chi_TM = np.sqrt((n*k)/gamma)
    E_TE_x = beta/norm
    E_TE_y = -1*alpha/norm
    E_TE_z = np.array(np.zeros(np.size(E_TE_x)))
    E_TM_x = alpha/norm
    E_TM_y = beta/norm
    E_TM_z = -1*norm/gamma
    eta_TE_x = (Tnet_TE*E_TE_x)/chi_TE
    eta_TE_y = (Tnet_TE*E_TE_y)/chi_TE
    eta_TE_z = (Tnet_TE*E_TE_z)/chi_TE
    eta_TM_x = (Tnet_TM*E_TM_x)/chi_TM
    eta_TM_y = (Tnet_TM*E_TM_y)/chi_TM
    eta_TM_z = (Tnet_TM*E_TM_z)/chi_TM
    prop = 2*np.count_nonzero(np.real(gamma))
    evan = 2*PWordtot - prop
    s = range(2*PWordtot)

    E_TE_x_array = np.zeros((nu_calc_pts**2,np.size(zrange)), dtype = 'complex')
    E_TE_y_array = np.zeros((nu_calc_pts**2,np.size(zrange)), dtype = 'complex')
    E_TE_z_array = np.zeros((nu_calc_pts**2,np.size(zrange)), dtype = 'complex')
    E_TM_x_array = np.zeros((nu_calc_pts**2,np.size(zrange)), dtype = 'complex')
    E_TM_y_array = np.zeros((nu_calc_pts**2,np.size(zrange)), dtype = 'complex')
    E_TM_z_array = np.zeros((nu_calc_pts**2,np.size(zrange)), dtype = 'complex')
    
    xpos = []
    ypos = []

    for z in range(len(zrange)):
        xy = 0
        for y in range(nu_calc_pts):
            for x in range(nu_calc_pts):
                expo = np.exp(1j*(alpha*xrange[x]+beta*yrange[y]-gamma*zrange[z]))
                E_TE_x = np.sum(eta_TE_x*expo)
                E_TE_y = np.sum(eta_TE_y*expo)
                E_TE_z = np.sum(eta_TE_z*expo)
                E_TM_x = np.sum(eta_TM_x*expo)
                E_TM_y = np.sum(eta_TM_y*expo)
                E_TM_z = np.sum(eta_TM_z*expo)
                E_TE_x_array[xy,z] = E_TE_x 
                E_TE_y_array[xy,z] = E_TE_y
                E_TE_z_array[xy,z] = E_TE_z
                E_TM_x_array[xy,z] = E_TM_x 
                E_TM_y_array[xy,z] = E_TM_y
                E_TM_z_array[xy,z] = E_TM_z
                xpos.append(xrange[x])
                ypos.append(yrange[y])
                xy += 1
    E_x_array = E_TE_x_array + E_TM_x_array
    E_y_array = E_TE_y_array + E_TM_y_array
    E_z_array = E_TE_z_array + E_TM_z_array
    Etot = np.sqrt(E_x_array*np.conj(E_x_array) + E_y_array*np.conj(E_y_array) + E_z_array*np.conj(E_z_array))

    # np.savetxt('sub2_eta_re.txt',np.real(eta),fmt = '%18.11f')
    # np.savetxt('sub2_eta_imag.txt',np.imag(eta),fmt = '%18.11f')
    # np.savetxt('sub2_k.txt',np.column_stack((np.real(k_sort),np.imag(gamma))),fmt = '%18.11f')
    # np.savetxt('sub2_Tnet_re.txt',np.real(np.column_stack((Tnet_TE,Tnet_TM))),fmt = '%18.11f')
    # np.savetxt('sub2_Tnet_imag.txt',np.imag(np.column_stack((Tnet_TE,Tnet_TM))),fmt = '%18.11f')
    # np.savetxt('sub2_Ecoeff_re.txt', np.real(E),fmt = '%18.11f')
    # np.savetxt('sub2_Ecoeff_imag.txt', np.imag(E),fmt = '%18.11f')
    # np.savetxt('sub2_chi_re.txt',np.real(chi),fmt = '%18.11f')
    # np.savetxt('sub2_chi_imag.txt',np.imag(chi),fmt = '%18.11f')

    
    fig1 = plt.figure(num=None, figsize=(12,21), dpi=80, facecolor='w', edgecolor='k')
    ax1 = fig1.add_subplot(2,1,1)
    # xpos = np.tile(xrange,np.size(yrange))
    # ypos = np.repeat(yrange,np.size(xrange))
    Ex = np.real(E_x_array[:,0])
    x_min = xrange[0]
    x_max = xrange[-1]
    y_min = yrange[0]
    y_max = yrange[-1]
    int_x,int_y = np.mgrid[x_min:x_max:nu_calc_pts*1j, y_min:y_max:nu_calc_pts*1j]
    int_Ex = griddata(xpos, ypos, Ex, int_x, int_y)
    cmap = plt.get_cmap('jet')
    CS = plt.contourf(int_x,int_y,int_Ex,15,cmap=cmap)
    circle = plt.Circle((0.5,0.5), radius=r_a/period, fc='none', ec = 'w',fill = False)
    ax1.add_artist(circle)
    plt.axis([x_min,x_max,y_min,y_max])
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(r'Real E_x')
    ax1.set_xlabel('x (normalised period)')
    ax1.set_ylabel('y (normalised period)')
    ax1.axis('scaled')

    ax1 = fig1.add_subplot(2,1,2)
    # xpos = np.tile(xrange,np.size(yrange))
    # ypos = np.repeat(yrange,np.size(xrange))
    Ex = np.real(E_y_array[:,0])
    x_min = xrange[0]
    x_max = xrange[-1]
    y_min = yrange[0]
    y_max = yrange[-1]
    int_Ex = griddata(xpos, ypos, Ex, int_x, int_y)
    cmap = plt.get_cmap('jet')
    CS = plt.contourf(int_x,int_y,int_Ex,15,cmap=cmap)
    circle = plt.Circle((0.5,0.5), radius=r_a/period, fc='none', ec = 'w',fill = False)
    ax1.add_artist(circle)
    plt.axis([x_min,x_max,y_min,y_max])
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(r'Real E_y')
    ax1.set_xlabel('x (normalised period)')
    ax1.set_ylabel('y (normalised period)')
    ax1.axis('scaled')


    plt.suptitle('E = E_TE+E_TM \n wavelength = %(wl)s ,period = %(d)s, PW = %(pw)s, # BM = %(bm)s,' % {'wl' : wavelengths[0], 'd' : period, 'pw' : pw, 'bm' : bm} + '\n' + 
                 '# prop. ords = %(prop)s, # evan. ords = %(evan)s , refractive index = %(n)s' % {'evan' : evan, 'prop' : prop, 'n' : n})
    plt.savefig('sub2_E_interface%(wl)s_contour.pdf'% {'wl' : wavelengths[0]})

def fields_3d(pstack, wl):
    """
    """
    from fortran import EMUstack
    import subprocess

    nnodes=6
    dir_name = "3D_Fields"
    if os.path.exists(dir_name):
        subprocess.call('rm %s -r' %dir_name, shell = True)
    os.mkdir(dir_name)
    dir_name = "3D_Fields/Anim"
    os.mkdir(dir_name)

    for lay in range(len(pstack.layers)-2): #remove -2 once semi inf TFs can be plotted
        lay+=1
        layer_name = 'Lay' + zeros_int_str(lay)

        meat = pstack.layers[lay]
        gmsh_file_pos = meat.structure.mesh_file[0:-5]

        vec_coef = np.concatenate((pstack.vec_coef_down[1],pstack.vec_coef_up[1]))
        # vec_coef_up = np.zeros(shape=(np.shape(pstack.vec_coef_down[1])),dtype='complex128')
        # vec_coef = np.concatenate((pstack.vec_coef_down[1],vec_coef_up))
        h_normed = float(meat.structure.height_nm)/float(meat.structure.period)
        wl_normed = pstack.layers[lay].wl_norm()

        EMUstack.gmsh_plot_field_3d(wl_normed, h_normed, meat.num_BM,   
            meat.E_H_field, meat.n_msh_el, meat.n_msh_pts, 
            nnodes, meat.type_el, meat.nb_typ_el, meat.n_effs, meat.table_nod,
            meat.k_z, meat.sol1, vec_coef, meat.x_arr, gmsh_file_pos, layer_name)
#######################################################################################


####Fabry-Perot resonances#############################################################
def Fabry_Perot_res(stacks_list, freq_list, kx_list, f_0, k_0, lay_interests=1):
    """ Calculate the Fabry-Perot resonance condition for a resonances within a layer.
    'lay_interest' : specifies which layer in the stack to find F-P resonances of.
    """
    n_freqs = len(freq_list)
    n_kxs   = len(kx_list)
    height = stacks_list[-1].heights_nm()[lay_interest-1] # assumes all stacks have equal height
    period = stacks_list[-1].period
    num_BMs = stacks_list[-1].layers[lay_interest].num_BM
    I_mat = np.matrix(np.eye(num_BMs),dtype='D')

    plot_mat = np.zeros(shape=(n_freqs,n_kxs),dtype='complex128')
    for i in range(n_kxs):
        kx_slice = stacks_list[i*n_freqs:(i+1)*n_freqs]
        j = 0
        for stack in kx_slice:
            P   = stack.layers[lay_interest].prop_fwd(height/period)
            R21 = stack.layers[lay_interest].R21
            FP_term = np.linalg.det(I_mat - R21*P*R21*P)
            plot_mat[j,i] = FP_term
            j += 1
    image = np.log(np.abs(plot_mat))

    fig = plt.figure()
    ax1 = fig.add_subplot(1,1,1)
    cax = ax1.imshow(image,cmap=plt.cm.gray_r, interpolation='none')

    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    shape = np.shape(plot_mat)
    majorLocator   = MultipleLocator(shape[1]-1)
    ax1.xaxis.set_major_locator(majorLocator)
    majorLocator   = MultipleLocator(shape[0]-1)
    ax1.yaxis.set_major_locator(majorLocator)
    xlims = [kx_list[0]/k_0, kx_list[-1]/k_0]
    ax1.set_xticklabels(xlims)
    ylims = [freq_list[0]/f_0, freq_list[-1]/f_0]
    ax1.set_yticklabels(ylims)

    cbar = fig.colorbar(cax)
    cbar.set_label(r'$ln(|I - R_{21}PR_{21}P|)$',size=18)
    ax1.set_xlabel(r'$k_\perp/k_0$')
    ax1.set_ylabel(r'$f/f_0$')
    ax1.axis('image')
    plt.savefig('Fabry-Perot_resonances', bbox_inches='tight')
#######################################################################################