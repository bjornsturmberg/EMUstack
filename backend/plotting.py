"""
    plotting.py is a subroutine of EMUstack that contains numerous plotting
    routines.

    Copyright (C) 2015  Bjorn Sturmberg, Kokou Dossou, Felix Lawrence

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
from __future__ import print_function

import objects
import mode_calcs

import numpy as np
from scipy import sqrt
import subprocess
from matplotlib.mlab import griddata
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os

# font = {'family' : 'normal',
#         'weight' : 'bold',
#         'size'   : 18}

# font = {'size'   : 14}
# matplotlib.rc('font', **font)
linesstrength = 2.5
title_font = 10


#### Natural constants ########################################################
ASTM15_tot_I   = 900.084            # Integral ASTM 1.5 solar irradiance W/m**2
Plancks_h      = 6.62606957*1e-34   # Planck's constant
speed_c        = 299792458          # Speed of light in vacuum
charge_e       = 1.602176565*1e-19  # Charge of an electron
###############################################################################


#### Short utility functions ##################################################
def clear_previous():
    """ Delete all files of specified type as well as field directories. """

    devnull = open(os.devnull, 'wb')

    type_list = ['*.npz', '*.pdf', '*.txt', '*.gif', '*.png', '*.log', '*.svg',
    'fields_vertically -r', 'in_plane_fields -r', 'Bloch_fields -r',
    'field_values -r', '3d_fields -r']
    for typ in type_list:
        try:
            files_rm = 'rm %s'% typ
            subprocess.call(files_rm, shell = True, stderr=devnull)
        except:
            pass


def zeros_int_str(zero_int):
    """ Convert integer into string with '0' in place of ' '. """
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
    """ Find maximum refractive index n in stacks_list. """
    ns = []
    for s in stacks_list:
        for l in s.layers:
            if isinstance(l, objects.Anallo):
                ns.append(l.n())
            if isinstance(l, objects.Simmo):
                wl = l.light.wl_nm
                ns.append(l.structure.background.n(wl))
                ns.append(l.structure.inclusion_a.n(wl))
                ns.append(l.structure.inclusion_b.n(wl))
    return np.real(np.max(ns))


def gen_params_string(stack, layer=1):
    """ Generate the string of simulation info that is to be printed \
        at the top of plots.
    """
    param_layer = stack[0].layers[layer].structure
    # Plot t,r,a for each layer & total, then save each to text files
    if isinstance(param_layer, objects.NanoStruct):
        params_2_print = 'ff = %5.3f, '% param_layer.ff
        params_2_print += 'd = %(period)d, a1 = %(diameter)d, '% {
        'period'        : param_layer.period, 'diameter' : param_layer.diameter1,}
        if param_layer.diameter2 != 0: params_2_print += 'a2 = %(rad)d '% {'rad' : param_layer.diameter2,}
        if param_layer.periodicity == '2D_array':
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
            if param_layer.inc_shape == 'square': params_2_print += '\nSquare ins '
            if param_layer.inc_shape == 'ellipse': params_2_print += '\nEllipticity = %(rad)5.3f '% {'rad' : param_layer.ellipticity}
        elif param_layer.periodicity == '1D_array':
            params_2_print += ''
        params_2_print += '%(BMs)dBMs, PW_radius = %(PWs)d, '% \
        {'BMs' : stack[0].layers[layer].num_BMs, \
        'PWs' : stack[0].layers[layer].max_order_PWs}
    else:
        params_2_print = 'PW_radius = %(PWs)d, ' \
        % {'PWs' : stack[0].layers[layer].max_order_PWs, }

    # light = stack[0].layers[layer].light
    # k_pll = light.k_pll * param_layer.period
    # params_2_print += r'$k_\parallel d$ = '
    # tmp3 = '(%(kx)1.4f, %(ky)1.4f), '% {'kx' : k_pll[0], 'ky' : k_pll[1]}
    # params_2_print += r'$\theta$ = %(theta)6.2f, $\phi$ = %(phi)6.2f, '% {
    # 'theta' : light.theta,'phi' : light.phi, }

    return params_2_print
###############################################################################


#### Standard plotting of spectra #############################################
def t_r_a_plots(stacks_list, xvalues=None, xlabel='', params_layer=1,
    active_layer_nu=1, stack_label=1, ult_eta=False, J_sc=False,
    weight_spec=False, extinct=False, add_height=0, add_name='',
    save_pdf=True, save_txt=False, set_y_lim=True, label_eV=False):
    """ Plot t, r, a for each layer & in total.

        Args:
            stacks_list  (list): Stack objects containing data to plot.

        Keyword Args:
            xvalues  (list): The values stacks_list is to be plotted as a \
                function of.

            params_layer  (int): The index in stacks_list of the layer for \
                which the geometric parameters are put in the title of \
                the plots.

            active_layer_nu  (int): The index in stacks_list (from bottom) \
                of the layer for which the ult_eta and/or J_sc are calculated.

            stack_label  (int): Label to differentiate plots of different \
                :Stack:s.

            ult_eta  (bool): If True; calculate the 'ultimate efficiency'.

            J_sc  (bool): If True; calculate the idealised short circuit \
                current.

            weight_spec  (bool): If True; plot t, r, a spectra weighted by \
                the ASTM 1.5 solar spectrum.

            extinct  (bool): If True; calculate the extinction ratio in \
                transmission.

            add_height  (float): Print the heights of :Stack: layer in title.

            add_name  (str): Add add_name to title.

            save_pdf  (bool): If True; save spectra as pdf files. \
                True by default.

            save_txt  (bool): If True; save spectra data to text \
                files.

            set_y_lim  (bool): If True; set y limits to (0,1).

            label_eV  (bool): If True; add energy label in eV.
    """

    height_list = stacks_list[0].heights_nm()[::-1]
    params_2_print = gen_params_string(stacks_list, params_layer)
    params_2_print += '\n'r'$h_t,...,h_b$ = '
    params_2_print += ''.join('%4f, ' % num for num in height_list)

    if xvalues == None:
        if stacks_list[0].layers[0].light.wl_nm != stacks_list[-1].layers[0].light.wl_nm:
            xvalues = [s.layers[0].light.wl_nm for s in stacks_list]
            xlabel = r'$\lambda$ (nm)'
        elif set(stacks_list[0].layers[0].light.k_pll) != set(stacks_list[-1].layers[0].light.k_pll):
            xvalues = [np.sqrt(s.layers[0].light.k_pll[0]**2 + s.layers[0].light.k_pll[1]**2) for s in stacks_list]
            xlabel = r'$|k_\parallel|$'
        else:
            xvalues = [s.layers[0].light.wl_nm for s in stacks_list]
            xlabel = r'$\lambda$ (nm)'
            print("t_r_a_plots is guessing you have a single wavelength, else specify xvalues.")

    if add_height!=0: add_name += "_" + zeros_int_str(add_height)
    stack_label = zeros_int_str(stack_label)

    a_list = []
    t_list = []
    r_list = []
    for stack in stacks_list:
        a_list.extend(stack.a_list)
        t_list.extend(stack.t_list)
        r_list.extend(stack.r_list)

    layers_steps = len(stacks_list[0].layers) - 1
    a_tot = []
    t_tot = []
    r_tot = []
    for i in range(len(xvalues)):
        a_tot.append(float(a_list[layers_steps-1+(i*layers_steps)]))
        t_tot.append(float(t_list[layers_steps-1+(i*layers_steps)]))
        r_tot.append(float(r_list[i]))

    if ult_eta is True or J_sc is True:
        active_abs = []
        active_layer_nu = len(stacks_list[0].layers) - active_layer_nu - 1
        if not 0 < active_layer_nu < len(stacks_list[0].layers)-1:
            raise ValueError("active_layer_nu must refer to a finite layer.")
        for i in range(len(xvalues)):
            active_abs.append(float(a_list[active_layer_nu - 1 +
                              i*layers_steps]))
        out = []

    if ult_eta is True:
        Efficiency = ult_efficiency(active_abs, xvalues,
                                    params_2_print=params_2_print,
                                    stack_label=stack_label, add_name=add_name)
        params_2_print += r'$\eta$ = %(Efficiency)6.2f' % {'Efficiency':
                                                           Efficiency*100, }
        params_2_print += ' %'
        out.append(Efficiency)

    if J_sc is True:
        J = J_short_circuit(active_abs, xvalues,
                            params_2_print=params_2_print,
                            stack_label=stack_label, add_name=add_name)
        params_2_print += r'$J_{sc}$ = %(J)6.2f'% {'J' : J, }
        params_2_print += r' mA/cm$^2$'
        out.append(J)

    if save_txt is True or save_pdf is True:
        total_h = sum(stacks_list[0].heights_nm()) # look at first wl result to find h.
        # Plot t,r,a for each layer.
        layers_plot('Lay_Absorb', a_list, xvalues, xlabel, total_h,
                    params_2_print, stack_label, add_name, save_pdf, save_txt,
                    label_eV, set_y_lim)
        layers_plot('Lay_Trans',  t_list, xvalues, xlabel, total_h,
                    params_2_print, stack_label, add_name, save_pdf, save_txt,
                    label_eV, set_y_lim)
        layers_plot('Lay_Reflec', r_list, xvalues, xlabel, total_h,
                    params_2_print, stack_label, add_name, save_pdf, save_txt,
                    label_eV, set_y_lim)

    # Plot total t,r,a on a single plot.
    if save_pdf == True:
        plot_name = 'Total_Spectra'
        total_tra_plot(plot_name, a_tot, t_tot, r_tot, xvalues, xlabel,
                       params_2_print, stack_label, add_name, label_eV,
                       set_y_lim)

    if weight_spec == True:
        # Plot totals weighted by solar irradiance.
        Irrad_spec_file = '../backend/data/ASTMG173'
        i_data = np.loadtxt('%s.txt' % Irrad_spec_file)
        i_spec = np.interp(xvalues, i_data[:,0], i_data[:,3])
        bandgap_wl = xvalues[-1]
        weighting = i_spec/i_spec.max()*(xvalues/bandgap_wl)
        a_weighted = a_tot*weighting
        t_weighted = t_tot*weighting
        r_weighted = r_tot*weighting
        plot_name = 'Weighted_Total_Spectra'
        total_tra_plot(plot_name, a_weighted, t_weighted, r_weighted, xvalues,
                       xlabel, params_2_print, stack_label, add_name, label_eV,
                       set_y_lim)

    if extinct == True:
        extinction_plot(t_tot, xvalues, params_2_print, stack_label, add_name)

    if J_sc == True:
        return out
    elif ult_eta == True:
        return out
    else:
        return


def layers_plot(spectra_name, spec_list, xvalues, xlabel, total_h,
                params_2_print, stack_label, add_name, save_pdf, save_txt,
                label_eV, set_y_lim):
    """ Plots one type of spectrum across all layers.

        Is called from t_r_a_plots.
    """

    fig = plt.figure(num=None, figsize=(9, 12))
    nu_layers = int(len(spec_list)/len(xvalues))
    h_array = np.ones(len(xvalues))*total_h
    for i in range(nu_layers):
        layer_spec = []
        for wl in range(len(xvalues)):
            layer_spec = np.append(layer_spec, spec_list[wl*nu_layers + i])
        av_array = zip(xvalues, layer_spec, h_array)
        ax1 = fig.add_subplot(nu_layers, 1, i+1)
        ax1.plot(xvalues, layer_spec, linewidth=linesstrength)
        if label_eV == True:
            ax2 = ax1.twiny()
            new_tick_values = np.linspace(10, 0.5, 20)
            new_tick_locations = tick_function(new_tick_values)
            new_tick_labels = ["%.1f" % z for z in new_tick_values]
            ax2.set_xticks(new_tick_locations)
            ax2.set_xticklabels(new_tick_labels)
            ax2.set_xlim((xvalues[0], xvalues[-1]))
        if i == 0:
            if xlabel == r'$\lambda$ (nm)' and label_eV == True:
                ax2.set_xlabel('Energy (eV)')
        elif i == nu_layers-1:
            ax1.set_xlabel(xlabel)
            ax1.set_ylabel('Total')
        if i != 0 and label_eV == True:
            ax2.set_xticklabels( () )
        if i != nu_layers-1:
            ax1.set_xticklabels( () )
        if spectra_name == 'Lay_Absorb':
            if i == 0: ax1.set_ylabel('Top Layer')
            if i == nu_layers - 2: ax1.set_ylabel('Bottom Layer')
            suptitle_w_params = 'Absorptance in each layer' + add_name + \
            '\n' + params_2_print
            plt.suptitle(suptitle_w_params, fontsize = title_font)
            lay_spec_name = 'Lay_Absorb'
            if i == nu_layers-1:
                ax1.set_ylabel('Total')
                lay_spec_name = 'Absorptance'
                if set_y_lim==True:
                    plt.ylim((0, 1))
        elif spectra_name == 'Lay_Trans':
            if i == 0: ax1.set_ylabel('Top Layer')
            if i == nu_layers-2: ax1.set_ylabel('Bottom Layer')
            suptitle_w_params = 'Transmittance in each layer'+add_name+'\n'+params_2_print
            plt.suptitle(suptitle_w_params,fontsize=title_font)
            lay_spec_name = 'Lay_Trans'
            if i == nu_layers-1:
                ax1.set_ylabel('Total')
                lay_spec_name = 'Transmittance'
                if set_y_lim==True:
                    plt.ylim((0, 1))
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
                if set_y_lim==True:
                    plt.ylim((0, 1))
        ax1.set_xlim((xvalues[0], xvalues[-1]))

        if save_txt == True:
            if i != nu_layers-1:
                np.savetxt('%(s)s_%(i)i_stack%(bon)s%(add)s.txt'% \
                    {'s' : lay_spec_name, 'i' : i, 'bon' : stack_label, \
                    'add' : add_name}, av_array, fmt = '%18.11f')
            else:
                np.savetxt('%(s)s_stack%(bon)s%(add)s.txt'% \
                    {'s' : lay_spec_name, 'bon' : stack_label, \
                    'add' : add_name}, av_array, fmt = '%18.11f')
        if save_pdf == True:
            plt.savefig('%(s)s_stack%(bon)s%(add)s'% \
                {'s': spectra_name, 'bon': stack_label, 'add' : add_name})
    fig.clf()


def total_tra_plot(plot_name, a_spec, t_spec, r_spec, xvalues, xlabel,
    params_2_print, stack_label, add_name, label_eV, set_y_lim):
    """ Plots total t, r, a spectra on one plot.

        Is called from t_r_a_plots, t_r_a_plots_subs
    """

    fig = plt.figure(num=None, figsize=(9, 12))

    ax1 = fig.add_subplot(3, 1, 1)
    ax1.plot(xvalues, a_spec, linewidth=linesstrength)
    ax1.set_ylabel('Absorptance')
    ax1.set_xlim((xvalues[0], xvalues[-1]))
    if label_eV == True:
        ax2 = ax1.twiny()
        new_tick_values = np.linspace(10, 0.5, 20)
        new_tick_locations = tick_function(new_tick_values)
        new_tick_labels = ["%.1f" % z for z in new_tick_values]
        ax2.set_xticks(new_tick_locations)
        ax2.set_xticklabels(new_tick_labels)
        ax2.set_xlim((xvalues[0], xvalues[-1]))
        if xlabel == r'$\lambda$ (nm)':
            ax2.set_xlabel('Energy (eV)')
    if set_y_lim==True:
        plt.ylim((0, 1))

    ax1 = fig.add_subplot(3, 1, 2)
    ax1.plot(xvalues, t_spec, linewidth=linesstrength)
    ax1.set_ylabel('Transmittance')
    ax1.set_xlim((xvalues[0], xvalues[-1]))
    if label_eV == True:
        ax2 = ax1.twiny()
        ax2.set_xticklabels( () )
        ax2.set_xticks(new_tick_locations)
        ax2.set_xlim((xvalues[0], xvalues[-1]))
    if set_y_lim==True:
        plt.ylim((0, 1))

    ax1 = fig.add_subplot(3, 1, 3)
    ax1.plot(xvalues, r_spec, linewidth=linesstrength)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel('Reflectance')
    ax1.set_xlim((xvalues[0], xvalues[-1]))
    if label_eV == True:
        ax2 = ax1.twiny()
        ax2.set_xticklabels( () )
        ax2.set_xticks(new_tick_locations)
        ax2.set_xlim((xvalues[0], xvalues[-1]))
    if set_y_lim==True:
        plt.ylim((0, 1))

    plt.suptitle(params_2_print, fontsize=title_font)
    # plt.suptitle(plot_name+add_name+'\n'+params_2_print)
    plt.savefig('%(s)s_stack%(bon)s%(add)s'% \
                {'s': plot_name, 'bon': stack_label, 'add': add_name})
    fig.clf()
###############################################################################


#### Plot spectra indicating Wood anomalies in substrate ######################
def t_r_a_plots_subs(stacks_list, wavelengths, period, sub_n,
    params_layer=1, active_layer_nu=1, stack_label=1, ult_eta=False,
    J_sc=False, weight_spec=False, extinct=False, add_name=''):
    """ Plot t, r, a indicating Wood anomalies in substrate for each layer \
        & total.

        Args:
            stacks_list  (list): Stack objects containing data to plot.

            wavelengths  (list): The wavelengths corresponding to stacks_list.

            period  (float): Period of :Stack:s.

            sub_n  (float): Refractive index of the substrate in which Wood \
                anomalies are considered.

        Keyword Args:
            params_layer  (int): The index in stacks_list of the layer for \
                which the geometric parameters are put in the title of \
                the plots.

            active_layer_nu  (int): The index in stacks_list (from bottom) \
                of the layer for which the ult_eta and/or J_sc are calculated.

            stack_label  (int): Label to differentiate plots of different \
                :Stack:s.

            ult_eta  (bool): If True, calculate the 'ultimate efficiency'.

            J_sc  (bool): If True, calculate the idealised short circuit \
                current.

            weight_spec  (bool): If True, plot t, r, a spectra weighted by \
                the ASTM 1.5 solar spectrum.

            extinct  (bool): If True, calculate the extinction ratio in \
                transmission.

             add_name  (str): Add add_name to title.
    """

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
    a_tot      = []
    t_tot      = []
    r_tot      = []
    for i in range(len(wavelengths)):
        a_tot.append(float(a_list[layers_steps-1+(i*layers_steps)]))
        t_tot.append(float(t_list[layers_steps-1+(i*layers_steps)]))
        r_tot.append(float(r_list[i]))

    if ult_eta == True or J_sc == True:
        active_abs = []
        active_layer_nu = len(stacks_list[0].layers) - active_layer_nu - 1
        if not 0 < active_layer_nu < len(stacks_list[0].layers)-1:
            raise ValueError("active_layer_nu must refer to a finite layer.")
        for i in range(len(xvalues)):
            active_abs.append(float(a_list[active_layer_nu - 1 + \
                i*layers_steps]))


    if ult_eta == True:
        Efficiency = ult_efficiency(active_abs, wavelengths,
                                    params_2_print=params_2_print,
                                    stack_label=stack_label, add_name=add_name)
        params_2_print += r'$\eta$ = %(Efficiency)6.2f'% {'Efficiency' : Efficiency*100, }
        params_2_print += ' %'

    if J_sc == True:
        J = J_short_circuit(active_abs, wavelengths,
                            params_2_print=params_2_print,
                            stack_label=stack_label, add_name=add_name)
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
    if ult_eta == True or J_sc == True: del active_abs
    return


def total_tra_plot_subs(plot_name, a_spec, t_spec, r_spec, wavelengths,
    params_2_print, stack_label, add_name, period, sub_n):
    """ Plots total t, r, a spectra with lines at first 6 Wood anomalies.

        Is called from t_r_a_plots_subs
    """

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

    fig = plt.figure(num=None, figsize=(9, 12))
    ax1 = fig.add_subplot(3, 1, 1)
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
    new_tick_values = np.linspace(10, 0.5, 20)
    new_tick_locations = tick_function(new_tick_values)
    new_tick_labels = ["%.1f" % z for z in new_tick_values]
    ax2.set_xticks(new_tick_locations)
    ax2.set_xticklabels(new_tick_labels)
    ax2.set_xlim((wavelengths[0], wavelengths[-1]))
    ax2.set_xlabel('Energy (eV)')
    plt.ylim((0, 1))
    ax1 = fig.add_subplot(3, 1, 2)
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
    ax1 = fig.add_subplot(3, 1, 3)
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
    fig.clf()
###############################################################################


#### Save J_sc & ult efficiency w/o spectra ###################################
def J_sc_eta_NO_plots(stacks_list, wavelengths, params_layer=1,
    active_layer_nu=1, stack_label=1, add_name=''):
    """ Calculate J_sc & ultimate efficiency but do not save or plot spectra.

        Args:
            stacks_list  (list): Stack objects containing data to plot.

            wavelengths  (list): The wavelengths corresponding to stacks_list.

        Keyword Args:
            params_layer  (int): The index in stacks_list of the layer for \
                which the geometric parameters are put in the title of \
                the plots.

            active_layer_nu  (int): The index in stacks_list (from bottom) \
                of the layer for which the ult_eta and/or J_sc are calculated.

            stack_label  (int): Label to differentiate plots of different \
                :Stack:s.

            add_name  (str): Add add_name to title.
    """

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
    active_layer_nu = len(stacks_list[0].layers) - active_layer_nu - 1
    if not 0 < active_layer_nu < len(stacks_list[0].layers)-1:
        raise ValueError("active_layer_nu must refer to a finite layer.")
    for i in range(len(wavelengths)):
        active_abs.append(float(a_list[active_layer_nu - 1 + i*layers_steps]))

    Efficiency = ult_efficiency(active_abs, wavelengths,
                                params_2_print=params_2_print,
                                stack_label=stack_label, add_name=add_name)

    J = J_short_circuit(active_abs, wavelengths,
                        params_2_print=params_2_print,
                        stack_label=stack_label, add_name=add_name)
    return
###############################################################################


#### Saving spectra to files ##################################################
def t_r_a_write_files(stacks_list, wavelengths, stack_label=1,
    add_name=''):
    """ Save t, r, a for each layer & total in text files.

        Args:
            stacks_list  (list): Stack objects containing data to plot.

            wavelengths  (list): The wavelengths corresponding to stacks_list.

        Keyword Args:
            stack_label  (int): Label to differentiate plots of different \
                :Stack:s.

            add_name  (str): Add add_name to title.
    """
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

    total_h = sum(stacks_list[0].heights_nm()) # look at first wl result to find h.
    layers_print('Lay_Absorb', a_list, wavelengths, total_h, stack_label)
    layers_print('Lay_Trans',  t_list, wavelengths, total_h, stack_label)
    layers_print('Lay_Reflec', r_list, wavelengths, total_h, stack_label)


def layers_print(spectra_name, spec_list, wavelengths, total_h,
    stack_label=1, add_name=''):
    """ Save spectra to text files.

        Is called from t_r_a_write_files.
    """

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
###############################################################################


#### Plot spectra on other scales #############################################
def extinction_plot(t_spec, wavelengths, params_2_print, stack_label,
    add_name):
    """ Plot extinction ratio in transmission extinct = log_10(1/t). """

    fig = plt.figure(num=None, figsize=(9, 12), dpi=80, facecolor='w', edgecolor='k')
    ax1 = fig.add_subplot(1,1,1)
    extinciton = np.log10(1.0/np.array(t_spec))
    ax1.plot(wavelengths, extinciton, linewidth=linesstrength)
    ax1.set_xlabel('Wavelength (nm)')
    ax1.set_ylabel('Extinction')
    plot_name = 'Extinction'
    plt.suptitle(plot_name + add_name + '\n' + params_2_print, \
        fontsize = title_font)
    plt.savefig('%(s)s_stack%(bon)s_%(add)s'% \
        {'s' : plot_name, 'bon' : stack_label,'add' : add_name})
    fig.clf()


def EOT_plot(stacks_list, wavelengths, pol='TM', params_layer=1,
    add_name='', savetxt=False):
    """ Plot T_{00} as in Martin-Moreno PRL 86 2001.

        Args:
            pol  (str): Take the (0,0) component for TE->TE scattering
                or the (0,0) component for TM->TM scattering.
    """

    height_list = stacks_list[0].heights_nm()[::-1]
    params_2_print = gen_params_string(stacks_list, params_layer)
    params_2_print += r'$h_t,...,h_b$ = '
    params_2_print += ''.join('%4d, ' % num for num in height_list)

    if pol == 'TM':
        num_pw_per_pol = stacks_list[0].layers[-1].structure.num_pw_per_pol
    else:
        num_pw_per_pol = 0

    T_00 = []
    for i in range(len(wavelengths)):
        t00  = stacks_list[i].T_net[num_pw_per_pol,num_pw_per_pol]
        t00c = t00.conjugate()
        T_00.append(np.real(t00*t00c))

    fig = plt.figure()
    ax1 = fig.add_subplot(2, 1, 1)
    ax1.plot(wavelengths, T_00, linewidth = linesstrength)
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
    plt.savefig('%(s)s%(add)s'% {'s' : plot_name, 'add' : add_name})
    fig.clf()

    if savetxt == True:
        np.savetxt('Trans_%(s)s%(add)s.txt'% {'s' : plot_name, 'add' : add_name}, \
        zip(wavelengths, T_00), fmt = '%10.6f')
        np.savetxt('Refl_%(s)s%(add)s.txt'% {'s' : plot_name, 'add' : add_name}, \
        zip(wavelengths, R_00), fmt = '%10.6f')
###############################################################################


#### Calculate short circuit current and ultimate efficiency ##################
def J_short_circuit(active_abs, wavelengths, params_2_print='',
                    stack_label='', add_name=''):
    """ Calculate the short circuit current J_sc under ASTM 1.5 illumination.
        Assuming every absorbed photon produces a pair of charge carriers.
    """

    Irrad_spec_file = '../backend/data/ASTMG173'
    i_data = np.loadtxt('%s.txt' % Irrad_spec_file)
    i_spec = np.interp(wavelengths, i_data[:, 0], i_data[:, 3])
    expression = i_spec*active_abs*wavelengths
    integral_tmp = np.trapz(expression, x=wavelengths)
    J = (charge_e/(Plancks_h*speed_c)) * integral_tmp * 1e-10  # in mA/cm^2
    nums_2_print = params_2_print.split()
    np.savetxt('J_sc_stack%(bon)s%(add)s.txt' % {'bon': stack_label,
               'add': add_name}, np.array([J]), fmt='%10.6f')
    return J


def ult_efficiency(active_abs, wavelengths, bandgap_wl=None,
                   params_2_print='', stack_label='', add_name=''):
    """ Calculate the photovoltaic ultimate efficiency achieved in the specified active layer.

        For definition see `Sturmberg et al., Optics Express, Vol. 19, Issue S5, pp. A1067-A1081 (2011)\
         <http://dx.doi.org/10.1364/OE.19.0A1067>`_.

        Args:
            bandgap_wl  (float): allows you to set the wavelength equivalent
                to the bandgap. Else it is assumed to be the maximum wavelength
                simulated.
    """
    if bandgap_wl is None:
        bandgap_wl = np.max(wavelengths)
    Irrad_spec_file = '../backend/data/ASTMG173'
    i_data       = np.loadtxt('%s.txt' % Irrad_spec_file)
    i_spec       = np.interp(wavelengths, i_data[:,0], i_data[:,3])
    expression   = i_spec*active_abs*wavelengths
    integral_tmp = np.trapz(expression, x=wavelengths)
    Efficiency   = integral_tmp/(bandgap_wl*ASTM15_tot_I)
    nums_2_print = params_2_print.split()
    np.savetxt('Efficiency_stack%(bon)s%(add)s.txt'% {'bon' : stack_label,'add' : add_name},\
        np.array([Efficiency]), fmt = '%8.6f')
    return Efficiency
###############################################################################


#### Plot dispersion diagrams & field concentrations function of wavelength ###
def omega_plot(stacks_list, wavelengths, params_layer=1, stack_label=1):
    """ Plots the dispersion diagram of each layer in one plot. \
        k_z has units nm^-1.

        Args:
            stacks_list  (list): Stack objects containing data to plot.

            wavelengths  (list): The wavelengths corresponding to stacks_list.

        Keyword Args:
            params_layer  (int): The index in stacks_list of the layer for \
                which the geometric parameters are put in the title of \
                the plots.

            stack_label  (int): Label to differentiate plots of different \
                :Stack:s.
    """

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
        ax1.set_ylabel(r'Real $k_z$ (d)'), ax2.set_ylabel(r'Imaginary $k_z$ (d)')
        ax3.set_ylabel(r'Frequency ($\omega$d/2$\pi$c)'), ax4.set_ylabel(r'Frequency ($\omega$d/2$\pi$c)')
        if l == 0:
            ax1.set_ylabel('Bottom Layer'), ax2.set_ylabel('Bottom Layer')
            ax3.set_ylabel('Bottom Layer'), ax4.set_ylabel('Bottom Layer')
            ax3.set_xlabel(r'Real $k_z$ (d)'), ax4.set_xlabel(r'Imaginary $k_z$ (d)')
            ax1.set_xlabel('Wavelength (nm)'), ax2.set_xlabel('Wavelength (nm)')
        # else:
            # ax1.set_xticklabels( () )
            # ax2.set_xticklabels( () )
            # ax3.set_xticklabels( () )
            # ax4.set_xticklabels( () )
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


def E_conc_plot(stacks_list, which_layer, which_modes, wavelengths,
    params_layer=1, stack_label=1):
    """ Plots the energy concentration (epsilon E_cyl / epsilon E_cell) of given layer.

        Args:
            stacks_list  (list): Stack objects containing data to plot.

            which_layer  (int): The index in stacks_list of the layer for \
                which the energy concentration is to be calculated.

            which_modes  (list): Indices of Bloch modes for which to calculate \
                the energy concentration.

            wavelengths  (list): The wavelengths corresponding to stacks_list.

        Keyword Args:
            params_layer  (int): The index in stacks_list of the layer for \
                which the geometric parameters are put in the title of \
                the plots.

            stack_label  (int): Label to differentiate plots of different \
                :Stack:s.
    """

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
        fig1.suptitle('Energy Concentration = ' + r'$E_{cyl} / E_{cell}$' + '\n' + \
            params_2_print, fontsize=title_font)
        fig1.savefig('Energy_Concentration_stack%(bon)s'% {'bon' : stack_label}, \
            bbox_inches='tight')
    else:
        print("\nERROR: plotting.E_conc_plot; \n" + \
            "Can only calculate energy concentration in NanoStruct layers.")
        print(repr(stacks_list[0].layers[which_layer]))
###############################################################################


#### Visualise scattering Matrices ############################################
def vis_scat_mats(scat_mat, nu_prop_PWs=0, wl=None, add_name='',
    max_scale=None):
    """ Plot given scattering matrix as greyscale images.

        Args:
            scat_mat  (np.matrix): A scattering matrix, which is \
                organised as \
                | TE -> TE | TM -> TE | \
                | TE -> TM | TM -> TM |

        Keyword Args:
            nu_prop_PWs  (int): Number of propagating PWs.

            wl  (int): Index in case of calling in a loop.

            add_name  (str): Add add_name to title.

            max_scale  (float): Limit maximum amplitude shown.
    """

    fig = plt.figure(num=None, figsize=(14, 6), dpi=80, facecolor='w', edgecolor='k')
    if wl != None: add_name = '-wl%(wl)i-%(ti)s' % {'wl' : wl, 'ti' : add_name}

    for i in [1,2]:
        ax1 = fig.add_subplot(1,2,i)
        if i==1: image = abs(np.real(scat_mat))
        if i==2: image = abs(np.imag(scat_mat))
        if max_scale != None:
            mat = ax1.matshow(image, vmax = max_scale, cmap=plt.cm.gray)
        else:
            mat = ax1.matshow(image, cmap=plt.cm.gray)
        scat_mat_dim_x = np.shape(scat_mat)[0]
        scat_mat_dim_y = np.shape(scat_mat)[1]
        half_dim_x = scat_mat_dim_x/2-0.5
        half_dim_y = scat_mat_dim_y/2-0.5
        ax1.plot([-0.5, scat_mat_dim_x-0.5],[half_dim_y,half_dim_y],'w', linewidth=2)
        ax1.plot([half_dim_x,half_dim_x],[-0.5, scat_mat_dim_y-0.5],'w', linewidth=2)
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
        if i==1: cbar.ax.set_ylabel('|Real(matrix)|', fontsize=14)
        if i==2: cbar.ax.set_ylabel('|Imag(matrix)|', fontsize=14)
        # if i==1: ax1.set_title('|Re(matrix)|', fontsize=14)
        # if i==2: ax1.set_title('|Imag(matrix)|', fontsize=14)
        ax1.set_xticklabels('')
        ax1.set_yticklabels('')
        ax1.set_xlabel('TE        |        TM \n Incoming Orders', fontsize=14)
        ax1.set_ylabel('Outgoing Orders \nTM       |         TE', fontsize=14)


    plt.suptitle('Scattering Matrices' + add_name)
    plt.savefig('Scat_mat' + add_name)
    fig.clf()


def vis_matrix(scat_mat, add_name='', max_scale=None, only_real=False):
    """ Plot given matrix as a greyscale image.

        Args:
            scat_mat  (np.matrix): A matrix.

        Keyword Args:
            add_name  (str): Add add_name to title.

            max_scale  (float): Limit maximum amplitude shown.

            only_real  (bool): Only plot the real part of matrix.
    """

    fig = plt.figure(num=None, figsize=(10, 8), dpi=80, facecolor='w', edgecolor='k')

    if only_real == True: real_im = [1]
    else: real_im = [1,2]

    for i in real_im:
        ax1 = fig.add_subplot(1,len(real_im),i)
        if i==1: image = abs(np.real(scat_mat))
        if i==2: image = abs(np.imag(scat_mat))
        if max_scale != None:
            mat = ax1.matshow(image, vmax = max_scale, cmap=plt.cm.gray)
        else:
            mat = ax1.matshow(image, cmap=plt.cm.gray)
        cbar = fig.colorbar(mat,extend='neither')
        if i==1: cbar.ax.set_ylabel('|Re(matrix)|', fontsize=14)
        if i==2: cbar.ax.set_ylabel('|Im(matrix)|', fontsize=14)
        ax1.set_xticklabels('')
        ax1.set_yticklabels('')
        # ax1.set_xlabel('TE        |        TM \n Incoming Orders', fontsize=14)
        # ax1.set_ylabel('Outgoing Orders \nTM       |         TE', fontsize=14)

    plt.suptitle(add_name)
    plt.savefig('Matrix' + add_name)
    fig.clf()
###############################################################################


#### Plot PW amplitudes function k-vector #####################################
def t_func_k_plot_1D(stacks_list, lay_interest=0, pol='TE'):
    """ PW amplitudes in transmission as a function of their in-plane k-vector.

        Args:
            stacks_list  (list): Stack objects containing data to plot.

        Keyword Args:
            lay_interest  (int): The index in stacks_list of the layer in \
                which amplitudes are calculated.

            pol  (str): Include transmission in Which polarisation.

    """
    fig = plt.figure(num=None, dpi=80, facecolor='w', edgecolor='k')
    ax1 = fig.add_subplot(1,1,1)

    # vec_coef sorted from top of stack, everything else is sorted from bottom
    vec_index = len(stacks_list[0].layers) - lay_interest - 1

    ## Old code if selecting betas = beta0 + 0 value out of 2D PW basis.
    # # Create arrays of grating order indexes (-p, ..., p)
    # max_ords = stacks_list[0].layers[-1].max_order_PWs
    # pxs = np.arange(-max_ords, max_ords + 1)

    store_alphas = []
    store_k_trans = []
    for stack in stacks_list:
        n_PW_p_pols = stack.layers[0].structure.num_pw_per_pol
        sort_order = stack.layers[lay_interest].sort_order

        select_trans  = abs(stack.vec_coef_down[vec_index][0:n_PW_p_pols])
        store_alphas  = np.append(store_alphas,stack.layers[lay_interest].alphas[sort_order])
        store_k_trans = np.append(store_k_trans,select_trans)

        ## Old code if selecting betas = beta0 + 0 value out of 2D PW basis.
        # k0 = stack.layers[0].k()
        # # Calculate k_x that correspond to k_y = beta0 (in normalized units)
        # alpha0, beta0 = stack.layers[0].k_pll_norm()
        # alphas = alpha0 + pxs * 2 * np.pi
        # on_axis_kzs = sqrt(k0**2 - alphas**2 - beta0**2)
        # full_k_space = stack.layers[0].k_z
        # # consider singular polarisation
        # n_PW_p_pols = stack.layers[0].structure.num_pw_per_pol
        # one_pol_k_space = full_k_space[0:n_PW_p_pols]

        # axis_indices = []
        # for a in on_axis_kzs:
        #     ix = np.in1d(one_pol_k_space.ravel(), a).reshape(one_pol_k_space.shape)
        #     axis_indices = np.append(axis_indices, np.ravel(np.array(np.where(ix))))
        # axis_indices = axis_indices.astype(int)

        # # Outgoing TE polarisation
        # if pol=='TE': trans_k = np.abs(stack.vec_coef_down[vec_index][0:n_PW_p_pols]).reshape(-1,)
        # # Outgoing TM polarisation
        # if pol=='TM': trans_k = np.abs(stack.vec_coef_down[vec_index][n_PW_p_pols-1:-1]).reshape(-1,)
        # trans_k_array = np.array(trans_k).reshape(-1,)

        # select_trans = trans_k_array[axis_indices]
        # store_alphas = np.append(store_alphas,alphas)
        # store_k_trans = np.append(store_k_trans,select_trans)

    sort_indices = np.argsort(store_alphas)
    plot_alphas = store_alphas[sort_indices]
    plot_k_trans = store_k_trans[sort_indices]

    ax1.plot(plot_alphas,plot_k_trans, linewidth=1.5)

    min_k_label = np.max(np.abs(plot_alphas))
    k0 = abs(stack.layers[lay_interest].k())
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
    fig.clf()
###############################################################################


#### Plot amplitudes of modes #################################################
def BM_amplitudes(stacks_list, xvalues=None, chosen_BMs=None,
    lay_interest=1, up_and_down=True, add_height=None, add_name='',
    save_pdf=True, save_npz=False, save_txt=False):
    """ Plot the amplitudes of Bloch modes in selected layer.

        Args:
            stacks_list  (list): Stack objects containing data to plot.

        Keyword Args:
            xvalues  (list): The values stacks_list is to be plotted as a \
                function of.

            chosen_BMs  (list): Bloch Modes to include, identified by their \
                indices in the scattering matrices (order most propagating \
                to most evanescent) eg. [0,2,4].

            lay_interest  (int): The index in stacks_list of the layer in \
                which amplitudes are calculated.

            up_and_down  (bool): Average the amplitudes of up & downward \
                propagating modes. Else include only downward in all layers\
                except for the superstrate where include only upward.

            add_height  (float): Print the heights of :Stack: layer in title.

            add_name  (str): Add add_name to title.

            save_pdf  (bool): If True save spectra as pdf files. \
                True by default.

            save_npz  (bool): If True, saves lists of BM amplitudes to file.

            save_txt  (bool): If True, saves lists of BM amps to txt file.
    """

    fig = plt.figure(num=None, dpi=80, facecolor='w', edgecolor='k')
    ax1 = fig.add_subplot(1,1,1)
    # vec_coef sorted from top of stack, everything else is sorted from bottom
    vec_index = len(stacks_list[-1].layers) - lay_interest - 1

    xlabel = 'xvalues'
    if xvalues==None:
        if stacks_list[0].layers[0].light.wl_nm != stacks_list[-1].layers[0].light.wl_nm:
            xvalues = [s.layers[0].light.wl_nm for s in stacks_list]
            xlabel = r'$\lambda$ (nm)'
        elif set(stacks_list[0].layers[0].light.k_pll) != set(stacks_list[-1].layers[0].light.k_pll):
            xvalues = [np.sqrt(s.layers[0].light.k_pll[0]**2 + s.layers[0].light.k_pll[1]**2) for s in stacks_list]
            xlabel = r'$|k_\parallel|$'
        else:
            xvalues = [s.layers[0].light.wl_nm for s in stacks_list]
            xlabel = r'$\lambda$ (nm)'
            print("BM_amplitudes is guessing you have a single wavelength, else specify xvalues.")

    if chosen_BMs == None: chosen_BMs = range(stacks_list[-1].layers[lay_interest].num_BMs)
    try:
        save_trans = []
        for BM in chosen_BMs:
            store_trans = []
            for stack in stacks_list:
                if not isinstance(stack.layers[lay_interest],objects.Simmo):
                    raise ValueError

                trans = np.abs(stack.vec_coef_down[vec_index][BM])
                if up_and_down == True: # Take average of up & downward propagating modes.
                    trans += np.abs(stack.vec_coef_up[vec_index][BM])
                    store_trans = np.append(store_trans,trans/2)
                else:
                    store_trans = np.append(store_trans,trans)
            if save_pdf == True:
                ax1.plot(xvalues,store_trans, label="BM %i" % BM)
            if save_npz == True or save_txt == True:
                save_trans.append(store_trans)

        if save_pdf == True or save_npz == True or save_txt == True:
            if add_height!= None: add_name += '_' + zeros_int_str(add_height)
            add_name = str(lay_interest) + add_name

        if save_pdf == True:
            handles, labels = ax1.get_legend_handles_labels()
            lgd = ax1.legend(handles, labels, loc='center left', bbox_to_anchor=(1.0,0.5))
            ax1.set_ylabel('BM Amplitude')
            ax1.set_xlabel(xlabel)
            plt.suptitle(add_name)
            plt.savefig('BM_amplitudes-lay_%s' % add_name, \
                fontsize=title_font, bbox_extra_artists=(lgd,), bbox_inches='tight')

        if save_npz == True:
            np.savez('BM_amplitudes-lay_%s' % add_name, save_trans=save_trans)

        if save_txt == True:
            np.savetxt('BM_amplitudes-lay_%s.txt' % add_name, \
                save_trans, fmt = '%18.11f')

    except ValueError:
        print("BM_amplitudes only works in NanoStruct layers."\
        "\nPlease select lay_interest !=%i.\n" % lay_interest)

    fig.clf()

def PW_amplitudes(stacks_list, xvalues=None, chosen_PWs=None,
    lay_interest=0, up_and_down=True, add_height=None, add_name='',
    save_pdf=True, save_npz=False, save_txt=False):
    """ Plot the amplitudes of plane wave orders in selected layer.

        Assumes dealing with 1D grating and only have 1D diffraction orders.
        Takes the average of up & downward propagating modes.

        Args:
            stacks_list  (list): Stack objects containing data to plot.

        Keyword Args:
            xvalues  (list): The values stacks_list is to be plotted as a \
                function of.

            chosen_PWs  (list): PW diffraction orders to include. \
                eg. [-1,0,2]. If 'None' are given all are plotted.

            lay_interest  (int): The index in stacks_list of the layer in \
                which amplitudes are calculated.

            up_and_down  (bool): Average the amplitudes of up & downward \
                propagating modes. Else include only downward in all layers\
                except for the superstrate where include only upward.

            add_height  (float): Print the heights of :Stack: layer in title.

            add_name  (str): Add add_name to title.

            save_pdf  (bool): If True save spectra as pdf files. \
                True by default.

            save_npz  (bool): If True, saves lists of PW amplitudes to file.

            save_txt  (bool): If True, saves lists of PW amps to txt file.
    """

    fig = plt.figure(num=None, dpi=80, facecolor='w', edgecolor='k')
    ax1 = fig.add_subplot(1,1,1)
    # vec_coef sorted from top of stack, everything else is sorted from bottom
    vec_index = len(stacks_list[-1].layers) - lay_interest - 1
    if chosen_PWs == None:
        # Create arrays of grating order indexes (-p, ..., p)
        max_ords = stacks_list[0].layers[-1].max_order_PWs
        chosen_PWs = np.arange(-max_ords, max_ords + 1)

    xlabel = 'xvalues'
    if xvalues==None:
        if stacks_list[0].layers[0].light.wl_nm != stacks_list[-1].layers[0].light.wl_nm:
            xvalues = [s.layers[0].light.wl_nm for s in stacks_list]
            xlabel = r'$\lambda$ (nm)'
        elif set(stacks_list[0].layers[0].light.k_pll) != set(stacks_list[-1].layers[0].light.k_pll):
            xvalues = [np.sqrt(s.layers[0].light.k_pll[0]**2 + s.layers[0].light.k_pll[1]**2) for s in stacks_list]
            xlabel = r'$|k_\parallel|$'
        else:
            xvalues = [s.layers[0].light.wl_nm for s in stacks_list]
            xlabel = r'$\lambda$ (nm)'
            print("PW_amplitudes is guessing you have a single wavelength, else specify xvalues.")

    try:
        save_trans = []
        for pxs in chosen_PWs:
            store_trans = []
            for stack in stacks_list:
                if not isinstance(stack.layers[lay_interest],objects.Anallo):
                    raise ValueError

                k0 = stack.layers[0].k()
                n_PW_p_pols = stack.layers[0].structure.num_pw_per_pol
                # Calculate k_x that correspond to k_y = beta0 = 0 (in normalized units)
                alpha0, beta0 = stack.layers[0].k_pll_norm()
                alphas = alpha0 + pxs * 2 * np.pi
                on_axis_kzs = sqrt(k0**2 - alphas**2 - beta0**2)
                full_k_space = stack.layers[0].k_z
                # Consider only transmission into singular polarization.
                one_pol_k_space = full_k_space[0:n_PW_p_pols]

                ix = np.in1d(one_pol_k_space.ravel(), on_axis_kzs).reshape(one_pol_k_space.shape)
                axis_indices = np.ravel(np.array(np.where(ix))).astype(int)
                # Substrate - only ever take downward propagating.
                if vec_index == len(stacks_list[-1].layers) - 1:
                    # Outgoing TE polarisation
                    trans = np.abs(stack.vec_coef_down[vec_index][axis_indices]).reshape(-1,)
                    # Outgoing TM polarisation
                    trans += np.abs(stack.vec_coef_down[vec_index][n_PW_p_pols+axis_indices]).reshape(-1,)
                    store_trans = np.append(store_trans,trans)
                # Superstrate - if not up & down then take only upwards propagating.
                elif vec_index == 0:
                    trans  = np.abs(stack.vec_coef_up[vec_index][axis_indices]).reshape(-1,)
                    trans += np.abs(stack.vec_coef_up[vec_index][n_PW_p_pols+axis_indices]).reshape(-1,)
                    if up_and_down == True: # Take average of up & downward propagating modes.
                        trans += np.abs(stack.vec_coef_down[vec_index][axis_indices]).reshape(-1,)
                        trans += np.abs(stack.vec_coef_down[vec_index][n_PW_p_pols+axis_indices]).reshape(-1,)
                        store_trans = np.append(store_trans,trans/2)
                    else:
                        store_trans = np.append(store_trans,trans)
                # Finite layer
                else:
                    trans = np.abs(stack.vec_coef_down[vec_index][axis_indices]).reshape(-1,)
                    trans += np.abs(stack.vec_coef_down[vec_index][n_PW_p_pols+axis_indices]).reshape(-1,)
                    if up_and_down == True: # Take average of up & downward propagating modes.
                        trans += np.abs(stack.vec_coef_up[vec_index][axis_indices]).reshape(-1,)
                        trans += np.abs(stack.vec_coef_up[vec_index][n_PW_p_pols+axis_indices]).reshape(-1,)
                        store_trans = np.append(store_trans,trans/2)
                    else:
                        store_trans = np.append(store_trans,trans)
            if save_pdf == True:
                ax1.plot(xvalues,store_trans, label="m = %i" % pxs)
            if save_npz == True or save_txt == True:
                save_trans.append(store_trans)

        if save_pdf == True or save_npz == True:
            if add_height!= None: add_name += '_' + zeros_int_str(add_height)
            add_name = str(lay_interest) + add_name

        if save_pdf == True:
            handles, labels = ax1.get_legend_handles_labels()
            lgd = ax1.legend(handles, labels, loc='center left', bbox_to_anchor=(1.0,0.5))
            ax1.set_ylabel(r'$|E|$')
            ax1.set_xlabel(xlabel)
            plt.suptitle(add_name)
            plt.savefig('PW_amplitudes-lay_%s' % add_name, \
                fontsize=title_font, bbox_extra_artists=(lgd,), bbox_inches='tight')

        if save_npz == True:
            np.savez('PW_amplitudes-lay_%s' % add_name, save_trans=save_trans)

        if save_txt == True:
            np.savetxt('PW_amplitudes-lay_%s.txt' % add_name, \
                save_trans, fmt = '%18.11f')

    except ValueError:
        print("PW_amplitudes only works in ThinFilm layers."\
        "\nPlease select lay_interest !=%i.\n" % lay_interest)

    fig.clf()


def evanescent_merit(stacks_list, xvalues=None, chosen_PWs=None,
    lay_interest=0, add_height=None, add_name='',
    save_pdf=True, save_txt=False):
    """ Plot a figure of merit for the 'evanescent-ness' of excited fields.

        Assumes dealing with 1D grating and only have 1D diffraction orders.

        Args:
            stacks_list  (list): Stack objects containing data to plot.

        Keyword Args:
            xvalues  (list): The values stacks_list is to be plotted as a \
                function of.

            chosen_PWs  (list): PW diffraction orders to include. \
                eg. [-1,0,2].

            lay_interest  (int): The index in stacks_list of the layer in \
                which amplitudes are calculated.

            add_height  (float): Print the heights of :Stack: layer in title.

            add_name  (str): Add add_name to title.

            save_pdf  (bool): If True save spectra as pdf files. \
                True by default.

            save_txt  (bool): If True, saves average value of \
                mean PW order to file.
    """
    # vec_coef sorted from top of stack, everything else is sorted from bottom
    vec_index = len(stacks_list[-1].layers) - lay_interest - 1
    if chosen_PWs == None:
        # Create arrays of grating order indexes (-p, ..., p)
        max_ords = stacks_list[0].layers[-1].max_order_PWs
        chosen_PWs = np.arange(-max_ords, max_ords + 1)
    n_H = max_n(stacks_list)

    xlabel = 'xvalues'
    if xvalues==None:
        if stacks_list[0].layers[0].light.wl_nm != stacks_list[-1].layers[0].light.wl_nm:
            xvalues = [s.layers[0].light.wl_nm for s in stacks_list]
            xlabel = r'$\lambda$ (nm)'
        elif set(stacks_list[0].layers[0].light.k_pll) != set(stacks_list[-1].layers[0].light.k_pll):
            xvalues = [np.sqrt(s.layers[0].light.k_pll[0]**2 + s.layers[0].light.k_pll[1]**2) for s in stacks_list]
            xlabel = r'$|k_\parallel|$'
        else:
            xvalues = [s.layers[0].light.wl_nm for s in stacks_list]
            xlabel = r'$\lambda$ (nm)'
            print("evanescent_merit is guessing you have a single wavelength, else specify xvalues.")

    store_m_p      = []
    store_m_ne     = []
    store_m_fe     = []
    store_x_p      = []
    store_x_ne     = []
    store_x_fe     = []
    store_tot_amps = []
    store_mean_ev  = []
    s = 0
    try:
        for stack in stacks_list:
            if not isinstance(stack.layers[lay_interest],objects.Anallo):
                raise ValueError

            merit_prop    = 0.0
            merit_near_ev = 0.0
            merit_far_ev  = 0.0
            sum_p_amps    = 0.0
            sum_amps      = 0.0
            for pxs in chosen_PWs:
                k0 = stack.layers[-1].k() # Incident k0
                k = stack.layers[lay_interest].k() # k in film
                n_PW_p_pols = stack.layers[lay_interest].structure.num_pw_per_pol
                # Calculate k_x that correspond to k_y = beta0 = 0 (in normalized units)
                alpha0, beta0 = stack.layers[lay_interest].k_pll_norm()
                alphas = alpha0 + pxs * 2 * np.pi
                this_k_pll2 = alphas**2 + beta0**2
                on_axis_kzs = sqrt(k**2 - alphas**2 - beta0**2)
                full_k_space = stack.layers[lay_interest].k_z
                # Consider only transmission into singular polarization
                one_pol_k_space = full_k_space[0:n_PW_p_pols]

                ix = np.in1d(one_pol_k_space.ravel(), on_axis_kzs).reshape(one_pol_k_space.shape)
                axis_indices = np.ravel(np.array(np.where(ix))).astype(int)

                # Outgoing TE polarisation
                trans = np.abs(stack.vec_coef_down[vec_index][axis_indices])
                # Outgoing TM polarisation
                trans += np.abs(stack.vec_coef_down[vec_index][n_PW_p_pols +axis_indices])

                if np.abs(this_k_pll2) < np.abs(k**2): merit_prop += trans
                if np.abs(k**2) < np.abs(this_k_pll2) < np.abs((n_H*k0)**2):
                    merit_near_ev += trans
                if np.abs(this_k_pll2) > np.abs((n_H*k0)**2): merit_far_ev += trans
                sum_p_amps += np.abs(pxs) * trans
                sum_amps += np.abs(stack.vec_coef_down[vec_index][axis_indices])**2 + np.abs(stack.vec_coef_down[vec_index][n_PW_p_pols +axis_indices])**2

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
            store_tot_amps = np.append(store_tot_amps,sum_amps)
            store_mean_ev  = np.append(store_mean_ev,sum_p_amps/sum_amps)

        if add_height!= None: add_name += '_'+zeros_int_str(add_height)
        add_name = str(lay_interest) + add_name
        if save_pdf == True:
            fig = plt.figure(num=None, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
            ax1 = fig.add_subplot(2,1,1)
            if len(store_m_p)  != 0: ax1.plot(store_x_p,store_m_p, 'b', label="prop")
            if len(store_m_ne) != 0: ax1.plot(store_x_ne,store_m_ne, 'g', label="near ev")
            if len(store_m_fe) != 0: ax1.plot(store_x_fe,store_m_fe, 'r', label="far ev")
            ax1.plot(xvalues, store_mean_ev, 'k', label=r'$\Sigma|p| |a_p| / \Sigma|a_p|$')
            ax2 = fig.add_subplot(2, 1, 2)
            ax2.plot(xvalues, store_tot_amps, 'k', label="prop")

            handles, labels = ax1.get_legend_handles_labels()
            lgd = ax1.legend(handles, labels, ncol=4, loc='upper center',
                             bbox_to_anchor=(0.5, 1.6))
            ax1.set_ylabel('Ev FoM')
            ax1.set_xticklabels(())
            ax2.set_ylabel(r'$\Sigma|a_p|^2$')
            ax2.set_xlabel(xlabel)
            plt.suptitle(add_name)
            plt.savefig('evanescent_merit-lay_%s' % add_name,
                        fontsize=title_font, bbox_extra_artists=(lgd,),
                        bbox_inches='tight')

        if save_txt is True:
            av_diff = [np.mean(store_m_p)]
            np.savetxt('prop_diff_order-lay_%s.txt' % add_name,
                       av_diff, fmt='%18.11f')
            av_diff = [np.mean(store_m_ne)]
            np.savetxt('near_ev_diff_order-lay_%s.txt' % add_name,
                       av_diff, fmt='%18.11f')
            av_diff = [np.mean(store_m_fe)]
            np.savetxt('far_ev_diff_order-lay_%s.txt' % add_name,
                       av_diff, fmt='%18.11f')
            av_diff = [np.mean(store_mean_ev)]
            np.savetxt('average_diff_order-lay_%s.txt' % add_name,
                       av_diff, fmt='%18.11f')

    except ValueError:
        print("evanescent_merit only works in ThinFilm layers."\
        "\nPlease select lay_interest !=%i.\n" % lay_interest)

    fig.clf()
###############################################################################


#### Field plotting routines ##################################################
def fields_in_plane(stacks_list, lay_interest=1, z_values=[0.1, 3.6],
                    nu_calc_pts=51):
    """
    Plot fields in the x-y plane at chosen values of z.

        Args:
            stacks_list  (list): Stack objects containing data to plot.

        Keyword Args:
            lay_interest  (int): the index of the layer considered within \
                the stack.

            z_values  (float): distance in nm from bottom surface of layer \
                at which to calculate fields. If layer is semi-inf substrate \
                then z_value is distance from top of this layer (i.e. bottom \
                interface of stack).

            nu_calc_pts  (int): fields are calculated over a mesh of \
                nu_calc_pts * nu_calc_pts points.
    """
    from fortran import EMUstack

    dir_name = "in_plane_fields"
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)

    # always make odd
    if nu_calc_pts % 2 == 0:
        nu_calc_pts += 1

    stack_num = 0
    for pstack in stacks_list:
        z_values = np.array(z_values)/float(pstack.layers[lay_interest].structure.period)
        num_lays = len(pstack.layers)
        # If selected z location is within a NanoStruct layer
        # plot fields in Bloch Mode Basis using fortran routine.
        if isinstance(pstack.layers[lay_interest], mode_calcs.Simmo):
            meat = pstack.layers[lay_interest]
            try:
                if not meat.structure.periodicity == '2D_array':
                    raise ValueError
                gmsh_file_pos = meat.structure.mesh_file[0:-5]
                eps = meat.n_effs**2
                h_normed = float(meat.structure.height_nm)/float(meat.structure.period)
                for z_value in z_values:
                    # Fortran routine naturally measure z top down
                    z_value = h_normed - z_value
                    wl_normed = pstack.layers[lay_interest].wl_norm()

                    nnodes = 6
                    if meat.E_H_field == 1:
                        EH_name = "E_"
                    else:
                        EH_name = "H_"
                    extra_name = EH_name + "Lay" + zeros_int_str(lay_interest)

                    # vec_coef sorted from top of stack, everything else is sorted from bottom
                    vec_index = num_lays - lay_interest - 1
                    vec_coef = np.concatenate((pstack.vec_coef_down[vec_index],
                                              pstack.vec_coef_up[vec_index]))

                    EMUstack.gmsh_plot_field(meat.num_BMs,
                                             meat.n_msh_el, meat.n_msh_pts,
                                             nnodes, meat.structure.nb_typ_el,
                                             meat.table_nod, meat.type_el, eps,
                                             meat.x_arr, meat.k_z, meat.sol1,
                                             vec_coef, h_normed, z_value,
                                             gmsh_file_pos,
                                             meat.structure.plot_real,
                                             meat.structure.plot_imag,
                                             meat.structure.plot_abs,
                                             extra_name)

            except ValueError:
                print("fields_in_plane cannot plot fields in 1D-arrays."\
                "\nPlease select a different lay_interest.\n")


        # If selected z location is within a ThinFilm layer
        # plot fields in Plane Wave Basis using python routine.
        else:
            wl = np.around(pstack.layers[-1].light.wl_nm,decimals=2)
            pw = pstack.layers[-1].max_order_PWs
            period = pstack.layers[-1].structure.period
            heights_list = []
            name_lay = ''

            for i in xrange(num_lays):
                if i == 0 or i == num_lays-1:pass
                else: heights_list.append(pstack.layers[i].structure.height_nm)
                try:
                    pstack.layers[i].structure.diameter1
                    diameter = pstack.layers[i].structure.diameter1
                except:pass

            if lay_interest == 0: name_lay = "0_Substrate"
            elif lay_interest == num_lays-1: name_lay = "%(i)s_Superstrate" %{'i':num_lays-1}
            else: name_lay = "Thin_Film_%(i)s" %{'i':lay_interest}

            x_range = np.linspace(0.0,1.0,nu_calc_pts)
            y_range = np.linspace(0.0,1.0,nu_calc_pts)
            z_range = np.array(z_values)

            s = pstack.layers[lay_interest].sort_order
            alpha_unsrt = np.array(pstack.layers[lay_interest].alphas)
            beta_unsrt = np.array(pstack.layers[lay_interest].betas)
            alpha = alpha_unsrt[s]
            if pstack.layers[lay_interest].structure.world_1d == True:
                beta = beta_unsrt
            else:
                beta = beta_unsrt[s]
            gamma = np.array(pstack.layers[lay_interest].calc_kz())
            n = pstack.layers[lay_interest].n()
            PWordtot = pstack.layers[lay_interest].structure.num_pw_per_pol
            prop = 2*(gamma.imag == 0).sum()
            evan = 2*PWordtot - prop

            vec_coef_down = np.array(pstack.vec_coef_down[num_lays-1-lay_interest]).flatten()
            vec_coef_down_TE = vec_coef_down[0:PWordtot]
            vec_coef_down_TM = vec_coef_down[PWordtot::]
            if lay_interest == 0:
                vec_coef_up = np.zeros((2*PWordtot), dtype = 'complex')
            else:
                vec_coef_up = np.array(pstack.vec_coef_up[num_lays-1-lay_interest]).flatten()

            vec_coef_up_TE = vec_coef_up[0:PWordtot]
            vec_coef_up_TM = vec_coef_up[PWordtot::]

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
            k = np.around(k,decimals=4)
            n = np.around(n,decimals=4)

            eta_TE_x_down = (vec_coef_down_TE*E_TE_x)/chi_TE
            eta_TE_y_down = (vec_coef_down_TE*E_TE_y)/chi_TE
            eta_TE_z_down = (vec_coef_down_TE*E_TE_z)/chi_TE
            eta_TM_x_down = (vec_coef_down_TM*E_TM_x)/chi_TM
            eta_TM_y_down = (vec_coef_down_TM*E_TM_y)/chi_TM
            eta_TM_z_down = (vec_coef_down_TM*E_TM_z)/chi_TM

            eta_TE_x_up = (vec_coef_up_TE*E_TE_x)/chi_TE
            eta_TE_y_up = (vec_coef_up_TE*E_TE_y)/chi_TE
            eta_TE_z_up = (vec_coef_up_TE*E_TE_z)/chi_TE
            eta_TM_x_up = (vec_coef_up_TM*E_TM_x)/chi_TM
            eta_TM_y_up = (vec_coef_up_TM*E_TM_y)/chi_TM
            eta_TM_z_up = (vec_coef_up_TM*E_TM_z)/chi_TM

            x_axis = np.zeros((nu_calc_pts,nu_calc_pts), dtype = 'float')
            y_axis = np.zeros((nu_calc_pts,nu_calc_pts), dtype = 'float')

            if lay_interest == 0: z_range = -1*z_range
            else: z_range = np.abs(z_range)

            x1 = x_range
            y1 = y_range
            z1 = z_range
            (y_axis,x_axis) = np.meshgrid(y1,x1)

            E_TE_x_array = np.zeros((np.size(x1),np.size(y1),np.size(z1)), dtype = 'complex')
            E_TE_y_array = np.zeros((np.size(x1),np.size(y1),np.size(z1)), dtype = 'complex')
            E_TE_z_array = np.zeros((np.size(x1),np.size(y1),np.size(z1)), dtype = 'complex')
            E_TM_x_array = np.zeros((np.size(x1),np.size(y1),np.size(z1)), dtype = 'complex')
            E_TM_y_array = np.zeros((np.size(x1),np.size(y1),np.size(z1)), dtype = 'complex')
            E_TM_z_array = np.zeros((np.size(x1),np.size(y1),np.size(z1)), dtype = 'complex')

            for z in xrange(np.size(z1)):
                for y in xrange(np.size(y1)):
                    for x in xrange(np.size(x1)):
                        if pstack.layers[lay_interest].structure.height_nm == 'semi_inf':
                            expo_down = np.exp(1j*(alpha*x1[x]+beta*y1[y]-gamma*z1[z]))
                        elif z1[z] <= float(pstack.layers[lay_interest].structure.height_nm)/period:
                            expo_down = np.exp(1j*(alpha*x1[x]+beta*y1[y]-gamma*(z1[z]-float(pstack.layers[lay_interest].structure.height_nm)/period)))
                        else:
                            raise ValueError("fields_in_plane: z_value exceeds thickness of layer, reduce it. ")
                        expo_up = np.exp(1j*(alpha*x1[x]+beta*y1[y]+gamma*z1[z]))

                        E_TE_x = np.sum(eta_TE_x_down*expo_down + eta_TE_x_up*expo_up)
                        E_TE_y = np.sum(eta_TE_y_down*expo_down + eta_TE_y_up*expo_up)
                        E_TE_z = np.sum(eta_TE_z_down*expo_down + eta_TE_z_up*expo_up)
                        E_TM_x = np.sum(eta_TM_x_down*expo_down + eta_TM_x_up*expo_up)
                        E_TM_y = np.sum(eta_TM_y_down*expo_down + eta_TM_y_up*expo_up)
                        E_TM_z = np.sum(eta_TM_z_down*expo_down + eta_TM_z_up*expo_up)
                        E_TE_x_array[x,y,z] = E_TE_x
                        E_TE_y_array[x,y,z] = E_TE_y
                        E_TE_z_array[x,y,z] = E_TE_z
                        E_TM_x_array[x,y,z] = E_TM_x
                        E_TM_y_array[x,y,z] = E_TM_y
                        E_TM_z_array[x,y,z] = E_TM_z
            E_x_array = E_TE_x_array + E_TM_x_array
            E_y_array = E_TE_y_array + E_TM_y_array
            E_z_array = E_TE_z_array + E_TM_z_array
            E_tot_array = np.sqrt(E_x_array*np.conj(E_x_array) + E_y_array*np.conj(E_y_array) + E_z_array*np.conj(E_z_array))

            E_super = [E_x_array, E_y_array, E_z_array, E_tot_array]
            E_sup_labels = [' E_x', ' E_y', ' E_z', ' |E|',]

            # Save figures
            for re_im in ['real','imag']:
                for z_of_xy in xrange(np.size(z1)):
                    fig1 = plt.figure(num=None, figsize=(12,21), dpi=80, facecolor='w', edgecolor='k')
                    for i in range(len(E_super)):
                        ax1 = fig1.add_subplot(4,1,i+1)
                        if re_im == 'real':
                            E_slice = np.real(E_super[i][:,:,z_of_xy])
                        elif re_im == 'imag':
                            E_slice = np.imag(E_super[i][:,:,z_of_xy])
                        x_min = x1[0]
                        x_max = x1[-1]
                        y_min = y1[0]
                        y_max = y1[-1]
                        cmap = plt.get_cmap('jet')
                        CS = plt.contourf(x_axis,y_axis,E_slice,35,cmap=cmap)
                        plt.axis([x_min,x_max,y_min,y_max])
                        cbar = plt.colorbar()
                        cbar.ax.set_ylabel(re_im + E_sup_labels[i])
                        ax1.set_xlabel('x (d)')
                        ax1.set_ylabel('y (d)')
                        ax1.axis('scaled')
                        ax1.set_ylim((y_min,y_max))

                    plt.suptitle('%(name)s \n %(p)s_E_xy_slice, z = %(z_pos)s, heights = %(h)s \n \
                        $\lambda$ = %(wl)f nm, period = %(d)f, PW = %(pw)i,' % \
                        {'name' : name_lay, 'h':heights_list, 'p' : re_im, 'z_pos' : z1[z_of_xy],'wl' : wl, \
                        'd' : period, 'pw' : pw} + '\n'
                        + '# prop. ords = %(prop)s, # evan. ords = %(evan)s, n = %(n)s,k = %(k)s'\
                        % {'evan' : evan, 'prop' : prop, 'n' : n, 'k' : k[0]})
                    plt.savefig('%(dir_name)s/stack_%(stack_num)s_%(p)s_E_%(name)s_slice=%(z_pos)s_wl=%(wl)s.pdf'% \
                        {'dir_name' : dir_name, 'wl' : wl, 'p' : re_im, 'z_pos' : z1[z_of_xy], \
                        'name' : name_lay,'stack_num':stack_num})

                    fig.clf()

        stack_num += 1

            # ## Old fortran plane wave plotting
            # extra_name = EH_name + "Lay" + zeros_int_str(TF_lay)
            # select_h = 0.0
            # lat_vec = [[1.0, 0.0], [0.0, 1.0]]
            # bun = pstack.layers[TF_lay] # superstrate # substrate
            # n_eff_0 = bun.n()
            # neq_PW = bun.structure.num_pw_per_pol
            # bloch_vec = bun.k_pll_norm()
            # ordre_ls = bun.max_order_PWs
            # index_pw = bun.sort_order
            # index_pw_inv = np.zeros(shape=(np.shape(index_pw)),dtype='int')
            # for s in range(len(index_pw)):
            #     s2 = index_pw[s]
            #     index_pw_inv[s2] = s
            # index_pw_inv += 1 # convert to fortran indices (starting from 1)

            # # vec_coef sorted from top of stack, everything else is sorted from bottom
            # vec_index = len(pstack.layers) - TF_lay - 1
            # vec_coef_down = pstack.vec_coef_down[vec_index]
            # if TF_lay == 0: vec_coef_up = np.zeros(shape=(np.shape(vec_coef_down)),dtype='complex128')
            # else: vec_coef_up = pstack.vec_coef_up[vec_index]

            # EMUstack.gmsh_plot_pw(meat.n_msh_el, meat.n_msh_pts, nnodes, neq_PW,
            #     bloch_vec, meat.table_nod, meat.x_arr, lat_vec, wl_normed,
            #     n_eff_0, vec_coef_down, vec_coef_up,
            #     index_pw_inv, ordre_ls, select_h,
            #     gmsh_file_pos, meat.structure.plot_real, meat.structure.plot_imag,
            #     meat.structure.plot_abs, extra_name)

            # # # Semi-inf case
            # # vec_coef_up = meat.R12[:,0] # TE polarisation
            # # vec_coef_up = meat.R12[:,neq_PW] # TM polarisation

            # # vec_coef_down = np.zeros(shape=(np.shape(vec_coef_up)),dtype='complex128')
            # # vec_coef_down[neq_PW] = 1.0


# field interpolators for Nanostruct layers
def fields_interpolator_in_plane(pstack, lay_interest=1, z_value=0.1):
    """
    Returns linear interpolators in the x-y plane at chosen value of z for
    Re[Ex],Re[Ey],Re[Ez],Im[Ex],Im[Ey],Im[Ez],|E|.

    Requires matplotlib v1.3.0 or later.

        Args:
            pstack: Stack object (not a list!!!) containing data to plot.

        Keyword Args:
            lay_interest  (int): the index of the layer considered within \
                the stack.

            z_value  (float): distance in nm from bottom surface of layer \
                at which to calculate fields. If layer is semi-inf substrate \
                then z_value is distance from top of this layer (i.e. bottom \
                interface of stack).
    """
    from fortran import EMUstack

    num_lays = len(pstack.layers)
    # If selected z location is within a NanoStruct layer
    # plot fields in Bloch Mode Basis using fortran routine.
    if isinstance(pstack.layers[lay_interest], mode_calcs.Simmo):
        meat = pstack.layers[lay_interest]

        try:
            if not meat.structure.periodicity == '2D_array':
                raise ValueError
            eps = meat.n_effs**2
            h_normed = float(meat.structure.height_nm)/float(meat.structure.period)

            # fortran routine naturally measure z top down
            z_value = z_value/float(meat.structure.period)
            z_value = h_normed - z_value
            wl_normed = pstack.layers[lay_interest].wl_norm()

            nnodes = 6

            # vec_coef sorted from top of stack, everything else is sorted from bottom
            vec_index = num_lays - lay_interest - 1
            vec_coef = np.concatenate((pstack.vec_coef_down[vec_index], pstack.vec_coef_up[vec_index]))

            # piling up of all the bloch modes to get all the fields
            m_E = EMUstack.field_value_plane(meat.num_BMs,
                                             meat.n_msh_el, meat.n_msh_pts, nnodes,
                                             meat.structure.nb_typ_el,
                                             meat.table_nod, meat.type_el,
                                             eps, meat.x_arr, meat.k_z, meat.sol1,
                                             vec_coef, h_normed, z_value)

            # unrolling data for the interpolators
            table_nod = meat.table_nod.T
            x_arr = meat.x_arr.T

            # dense triangulation with multiple points
            v_x6p = np.zeros(6*meat.n_msh_el)
            v_y6p = np.zeros(6*meat.n_msh_el)
            v_Ex6p = np.zeros(6*meat.n_msh_el, dtype=np.complex128)
            v_Ey6p = np.zeros(6*meat.n_msh_el, dtype=np.complex128)
            v_Ez6p = np.zeros(6*meat.n_msh_el, dtype=np.complex128)
            v_triang6p = []

            i = 0
            for i_el in np.arange(meat.n_msh_el):

                # triangles
                idx = np.arange(6*i_el, 6*(i_el+1))
                triangles = [[idx[0], idx[3], idx[5]],
                             [idx[1], idx[4], idx[3]],
                             [idx[2], idx[5], idx[4]],
                             [idx[3], idx[4], idx[5]]]
                v_triang6p.extend(triangles)

                for i_node in np.arange(6):

                    # index for the coordinates
                    i_ex = table_nod[i_el, i_node]-1

                    # values
                    v_x6p[i] = x_arr[i_ex, 0]
                    v_y6p[i] = x_arr[i_ex, 1]
                    v_Ex6p[i] = m_E[i_el, i_node, 0]
                    v_Ey6p[i] = m_E[i_el, i_node, 1]
                    v_Ez6p[i] = m_E[i_el, i_node, 2]

                    i += 1

            v_E6p = np.sqrt(np.abs(v_Ex6p)**2 +
                            np.abs(v_Ey6p)**2 +
                            np.abs(v_Ez6p)**2)

            # dense triangulation with unique points
            v_triang1p = []
            for i_el in np.arange(meat.n_msh_el):

                # triangles
                triangles = [[table_nod[i_el,0]-1,table_nod[i_el,3]-1,table_nod[i_el,5]-1],
                             [table_nod[i_el,1]-1,table_nod[i_el,4]-1,table_nod[i_el,3]-1],
                             [table_nod[i_el,2]-1,table_nod[i_el,5]-1,table_nod[i_el,4]-1],
                             [table_nod[i_el,3]-1,table_nod[i_el,4]-1,table_nod[i_el,5]-1]]
                v_triang1p.extend(triangles)

            # triangulations
            triang6p = matplotlib.tri.Triangulation(v_x6p,v_y6p,v_triang6p)
            triang1p = matplotlib.tri.Triangulation(x_arr[:,0],x_arr[:,1],v_triang1p)

            # building interpolators: triang1p for the finder, triang6p for the values
            finder = matplotlib.tri.TrapezoidMapTriFinder(triang1p)
            ReEx = matplotlib.tri.LinearTriInterpolator(triang6p,v_Ex6p.real,trifinder=finder)
            ImEx = matplotlib.tri.LinearTriInterpolator(triang6p,v_Ex6p.imag,trifinder=finder)
            ReEy = matplotlib.tri.LinearTriInterpolator(triang6p,v_Ey6p.real,trifinder=finder)
            ImEy = matplotlib.tri.LinearTriInterpolator(triang6p,v_Ey6p.imag,trifinder=finder)
            ReEz = matplotlib.tri.LinearTriInterpolator(triang6p,v_Ez6p.real,trifinder=finder)
            ImEz = matplotlib.tri.LinearTriInterpolator(triang6p,v_Ez6p.imag,trifinder=finder)
            AbsE = matplotlib.tri.LinearTriInterpolator(triang6p,v_E6p,trifinder=finder)

            return ReEx, ImEx, ReEy, ImEy, ReEz, ImEz, AbsE

        except ValueError as e:
            print(e)
            print("fields_in_plane cannot plot fields in 1D-arrays."\
                  "\nPlease select a different lay_interest.\n")

    # If selected z location is within a ThinFilm layer
    # plot fields in Plane Wave Basis using python routine.
    else:
        raise Exception("Not a NanoStruct Layer")


def fields_vertically(stacks_list, factor_pts_vert=31, nu_pts_hori=41,
                      semi_inf_height=1.0, gradient=None, no_incoming=False,
                      add_name='', force_eq_ratio=False, colour_res=30,
                      plot_boundaries=True):
    """
    Plot fields in the x-y plane at chosen values of z, where z is \
    calculated from the bottom of chosen layer.

        Args:
            stacks_list  (list): Stack objects containing data to plot.

        Keyword Args:
            factor_pts_vert  (int): sampling factor for fields vertically. \
            Calculated as factor_pts_vert * (epsilon*h/wl).

            nu_pts_hori  (int): in-plane fields are calculated over a mesh of \
                nu_pts_hori * nu_pts_hori points.

            semi_inf_height  (float): distance to which fields are plotting \
                in semi-infinite (sub)superstrates.

            gradient  (float): further slices calculated with given gradient \
                and -gradient. It is entitled 'specified_diagonal_slice'.\
                These slices are only calculated for ThinFilm layers.

            no_incoming  (bool): if True, plots fields in superstrate in the \
                absence of the incident driving field (i.e. only showing \
                upward propagating scattered field).

            add_name  (str): concatenate add_name to title.

            force_eq_ratio  (bool): each layer plotted on equal space.

            colour_res  (int): number of colour intervals to use.

            plot_boundaries  (bool): plot boundaries of inclusions.
    """

    from fortran import EMUstack

    dir_name = "fields_vertically"
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)

    # always make odd
    if nu_pts_hori % 2 == 0: nu_pts_hori += 1

    stack_num = 0
    for pstack in stacks_list:
        num_lays  = len(pstack.layers)
        period    = pstack.layers[0].structure.period
        pw        = pstack.layers[0].max_order_PWs
        wl_normed = pstack.layers[0].wl_norm()

        # Stack layers on top of one another to find heights of interfaces.
        ind_h_list = pstack.heights_norm()
        h_list = [-semi_inf_height, 0]
        for i in range(len(ind_h_list)):
            h_list.append(h_list[-1]+ind_h_list[i])
        h_list.append(h_list[-1]+semi_inf_height)
        ind_h_list.append(semi_inf_height)
        ind_h_list.insert(0, semi_inf_height)
        min_h = np.min(ind_h_list)
        h_ratio = [h/min_h for h in ind_h_list[::-1]]

        layer = pstack.layers[0]
        if isinstance(layer,mode_calcs.Anallo):
            if layer.structure.world_1d == True:
                slice_types = ['xz']
            else:
                slice_types = ['xz']#,'yz','diag+','diag-'] #rest not implemented
                if gradient != None: slice_types.append(('special+','special-'))
        else:
            if layer.structure.periodicity == '1D_array':
                slice_types = ['xz']
            else:
                slice_types = ['xz']#,'yz','diag+','diag-'] #rest not implemented
                if gradient != None: slice_types.append(('special+','special-'))

        for sli in slice_types:
            E_fields = ['Re(E_x)', 'Im(E_x)', 'Re(E_y)', 'Im(E_y)', 'Re(E_z)',
                        'Im(E_z)', 'Re(E)',  'eps_abs(E)']
            E_fields_tot = []
            epsE_fields_tot = []
            for E in E_fields:
                fig = plt.figure(figsize=(9, 12))
                if force_eq_ratio == True:
                    h_ratio = [1 for h in ind_h_list[::-1]]
                    gs = gridspec.GridSpec(num_lays, 1, wspace=0.0, hspace=0.0,
                                           width_ratios=[1, 1],
                                           height_ratios=h_ratio)
                else:
                    gs = gridspec.GridSpec(num_lays, 1, wspace=0.0, hspace=0.0,
                                           width_ratios=[1, 1],
                                           height_ratios=h_ratio)
                # Find maximum field value to normalise plots.
                max_E = None
                min_E = None
                for max_E_field in [0, 1]:
                    # Iterate through layers
                    for lay in range(num_lays):
                        layer = pstack.layers[lay]
                        ax1 = plt.subplot(gs[num_lays-lay-1])
                        if isinstance(layer, mode_calcs.Simmo):
                            eps = layer.n_effs**2
                        else:
                            eps = layer.n()**2
                        nu_pts_vert = np.round(factor_pts_vert*(np.real(np.max(eps))*ind_h_list[lay] / wl_normed))
                        if nu_pts_vert < 11: nu_pts_vert = 11
                        if lay == 0:
                            z_range = np.linspace(h_list[lay],0.0,nu_pts_vert)
                        else:
                            z_range = np.linspace(0.0, ind_h_list[lay], nu_pts_vert)
                        # vec_coef sorted from top, everything else sorted from bottom
                        vec_index = num_lays - lay - 1

                        # If NanoStruct layer
                        if isinstance(layer, mode_calcs.Simmo):
                            name_lay = "%i_NanoStruct"% lay
                            gmsh_file_pos = 'stack_%(stack_num)s_lay_%(name)s_'% \
                                            {'name': name_lay, 'stack_num': stack_num}

                            vec_coef_fem = np.concatenate((pstack.vec_coef_down[vec_index],
                                                          pstack.vec_coef_up[vec_index]))

                            # If 2D_array plot fields using fortran routine.
                            if layer.structure.periodicity == '2D_array':
                                if max_E_field == 0 and E == 'Re(E_x)':
                                    E_fields_tot.append(None)
                                    epsE_fields_tot.append(None)
                                elif max_E_field == 1:
                                    if not os.path.exists(dir_name+"/gmsh_BMs"):
                                        os.mkdir(dir_name+"/gmsh_BMs")
                                    if not os.path.exists(dir_name+"/gmsh_BMs/anim"):
                                        os.mkdir(dir_name+"/gmsh_BMs/anim")

                                    # TODO fix the scaling issue in 2D FEM plotting!
                                    # Then remove these arbitrary scaling params.
                                    scale_plot = 2.0
                                    shift_x_plot = -.5
                                    h_normed = float(layer.structure.height_nm)/float(layer.structure.period)
                                    shift_v_plot = h_normed*0.75
                                    nnodes = 6
                                    EMUstack.gmsh_plot_slice(layer.E_H_field, layer.num_BMs,
                                        layer.n_msh_el, layer.n_msh_pts, nnodes, layer.type_el,
                                        layer.structure.nb_typ_el, eps, layer.table_nod,
                                        layer.x_arr, layer.k_z, layer.sol1, vec_coef_fem,
                                        h_normed, wl_normed, gmsh_file_pos,
                                        scale_plot, shift_v_plot, shift_x_plot)

                                    z_plot = np.linspace(h_list[lay], h_list[lay+1],2)
                                    (y_axis_plot, x_axis) = np.meshgrid(z_plot,[0,1])

                                    CS = plt.contourf(x_axis, y_axis_plot,[[0,0],[0,0]],
                                        cmap=plt.cm.binary, vmin=0, vmax=1)
                                    CS = plt.contourf(x_axis, y_axis_plot,[[0,0],[0,0]],
                                        cmap=plt.cm.binary, vmin=0, vmax=1)
                                    CS.set_clim(0, 1)
                                    if abs(ind_h_list[lay]) < 0.05 * np.sum(ind_h_list):
                                        tick_half = (h_list[lay]+h_list[lay+1])/2.
                                        ax1.yaxis.set_ticks([tick_half])
                                        ticktitles = ['%3.2f' % tick_half]
                                        ax1.yaxis.set_ticklabels(ticktitles)
                                    elif abs(ind_h_list[lay]) < 0.25 * np.sum(ind_h_list):
                                        tick_half = (h_list[lay]+h_list[lay+1])/2.
                                        tick_int = ind_h_list[lay]/4.
                                        ticks = [tick_half-tick_int, tick_half, tick_half+tick_int]
                                        ax1.yaxis.set_ticks(ticks)
                                        ticktitles = []
                                        [ticktitles.append('%3.2f' % tick) for tick in ticks]
                                        ax1.yaxis.set_ticklabels(ticktitles)

                            # If 1D_array plot fields as follows.
                            else:
                                struct = layer.structure
                                E_slice = np.zeros((struct.n_msh_pts,nu_pts_vert), dtype = 'complex128')

                                boundary = []
                                for i in range(len(struct.type_el) - 1):
                                    if struct.type_el[i] != struct.type_el[i+1]:
                                        boundary.append(struct.x_arr[2*(i+1)])

                                if E[-3:] != '(E)':
                                    if E[5] == 'x': comp = 0
                                    if E[5] == 'y': comp = 1
                                    if E[5] == 'z': comp = 2
                                    ### sol1 = sol_P2([Ex,Ey,Ez],P2_interpolation_points,nval,nel)
                                    for BM in range(layer.num_BMs):
                                        BM_sol = layer.sol1[comp,:,BM,:]
                                        beta = layer.k_z[BM]
                                        for h in range(len(z_range)):
                                            hz = z_range[h]
                                            P_down = np.exp(1j*beta*(ind_h_list[lay]-hz))   # Introduce Propagation in -z
                                            P_up = np.exp(1j*beta*hz) # Introduce Propagation in +z

                                            coef_down = vec_coef_fem[BM] * P_down
                                            coef_up = vec_coef_fem[BM+layer.num_BMs] * P_up
                                            if E[5] == 'z':
                                                coef_tot = (coef_up - coef_down)/beta # Taking into account the change of variable for Ez
                                            else:
                                                coef_tot = coef_up + coef_down
                                            coef_tot = coef_tot[0,0]
                                            E_slice[0,h] += BM_sol[0,0] * coef_tot
                                            for x in range(struct.n_msh_el - 1):
                                                E_slice[2*x+1,h] += BM_sol[1,x] * coef_tot
                                                # E_slice[2*x+2,h] += (BM_sol[2,x] + BM_sol[0,x+1] / 2.) * coef_tot
                                                E_slice[2*x+2,h] += (BM_sol[2,x]) * coef_tot
                                            E_slice[2*x+3,h] += BM_sol[1,-1] * coef_tot
                                            E_slice[2*x+4,h] += BM_sol[2,-1] * coef_tot

                                    if max_E_field == 0 and E[0] == 'R':
                                        if E[5] == 'x':
                                            E_fields_tot.append(E_slice*np.conj(E_slice))
                                            epsE_fields_tot.append(np.zeros(np.shape(E_slice)))
                                        elif E[5] == 'y' or E[5] == 'z':
                                            E_fields_tot[lay] += E_slice*np.conj(E_slice)
                                elif E == 'eps_abs(E)':
                                    if max_E_field == 0:
                                        type_el = np.vstack((struct.type_el, struct.type_el)).reshape((-1,),order='F')
                                        # type_el = np.append(type_el, type_el[-1])
                                        type_el = np.append(type_el[-1], type_el)
                                        type_el[type_el == 1] = np.abs(eps[0])
                                        type_el[type_el == 2] = np.abs(eps[1])
                                        type_el = np.diag(type_el)
                                        epsE_fields_tot[lay] = np.dot(type_el, (E_fields_tot[lay]))
                                    else:
                                        E_slice = epsE_fields_tot[lay]
                                elif E[-3:] == '(E)' and max_E_field == 1:
                                    E_slice = np.sqrt(E_fields_tot[lay])

                                if E[0] == 'R' or E[0] == 'e':
                                    E_slice = np.real(E_slice)
                                elif E[0] == 'I':
                                    E_slice = np.imag(E_slice)

                                if E[-3:] != '(E)':  # Re(x,y,z), Im(x,y,z)
                                    max_E_lay = np.max(E_slice)
                                    if max_E_lay > max_E or max_E == None:
                                        max_E = max_E_lay
                                    min_E_lay = np.min(E_slice)
                                    if min_E_lay < min_E or min_E == None:
                                        min_E = min_E_lay
                                elif E == 'Re(E)':
                                    max_E_lay = np.max(np.real(np.sqrt(E_fields_tot[lay])))
                                    if max_E_lay > max_E or max_E == None:
                                        max_E = max_E_lay
                                    min_E_lay = np.min(np.real(np.sqrt(E_fields_tot[lay])))
                                    if min_E_lay < min_E or min_E == None:
                                        min_E = min_E_lay
                                else:
                                    max_E_lay = np.max(np.abs(epsE_fields_tot[lay]))
                                    if max_E_lay > max_E or max_E == None:
                                        max_E = max_E_lay
                                    min_E_lay = np.min(np.abs(epsE_fields_tot[lay]))
                                    if min_E_lay < min_E or min_E == None:
                                        min_E = min_E_lay

                                if lay == 0:
                                    z_plot = np.linspace(h_list[lay], 0.0, nu_pts_vert)
                                else:
                                    z_plot = np.linspace(h_list[lay], h_list[lay+1],nu_pts_vert)
                                (y_axis_plot, x_axis) = np.meshgrid(z_plot,struct.x_arr)

                                if max_E_field == 1:
                                    if E == 'eps_abs(E)' or E == 'Re(E)':
                                        choice_cmap = plt.cm.hot
                                    else:
                                        choice_cmap = plt.cm.jet
                                    CS = plt.contourf(x_axis, y_axis_plot,
                                                      E_slice, colour_res,
                                                      cmap=choice_cmap,
                                                      vmin=min_E, vmax=max_E)
                                    CS = plt.contourf(x_axis, y_axis_plot,
                                                      E_slice, colour_res,
                                                      cmap=choice_cmap,
                                                      vmin=min_E, vmax=max_E)
                                    CS.set_clim(min_E, max_E)
                                    ax1.set_xlim((x_range[0], x_range[-1]))
                                    if abs(ind_h_list[lay]) < 0.05 * np.sum(ind_h_list):
                                        tick_half = (h_list[lay]+h_list[lay+1])/2.
                                        ax1.yaxis.set_ticks([tick_half])
                                        ticktitles = ['%3.2f' % tick_half]
                                        ax1.yaxis.set_ticklabels(ticktitles)
                                    elif abs(ind_h_list[lay]) < 0.25 * np.sum(ind_h_list):
                                        tick_half = (h_list[lay]+h_list[lay+1])/2.
                                        tick_int = ind_h_list[lay]/4.
                                        ticks = [tick_half-tick_int, tick_half, tick_half+tick_int]
                                        ax1.yaxis.set_ticks(ticks)
                                        ticktitles = []
                                        [ticktitles.append('%3.2f' % tick) for tick in ticks]
                                        ax1.yaxis.set_ticklabels(ticktitles)

                                    if plot_boundaries == True:
                                        start, end = ax1.get_ylim()
                                        for b in boundary:
                                            ax1.plot([b, b], [start, end], 'w',
                                                     linewidth=2)

                                # else:
                                # scale_plot = 2.0
                                # shift_x_plot = -.5
                                # shift_v_plot = h_normed*0.75
                                # EMUstack.gmsh_plot_slice_1d(layer.E_H_field, layer.num_BMs,
                                #     struct.n_msh_el, struct.n_msh_pts, struct.type_el,
                                #     struct.nb_typ_el, n_eff, struct.table_nod,
                                #     struct.x_arr, layer.k_z, layer.sol1, vec_coef,
                                #     h_normed, wl_normed, gmsh_file_pos,
                                #     scale_plot, shift_v_plot, shift_x_plot)

                        # ThinFilm layer, using Plane Wave basis.
                        else:
                            x_range = np.linspace(0.0, 1.0, nu_pts_hori)
                            y_range = np.linspace(0.0, 1.0, nu_pts_hori)

                            s = layer.sort_order
                            alpha_unsrt = np.array(layer.alphas)
                            beta_unsrt = np.array(layer.betas)
                            alpha = alpha_unsrt[s]
                            if layer.structure.world_1d is True:
                                beta = beta_unsrt
                            else:
                                beta = beta_unsrt[s]
                            gamma = np.array(layer.calc_kz())
                            n = layer.n()
                            num_pws = layer.structure.num_pw_per_pol

                            if lay == num_lays-1 and no_incoming == True:
                                vec_coef_down = np.zeros((2*num_pws), dtype='complex')
                            else:
                                vec_coef_down = np.array(pstack.vec_coef_down[vec_index]).flatten()
                            vec_coef_down_TE = vec_coef_down[0:num_pws]
                            vec_coef_down_TM = vec_coef_down[num_pws::]

                            if lay == 0:
                                vec_coef_up = np.zeros((2*num_pws), dtype='complex')
                            else:
                                vec_coef_up = np.array(pstack.vec_coef_up[vec_index]).flatten()
                            vec_coef_up_TE = vec_coef_up[0:num_pws]
                            vec_coef_up_TM = vec_coef_up[num_pws::]

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

                            if E[-3:] != '(E)':
                                if E[5] == 'x':
                                    TE_coef = E_TE_x
                                    TM_coef = E_TM_x
                                elif E[5] == 'y':
                                    TE_coef = E_TE_y
                                    TM_coef = E_TM_y
                                elif E[5] == 'z':
                                    TE_coef = E_TE_z
                                    TM_coef = E_TM_z

                            eta_TE_down = (vec_coef_down_TE*TE_coef)/chi_TE
                            eta_TM_down = (vec_coef_down_TM*TM_coef)/chi_TM
                            eta_TE_up = (vec_coef_up_TE*TE_coef)/chi_TE
                            eta_TM_up = (vec_coef_up_TM*TM_coef)/chi_TM

                            x_axis = np.zeros((nu_pts_hori, nu_pts_hori), dtype='float')
                            y_axis = np.zeros((nu_pts_hori, nu_pts_hori), dtype='float')

                            if sli == 'xz':
                                x1 = x_range
                                if layer.structure.world_1d == True:
                                    y1 = np.array([0])
                                else:
                                    y1 = np.array([0, 0.5])
                                (y_axis, x_axis) = np.meshgrid(z_range, x1)
                            elif sli == 'yz':
                                x1 = np.array([0, 0.5])
                                y1 = y_range
                                (y_axis, x_axis) = np.meshgrid(z_range, y1)
                            elif sli == 'diag+':
                                x1 = x_range
                                y1 = np.array([0])
                                y2 = x_range
                                (y_axis, x_axis) = np.meshgrid(z_range, sqrt(2)*x1)
                            elif sli == 'diag-':
                                x1 = x_range[::-1]
                                y1 = np.array([0])
                                y2 = x_range
                                (y_axis, x_axis) = np.meshgrid(z_range, sqrt(2)*x1[::-1])
                            elif sli == 'special+':
                                x1 = x_range
                                y1 = np.array([0])
                                y2 = gradient*x_range
                                for i in xrange(np.size(y2)):
                                    if y2[i] > 1:
                                        y2 = np.resize(y2, (i,))
                                        x1 = np.resize(x1, (i,))
                                        break
                                (y_axis, x_axis) = np.meshgrid(z_range,sqrt(1+gradient**2)*x1)
                            elif sli == 'special-':
                                x1 = x_range[::-1]
                                y1 = np.array([0])
                                y2 = gradient*x_range
                                for i in xrange(np.size(y2)):
                                    if y2[i] > 1:
                                        y2 = np.resize(y2,(i,))
                                        x1 = np.resize(x1,(i,))
                                        break
                                (y_axis, x_axis) = np.meshgrid(z_range, sqrt(1+gradient**2)*(x1[::-1]-x1[-1]))

                            E_field = np.zeros((np.size(x1), np.size(y1),
                                               np.size(z_range)), dtype='complex')

                            for z in xrange(np.size(z_range)):
                                for y in xrange(np.size(y1)):
                                    for x in xrange(np.size(x1)):
                                        if sli == 'diag+' or sli == 'diag-' or sli == 'special+' or sli == 'special-':
                                            if layer.structure.height_nm == 'semi_inf':
                                                expo_down = np.exp(1j*(alpha*x1[x]+beta*y2[x]-gamma*z_range[z]))
                                            else:
                                                expo_down = np.exp(1j*(alpha*x1[x]+beta*y2[x]-gamma*(z_range[z]-ind_h_list[lay])))
                                            expo_up = np.exp(1j*(alpha*x1[x]+beta*y2[x]+gamma*z_range[z]))
                                        else:
                                            if lay == 0 or lay == num_lays-1:
                                                expo_down = np.exp(1j*(alpha*x1[x]+beta*y1[y]-gamma*z_range[z]))
                                            else:
                                                expo_down = np.exp(1j*(alpha*x1[x]+beta*y1[y]-gamma*(z_range[z]-ind_h_list[lay])))
                                            expo_up = np.exp(1j*(alpha*x1[x]+beta*y1[y]+gamma*z_range[z]))

                                        E_field[x,y,z] = np.sum(eta_TE_down*expo_down + eta_TE_up*expo_up) + np.sum(eta_TM_down*expo_down + eta_TM_up*expo_up)

                            if lay == 0:
                                z_plot = np.linspace(h_list[lay],0.0,nu_pts_vert)
                            else:
                                z_plot = np.linspace(h_list[lay],h_list[lay+1],nu_pts_vert)
                            (y_axis_plot, x_axis) = np.meshgrid(z_plot, x_range)

                            if sli == 'xz':
                                E_slice = E_field[:, 0, :]
                                if max_E_field == 0:
                                    if E[-3:] != '(E)' and E[0] == 'R':
                                        if E[5] == 'x':
                                            E_fields_tot.append(E_slice*np.conj(E_slice))
                                            epsE_fields_tot.append(np.zeros(np.shape(E_slice)))
                                        elif E[5] == 'y' or E[5] == 'z':
                                            E_fields_tot[lay] += E_slice*np.conj(E_slice)
                                    elif E == 'eps_abs(E)':
                                        epsE_fields_tot[lay] = np.real(eps) * (E_fields_tot[lay])
                                elif max_E_field == 1 and E == 'eps_abs(E)':
                                    E_slice = epsE_fields_tot[lay]
                                elif max_E_field == 1 and E[-3:] == '(E)':
                                    E_slice = np.sqrt(E_fields_tot[lay])

                                if E[0] == 'R' or E[0] == 'e':
                                    E_slice = np.real(E_slice)
                                elif E[0] == 'I':
                                    E_slice = np.imag(E_slice)

                                if E[-3:] != '(E)':
                                    max_E_lay = np.max(E_slice)
                                    if max_E_lay > max_E or max_E == None:
                                        max_E = max_E_lay
                                    min_E_lay = np.min(E_slice)
                                    if min_E_lay < min_E or min_E == None:
                                        min_E = min_E_lay
                                elif E == 'Re(E)':
                                    max_E_lay = np.max(np.real(np.sqrt(E_fields_tot[lay])))
                                    if max_E_lay > max_E or max_E == None:
                                        max_E = max_E_lay
                                    min_E_lay = np.min(np.real(np.sqrt(E_fields_tot[lay])))
                                    if min_E_lay < min_E or min_E == None:
                                        min_E = min_E_lay
                                else:
                                    max_E_lay = np.max(np.real(epsE_fields_tot[lay]))
                                    if max_E_lay > max_E or max_E == None:
                                        max_E = max_E_lay
                                    min_E_lay = np.min(np.real(epsE_fields_tot[lay]))
                                    if min_E_lay < min_E or min_E == None:
                                        min_E = min_E_lay

                                if max_E_field == 1:
                                    if E == 'eps_abs(E)' or E == 'Re(E)':
                                        choice_cmap = plt.cm.hot
                                    else:
                                        choice_cmap = plt.cm.jet
                                    CS = plt.contourf(x_axis,y_axis_plot,E_slice,
                                        colour_res, cmap=choice_cmap, vmin=min_E, vmax=max_E)
                                    # hack to remove white contour lines
                                    CS = plt.contourf(x_axis,y_axis_plot,E_slice,
                                        colour_res, cmap=choice_cmap, vmin=min_E, vmax=max_E)
                                    CS.set_clim(min_E,max_E)
                                    ax1.set_xlim((x_range[0],x_range[-1]))
                                    if abs(ind_h_list[lay]) < 0.05 * np.sum(ind_h_list):
                                        tick_half = (h_list[lay]+h_list[lay+1])/2.
                                        ax1.yaxis.set_ticks([tick_half])
                                        ticktitles = ['%3.2f' % tick_half]
                                        ax1.yaxis.set_ticklabels(ticktitles)
                                    elif abs(ind_h_list[lay]) < 0.25 * np.sum(ind_h_list):
                                        tick_half = (h_list[lay]+h_list[lay+1])/2.
                                        tick_int = ind_h_list[lay]/4.
                                        ticks = [tick_half-tick_int, tick_half, tick_half+tick_int]
                                        ax1.yaxis.set_ticks(ticks)
                                        ticktitles = []
                                        [ticktitles.append('%3.2f' % tick) for tick in ticks]
                                        ax1.yaxis.set_ticklabels(ticktitles)

                            # elif sli == 'yz':
                            #     E_slice = E_field[0, :, :]
                            #     if max_E_field == 0:
                            #         if E[-3:] != '(E)' and E[0] == 'R':
                            #             if E[5] == 'x':
                            #                 E_fields_tot.append(E_slice*np.conj(E_slice))
                            #                 epsE_fields_tot.append(np.zeros(np.shape(E_slice)))
                            #             elif E[5] == 'y' or E[5] == 'z':
                            #                 E_fields_tot[lay] += E_slice*np.conj(E_slice)
                            #         elif E == 'eps_abs(E)':
                            #             epsE_fields_tot[lay] = np.real(eps) * (E_fields_tot[lay])
                            #     elif max_E_field == 1 and E == 'eps_abs(E)':
                            #         E_slice = epsE_fields_tot[lay]
                            #     elif max_E_field == 1 and E[-3:] == '(E)':
                            #         E_slice = np.sqrt(E_fields_tot[lay])

                            #     if E[0] == 'R' or E[0] == 'e':
                            #         E_slice = np.real(E_slice)
                            #     elif E[0] == 'I':
                            #         E_slice = np.imag(E_slice)

                            #     if E[-3:] != '(E)':
                            #         max_E_lay = np.max(E_slice)
                            #         if max_E_lay > max_E or max_E == None:
                            #             max_E = max_E_lay
                            #         min_E_lay = np.min(E_slice)
                            #         if min_E_lay < min_E or min_E == None:
                            #             min_E = min_E_lay
                            #     elif E == 'Re(E)':
                            #         max_E_lay = np.max(np.real(np.sqrt(E_fields_tot[lay])))
                            #         if max_E_lay > max_E or max_E == None:
                            #             max_E = max_E_lay
                            #         min_E_lay = np.min(np.real(np.sqrt(E_fields_tot[lay])))
                            #         if min_E_lay < min_E or min_E == None:
                            #             min_E = min_E_lay
                            #     else:
                            #         max_E_lay = np.max(np.real(epsE_fields_tot[lay]))
                            #         if max_E_lay > max_E or max_E == None:
                            #             max_E = max_E_lay
                            #         min_E_lay = np.min(np.real(epsE_fields_tot[lay]))
                            #         if min_E_lay < min_E or min_E == None:
                            #             min_E = min_E_lay

                            #     if max_E_field == 1:
                            #         if E == 'eps_abs(E)' or E == 'Re(E)':
                            #             choice_cmap = plt.cm.hot
                            #         else:
                            #             choice_cmap = plt.cm.jet
                            #         CS = plt.contourf(x_axis, y_axis_plot, E_slice,
                            #             colour_res, cmap=choice_cmap, vmin=min_E, vmax=max_E)
                            #         CS.set_clim(min_E, max_E)
                            #         ax1.set_xlim((x_range[0], x_range[-1]))
                            #         if abs(ind_h_list[lay]) < 0.05 * np.sum(ind_h_list):
                            #             tick_half = (h_list[lay]+h_list[lay+1])/2.
                            #             ax1.yaxis.set_ticks([tick_half])
                            #             ticktitles = ['%3.2f' % tick_half]
                            #             ax1.yaxis.set_ticklabels(ticktitles)
                            #         elif abs(ind_h_list[lay]) < 0.25 * np.sum(ind_h_list):
                            #             tick_half = (h_list[lay]+h_list[lay+1])/2.
                            #             tick_int = ind_h_list[lay]/4.
                            #             ticks = [tick_half-tick_int, tick_half, tick_half+tick_int]
                            #             ax1.yaxis.set_ticks(ticks)
                            #             ticktitles = []
                            #             [ticktitles.append('%3.2f' % tick) for tick in ticks]
                            #             ax1.yaxis.set_ticklabels(ticktitles)

                        if lay == 0 and sli == 'xz':
                            ax1.set_xlabel('x (d)')
                        elif lay == 0 and sli == 'yz':
                            ax1.set_xlabel('y (d)')
                        else:
                            ax1.set_xticklabels(())

                cax = fig.add_axes([0.6, 0.1, 0.03, 0.8])
                cb = fig.colorbar(plt.contourf([0, 1], [0, 1], [[min_E, min_E],
                                  [max_E, max_E]], colour_res,
                                  cmap=choice_cmap, vmin=min_E, vmax=max_E),
                                  cax=cax)
                cb.set_clim(min_E, max_E)
                if E == 'eps_abs(E)':
                    cax.set_ylabel(r'$|\epsilon|$ |E|$^2$')
                elif E == 'Re(E)':
                    cax.set_ylabel(r'|E|')
                else:
                    cax.set_ylabel(E)

                plt.savefig('%(dir_name)s/stack_%(stack_num)s_E_%(slice)s_slice_%(comp)s%(add)s.pdf'% \
                    {'dir_name': dir_name, 'comp': E, 'slice': sli,\
                     'stack_num': stack_num, 'add': add_name},
                        bbox_inches='tight')

                fig.clf()

                            # elif sli == 'yz':
                            #     for x_of_yz in xrange(np.size(x1)):
                            #         fig1 = plt.figure(num=None, figsize=(12, 21), dpi=80, facecolor='w', edgecolor='k')
                            #         for i in range(len(E_super)):
                            #             ax1 = fig1.add_subplot(4, 1, i+1)
                            #             if re_im == 'real':
                            #                 E_slice = np.real(E_super[i][x_of_yz,:,:])
                            #             elif re_im == 'imag':
                            #                 E_slice = np.imag(E_super[i][x_of_yz,:,:])
                            #             y_min = y_range[0]
                            #             y_max = y_range[-1]
                            #             z_min = z_range[0]
                            #             z_max = z_range[-1]
                            #             cmap = plt.get_cmap('jet')
                            #             CS = plt.contourf(x_axis,y_axis,E_slice,15,cmap=cmap)
                            #             plt.axis([y_min,y_max,z_min,z_max])
                            #             cbar = plt.colorbar()
                            #             cbar.ax.set_ylabel(r'%(p)s E_x'%{'p':p})
                            #             ax1.set_xlabel('y (d)')
                            #             ax1.set_ylabel('z (d)')
                            #             if scale_axis == True: ax1.axis('scaled')
                            #             ax1.xaxis.set_ticks([x_min,x_max])
                            #             ax1.set_ylim((z_min,z_max))
                            #             if np.abs(z_max-z_min) < x_max: ax1.yaxis.set_ticks([z_min,z_max])

                                    # plt.suptitle('%(name)s \n E_yz_slice_%(p)s, x = %(x_pos)s, heights = %(h)s \n \
                                    #     $\lambda$ = %(wl)snm, period = %(d)f, PW = %(pw)i, %(add)s' % \
                                    #     {'name' : name_lay, 'h':heights_list, 'p' : p, 'x_pos' : x1[x_of_yz],'wl' : wl, \
                                    #     'd' : period, 'pw' : pw, 'add' : add_name} + '\n'
                                    #     + '# prop. ords = %(prop)s, # evan. ords = %(evan)s , \
                                    #     n = %(n)s, k = %(k)s' % {'evan' : evan, 'prop' : prop, 'n' : n, 'k' : k[0]})
                                    # plt.savefig('%(dir_name)s/stack_%(stack_num)s_lay_%(name)s_E_yz_slice=%(x_pos)s_wl=%(wl)s_%(p)s%(add)s.pdf'% \
                                    #     {'dir_name' : dir_name, 'p':p, 'wl' : wl, 'x_pos' : x1[x_of_yz],\
                                    #     'name' : name_lay,'stack_num':stack_num, 'add' : add_name})

                            # elif sli == 'diag+':
                            #     diag = 1
                            #     fig1 = plt.figure(num=None, figsize=(12,21), dpi=80, facecolor='w', edgecolor='k')
                            #     for i in range(len(E_super)):
                            #         ax1 = fig1.add_subplot(4,1,i+1)
                            #         if re_im == 'real':
                            #             E_slice = np.real(E_super[i][:,0,:])
                            #         elif re_im == 'imag':
                            #             E_slice = np.imag(E_super[i][:,0,:])
                            #         y_min = y_range[0]
                            #         y_max = np.around(np.sqrt(2)*y_range[-1],decimals=2)
                            #         z_min = z1[0]
                            #         z_max = z1[-1]
                            #         cmap = plt.get_cmap('jet')
                            #         CS = plt.contourf(x_axis,y_axis,E_slice,15,cmap=cmap)
                            #         plt.axis([y_min,y_max,z_min,z_max])
                            #         cbar = plt.colorbar()
                            #         cbar.ax.set_ylabel(r'%(p)s E_x'%{'p':p})
                            #         ax1.set_xlabel('x=%(diag)sy (d)'%{'diag':diag})
                            #         ax1.set_ylabel('z (d)')
                            #         if scale_axis == True: ax1.axis('scaled')
                            #         ax1.xaxis.set_ticks([y_min,y_max])
                            #         ax1.set_xlim((y_min,y_max))
                            #         ax1.set_ylim((z_min,z_max))
                            #         if np.abs(z_max-z_min) < x_max: ax1.yaxis.set_ticks([z_min,z_max])

                            #     # plt.suptitle('%(name)s \n E_diagonal_slice_%(p)s, y = %(diag)sx, heights = %(h)s \n\
                            #     #     $\lambda$ = %(wl)f, period = %(d)f, PW = %(pw)i, %(add)s' % \
                            #     #     {'name' : name_lay, 'h':heights_list, 'p':p,'diag' : diag,'wl' : wl, 'd' : period, \
                            #     #     'pw' : pw, 'add' : add_name} + '\n'
                            #     #     + '# prop. ords = %(prop)s, # evan. ords = %(evan)s , n = %(n)s, k = %(k)s' % \
                            #     #     {'evan' : evan, 'prop' : prop, 'n' : n, 'k' : k[0]})
                            #     # plt.savefig('%(dir_name)s/stack_%(stack_num)s_lay_%(name)s_E_diag_slice_y=%(diag)sx_wl=%(wl)s_%(p)s%(add)s.pdf'% \
                            #     #     {'dir_name' : dir_name, 'p':p,'wl' : wl, 'diag' : diag, \
                            #     #     'name' : name_lay,'stack_num':stack_num, 'add' : add_name})

                            # elif sli == 'diag-':
                            #     diag = -1
                            #     fig1 = plt.figure(num=None, figsize=(12,21), dpi=80, facecolor='w', edgecolor='k')
                            #     for i in range(len(E_super)):
                            #         ax1 = fig1.add_subplot(4,1,i+1)
                            #         if re_im == 'real':
                            #             E_slice = np.real(E_super[i][:,0,:])
                            #         elif re_im == 'imag':
                            #             E_slice = np.imag(E_super[i][:,0,:])
                            #         y_min = y_range[0]
                            #         y_max = np.around(np.sqrt(2)*y_range[-1],decimals=2)
                            #         z_min = z1[0]
                            #         z_max = z1[-1]
                            #         cmap = plt.get_cmap('jet')
                            #         CS = plt.contourf(x_axis,y_axis,E_slice,15,cmap=cmap)
                            #         plt.axis([y_min,y_max,z_min,z_max])
                            #         cbar = plt.colorbar()
                            #         cbar.ax.set_ylabel(r'%(p)s E_x'%{'p':p})
                            #         ax1.set_xlabel('x=%(diag)sy (d)'%{'diag':diag})
                            #         ax1.set_ylabel('z (d)')
                            #         if scale_axis == True: ax1.axis('scaled')
                            #         ax1.xaxis.set_ticks([y_min,y_max])
                            #         ax1.set_xlim((y_min,y_max))
                            #         ax1.set_ylim((z_min,z_max))
                            #         if np.abs(z_max-z_min) < x_max: ax1.yaxis.set_ticks([z_min,z_max])

                            #     # plt.suptitle('%(name)s \n E_diagonal_slice_%(p)s, y = %(diag)sx, heights = %(h)s \n\
                            #     #     $\lambda$ = %(wl)f nm, period = %(d)f, PW = %(pw)i, %(add)s' % \
                            #     #         {'name' : name_lay, 'h':heights_list,'diag' : diag, 'p' : p,'wl' : wl, \
                            #     #         'd' : period, 'pw' : pw, 'add' : add_name} + '\n'  +
                            #     #     '# prop. ords = %(prop)s, # evan. ords = %(evan)s , n = %(n)s, k = %(k)s' % \
                            #     #     {'evan' : evan, 'prop' : prop, 'n' : n, 'k' : k[0]})
                            #     # plt.savefig('%(dir_name)s/stack_%(stack_num)s_lay_%(name)s_E_diag_slice_y=%(diag)sx_wl=%(wl)s_%(p)s%(add)s.pdf'% \
                            #     #     {'dir_name' : dir_name, 'p':p,'wl' : wl, 'diag' : diag, \
                            #     #     'name' : name_lay,'stack_num':stack_num, 'add' : add_name})

                            # elif sli == 'special+':
                            #     diag = 1
                            #     fig1 = plt.figure(num=None, figsize=(12,21), dpi=80, facecolor='w', edgecolor='k')
                            #     for i in range(len(E_super)):
                            #         ax1 = fig1.add_subplot(4,1,i+1)
                            #         if re_im == 'real':
                            #             E_slice = np.real(E_super[i][:,0,:])
                            #         elif re_im == 'imag':
                            #             E_slice = np.imag(E_super[i][:,0,:])
                            #         y_min = y_range[0]
                            #         y_max = np.around(sqrt(1+gradient**2)*x1[-1],decimals=2)
                            #         z_min = z1[0]
                            #         z_max = z1[-1]
                            #         cmap = plt.get_cmap('jet')
                            #         CS = plt.contourf(x_axis,y_axis,E_slice,15,cmap=cmap)
                            #         plt.axis([y_min,y_max,z_min,z_max])
                            #         cbar = plt.colorbar()
                            #         cbar.ax.set_ylabel(r'%(p)s E_x'%{'p':p})
                            #         ax1.set_xlabel('y=%(diag)sx (d)'%{'diag':diag*gradient})
                            #         ax1.set_ylabel('z (d)')
                            #         if scale_axis == True: ax1.axis('scaled')
                            #         ax1.xaxis.set_ticks([y_min,y_max])
                            #         ax1.set_ylim((z_min,z_max))
                            #         ax1.set_xlim((y_min,y_max))
                            #         if np.abs(z_max-z_min) < x_max: ax1.yaxis.set_ticks([z_min,z_max])

                            #     # plt.suptitle('%(name)s \n E_specified_diagonal_slice_%(p)s, y = %(diag*gradient)sx, (x,y):(0,0) to (%(x)s,1), heights = %(h)s \n\
                            #     #     $\lambda$ = %(wl)f, period = %(d)f, PW = %(pw)i, %(add)s' % \
                            #     #     {'name' : name_lay, 'h':heights_list,'diag*gradient':diag*gradient, 'p':p,'wl' : wl, 'd' : period, \
                            #     #     'pw' : pw, 'x':x1[-1], 'add' : add_name} + '\n' +
                            #     #     '# prop. ords = %(prop)s, # evan. ords = %(evan)s , n = %(n)s, k = %(k)s' % \
                            #     #     {'evan' : evan, 'prop' : prop, 'n' : n, 'k' : k[0]})
                            #     # plt.savefig('%(dir_name)s/stack_%(stack_num)s_lay_%(name)s_E_specified_diagonal_slice_y=%(diag*gradient)sx_wl=%(wl)s_%(p)s%(add)s.pdf'% \
                            #     #     {'dir_name' : dir_name, 'diag*gradient':diag*gradient, 'p':p,'wl' : wl,\
                            #     #     'name' : name_lay,'stack_num':stack_num, 'add' : add_name})

                            # elif sli == 'special-':
                            #     diag = -1
                            #     fig1 = plt.figure(num=None, figsize=(12,21), dpi=80, facecolor='w', edgecolor='k')
                            #     for i in range(len(E_super)):
                            #         ax1 = fig1.add_subplot(4,1,i+1)
                            #         if re_im == 'real':
                            #             E_slice = np.real(E_super[i][:,0,:])
                            #         elif re_im == 'imag':
                            #             E_slice = np.imag(E_super[i][:,0,:])
                            #         y_min = y_range[0]
                            #         y_max = np.around(sqrt(1+gradient**2)*(x1[0]-x1[-1]),decimals=2)
                            #         z_min = z1[0]
                            #         z_max = z1[-1]
                            #         cmap = plt.get_cmap('jet')
                            #         CS = plt.contourf(x_axis,y_axis,E_slice,15,cmap=cmap)
                            #         plt.axis([y_min,y_max,z_min,z_max])
                            #         cbar = plt.colorbar()
                            #         cbar.ax.set_ylabel(r'%(p)s E_x'%{'p':p})
                            #         ax1.set_xlabel('y=%(diag)sx (d)'%{'diag':diag*gradient})
                            #         ax1.set_ylabel('z (d)')
                            #         if scale_axis == True: ax1.axis('scaled')
                            #         ax1.xaxis.set_ticks([y_min,y_max])
                            #         ax1.set_ylim((z_min,z_max))
                            #         ax1.set_xlim((y_min,y_max))
                            #         if np.abs(z_max-z_min) < x_max: ax1.yaxis.set_ticks([z_min,z_max])

                            #     # plt.suptitle('%(name)s \n E_specified_diagonal_slice_%(p)s, y = %(diag*gradient)sx, (x,y):(1,0) to (%(x)s,1), heights = %(h)s \n\
                            #     #     $\lambda$ = %(wl)f, period = %(d)f, PW = %(pw)i, %(add)s' % \
                            #     #     {'name' : name_lay, 'h':heights_list,'diag*gradient':diag*gradient, 'p':p,'wl' : wl, 'd' : period, \
                            #     #     'pw' : pw,'x':x1[-1], 'add' : add_name} + '\n' +
                            #     #     '# prop. ords = %(prop)s, # evan. ords = %(evan)s , n = %(n)s, k = %(k)s' % \
                            #     #     {'evan' : evan, 'prop' : prop, 'n' : n, 'k' : k[0]})
                            #     # plt.savefig('%(dir_name)s/stack_%(stack_num)s_lay_%(name)s_E_specified_diagonal_slice_y=%(diag*gradient)sx_wl=%(wl)s_%(p)s%(add)s.pdf'% \
                            #     #     {'dir_name' : dir_name, 'diag*gradient':diag*gradient, 'p':p,'wl' : wl,\
                            #     #     'name' : name_lay,'stack_num':stack_num, 'add' : add_name})

        stack_num += 1


def field_values(stacks_list, lay_interest=0, xyz_values=[(0.1,0.1,0.1)]):
    """
    Save electric field values at given x-y-z points. Points must be within \
    a ThinFilm layer. In txt file fields are given as \
    Re(Ex) Im(Ex) Re(Ey) Im(Ey) Re(Ez) Im(Ez) |E|

        Args:
            stacks_list  (list): Stack objects containing data to plot.

        Keyword Args:
            lay_interest  (int): the index of the layer considered within \
                the stack. Must be a ThinFilm layer.

            xyz_values  (list): list of distances in normalised units of (d) \
                from top surface of layer at which to calculate fields. \
                For semi-inf substrate then z_value is distance from \
                top of this layer (i.e. bottom interface of stack).
    """

    dir_name = 'field_values'
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)

    stack_num = 0
    super_points = []
    for pstack in stacks_list:
        try:
            if not isinstance(pstack.layers[lay_interest],mode_calcs.Anallo):
                raise ValueError

            num_lays = len(pstack.layers)
            wl = np.around(pstack.layers[-1].light.wl_nm,decimals=2)

            if lay_interest == 0: name_lay = "0_Substrate"
            elif lay_interest == num_lays-1: name_lay = "%i_Superstrate" % (num_lays-1)
            else: name_lay = "%i_Thin_Film" % lay_interest

            n = pstack.layers[lay_interest].n()
            period = pstack.layers[-1].structure.period
            PWordtot = pstack.layers[lay_interest].structure.num_pw_per_pol
            s = pstack.layers[lay_interest].sort_order
            alpha_unsrt = np.array(pstack.layers[lay_interest].alphas)
            beta_unsrt = np.array(pstack.layers[lay_interest].betas)
            alpha = alpha_unsrt[s]
            if pstack.layers[lay_interest].structure.world_1d is True:
                beta = beta_unsrt
            else:
                beta = beta_unsrt[s]
            gamma = np.array(pstack.layers[lay_interest].calc_kz())

            vec_coef_down = np.array(pstack.vec_coef_down[num_lays-1-lay_interest]).flatten()
            vec_coef_down_TE = vec_coef_down[0:PWordtot]
            vec_coef_down_TM = vec_coef_down[PWordtot::]
            if lay_interest == 0:
                vec_coef_up = np.zeros((2*PWordtot), dtype = 'complex')
            else:
                vec_coef_up = np.array(pstack.vec_coef_up[num_lays-1-lay_interest]).flatten()

            vec_coef_up_TE = vec_coef_up[0:PWordtot]
            vec_coef_up_TM = vec_coef_up[PWordtot::]

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

            eta_TE_x_down = (vec_coef_down_TE*E_TE_x)/chi_TE
            eta_TE_y_down = (vec_coef_down_TE*E_TE_y)/chi_TE
            eta_TE_z_down = (vec_coef_down_TE*E_TE_z)/chi_TE
            eta_TM_x_down = (vec_coef_down_TM*E_TM_x)/chi_TM
            eta_TM_y_down = (vec_coef_down_TM*E_TM_y)/chi_TM
            eta_TM_z_down = (vec_coef_down_TM*E_TM_z)/chi_TM

            eta_TE_x_up = (vec_coef_up_TE*E_TE_x)/chi_TE
            eta_TE_y_up = (vec_coef_up_TE*E_TE_y)/chi_TE
            eta_TE_z_up = (vec_coef_up_TE*E_TE_z)/chi_TE
            eta_TM_x_up = (vec_coef_up_TM*E_TM_x)/chi_TM
            eta_TM_y_up = (vec_coef_up_TM*E_TM_y)/chi_TM
            eta_TM_z_up = (vec_coef_up_TM*E_TM_z)/chi_TM

            calc_E_TE_x_array = np.zeros(len(xyz_values),dtype='complex')
            calc_E_TE_y_array = np.zeros(len(xyz_values),dtype='complex')
            calc_E_TE_z_array = np.zeros(len(xyz_values),dtype='complex')
            calc_E_TM_x_array = np.zeros(len(xyz_values),dtype='complex')
            calc_E_TM_y_array = np.zeros(len(xyz_values),dtype='complex')
            calc_E_TM_z_array = np.zeros(len(xyz_values),dtype='complex')

            for i in xrange(len(xyz_values)):

                (x1,y1,z1) = np.array(xyz_values[i])#/float(pstack.layers[lay_interest].structure.period)

                if pstack.layers[lay_interest].structure.world_1d == True: y1 = 0

                if lay_interest == 0: z1 = -1*z1
                else: z1 = np.abs(z1)

                if pstack.layers[lay_interest].structure.height_nm == 'semi_inf':
                    calc_expo_down = np.exp(1j*(alpha*x1+beta*y1-gamma*z1))
                else:
                    calc_expo_down = np.exp(1j*(alpha*x1+beta*y1-gamma*(z1-float(pstack.layers[lay_interest].structure.height_nm)/period)))
                calc_expo_up = np.exp(1j*(alpha*x1+beta*y1+gamma*z1))

                calc_E_TE_x = np.sum(eta_TE_x_down*calc_expo_down + eta_TE_x_up*calc_expo_up)
                calc_E_TE_y = np.sum(eta_TE_y_down*calc_expo_down + eta_TE_y_up*calc_expo_up)
                calc_E_TE_z = np.sum(eta_TE_z_down*calc_expo_down + eta_TE_z_up*calc_expo_up)
                calc_E_TM_x = np.sum(eta_TM_x_down*calc_expo_down + eta_TM_x_up*calc_expo_up)
                calc_E_TM_y = np.sum(eta_TM_y_down*calc_expo_down + eta_TM_y_up*calc_expo_up)
                calc_E_TM_z = np.sum(eta_TM_z_down*calc_expo_down + eta_TM_z_up*calc_expo_up)
                calc_E_TE_x_array[i] = calc_E_TE_x
                calc_E_TE_y_array[i] = calc_E_TE_y
                calc_E_TE_z_array[i] = calc_E_TE_z
                calc_E_TM_x_array[i] = calc_E_TM_x
                calc_E_TM_y_array[i] = calc_E_TM_y
                calc_E_TM_z_array[i] = calc_E_TM_z
            calc_E_x_array = calc_E_TE_x_array + calc_E_TM_x_array
            calc_E_y_array = calc_E_TE_y_array + calc_E_TM_y_array
            calc_E_z_array = calc_E_TE_z_array + calc_E_TM_z_array
            calc_E_tot_array = np.sqrt(calc_E_x_array*np.conj(calc_E_x_array) \
                                +calc_E_y_array*np.conj(calc_E_y_array)+calc_E_z_array*np.conj(calc_E_z_array))

            np.savez('%(dir_name)s/%(stack_num)s_E_calc_points_%(name)s_wl_%(wl)s'% \
                                {'dir_name':dir_name, 'wl':wl,'name' : name_lay,'stack_num':stack_num},\
                                calc_E_x_array=calc_E_x_array,calc_E_y_array=calc_E_y_array,\
                                calc_E_z_array=calc_E_z_array,calc_E_tot_array=calc_E_tot_array)

            np.savetxt('%(dir_name)s/%(stack_num)s_E_calc_points_%(name)s_wl_%(wl)s.txt'% \
                                {'dir_name':dir_name, 'wl':wl,'name' : name_lay,'stack_num':stack_num},\
                                np.array([np.real(calc_E_x_array), np.imag(calc_E_x_array), np.real(calc_E_y_array),\
                                np.imag(calc_E_y_array), np.real(calc_E_z_array), np.imag(calc_E_z_array), np.real(calc_E_tot_array)]))
            stack_num += 1
            super_points.append([calc_E_x_array,calc_E_y_array,calc_E_z_array,calc_E_tot_array])

        except ValueError:
            print("field_values can only calculate field values in ThinFilms."\
            "\nPlease select a different lay_interest.\n")

        return super_points


def fields_3d(stacks_list, lay_interest=1):
    """
    Plot fields in 3D using gmsh.

        Args:
            stacks_list  (list): Stack objects containing data to plot.

        Keyword Args:
            lay_interest  (int): the index of the layer considered within \
            the stack.
    """
    from fortran import EMUstack
    import subprocess

    nnodes=6
    dir_name = "3d_fields"
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)
    dir_name = "3d_fields/anim"
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)

    stack_num = 0
    for pstack in stacks_list:
        try:
            if not isinstance(pstack.layers[lay_interest],mode_calcs.Simmo):
                raise ValueError

            meat = pstack.layers[lay_interest]
            if not meat.structure.periodicity == '2D_array':
                raise ValueError

            gmsh_file_pos = meat.structure.mesh_file[0:-5]

            # vec_coef sorted from top of stack, everything else is sorted from bottom
            vec_index = len(pstack.layers) - lay_interest - 1

            vec_coef = np.concatenate((pstack.vec_coef_down[vec_index],pstack.vec_coef_up[vec_index]))
            # vec_coef_up = np.zeros(shape=(np.shape(pstack.vec_coef_down[vec_index])),dtype='complex128')
            # vec_coef = np.concatenate((pstack.vec_coef_down[vec_index],vec_coef_up))
            h_normed = float(meat.structure.height_nm)/float(meat.structure.period)
            wl_normed = pstack.layers[lay_interest].wl_norm()

            layer_name = 'Lay_' + zeros_int_str(lay_interest) + 'Stack_' + str(stack_num)

            EMUstack.gmsh_plot_field_3d(wl_normed, h_normed, meat.num_BMs,
                meat.E_H_field, meat.n_msh_el, meat.n_msh_pts,
                nnodes, meat.type_el, meat.structure.nb_typ_el, meat.table_nod,
                meat.k_z, meat.sol1, vec_coef, meat.x_arr, gmsh_file_pos, layer_name)

            stack_num += 1
        except ValueError:
            print("fields_3d can only plot 3D fields within 2D_array "\
            "Nanostruct layers. \nPlease select a different lay_interest.\n")


def Bloch_fields_1d(stacks_list, lay_interest=None):
    """
    Plot Bloch mode fields along the x axis.

        Args:
            stacks_list  (list): Stack objects containing data to plot.

        Keyword Args:
            lay_interest  (int): the index of the layer considered within \
                the stack. Must be a 1D_array NanoStruct layer. By default \
                routine finds all such layers.
    """

    dir_name = "Bloch_fields"
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)

    stack_num = 0
    for pstack in stacks_list:
        # by default find all 1D layers
        if lay_interest == None:
            lay_interest = []
            for l in range(len(pstack.layers)):
                if isinstance(pstack.layers[l],mode_calcs.Simmo):
                    if pstack.layers[l].structure.periodicity == '1D_array':
                        lay_interest.append(l)
        else: lay_interest = [lay_interest]

        num_lays = len(pstack.layers)
        for lay in lay_interest:
            try:
                meat = pstack.layers[lay]
                if not isinstance(meat,mode_calcs.Simmo):
                    raise ValueError
                if not meat.structure.periodicity == '1D_array':
                    raise ValueError

                eps = meat.n_effs**2

                struct = meat.structure
                boundary = []
                for i in range(len(struct.type_el) - 1):
                    if struct.type_el[i] != struct.type_el[i+1]:
                        boundary.append(struct.x_arr[2*(i+1)])

                fields = ['Re(E_x)', 'Im(E_x)', 'Re(E_y)', 'Im(E_y)', 'Re(E_z)', 'Im(E_z)']
                for BM in range(meat.num_BMs):
                    fig = plt.figure(figsize=(12, 12))
                    plot_sol_abs = np.zeros(struct.n_msh_pts)
                    plot_sol_eps_abs = np.zeros(struct.n_msh_pts)
                    # plot_sol_2_FT    = np.zeros(struct.n_msh_pts)
                    for i in range(len(fields)):
                        E = fields[i]
                        # sol_P2([Ex,Ey,Ez],P2_interpolation_points,nval,nel)
                        BM_sol = meat.sol1[i/2,:,BM,:]
                        if E[0] == 'R': plot_sol = [np.real(BM_sol[0,0])]
                        else: plot_sol = [np.imag(BM_sol[0,0])]
                        plot_sol_abs[0] += np.real(np.sqrt((BM_sol[0,0])* np.conj(BM_sol[0,0])))
                        plot_sol_eps_abs[0] += np.real(eps[struct.type_el[0]-1] * np.sqrt((BM_sol[0,0])* np.conj(BM_sol[0,0])))
                        for x in range(struct.n_msh_el - 1):
                            if E[0] == 'R':
                                plot_sol.append(np.real(BM_sol[1,x]))
                                plot_sol.append(np.real(BM_sol[2,x] + BM_sol[0,x+1]) / 2)
                            else:
                                plot_sol.append(np.imag(BM_sol[1,x]))
                                plot_sol.append(np.imag(BM_sol[2,x] + BM_sol[0,x+1]) / 2)
                            plot_sol_abs[2*x+1] += np.real(np.sqrt((BM_sol[1,x])* np.conj(BM_sol[1,x])))
                            plot_sol_abs[2*x+2] += np.real(np.sqrt(((BM_sol[2,x] + BM_sol[0,x+1]) / 2)* np.conj((BM_sol[2,x] + BM_sol[0,x+1]) / 2)))
                            plot_sol_eps_abs[2*x+1] += np.real(eps[struct.type_el[x]-1] * np.sqrt((BM_sol[1,x])* np.conj(BM_sol[1,x])))
                            plot_sol_eps_abs[2*x+2] += np.real(eps[struct.type_el[x]-1] * np.sqrt(((BM_sol[2,x] + BM_sol[0,x+1]) / 2)* np.conj((BM_sol[2,x] + BM_sol[0,x+1]) / 2)))

                        if E[0] == 'R':
                            plot_sol.append(np.real(BM_sol[1,-1]))
                            plot_sol.append(np.real(BM_sol[2,-1]))
                        else:
                            plot_sol.append(np.imag(BM_sol[1,-1]))
                            plot_sol.append(np.imag(BM_sol[2,-1]))
                        plot_sol_abs[-2] += np.real(np.sqrt((BM_sol[1,x])* np.conj(BM_sol[1,x])))
                        plot_sol_abs[-1] += np.real(np.sqrt((BM_sol[2,x])* np.conj(BM_sol[2,x])))
                        plot_sol_eps_abs[-2] += np.real(eps[struct.type_el[-1]-1] * np.sqrt((BM_sol[1,x])* np.conj(BM_sol[1,x])))
                        plot_sol_eps_abs[-1] += np.real(eps[struct.type_el[-1]-1] * np.sqrt((BM_sol[2,x])* np.conj(BM_sol[2,x])))

                        ax1 = fig.add_subplot(6, 2, i+1)
                        ax1.plot(struct.x_arr, plot_sol)
                        ax1.set_ylabel(E)
                        if E == 'Re(E_x)' or E == 'Im(E_x)':
                            ax1.xaxis.set_label_position("top")
                            ax1.xaxis.tick_top()
                            ax1.xaxis.set_ticks_position('both')
                            ax1.set_xlabel('x (d)')
                        else:
                            ax1.set_xticklabels( () )
                        start, end = ax1.get_ylim()
                        ax1.yaxis.set_ticks(np.linspace(start, end, 3))
                        for b in boundary:
                            ax1.plot([b,b],[start,end],'k')
                        if i%2 == 1:
                            ax1.yaxis.set_label_position("right")
                            ax1.yaxis.tick_right()
                            ax1.yaxis.set_ticks_position('both')

                        if E == 'Re(E_x)': plot_sol_2_FTx = plot_sol
                        if E == 'Im(E_x)': plot_sol_2_FT2x = plot_sol
                        if E == 'Re(E_y)': plot_sol_2_FT = plot_sol
                        if E == 'Im(E_y)': plot_sol_2_FT2 = plot_sol
                    ax1 = fig.add_subplot(6, 2, 7)
                    ax1.plot(struct.x_arr, plot_sol_abs)
                    ax1.set_xticklabels( () )
                    start, end = ax1.get_ylim()
                    ax1.yaxis.set_ticks(np.linspace(start, end, 3))
                    for b in boundary:
                        ax1.plot([b,b],[start,end],'k')
                    ax1.set_ylabel(r'|E|')
                    ax1 = fig.add_subplot(6, 2, 8)
                    ax1.plot(struct.x_arr, plot_sol_eps_abs)
                    ax1.set_xticklabels( () )
                    start, end = ax1.get_ylim()
                    ax1.yaxis.set_ticks(np.linspace(start, end, 3))
                    for b in boundary:
                        ax1.plot([b,b],[start,end],'k')
                    ax1.yaxis.set_label_position("right")
                    ax1.yaxis.tick_right()
                    ax1.yaxis.set_ticks_position('both')
                    ax1.set_ylabel(r'Re($\epsilon$) |E|')

                    ax1 = fig.add_subplot(6, 2, 9)
                    xshape = struct.x_arr.size
                    # plot_sol_2_FT = np.sin(3*np.linspace(0,2*np.pi,xshape))
                    FT = np.fft.fft(plot_sol_2_FTx)
                    freqs = np.fft.fftfreq(xshape,struct.x_arr[1])
                    freqs = np.fft.fftshift(freqs)
                    FT = np.fft.fftshift(FT)
                    ax1.plot(freqs, np.abs(FT),'ro', linewidth=3.0)
                    start, end = ax1.get_ylim()
                    ax1.yaxis.set_ticks(np.linspace(start, end, 3))
                    ax1.set_ylabel(r'FT Re(E$_x$)')
                    ax1.set_xlim((0,10))

                    ax1 = fig.add_subplot(6, 2, 10)
                    xshape = struct.x_arr.size
                    FT = np.fft.fft(plot_sol_2_FT2x)
                    freqs = np.fft.fftfreq(xshape,struct.x_arr[1])
                    freqs = np.fft.fftshift(freqs)
                    FT = np.fft.fftshift(FT)
                    ax1.plot(freqs, np.abs(FT),'ro', linewidth=3.0)
                    start, end = ax1.get_ylim()
                    ax1.yaxis.set_ticks(np.linspace(start, end, 3))
                    ax1.set_ylabel(r'FT Im(E$_x$)')
                    ax1.set_xlim((0,10))
                    ax1.yaxis.set_label_position("right")
                    ax1.yaxis.tick_right()
                    ax1.yaxis.set_ticks_position('both')

                    ax1 = fig.add_subplot(6, 2, 11)
                    xshape = struct.x_arr.size
                    FT = np.fft.fft(plot_sol_2_FT)
                    freqs = np.fft.fftfreq(xshape,struct.x_arr[1])
                    freqs = np.fft.fftshift(freqs)
                    FT = np.fft.fftshift(FT)
                    ax1.plot(freqs, np.abs(FT),'ro', linewidth=3.0)
                    start, end = ax1.get_ylim()
                    ax1.yaxis.set_ticks(np.linspace(start, end, 3))
                    ax1.set_ylabel(r'FT Re(E$_y$)')
                    ax1.set_xlim((0,10))
                    ax1.set_xlabel('Diffraction order')

                    ax1 = fig.add_subplot(6, 2, 12)
                    xshape = struct.x_arr.size
                    FT = np.fft.fft(plot_sol_2_FT2)
                    freqs = np.fft.fftfreq(xshape,struct.x_arr[1])
                    freqs = np.fft.fftshift(freqs)
                    FT = np.fft.fftshift(FT)
                    ax1.plot(freqs, np.abs(FT),'ro', linewidth=3.0)
                    start, end = ax1.get_ylim()
                    ax1.yaxis.set_ticks(np.linspace(start, end, 3))
                    ax1.set_ylabel(r'FT Im(E$_y$)')
                    ax1.set_xlim((0,10))
                    ax1.set_xlabel('Diffraction order')
                    ax1.yaxis.set_label_position("right")
                    ax1.yaxis.tick_right()
                    ax1.yaxis.set_ticks_position('both')

                    name_lay = "layer-%(lay)i-BM_%(BM)s"% \
                        {'lay' : lay, 'BM' : zeros_int_str(BM)}

                    plt.suptitle(r'BM %(BM)i, k$_z$ = %(re)7.3f + %(im)7.3fi (d)'% \
                        {'re' : np.real(meat.k_z[BM]), 'im' : np.imag(meat.k_z[BM]), 'BM' : BM})
                    plt.savefig(dir_name + '/' + name_lay, bbox_inches = 'tight')
                    fig.clf()

            except ValueError:
                print("fields_1d can only plot 1D fields of 1D_array "\
                "NanoStruct layers. \nPlease select different lay_interest.\n")

###############################################################################

#### Fabry-Perot resonances ###################################################
def Fabry_Perot_res(stacks_list, freq_list, kx_list, f_0, k_0, lay_interest=1):
    """ Calculate the Fabry-Perot resonance condition for a resonances within a layer.

        This is equivalent to finding the slab waveguide modes of the layer.

        Args:
            stacks_list  (list): Stack objects containing data to plot.

            freq_list  (list): Frequencies included.

            kx_list  (list): In-plane wavenumbers included.

            f_0  (list): Frequency w.r.t. which axis is normalised.

            k_0  (list): In-plane wavenumber w.r.t. which axis is normalised.

        Keyword Args:
            lay_interest  (int): The index in stacks_list of the layer of\
            which F-P resonances are calculated.
    """

    n_freqs = len(freq_list)
    n_kxs = len(kx_list)
    # assuming all stacks have equal height
    height = stacks_list[-1].heights_nm()[lay_interest-1]
    period = stacks_list[-1].period
    num_BMss = stacks_list[-1].layers[lay_interest].num_BMs
    I_mat = np.matrix(np.eye(num_BMss),dtype='D')

    plot_mat = np.zeros(shape=(n_freqs,n_kxs),dtype='complex128')
    for i in range(n_kxs):
        kx_slice = stacks_list[i*n_freqs:(i+1)*n_freqs]
        j = 0
        for stack in kx_slice:
            P = stack.layers[lay_interest].prop_fwd(height/period)
            R21 = stack.layers[lay_interest].R21
            FP_term = np.linalg.det(I_mat - R21*P*R21*P)
            plot_mat[j, i] = FP_term
            j += 1
    image = np.log(np.abs(plot_mat))

    fig = plt.figure()
    ax1 = fig.add_subplot(1, 1, 1)
    cax = ax1.imshow(image, cmap=plt.cm.gray_r, interpolation='none')

    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    shape = np.shape(plot_mat)
    majorLocator = MultipleLocator(shape[1]-1)
    ax1.xaxis.set_major_locator(majorLocator)
    majorLocator = MultipleLocator(shape[0]-1)
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
    plt.savefig('Fabry-Perot_resonances', bbox_inches = 'tight')
    fig.clf()
###############################################################################
