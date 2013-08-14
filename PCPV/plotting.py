#
import objects
import mode_calcs

import numpy as np
from matplotlib.mlab import griddata
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt




#######################################################################################
def gen_params_string(param_layer, light, max_num_BMs=0):
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
        # params_2_print = 'Homogeneous Film'
        params_2_print = 'max_PW_order = %(max_order_PWs)d, '% {'max_order_PWs' : light.max_order_PWs, }


    # k_pll = light.k_pll * param_layer.period
    # params_2_print += r'$k_\parallel d$ = '
    # tmp3 = '(%(kx)1.4f, %(ky)1.4f), '% {'kx' : k_pll[0], 'ky' : k_pll[1]}
    # params_2_print += r'$\theta$ = %(theta)6.2f, $\phi$ = %(phi)6.2f, '% {
    # 'theta' : light.theta,'phi' : light.phi, } 

    return params_2_print
#######################################################################################


#######################################################################################
def t_r_a_plots(stack_wl_list, wavelengths, params_2_print, active_layer_nu=0, stack_label=1, add_name=''):
    # Plot t,r,a for each layer & total, then save each to text files

    height_list = stack_wl_list[0].heights_nm()[::-1]
    params_2_print += '\n'r'$h_t,...,h_b$ = '
    params_2_print += ''.join('%4d, ' % num for num in height_list)

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
        t_tot.append(float(a_list[layers_steps-1+(i*layers_steps)]))
        r_tot.append(float(a_list[layers_steps-1+(i*layers_steps)]))
    Efficiency, Irradiance = ult_efficiency(active_abs, wavelengths)
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
    weighting = Irradiance/Irradiance.max()
    a_weighted = a_tot*weighting
    t_weighted = t_tot*weighting
    r_weighted = r_tot*weighting
    plot_name = 'Weighted_Total_Spectra'
    total_tra_plot(plot_name, a_weighted, t_weighted, r_weighted, wavelengths, params_2_print, stack_label, add_name)
    return Efficiency


def ult_efficiency(active_abs, wavelengths):

    # TODO make E_g a property of material, not just longst wl included.

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
    # np.savetxt('Efficiency.txt', [Efficiency, diameter1, diameter2,
    # period, ff], fmt = '%12.8f')
    return Efficiency, i_spec


def layers_plot(spectra_name, spec_list, wavelengths, total_h, params_2_print, stack_label, add_name):
    fig = plt.figure(num=None, figsize=(9, 12), dpi=80, facecolor='w', edgecolor='k')
    nu_layers = len(spec_list)/len(wavelengths)
    h_array = np.ones(len(wavelengths))*total_h
    h  = 6.626068e-34;
    c  = 299792458;
    eV = 1.60217646e-19;
    energies  = (h*c/(wavelengths*1e-9))/eV
    for i in range(nu_layers):
        layer_spec = []
        for wl in range(len(wavelengths)):
            layer_spec = np.append(layer_spec,spec_list[wl*nu_layers + i])
        ax1 = fig.add_subplot(nu_layers,1,i+1)
        ax1.plot(wavelengths, layer_spec)
        ax2 = ax1.twiny()
        ax2.plot(energies, layer_spec, alpha=0)
        if i == nu_layers-1:
            ax1.set_xlabel('Wavelength (nm)')
            ax2.set_xlabel('Energy (eV)')
            ax1.set_ylabel('Total')
        else:
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
        ax2.set_xlim((energies[0], energies[-1]))
        plt.ylim((0, 1))

        if i != nu_layers-1: 
            np.savetxt('%(s)s_%(i)i_stack%(bon)i%(add)s.txt'% {'s' : lay_spec_name, 'i' : i, 
                'bon' : stack_label,'add' : add_name}, av_array, fmt = '%18.11f')
        else:
            np.savetxt('%(s)s_stack%(bon)i%(add)s.txt'% {'s' : lay_spec_name, 
                'bon' : stack_label,'add' : add_name}, av_array, fmt = '%18.11f')

        plt.savefig('%(s)s_stack%(bon)i%(add)s'% {'s' : spectra_name, 'bon' : stack_label,'add' : add_name})

def total_tra_plot(plot_name, a_spec, t_spec, r_spec, wavelengths, params_2_print, stack_label, add_name):
    h  = 6.626068e-34;
    c  = 299792458;
    eV = 1.60217646e-19;
    energies  = (h*c/(wavelengths*1e-9))/eV
    # num_e_ticks = 4
    # e_ticks_tmp = np.linspace(energies.min(), energies.max(), num_e_ticks)
    # e_ticks = [ round(elem, 2) for elem in e_ticks_tmp ]
    fig = plt.figure(num=None, figsize=(9, 12), dpi=80, facecolor='w', edgecolor='k')
    ax1 = fig.add_subplot(3,1,1)
    ax1.plot(wavelengths, a_spec)
    ax1.set_xlabel('Wavelength (nm)')
    ax1.set_ylabel('Absorptance')
    ax1.set_xlim((wavelengths[0], wavelengths[-1]))
    ax2 = ax1.twiny()
    ax2.plot(energies, a_spec, alpha=0)
    # ax2.set_xticks(e_ticks)
    ax2.set_xlabel('Energy (eV)')
    ax2.set_xlim((energies[0], energies[-1]))
    plt.ylim((0, 1))
    ax1 = fig.add_subplot(3,1,2)
    ax1.plot(wavelengths, t_spec)
    ax1.set_xlabel('Wavelength (nm)')
    ax1.set_ylabel('Transmittance')
    ax1.set_xlim((wavelengths[0], wavelengths[-1]))
    ax2 = ax1.twiny()
    ax2.plot(energies, t_spec, alpha=0)
    ax2.set_xticklabels( () )
    ax2.set_xlim((energies[0], energies[-1]))
    plt.ylim((0, 1))
    ax1 = fig.add_subplot(3,1,3)
    ax1.plot(wavelengths, r_spec)
    ax1.set_xlabel('Wavelength (nm)')
    ax1.set_ylabel('Reflectance')
    ax1.set_xlim((wavelengths[0], wavelengths[-1]))
    ax2 = ax1.twiny()
    ax2.plot(energies, r_spec, alpha=0)
    ax2.set_xticklabels( () )
    ax2.set_xlim((energies[0], energies[-1]))
    plt.ylim((0, 1))
    plt.suptitle(plot_name+add_name+'\n'+params_2_print)
    plt.savefig('%(s)s_stack%(bon)i%(add)s'% {'s' : plot_name, 'bon' : stack_label,'add' : add_name})
#######################################################################################



#######################################################################################
def omega_plot(stack_wl_list, wavelengths, params_2_print, stack_label=1):
    num_layers = len(stack_wl_list[0].layers)
    fig1 = plt.figure(num=None, figsize=(9, 12), dpi=80, facecolor='w', edgecolor='k')
    fig2 = plt.figure(num=None, figsize=(9, 12), dpi=80, facecolor='w', edgecolor='k')
    for l in range(num_layers):
        ax1 = fig1.add_subplot(num_layers,1,l+1)
        ax2 = fig2.add_subplot(num_layers,1,l+1)
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
        ax1.set_ylabel(r'Real $k_z$')
        ax2.set_ylabel(r'Imaginary $k_z$')
        if l == 0: ax1.set_ylabel('Top Layer'), ax2.set_ylabel('Top Layer')
        if l != num_layers-1:
            ax1.set_xticklabels( () )
            ax2.set_xticklabels( () )
        else:
            ax1.set_ylabel('Bottom Layer'), ax2.set_ylabel('Bottom Layer')
            ax1.set_xlabel('Wavelength (nm)'), ax2.set_xlabel('Wavelength (nm)')
        plt.xlim((wavelengths[0], wavelengths[-1]))
    fig1.suptitle(r'Real $k_z$'+params_2_print+'\n')
    fig2.suptitle(r'Imaginary $k_z$'+params_2_print+'\n')
    fig1.savefig('Disp_Diagram_Re_stack%(bon)i'% {'bon' : stack_label})
    fig2.savefig('Disp_Diagram_Im_stack%(bon)i'% {'bon' : stack_label})
    # np.savetxt('Disp_Data_stack%(bon)i.txt'% {'bon' : stack_label}, av_array, fmt = '%18.11f')

def E_conc_plot(stack_wl_list, which_layer, which_modes, wavelengths, params_2_print, stack_label=1):
    if isinstance(stack_wl_list[0].layers[which_layer], mode_calcs.Simmo):
        num_layers = len(stack_wl_list[0].layers)
        fig1 = plt.figure(num=None, figsize=(9, 4), dpi=80, facecolor='w', edgecolor='k')
        ax1 = fig1.add_subplot(1,1,1)
        # ax2 = fig1.add_subplot(2,1,2)
        # E_conc = []
        for i in range(len(wavelengths)):
            for mode in which_modes:
                E_conc_tmp = np.real(stack_wl_list[i].layers[which_layer].mode_pol[3,mode])
                # E_conc.append(E_conc_tmp)
                ax1.plot(wavelengths[i],E_conc_tmp,'bo')
        plt.xlim((wavelengths[0], wavelengths[-1]))
        ax1.set_xlabel('Wavelength (nm)')
        ax1.set_ylabel(r'$E_{cyl} / E_{cell}$')
        # ax2.plot(wavelengths,E_conc,'k')
        # ax2.set_xlabel('Wavelength (nm)')
        # ax2.set_ylabel(r'$E_{cyl} / E_{cell}$')
        fig1.suptitle('Energy Concentration = '+r'$E_{cyl} / E_{cell}$'+'\n'+params_2_print)
        fig1.savefig('Energy_Concentration_stack%(bon)i'% {'bon' : stack_label})
    else:
        print "\n ERROR: plotting.E_conc_plot; \n Can only calculate energy concentration in NanoStruct layers."
        print repr(stack_wl_list[0].layers[which_layer])

# np.genfromtxt can deal with incomplete data!
#     num_BMs     = np.genfromtxt(file_name, usecols=(1))
#######################################################################################






#######################################################################################
# plot scattering matrices as grayscale images
def vis_scat_mats(scat_mat,wl=0,extra_title=''):
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