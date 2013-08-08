#
import objects

import numpy as np
from matplotlib.mlab import griddata
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import os




#######################################################################################
def t_r_a_plots(stack_wl_list, param_layer, light, max_order_PWs, max_num_BMs=0,
    active_layer_nu=1, stack_label=1, additional_info=''):
    # Plot t,r,a for each layer & total, then save each to text files
    if isinstance(param_layer,objects.NanoStruct):
        if param_layer.geometry == 'NW_array':
            params_2_print = 'd = %(period)d, a1 = %(radius)d, '% {
            'period'        : param_layer.period, 'radius' : param_layer.radius1,}
            if param_layer.radius2 != 0: params_2_print += 'a2 = %(rad)d '% {'rad' : param_layer.radius2,}
            if param_layer.radius3 != 0: params_2_print += 'a3 = %(rad)d '% {'rad' : param_layer.radius3,}
            if param_layer.radius4 != 0: params_2_print += 'a4 = %(rad)d '% {'rad' : param_layer.radius4,}
            if param_layer.radius5 != 0: params_2_print += '\na5 = %(rad)d '% {'rad' : param_layer.radius5,}
            if param_layer.radius6 != 0: params_2_print += 'a6 = %(rad)d '% {'rad' : param_layer.radius6,}
            if param_layer.radius7 != 0: params_2_print += 'a7 = %(rad)d '% {'rad' : param_layer.radius7,}
            if param_layer.radius8 != 0: params_2_print += 'a8 = %(rad)d '% {'rad' : param_layer.radius8,}
            if param_layer.radius9 != 0: params_2_print += 'a9 = %(rad)d \n'% {'rad' : param_layer.radius9,}
            if param_layer.radius10 != 0: params_2_print += 'a10 = %(rad)d '% {'rad' : param_layer.radius10,}
            if param_layer.radius11 != 0: params_2_print += 'a11 = %(rad)d '% {'rad' : param_layer.radius11,}
            if param_layer.radius12 != 0: params_2_print += 'a12 = %(rad)d '% {'rad' : param_layer.radius12,}
            if param_layer.radius13 != 0: params_2_print += 'a13 = %(rad)d '% {'rad' : param_layer.radius13,}
            if param_layer.radius14 != 0: params_2_print += 'a14 = %(rad)d '% {'rad' : param_layer.radius14,}
            if param_layer.radius15 != 0: params_2_print += 'a15 = %(rad)d '% {'rad' : param_layer.radius15,}
            if param_layer.radius16 != 0: params_2_print += 'a16 = %(rad)d \n'% {'rad' : param_layer.radius16,}
            params_2_print += 'ff = %5.3f, '% param_layer.ff
            if param_layer.square == True: params_2_print += '\nSquare NWs '
            if param_layer.ellipticity == True: params_2_print += '\nEllipticity = %(rad)5.3f '% {'rad' : param_layer.ellipticity,}
            # k_pll = light.k_pll * param_layer.period
            # params_2_print += r'$k_\parallel d$ = '
            # tmp3 = '(%(kx)1.4f, %(ky)1.4f), '% {'kx' : k_pll[0], 'ky' : k_pll[1]}
            # params_2_print += r'$\theta$ = %(theta)6.2f, $\phi$ = %(phi)6.2f, '% {
            # 'theta' : light.theta,'phi' : light.phi, } 
            params_2_print += '\nmax_BMs = %(max_num_BMs)d, max_PW_order = %(max_order_PWs)d, '% {
            'max_num_BMs'   : max_num_BMs,'max_order_PWs' : max_order_PWs, }


        elif param_layer.geometry == '1D_grating':
            params_2_print = ''
    else:
        params_2_print = 'Homogeneous Film\n'
        params_2_print += 'h = %10.2f, '% param_layer.height_nm
        params_2_print += 'max_PW_order = %(max_order_PWs)d, '% {'max_order_PWs' : max_order_PWs, }

    params_2_print += '\nh = %10.2f, '% param_layer.height_nm




    wavelengths = np.array([s.layers[0].light.wl_nm for s in stack_wl_list]) #look at first layer to find wls.
    a_list = []
    t_list = []
    r_list = []
    for stack in stack_wl_list:
        a_list.extend(stack.a_list)
        t_list.extend(stack.t_list)
        r_list.extend(stack.r_list)

    num_layers = len(stack_wl_list[0].layers) - 1
    active_abs = []
    for i in range(len(wavelengths)): 
        active_abs.append(float(a_list[active_layer_nu + i*num_layers]))
    Efficiency = ult_efficiency(active_abs, wavelengths)
    params_2_print += '\n' r'$\eta$ = %(Efficiency)6.2f'% {'Efficiency' : Efficiency*100, }
    params_2_print += ' %'
    print params_2_print

    total_h = sum(stack_wl_list[0].heights_nm()) #look at first wl result to find h.
    layers_plot('Lay_Absorb', a_list, wavelengths, total_h, params_2_print, stack_label)
    layers_plot('Lay_Trans',  t_list, wavelengths, total_h, params_2_print, stack_label)
    layers_plot('Lay_Reflec', r_list, wavelengths, total_h, params_2_print, stack_label)

    return Efficiency


def ult_efficiency(active_abs, wavelengths):
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
    # np.savetxt('Efficiency.txt', [Efficiency, radius1, radius2,
    # period, ff], fmt = '%12.8f')
    return Efficiency



def layers_plot(spectra_name, spec_list, wavelengths, total_h, params_2_print, stack_label):
    fig = plt.figure(num=None, figsize=(9, 12), dpi=80, facecolor='w', edgecolor='k')
    nu_layers = len(spec_list)/len(wavelengths)
    h_array = np.ones(len(wavelengths))*total_h
    for i in range(nu_layers):
        layer_data = []
        for wl in range(len(wavelengths)):
            layer_data = np.append(layer_data,spec_list[wl*nu_layers + i])

        ax1 = fig.add_subplot(nu_layers,1,i+1)
        ax1.plot(wavelengths, layer_data)
        if i == nu_layers-1:
            ax1.set_xlabel('Wavelength (nm)')
            ax1.set_ylabel('Total')
        else:
            ax1.set_xticklabels( () )
        if spectra_name == 'Lay_Absorb':
            if i == 0: ax1.set_ylabel('Top Layer')
            if i == nu_layers-2: ax1.set_ylabel('Bottom Layer')
            suptitle_w_params = 'Absorptance in each layer\n'+params_2_print
            plt.suptitle(suptitle_w_params)
            lay_spec_name = 'Lay_Absorb'
            if i == nu_layers-1: 
                ax1.set_ylabel('Total')
                lay_spec_name = 'Absorptance'
        elif spectra_name == 'Lay_Trans':
            if i == 0: ax1.set_ylabel('Top Layer')
            if i == nu_layers-2: ax1.set_ylabel('Bottom Layer')
            suptitle_w_params = 'Transmittance in each layer\n'+params_2_print
            plt.suptitle(suptitle_w_params)
            lay_spec_name = 'Lay_Trans'
            if i == nu_layers-1: 
                ax1.set_ylabel('Total')
                lay_spec_name = 'Transmittance'
        elif spectra_name == 'Lay_Reflec':
            if i == 0: ax1.set_ylabel('Top Layer')
            if i == nu_layers-3: ax1.set_ylabel('Bottom Layer')
            if i == nu_layers-2: ax1.set_ylabel('Substrate')
            suptitle_w_params = 'Reflectance in each layer\n'+params_2_print
            plt.suptitle(suptitle_w_params)
            lay_spec_name = 'Lay_Reflec'
            if i == nu_layers-1: 
                ax1.set_ylabel('Total')
                lay_spec_name = 'Reflectance'
        av_array = zip(wavelengths, layer_data, h_array)
        plt.xlim((wavelengths[0], wavelengths[-1]))
        plt.ylim((0, 1))
        if i != nu_layers-1: 
            np.savetxt('%(s)s_%(i)i_stack%(bon)i.txt'% {'s' : lay_spec_name, 'i' : i, 'bon' : stack_label}, av_array, fmt = '%18.11f')
        else:
            np.savetxt('%(s)s_stack%(bon)i.txt'% {'s' : lay_spec_name, 'bon' : stack_label}, av_array, fmt = '%18.11f')

        plt.savefig('%(s)s_stack%(bon)i'% {'s' : spectra_name, 'bon' : stack_label})

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











# def omega_plot_concentration(st, BM_min,BM_max):
#     fig = plt.figure(num=None, figsize=(9, 4), dpi=80, facecolor='w', edgecolor='k')
#     ax1 = fig.add_subplot(1,1,1)
#     format_st     = '%04d' % st

#     file_name   = '../000-simo_template/omega_pol.txt'#_st%s.txt' % format_st
#     # file_name   = '../000-simo_template/omega_pol_st%s.txt' % format_st
#     wavelengths = np.genfromtxt(file_name, usecols=(2))
#     num_BMs     = np.genfromtxt(file_name, usecols=(1))
#     start_col   = 5

#     for i in range(BM_min,BM_max+1):
#         e_conc = np.genfromtxt(file_name, usecols=(start_col+i),   invalid_raise=False)
#         # print e_conc

#         trim_wls = wavelengths[0:len(e_conc)]
#         ax1.plot(trim_wls, e_conc, 'ro')
#         ax1.set_xlabel('Wavelength (nm)')
#         ax1.set_ylabel(r'$E_{cyl} / E_{cell}$')
#         ax1.set_xlim((wavelengths[0], wavelengths[-1]))
#         # ax1.set_xlim((400, 700))

#     plt.savefig('../000-simo_template/Energy_c')

# if __name__ == "__main__":
#     import sys
#     omega_plot_concentration(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]))


# def average_spec(spec_name,av_spec_name,num_wl,num_h):
#     data      = np.loadtxt('%s.txt' % spec_name)
#     av_wl     = []
#     av_spec   = []
#     av_h      = []
#     for i in np.linspace(0,num_wl-1,num_wl):
#         av_wl   = np.append(av_wl,data[i*num_h,0])
#         av_tmp  = np.mean(data[i*num_h:(i+1)*num_h,1])
#         av_spec = np.append(av_spec,av_tmp)
#         av_tmp  = np.mean(data[i*num_h:(i+1)*num_h,2])
#         av_h    = np.append(av_h,av_tmp)
#     # Save averages to file
#     av_array = zip(av_wl, av_spec, av_h)
#     np.savetxt('%s.txt' % av_spec_name, av_array, fmt = '%18.12f')


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


# def spectra_h(wavelengths,h_wl, h, eta_calc,i):
#     fig = plt.figure(num=None, figsize=(9, 4), dpi=80, facecolor='w', edgecolor='k')
#     ax1 = fig.add_subplot(1,1,1)
#     ax1.plot(wavelengths, h_wl)
#     ax1.set_xlabel('Wavelength (nm)')
#     ax1.set_ylabel('Absorptance')
#     plt.axis([wavelengths[0], wavelengths[-1], 0, 1])
#     tmp1 = 'h = %(h)d (nm), '% {
#     'h'         : h, }
#     tmp7 = r'$\eta$ = %(eta_calc)6.2f'% {
#     'eta_calc'  : eta_calc*100, }
#     tmp8 = ' %'

#     imp_facts = tmp1 + tmp7 + tmp8
#     plt.suptitle(imp_facts)
#     if i < 10:
#         plt.savefig('Animated/000%d' % i)
#     elif i < 100:
#         plt.savefig('Animated/00%d' % i)
#     elif i < 1000:
#         plt.savefig('Animated/0%d' % i)
#     else:
#         plt.savefig('Animated/%d' % i)


# def efficiency(wavelengths,i_spec,spec):
#     #  Total solar irradiance - integral of I(lambda) from 310nm-4000nm
#     #  intergral done in Mathematica (OtherCode/Silicon_ASTM/ASTMG173.nb)
#     tot_irradiance = 900.084
#     bandgap_wl     = wavelengths[-1] #have as property of material.
#     expression     = i_spec*spec*wavelengths
#     integral_tmp   = np.trapz(expression, x=wavelengths)
#     Efficiency     = integral_tmp/(bandgap_wl*tot_irradiance)
#     return Efficiency

# def irradiance(Absorb, W_Absorb, Trans, W_Trans, Reflec, W_Reflec, radius1, radius2=0,
#         period=0, ff=0):
#     data         = np.loadtxt('%s.txt' % Absorb)
#     wavelengths  = data[:,0]
#     spec         = data[:,1]
#     h_av         = data[:,2] 
#     i_data       = np.loadtxt('%s.txt' % Irrad_spec_file)
#     i_spec       = np.interp(wavelengths, i_data[:,0], i_data[:,2])

#     # call efficiency function 
#     Efficiency = efficiency(wavelengths,i_spec,spec)
#     np.savetxt('Efficiency.txt', [Efficiency, radius1, radius2,
#         period, ff], fmt = '%12.8f')

#     # weighted absorptance
#     weighting(i_spec, wavelengths, spec, h_av, W_Absorb)
#     # weighted transmittance
#     data         = np.loadtxt('%s.txt' % Trans)
#     spec         = data[:,1] 
#     weighting(i_spec, wavelengths, spec, h_av, W_Trans)
#     # weighted reflectance
#     data         = np.loadtxt('%s.txt' % Reflec)
#     spec         = data[:,1] 
#     weighting(i_spec, wavelengths, spec, h_av, W_Reflec)
#     return Efficiency


# def weighting(weight_by, wavelengths, spectrum, h_av, w_spectum):
#     # spec_data      = np.loadtxt('%s.txt' % spectrum)
#     # spec           = spec_data[:,1]
#     weighted       = spectrum*weight_by/weight_by.max()
#     weigthed_array = zip(wavelengths, weighted, h_av)
#     np.savetxt('%s.txt' % w_spectum, weigthed_array, fmt = '%18.11f')


# def tra_plot(spectra_name, spec_list, layer, light, max_num_BMs, max_order_PWs, Efficiency):
#     fig = plt.figure(num=None, figsize=(9, 12), dpi=80, facecolor='w', edgecolor='k')
#     for i in range(len(spec_list)):
#         ax1 = fig.add_subplot(3,1,i+1)
#         spec_name = spec_list.pop(0)
#         s_data  = np.loadtxt('%s.txt' % spec_name)
#         wavelengths = s_data[:,0]
#         spectrum = s_data[:,1]
#         ax1.plot(wavelengths, spectrum)
#         ax1.set_xlabel('Wavelength (nm)')
#         ax1.set_ylabel(spec_name)
#         plt.axis([wavelengths[0], wavelengths[-1], 0.0, 1])
#     tmp1 = 'd = %(period)d, a1 = %(radius)d '% {
#     'period'        : layer.period, 'radius' : layer.radius1,}
#     tmp10 ='ff = %(ff)4.2f, \nh_1 = %(h_one)d, '% {
#     'ff'            : layer.ff, 
#     'h_one'         : layer.height_nm, }
#     tmp2 = r'$k_\parallel d$ = '
#     k_pll = light.k_pll * layer.period
#     tmp3 = '(%(kx)1.4f, %(ky)1.4f), '% {'kx' : k_pll[0], 'ky' : k_pll[1]}
#     tmp6 = 'max_BMs = %(max_num_BMs)d, max_PW_order = %(max_order_PWs)d'% {
#     'max_num_BMs'   : max_num_BMs,
#     'max_order_PWs' : max_order_PWs, }
#     tmp7 = '\n' r'$\eta$ = %(Efficiency)6.2f'% {'Efficiency' : Efficiency*100, }
#     tmp8 = ' %'
#     tmp_end = tmp10 + tmp2 + tmp3 + tmp6 + tmp7 + tmp8

#     if layer.radius2 + layer.radius3 + layer.radius4 == 0:
#         imp_facts = tmp1 + tmp_end
#     elif layer.radius3 + layer.radius4 == 0:
#         tmp11 = 'a2 = %(radius)d '% {'radius' : layer.radius2,}
#         imp_facts = tmp1 +  tmp11 + tmp_end
#     elif layer.radius4 == 0:
#         tmp11 = 'a2 = %(radius)d '% {'radius' : layer.radius2,}
#         tmp12 = 'a3 = %(radius)d '% {'radius' : layer.radius3,}
#         imp_facts = tmp1 + tmp11 + tmp12 + tmp_end
#     else:
#         tmp11 = 'a2 = %(radius)d '% {'radius' : layer.radius2,}
#         tmp12 = 'a3 = %(radius)d '% {'radius' : layer.radius3,}
#         tmp13 = 'a4 = %(radius)d '% {'radius' : layer.radius4,}
#         imp_facts = tmp1 + tmp11 + tmp12 + tmp13 + tmp_end

#     plt.suptitle(imp_facts)
#     plt.savefig(spectra_name)


# def overlay_plot(spectra_name, spec_list, layer, light, max_num_BMs, max_order_PWs, Efficiency):
#     fig = plt.figure(num=None, figsize=(9, 12), dpi=80, facecolor='w', edgecolor='k')
#     ax1 = fig.add_subplot(3,1,2)
#     line_list = ['k-','r-.','b--','g:','m,']
#     for i in range(len(spec_list)):
#         spec_name   = spec_list.pop(0)
#         line_name   = line_list.pop(0)
#         s_data      = np.loadtxt('%s.txt' % spec_name)
#         wavelengths = s_data[:,0]
#         spectrum    = s_data[:,1]
#         ax1.plot(wavelengths, spectrum, line_name)
#         ax1.set_xlabel('Wavelength (nm)')
#         ax1.set_ylabel(spec_name)
#         plt.axis([wavelengths[0], wavelengths[-1], 0, 1])
#     plt.legend( ('posxy 100','posxy 95','posxy 98') )
#     tmp1 = 'd = %(period)d, a1 = %(radius)d '% {
#     'period'        : layer.period, 'radius' : layer.radius1,}
#     tmp10 ='ff = %(ff)4.2f, \nh_1 = %(h_one)d, '% {
#     'ff'            : layer.ff, 
#     'h_one'         : layer.height_nm, }
#     tmp2 = r'$k_\parallel d$ = '
#     k_pll = light.k_pll * layer.period
#     tmp3 = '(%(kx)1.4f, %(ky)1.4f), '% {'kx' : k_pll[0], 'ky' : k_pll[1]}
#     tmp4 = r'$\phi$ = '
#     tmp5 = '%(phi)6.2f, '% {'phi' : light.phi, }
#     tmp6 = 'max_BMs = %(max_num_BMs)d, max_PW_order = %(max_order_PWs)d'% {
#     'max_num_BMs'   : max_num_BMs,
#     'max_order_PWs' : max_order_PWs, }
#     tmp7 = '\n' r'$\eta$ = %(Efficiency)6.2f'% {'Efficiency' : Efficiency*100, }
#     tmp8 = ' %'
#     tmp_end = tmp10 + tmp2 + tmp3 + tmp6 + tmp7 + tmp8

#     if layer.radius2 + layer.radius3 + layer.radius4 == 0:
#         imp_facts = tmp1 + tmp_end
#     elif layer.radius3 + layer.radius4 == 0:
#         tmp11 = 'a2 = %(radius)d '% {'radius' : layer.radius2,}
#         imp_facts = tmp1 +  tmp11 + tmp_end
#     elif layer.radius4 == 0:
#         tmp11 = 'a2 = %(radius)d '% {'radius' : layer.radius2,}
#         tmp12 = 'a3 = %(radius)d '% {'radius' : layer.radius3,}
#         imp_facts = tmp1 + tmp11 + tmp12 + tmp_end
#     else:
#         tmp11 = 'a2 = %(radius)d '% {'radius' : layer.radius2,}
#         tmp12 = 'a3 = %(radius)d '% {'radius' : layer.radius3,}
#         tmp13 = 'a4 = %(radius)d '% {'radius' : layer.radius4,}
#         imp_facts = tmp1 + tmp11 + tmp12 + tmp13 + tmp_end

#     plt.suptitle(imp_facts)
#     plt.savefig(spectra_name)


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

#     tmp1 = 'd = %(period)d, a1 = %(radius)d '% {
#     'period'        : layer.period, 'radius' : layer.radius1,}
#     tmp10 ='ff = %(ff)4.2f, \nh_1 = %(h_one)d, '% {
#     'ff'            : layer.ff, 
#     'h_one'         : layer.height_nm, }
#     tmp2 = r'$k_\parallel d$ = '
#     k_pll = light.k_pll * layer.period
#     tmp3 = '(%(kx)1.4f, %(ky)1.4f), '% {'kx' : k_pll[0], 'ky' : k_pll[1]}
#     tmp6 = 'max_BMs = %(max_num_BMs)d, max_PW_order = %(max_order_PWs)d'% {
#     'max_num_BMs'   : max_num_BMs,
#     'max_order_PWs' : max_order_PWs, }
#     tmp7 = '\n' r'$\eta$ = %(Efficiency)6.2f'% {'Efficiency' : Efficiency*100, }
#     tmp8 = ' %'
#     tmp_end = tmp10 + tmp2 + tmp3 + tmp6 + tmp7 + tmp8

#     if layer.radius2 + layer.radius3 + layer.radius4 == 0:
#         imp_facts = tmp1 + tmp_end
#     elif layer.radius3 + layer.radius4 == 0:
#         tmp11 = 'a2 = %(radius)d '% {'radius' : layer.radius2,}
#         imp_facts = tmp1 +  tmp11 + tmp_end
#     elif layer.radius4 == 0:
#         tmp11 = 'a2 = %(radius)d '% {'radius' : layer.radius2,}
#         tmp12 = 'a3 = %(radius)d '% {'radius' : layer.radius3,}
#         imp_facts = tmp1 + tmp11 + tmp12 + tmp_end
#     else:
#         tmp11 = 'a2 = %(radius)d '% {'radius' : layer.radius2,}
#         tmp12 = 'a3 = %(radius)d '% {'radius' : layer.radius3,}
#         tmp13 = 'a4 = %(radius)d '% {'radius' : layer.radius4,}
#         imp_facts = tmp1 + tmp11 + tmp12 + tmp13 + tmp_end

#     plt.suptitle(imp_facts)
#     plt.savefig(name_out)

# # Plot vs Energy
#     fig  = plt.figure(num=None, figsize=(12, 9), dpi=80, facecolor='w', edgecolor='k')
#     ax1  = fig.add_subplot(111)
#     h  = 6.626068e-34;
#     c  = 299792458;
#     eV = 1.60217646e-19;
#     e_data  = (h*c/(wl_data*1e-9))/eV
#     e_1    = e_data[0]
#     e_2    = e_data[-1]
#     d_e    = (e_2-e_1)/(500-1)

#     # Interpolate the data
#     int_e,int_h = np.mgrid[e_1:e_2+d_e:d_e, h_1:h_2+h_int:h_int]
#     int_Eff_e   = griddata(e_data, h_data, s_data, int_e, int_h)
#     int_Eff_e_t = np.fliplr(int_Eff_e.T)

#     CS = plt.imshow(int_Eff_e_t, cmap=plt.cm.hot, norm=None, aspect='auto', interpolation=None,
#       alpha=None, vmin=0, vmax=s_data.max(), origin='lower', extent=(e_2, e_1, h_1, h_2))
#     ax1.set_xlabel('Energy (eV)')
#     ax1.set_ylabel('Height (nm)')
#     cbar = plt.colorbar(extend='neither',ticks=tick_array)
#     cbar.set_ticklabels(tick_array)
#     cbar.ax.set_ylabel('Absorptance')
#     plt.suptitle(imp_facts)
#     plt.savefig('%s_energy' % name_out)

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


# def omega_plot(complete_st, layer, light, max_num_BMs, max_order_PWs, Efficiency):
#     fig = plt.figure(num=None, figsize=(9, 12), dpi=80, facecolor='w', edgecolor='k')
#     nu_layers = len(complete_st)
#     for i in range(0, nu_layers):
#         #reverse order so top layer gets plotted up top
#         st = nu_layers-1-i
#         format_st     = '%04d' % st
#         ax1 = fig.add_subplot(nu_layers,1,i+1)
#         # ax2 = fig.add_subplot(nu_layers,1,i+1)
#         if isinstance(complete_st[nu_layers-1-i],objects.ThinFilm):
#             file_name   = 'beta_st%s.txt' % format_st
#             wavelengths = np.genfromtxt(file_name, usecols=(1))
#             num_BMs     = np.genfromtxt(file_name, usecols=(0))
#             start_col   = 2
#         else:
#             file_name   = 'omega_st%s.txt' % format_st
#             wavelengths = np.genfromtxt(file_name, usecols=(2))
#             num_BMs     = np.genfromtxt(file_name, usecols=(1))
#             start_col   = 5

#         count = 1
#         for i in range(len(num_BMs)):
#             prop = []
#             prop_im = []
#             re = np.genfromtxt(file_name, usecols=(start_col+2*i),   invalid_raise=False)
#             im = np.genfromtxt(file_name, usecols=(start_col+2*i+1), invalid_raise=False)
#             nu_real_modes = len(re)
#             for j in range(nu_real_modes):
#                 if re[j] > im[j]:
#                     prop.append(re[j])
#                     prop_im.append(im[j])
#             count +=1
#             if len(prop) == 0:
#                 for j in range(nu_real_modes):
#                     prop.append(re[j])
#                     prop_im.append(im[j])
#                 count +=1

#             trim_wls = wavelengths[0:len(prop)]
#             ax1.plot(trim_wls, prop, 'ro')
#             ax1.set_xlabel('Wavelength (nm)')
#             ax1.set_ylabel(r'Re(k$_z$)')
#             # plt.xlim((wavelengths[0], wavelengths[-1]))
#             # ax2.plot(trim_wls, prop_im, 'bo')
#             # ax2.set_xlabel('Wavelength (nm)')
#             # ax2.set_ylabel(r'Im(k$_z$)')
#         ax1.set_xlim((wavelengths[0], wavelengths[-1]))
#         # ax2.set_xlim((wavelengths[0], wavelengths[-1]))

#     tmp1 = 'd = %(period)d, a1 = %(radius)d '% {
#     'period'        : layer.period, 'radius' : layer.radius1,}
#     tmp10 ='ff = %(ff)4.2f, \nh_1 = %(h_one)d,  '% {
#     'ff'            : layer.ff, 
#     'h_one'         : layer.height_nm, }
#     tmp2 = r'$k_\parallel d$ = '
#     k_pll = light.k_pll * layer.period
#     tmp3 = '(%(kx)1.4f, %(ky)1.4f), '% {'kx' : k_pll[0], 'ky' : k_pll[1]}
#     tmp6 = 'max_BMs = %(max_num_BMs)d, max_PW_order = %(max_order_PWs)d'% {
#     'max_num_BMs'   : max_num_BMs,
#     'max_order_PWs' : max_order_PWs, }
#     tmp7 = '\n' r'$\eta$ = %(Efficiency)6.2f'% {'Efficiency' : Efficiency*100, }
#     tmp8 = ' %'
#     tmp_end = tmp10 + tmp2 + tmp3 + tmp6 + tmp7 + tmp8

#     if layer.radius2 + layer.radius3 + layer.radius4 == 0:
#         imp_facts = tmp1 + tmp_end
#     elif layer.radius3 + layer.radius4 == 0:
#         tmp11 = 'a2 = %(radius)d '% {'radius' : layer.radius2,}
#         imp_facts = tmp1 +  tmp11 + tmp_end
#     elif layer.radius4 == 0:
#         tmp11 = 'a2 = %(radius)d '% {'radius' : layer.radius2,}
#         tmp12 = 'a3 = %(radius)d '% {'radius' : layer.radius3,}
#         imp_facts = tmp1 + tmp11 + tmp12 + tmp_end
#     else:
#         tmp11 = 'a2 = %(radius)d '% {'radius' : layer.radius2,}
#         tmp12 = 'a3 = %(radius)d '% {'radius' : layer.radius3,}
#         tmp13 = 'a4 = %(radius)d '% {'radius' : layer.radius4,}
#         imp_facts = tmp1 + tmp11 + tmp12 + tmp13 + tmp_end

#     plt.suptitle(imp_facts)
#     plt.savefig('Disp_Diagram')







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