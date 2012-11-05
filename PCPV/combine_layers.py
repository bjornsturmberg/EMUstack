

import numpy as np
from plotting      import layers_plot

def load_scat_mat(name, st, p):
    # reshape matrices to be consistent with pcpv.exe output
    format_title = '%04d' % st
    format_p     = '%04d' % p

    file_name = "st%(st)s_wl%(wl)s_%(mat_name)s.txt" % {
        'st' : format_title, 'wl' : format_p, 'mat_name' : name }
    data   = np.loadtxt(file_name)
    num_1  = max(data[:,0])
    num_2  = max(data[:,1])
    matrix = np.mat(data[:,2] + data[:,3]*(0+1j))
    matrix = np.reshape(matrix, (num_2, num_1))
    return matrix



def calc_tra_layers(solar_cell, rnet_list, inverted_t21_list, P_list,
            f3_minus, t21s, r12s, nu_intfaces, wavelength, p, name=''):
    """ 
    f3_minus is downward incident field
    f2_minus is downward field in layer
    f1_minus is downward outgoing field
    f3_plus  is upward outgoing field   
    """
    t_list = []
    r_list = []
    a_list = []

    # R, T matrices, r, t scalar reflectance, transmittance. i is absorptive layer
    # note that rnet and other lists are out of sync by 1 as rnet has extra entry.
    for layer_nu in range(nu_intfaces-1,0,-1): # layer_nu is position layer in solar_cell
        print layer_nu
        layer_in = layer_nu + 1
        if isinstance(solar_cell[layer_in],objects.NanoStruct):
            prop_in = 2
            print 'in', solar_cell[layer_in].label_nu
        elif isinstance(solar_cell[layer_in],objects.ThinFilm):
            prop_in = solar_cell[layer_in].nu_prop_ords
            print 'in_TF', solar_cell[layer_in].label_nu
            print 'prop_in', prop_in
        layer_out = layer_nu - 1
        if isinstance(solar_cell[layer_out],objects.NanoStruct):
            prop_out = 2
            print 'out', solar_cell[layer_out].label_nu
        elif isinstance(solar_cell[layer_out],objects.ThinFilm):
            prop_out = solar_cell[layer_out].nu_prop_ords
            print 'out_TF', solar_cell[layer_out].label_nu
            print 'prop_out', prop_out

# CONJ OR ABS()**2?????????????
        # d1 = abs(f3_minus[0:prop_in])
        d1 = f3_minus[0:prop_in]
        delta = np.linalg.norm(d1)**2
        if layer_nu == nu_intfaces-1: delta_in = delta # top layer

        rnet_top = rnet_list[layer_nu] #rnet at top of layer layer_nu
        f3_plus  = rnet_top*f3_minus
        P        = P_list[layer_nu-1]
        # inverted_t21_list has # el of layers. therefore python index layer_nu-1
        # I_mat          = np.matrix(np.eye(len(P)),dtype='D')
        # to_invert      = (I_mat)# - r12s[layer_nu]*P*rnet_list[layer_nu-1]*P)
        # inverted_t21   = np.linalg.solve(to_invert,t21s[layer_nu])
        # f2_minus = np.matrix(inverted_t21)*f3_minus 
        # f2_minus = t21s[layer_nu]*f3_minus
        f2_minus = np.matrix(inverted_t21_list[layer_nu-1])*f3_minus 

        f1 = f2_minus[0:2]
        f  = np.linalg.norm(f1)**2
        # print t21s[layer_nu]
        # print f3_minus
        # print f2_minus[0:10]

        f3_minus = P*f2_minus/prop_out
        f1_minus = t21s[layer_nu-1]*f3_minus

        # t1 = abs(f1_minus[0:prop_out])
        t1 = f1_minus[0:prop_out]
        # print f1_minus
        # print t1
        t  = np.linalg.norm(t1)**2
        # r1 = abs(f3_plus[0:prop_in])
        r1 = f3_plus[0:prop_in]
        r  = np.linalg.norm(r1)**2
        a  = delta - r - t
        t_list.append(t)
        r_list.append(r)
        a_list.append(a)
        print 'delta', delta
        print delta - r
        print 'f', f

    a_tot = delta_in - t_list[-1] - r_list[0]
    save_tra_layers(t_list, 'T_Layers'+ name, wavelength, p)
    save_tra_layers(r_list, 'R_Layers'+ name, wavelength, p)
    save_tra_layers(a_list, 'A_Layers'+ name, wavelength, p)
    save_tra_layers([t_list[-1]], 'T_Lambda'+ name, wavelength, p)
    save_tra_layers([r_list[0]],  'R_Lambda'+ name, wavelength, p)
    save_tra_layers([a_tot],      'A_Lambda'+ name, wavelength, p)

    print 't', t_list
    print 'r', r_list
    print 'a', a_list

    print 't', np.sum(t_list)
    print 'r', np.sum(r_list)
    print 'a', np.sum(a_list)

    print 't', [t_list[-1]]
    print 'r', [r_list[0]]
    print 'a', [a_tot]






def net_scat_mats(solar_cell, wavelengths):
    nu_intfaces    = 2*(len(solar_cell)-1)
    nu_TFs         = len(solar_cell)-2 # adujsted by 1 for python index start at 0
    # Transmittance = []
    # Reflectance   = []
    # Absorptance   = []
    # T_save        = []
    # R_save        = []
    # A_save        = []
    neq_PW         = solar_cell[0].nu_tot_ords # assumes incident from homogeneous film
    PW_pols        = 2*neq_PW
    num_prop_air   = solar_cell[0].nu_tot_ords#solar_cell[-1].num_prop_air
    num_prop_in    = solar_cell[-1].num_prop_TF
    num_prop_out   = solar_cell[0].num_prop_TF
    zero_mat       = np.matrix(np.zeros((PW_pols, 1),float))
    I_air          = np.matrix(np.eye(PW_pols),dtype='D')
    inc            = solar_cell[-1].set_ord_in 
    out            = solar_cell[0].set_ord_out

    t_list = []
    r_list = []
    a_list = []

    for p in range(len(wavelengths)):
        p += 1

        """ Calculate net scattering matrices starting at the bottom
            1 is infintesimal air layer
            2 is medium in layer (symetric as air on each side)
            (r)t12 and (r)tnet lists run from bottom to top!
        """
        r12_list = []
        r21_list = []
        t12_list = []
        t21_list = []
        P_list   = []
        for st1 in solar_cell:
            r12_list.append(load_scat_mat('R12', st1.label_nu, p).T)
            r21_list.append(load_scat_mat('R21', st1.label_nu, p).T)
            # potential to save on one transpose
            t12_list.append(load_scat_mat('T12', st1.label_nu, p).T)
            t21_list.append(t12_list[-1].T)
            # t21s.append(tmat_through_air(st2, st1, p))

# initiate (r)tnet as substrate top interface
        tnet      = t12_list[0]
        rnet      = r12_list[0]
        tnet_list = []
        rnet_list = []
        tnet_list.append(tnet)
        rnet_list.append(rnet)
        T_hat_list = []
        R_hat_list = []

        inv_t21_list   = []
        inv_t12_list   = []
        for i in range(nu_TFs):
            i = i + 1 # start in bottom air layer not substrate
# through air layer at bottom of TF
            to_invert      = (I_air - r12_list[i]*rnet)
            inverted_t21   = np.linalg.solve(to_invert,t21_list[i])
            tnet           = tnet*inverted_t21
            rnet           = r21_list[i] + t12_list[i]*rnet*inverted_t21
            inv_t21_list.append(inverted_t21)
            tnet_list.append(tnet)
            rnet_list.append(rnet)
# through TF layer
            P              = load_scat_mat('P', solar_cell[i].label_nu, p).T
            I_TF           = np.matrix(np.eye(len(P)),dtype='D')
            to_invert      = (I_TF - r21_list[i]*P*rnet*P)
            inverted_t12   = np.linalg.solve(to_invert,t12_list[i])
            P_inverted_t12 = P*inverted_t12
            tnet           = tnet*P_inverted_t12
            rnet           = r12_list[i] + t21_list[i]*P*rnet*P_inverted_t12

            to_invert_hat      = (I_TF - r21_list[i]*P*r21_list[i]*P)
            inverted_t12_hat   = np.linalg.solve(to_invert_hat,t12_list[i])
            P_inverted_t12_hat = P*inverted_t12_hat
            T_hat              = t21_list[i]*P_inverted_t12_hat
            R_hat              = r12_list[i] + t21_list[i]*P*r21_list[i]*P_inverted_t12_hat

            P_list.append(P)
            inv_t12_list.append(inverted_t12)
            tnet_list.append(tnet)
            rnet_list.append(rnet)
            T_hat_list.append(T_hat)
            R_hat_list.append(R_hat)

# into top semi-infinite medium
        to_invert = (I_air - r12_list[-1]*rnet)
        inverted_t21 = np.linalg.solve(to_invert,t21_list[-1])
        tnet = tnet*inverted_t21
        rnet = r21_list[-1] + t12_list[-1]*rnet*inverted_t21
        inv_t21_list.append(inverted_t21)
        tnet_list.append(tnet)
        rnet_list.append(rnet)


        """ Calculate field expansions for all layers (including air) starting at top
            These correctly describe energy fluxes in and out of each TF layer,
            from which layer absorptances can be calculated.
            Ordering is now top to bottom (inverse of above)! ie f1 is superstrate (top)
        """
        fminus_list = []
        fplus_list  = []
        fminus_eng  = []
        fplus_eng   = []

#   incoming from semi-inf
        d_minus  = zero_mat
        # for TM polarisation
        d_minus[inc,0] = float(1.0)
        # for TE polarisation
        # f1_minus[neq_PW+inc,0] = float(1.0)
        # for Right circular polarisation
        # calc_tra_layers(rnet_list, inverted_t21_list, P_list,
        #     f3_minus, t21s, nu_intfaces, wavelengths[p-1], p, '_R')
        # for Left circular polarisation
        # calc_tra_layers(rnet_list, inverted_t21_list, P_list,
        #     f3_minus, t21s, nu_intfaces, wavelengths[p-1], p, '_L')
        # for Circular Dichroism

        energy1 = np.linalg.norm(d_minus[0:num_prop_in])**2
        energy2 = np.linalg.norm(d_minus[neq_PW:neq_PW+num_prop_in])**2
        fminus_eng.append(energy1 + energy2)

#   up into semi-inf off top air gap
        d_plus  = rnet_list[-1]*d_minus
        fminus_list.append(d_minus)
        fplus_list.append(d_plus)
        energy1 = np.linalg.norm(d_plus[0:num_prop_in])**2
        energy2 = np.linalg.norm(d_plus[neq_PW:neq_PW+num_prop_in])**2
        fplus_eng.append(energy1 + energy2)

#   incoming from semi-inf into top air gap
        f1_minus = inv_t21_list[-1]*d_minus
        fminus_list.append(f1_minus)
        # energy1  = np.linalg.norm(f1_minus[0:num_prop_air])**2
        # energy2  = np.linalg.norm(f1_minus[neq_PW:neq_PW+num_prop_air])**2
        # fminus_eng.append(energy1 + energy2)
        fminus_eng.append(np.linalg.norm(f1_minus)**2)


        for i in range(nu_TFs):
#   d
            f1_plus = rnet_list[-2*i-2]*f1_minus
            fplus_list.append(f1_plus)
            # energy1 = np.linalg.norm(f1_plus[0:num_prop_air])**2
            # energy2 = np.linalg.norm(f1_plus[neq_PW:neq_PW+num_prop_air])**2
            # fplus_eng.append(energy1 + energy2)
            fplus_eng.append(np.linalg.norm(f1_plus)**2)

            f2_minus = inv_t12_list[-i-1]*f1_minus

            f2_plus  = rnet_list[-2*i-3]*P_list[-i-1]*f2_minus

            f1_minus = inv_t21_list[-i-2]*P_list[-i-1]*f2_minus
            fminus_list.append(f1_minus)
            # energy1  = np.linalg.norm(f1_minus[0:num_prop_air])**2
            # energy2  = np.linalg.norm(f1_minus[neq_PW:neq_PW+num_prop_air])**2
            # fminus_eng.append(energy1 + energy2)
            fminus_eng.append(np.linalg.norm(f1_minus)**2)
          
# bottom air to semi-inf substrate
        f1_plus = rnet_list[0]*f1_minus
        fplus_list.append(f1_plus)
        # energy1 = np.linalg.norm(f1_plus[0:num_prop_air])**2
        # energy2 = np.linalg.norm(f1_plus[neq_PW:neq_PW+num_prop_air])**2
        # fplus_eng.append(energy1 + energy2)
        fplus_eng.append(np.linalg.norm(f1_plus)**2)

        f2_minus = t12_list[0]*f1_minus
        fminus_list.append(f2_minus)
        # energy1 = np.linalg.norm(f2_minus[0:num_prop_out])**2
        # energy2 = np.linalg.norm(f2_minus[neq_PW:neq_PW+num_prop_out])**2
        # fminus_eng.append(energy1 + energy2)
        fminus_eng.append(np.linalg.norm(f2_minus)**2)



        t_layer = fminus_eng[1]/fminus_eng[0]
        r_layer = fminus_eng[0] - t_layer
        t_list.append(t_layer)
        r_list.append(r_layer)
        print 't', t_layer
        print 'r', r_layer
        print ''
        for i in range(1,nu_TFs+1):
            t_layer = fminus_eng[i+1]#/fminus_eng[i]
            r_layer = fplus_eng[i]/fminus_eng[i]
            # t_layer = abs(T_hat_list[i-1])**2# fminus_eng[i+1]/fminus_eng[i]
            # r_layer = abs(R_hat_list[i-1])**2#fplus_eng[i]/fminus_eng[i]
            a_layer = fminus_eng[i] + fplus_eng[i+1] - fminus_eng[i+1] - fplus_eng[i]
            t_list.append(t_layer)
            r_list.append(r_layer)
            a_list.append(a_layer)
            print 't', t_layer
            print 'r', r_layer
            print 'a', a_layer
        t_layer = fminus_eng[-1]#/fminus_eng[-2]
        r_layer = fplus_eng[-1]/fminus_eng[-1]
        t_list.append(t_layer)
        r_list.append(r_layer)
        print 't', t_layer
        print 'r', r_layer
        print ''

        print 'final t', fminus_eng[-1]/fminus_eng[0]
        print 'final r', fplus_eng[0]/fminus_eng[0]

        final_t  = fminus_eng[-1]/fminus_eng[0]
        final_r  = fplus_eng[0]/fminus_eng[0]
        abs_tot  = sum(a_list[-nu_TFs:-1]) + a_list[-1]
        t_list.append(final_t)
        r_list.append(final_r)
        a_list.append(abs_tot)
        abs_diff = (fminus_eng[0] - fminus_eng[-1] - fplus_eng[0]) - abs_tot
        a_list.append(abs_diff)
        print 'abs_tot', abs_tot
        print ''        
        print 'abs_diff', abs_diff


    layers_plot('Lay_Trans', t_list, wavelengths)
    layers_plot('Lay_Reflec', r_list, wavelengths)
    layers_plot('Lay_Absorb', a_list, wavelengths)










        # # select elements of Tnet, Rnet matrices to calculate absorption etc.
        # neq_PW   = solar_cell[0].nu_tot_ords 
        # # select_ord_in  = solar_cell[-1].set_ord_in
        # # select_ord_out = solar_cell[0].set_ord_out


        # inc      = 0#solar_cell[0].zero_ord
        # Lambda_t = 0
        # for i in range(num_prop_out):
        #     #TM
        #     # Lambda_t = Lambda_t + abs(tnet[i,neq_PW+inc])**2 + abs(tnet[neq_PW+i,neq_PW+inc])**2
        #     #TE
        #     Lambda_t = Lambda_t + abs(tnet[i,inc])**2 + abs(tnet[neq_PW+i,inc])**2


        # # inc      = solar_cell[-1].zero_ord
        # Lambda_r = 0
        # for i in range(num_prop_in):
        #     #TM
        #     # Lambda_r = Lambda_r + abs(rnet[i,neq_PW+inc])**2 + abs(rnet[neq_PW+i,neq_PW+inc])**2
        #     #TE
        #     Lambda_r = Lambda_r + abs(rnet[i,inc])**2 + abs(rnet[neq_PW+i,inc])**2

        # absorption = 1 - Lambda_r - Lambda_t
        # print 't', Lambda_t
        # print 'r', Lambda_r
        # print 'a', absorption



    # #     Transmittance.append(Lambda_t)
    # #     Reflectance.append(Lambda_r)
    # #     Absorptance.append(absorption)
    # # # T_save.append((wavelengths[p-1], Lambda_t))
    # # # R_save.append((wavelengths[p-1], Lambda_r))
    # # # A_save.append((wavelengths[p-1], absorption))

    # # # np.savetxt('Transmittance.txt', T_save, fmt=['%18.10f','%25.17G'], delimiter='')
    # # # np.savetxt('Reflectance.txt', R_save, fmt=['%18.10f','%25.17G'], delimiter='')
    # # # np.savetxt('Absorptance.txt', A_save, fmt=['%18.10f','%25.17G'], delimiter='')

    # # relative_wave_nu = 120/wavelengths

    # # tmp_plot(relative_wave_nu, Reflectance, 'R_Spec')
    # # tmp_plot(relative_wave_nu, Transmittance, 'T_Spec')
    # # tmp_plot(relative_wave_nu, Absorptance, 'A_Spec')