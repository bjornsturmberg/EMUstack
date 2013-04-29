

import numpy as np
import subprocess
from plotting      import layers_plot

import cat_n_clean
# import matplotlib
# matplotlib.use('pdf')
import matplotlib.pyplot as plt
import objects
import temporary_bullshit as bs



def net_scat_mats(solar_cell, wavelengths, simo_para):

# test that all structures have the same period
    for cell in solar_cell:
        if cell.period == solar_cell[0].period:
            pass
        else:
            raise ValueError, "All layers in a multilayer stack must have the same period!"

    nu_intfaces    = 2*(len(solar_cell)-1)
    nu_TFs         = len(solar_cell)-2 # adujsted by 1 for python index start at 0
    neq_PW         = solar_cell[0].nu_tot_ords # assumes incident from homogeneous film
    PW_pols        = 2*neq_PW
    num_p_air_lst  = solar_cell[-1].num_prop_air
    num_p_in_lst   = solar_cell[-1].num_prop_TF
    num_p_out_lst  = solar_cell[0].num_prop_TF
    zero_vec       = np.matrix(np.zeros((PW_pols, 1),complex))
    I_air          = np.matrix(np.eye(PW_pols),dtype='D')
    inc            = solar_cell[-1].set_ord_in 
    out            = solar_cell[0].set_ord_out

    t_list = []
    r_list = []
    a_list = []

    total_h = 0
    for st in solar_cell[1:-1]:
        total_h += st.height_1

    for p in range(len(wavelengths)):
        num_prop_air = num_p_air_lst[p]
        num_prop_in  = num_p_in_lst[p]
        num_prop_out = num_p_out_lst[p]
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
            r12_list.append(bs.load_scat_mat('R12', st1.label_nu, p))
            r21_list.append(bs.load_scat_mat('R21', st1.label_nu, p))
            # potential to save on one transpose
            t12_list.append(bs.load_scat_mat('T12', st1.label_nu, p))
            t21_list.append(bs.load_scat_mat('T21', st1.label_nu, p))
            #t21_list.append(t12_list[-1].T)

# initiate (r)tnet as substrate top interface
        tnet_list = []
        rnet_list = []
        tnet      = t12_list[0]
        rnet      = r12_list[0]
        tnet_list.append(tnet)
        rnet_list.append(rnet)

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
            P              = bs.load_scat_mat_old('P', solar_cell[i].label_nu, p).T
            I_TF           = np.matrix(np.eye(len(P)),dtype='D')
            to_invert      = (I_TF - r21_list[i]*P*rnet*P)
            inverted_t12   = np.linalg.solve(to_invert,t12_list[i])
            P_inverted_t12 = P*inverted_t12
            tnet           = tnet*P_inverted_t12
            rnet           = r12_list[i] + t21_list[i]*P*rnet*P_inverted_t12

            to_invert_hat      = (I_TF - r21_list[i]*P*r21_list[i]*P)
            inverted_t12_hat   = np.linalg.solve(to_invert_hat,t12_list[i])
            P_inverted_t12_hat = P*inverted_t12_hat

            P_list.append(P)
            inv_t12_list.append(inverted_t12)
            tnet_list.append(tnet)
            rnet_list.append(rnet)

# into top semi-infinite medium
        to_invert    = (I_air - r12_list[-1]*rnet)
        inverted_t21 = np.linalg.solve(to_invert,t21_list[-1])
        tnet         = tnet*inverted_t21
        rnet         = r21_list[-1] + t12_list[-1]*rnet*inverted_t21
        inv_t21_list.append(inverted_t21)
        tnet_list.append(tnet)
        rnet_list.append(rnet)

# print (r)tnet matrices
        # np.savetxt('blah.txt', rnet)

        """ Calculate field expansions for all layers (including air) starting at top
            Ordering is now top to bottom (inverse of above)! ie f1 is superstrate (top)
            Calculate net downward energy flux in each infintesimal air layer & super/substrates
            (see appendix C in Dossou et al. JOSA 2012)
        """

        down_fluxes = []
        up_flux     = []

# Start by composing U matrix which is same for all air layers.
# diagonal with 1 for propagating, i for evanescent TE and -i for evanescent TM plane wave orders

        U_mat = np.matrix(np.zeros((2*PW_pols, 2*PW_pols),complex))
        for i in range(0,num_prop_air):
            U_mat[i,i]                               = 1.0
            U_mat[neq_PW+i,neq_PW+i]                 = 1.0
            U_mat[PW_pols+i,PW_pols+i]               = -1.0
            U_mat[PW_pols+neq_PW+i,PW_pols+neq_PW+i] = -1.0
        for i in range(num_prop_air,neq_PW):
            U_mat[i,PW_pols+i]                       = -1.0j
            U_mat[neq_PW+i,PW_pols+neq_PW+i]         = 1.0j
            U_mat[PW_pols+i,i]                       = 1.0j
            U_mat[PW_pols+neq_PW+i,neq_PW+i]         = -1.0j


#   incoming from semi-inf
        d_minus        = zero_vec
        # for TE polarisation
        d_minus[inc,0] = 1.0
        # for TM polarisation
        # f1_minus[neq_PW+inc,0] = float(1.0)
        # for Right circular polarisation
        # calc_tra_layers(rnet_list, inverted_t21_list, P_list,
        #     f3_minus, t21s, nu_intfaces, wavelengths[p-1], p, '_R')
        # for Left circular polarisation
        # calc_tra_layers(rnet_list, inverted_t21_list, P_list,
        #     f3_minus, t21s, nu_intfaces, wavelengths[p-1], p, '_L')
        # for Circular Dichroism


    # total incoming flux
        flux_TE = np.linalg.norm(d_minus[0:num_prop_in])**2
        flux_TM = np.linalg.norm(d_minus[neq_PW:neq_PW+num_prop_in])**2
        down_fluxes.append(flux_TE + flux_TM)

#   up into semi-inf off top air gap
        d_plus  = rnet_list[-1]*d_minus
    # total reflected flux
        flux_TE = np.linalg.norm(d_plus[0:num_prop_in])**2
        flux_TM = np.linalg.norm(d_plus[neq_PW:neq_PW+num_prop_in])**2
        up_flux.append(flux_TE + flux_TM)

#   incoming from semi-inf into top air gap
        f1_minus = inv_t21_list[-1]*d_minus

        for i in range(nu_TFs):
            f1_plus = rnet_list[-2*i-2]*f1_minus
    # net downward flux in infintesimal air layer
            f_mat   = np.matrix(np.concatenate((f1_minus,f1_plus)))
            f_mat_H = f_mat.conj().transpose()
            flux    = f_mat_H*U_mat*f_mat
            down_fluxes.append(flux)

            f2_minus = inv_t12_list[-i-1]*f1_minus
            f2_plus  = rnet_list[-2*i-3]*P_list[-i-1]*f2_minus

            f1_minus = inv_t21_list[-i-2]*P_list[-i-1]*f2_minus

# bottom air to semi-inf substrate
        f1_plus  = rnet_list[0]*f1_minus

        f2_minus = tnet_list[0]*f1_minus
        flux_TE  = np.linalg.norm(f2_minus[0:num_prop_out])**2
        flux_TM  = np.linalg.norm(f2_minus[neq_PW:neq_PW+num_prop_out])**2
        down_fluxes.append(flux_TE + flux_TM)

# calculate absorptance in each layer
        for i in range(1,len(down_fluxes)-1):
            a_layer = abs(abs(down_fluxes[i])-abs(down_fluxes[i+1]))
            a_list.append(a_layer)
        a_layer = abs(down_fluxes[0]-down_fluxes[-1]-up_flux[0])
        a_list.append(a_layer)
        # sum_abs = sum(a_list[-nu_TFs:])
        # a_list.append(a_layer - sum_abs)

# calculate reflectance in each layer
        for i in range(1,len(up_flux)-1):
            r_layer = abs(abs(up_flux[i])/abs(down_flux[i]))
            r_list.append(r_layer)
        r_layer = abs(up_flux[0]/down_fluxes[0])
        r_list.append(r_layer)

# calculate transmittance in each layer
        for i in range(0,len(down_fluxes)-2):
            t_layer = abs(abs(down_fluxes[i+2])/abs(down_fluxes[i]))
            t_list.append(t_layer)
        t_layer = abs(down_fluxes[-1]/down_fluxes[0])
        t_list.append(t_layer)


# plot scattering matrices as grayscale images
        if simo_para.plot_scat_mats == True:
            im = np.real(abs(rnet_list[-1]))
            plt.matshow(im,cmap=plt.cm.gray)
            cbar = plt.colorbar(extend='neither')
            plt.xlabel('Incoming Orders')
            plt.ylabel('Outgoing Orders')
            plt.suptitle('Net Reflection Scattering Matrix')
            plt.savefig('Rnet_wl%i' % p)
            # for i in range(len(rnet_list)):
            #     im = np.real(rnet_list[i])
            #     plt.matshow(im,cmap=plt.cm.gray)
            #     plt.savefig('0testmat%i' % i)

# plot t,r,a plots each containing results for each layer and total. Also save text t,r,a to files
    layers_plot('Lay_Absorb', a_list, wavelengths, total_h)
    layers_plot('Lay_Trans',  t_list, wavelengths, total_h)
    layers_plot('Lay_Reflec', r_list, wavelengths, total_h)

