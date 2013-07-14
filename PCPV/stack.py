import numpy as np
from plotting import layers_plot
from plotting import k_plot
from objects import Anallo, Simmo
from mode_calcs import r_t_mat
from scipy import sqrt

# import matplotlib
# matplotlib.use('pdf')
# import matplotlib.pyplot as plt


class Stack(object):
    """ Represents a stack of layers evaluated at one frequency.

        This includes the semi-infinite input and output layers.

        INPUTS:

          - `layers` : a tuple of :ThinFilm:s and :NanoStruct:s 
            ordered from top to bottom layer.

          - `heights_nm` : a tuple of the heights of the inside layers,
            i.e., all layers except for the top and bottom. This 
            overrides any heights specified in the :ThinFilm: or
            :NanoStruct: objects.

    """
    def __init__(self, layers, heights_nm = None):
        self.layers = tuple(layers)
        self._heights_nm = heights_nm
        self.period = float(layers[0].structure.period)
        self._check_periods_are_consistent()

    def heights_nm(self):
        if None != self._heights_nm:
            return self._heights_nm
        else:
            return [float(lay.structure.height_nm) for lay in self.layers[1:-1]]

    def heights_norm(self):
        return [h / self.period for h in self.heights_nm()]

    def total_height(self):
        return sum(self.heights())

    def structures(self):
        return (lay.structure for lay in self.layers)

    def calc_R_T_net(self, save_working = True):
        """ Calculate the scattering matrices for the stack as a whole.

            INPUTS:

            - `save_working` : If `True`, then store net reflection
                and transmission matrices at each part of the stack, in
                `self.R_net_list` and `self.T_net_list`, ordered with 
                the reflection and tranmsission of the first/topmost
                finitely thick layer first.

            OUTPUTS:

            - `R_net` : Net reflection matrix

            - `T_net` : Net transmission matrix.
        """
        self._check_periods_are_consistent()

        if save_working:
            self.R_net_list, self.T_net_list = [], []

        #TODO: swap order of layers
        # Reflection and transmission at bottom of structure
        R_net, T_net = r_t_mat(self.layers[1], self.layers[0])[:2]

        lays = self.layers
        for lay, lay_t, h in zip(lays[1:-1], lays[2:], self.heights_norm()):
            if save_working:
                self.R_net_list.insert(0, R_net)
                self.T_net_list.insert(0, T_net)

            # lay (2) is the layer we're in right now
            # lay_t (1) is the layer above it
            # tf = T in forwards direction (down into lay)
            R12, T12, R21, T21 = r_t_mat(lay_t, lay)
            P = lay.prop_fwd(h)
            idm = np.eye(len(P))

            # Matrix that maps vector of forward modes in medium 1
            # at 1-2 interface, to fwd modes in medium 2 at 2-3 interface
            # P * (I - R21 * P * R2net * P)^-1 * T12
            f12 = P * np.linalg.solve(idm - R21 * P * R_net * P, T12)
            T_net = T_net * f12
            R_net = R12 + T21 * P * R_net * f12

        self.R_net, self.T_net = R_net, T_net
        return self.R_net, self.T_net

    def calc_lay_amplitudes(self, incoming_amplitudes):
        """ Return the mode amplitudes at the bottom of each layer.

            OUTPUTS:

            - `f_down_list` : List of vectors of amplitudes of 
                downwards/forward modes

            - `f_up_list` : Amplitudes of upwards/backward modes

            Both lists start with the amplitudes in the first finitely-
            thick layer.

            N.B. this is numerically unstable when T_net is nearly singular.
            This could be overcome by looping through the interfaces once more
            iteratively to find f_down_list and f_up_list.
        """
        # Calculate the amplitudes of transmitted modes using T_net
        f_out = self.T_net * incoming_amplitudes

        # Now work backwards to find what incident field at each interface
        # leads to this superposition of transmitted modes.
        f_down_list = [np.linalg.solve(T_net_l, out) for T_net_l in self.T_net_list]
        
        # And from these, we can use each R_net to find the upward amplitudes
        f_up_list  = [R * f_d for R, f_d in zip(self.R_net_list, f_down_list)]
        
        return f_down_list, f_up_list

    def calc_scat(self, pol = 'TE', incoming_amplitudes = None):
        """ Calculate the transmission and reflection matrices of the stack"""
        # TODO: Switch to calc_R_T_net, which does not use infinitesimal air 
        # layers. This will require rewriting the parts that calculate fluxes
        # through each layer.

        self._check_periods_are_consistent()

        nu_intfaces     = 2*(len(self.layers)-1)
        neq_PW          = self.layers[0].structure.num_pw_per_pol # assumes incident from homogeneous film
        PW_pols         = 2*neq_PW
        I_air           = np.matrix(np.eye(PW_pols),dtype='D')

        """ Calculate net scattering matrices starting at the bottom
            1 is infintesimal air layer
            2 is medium in layer (symmetric as air on each side)
            (r)t12 and (r)tnet lists run from bottom to top!
        """
        r12_list = []
        r21_list = []
        t12_list = []
        t21_list = []
        P_list   = []
        for st1 in self.layers:
            R12, T12, R21, T21 = r_t_mat(st1.air_ref(), st1)
            r12_list.append(R12)
            r21_list.append(R21)
            t12_list.append(T12)
            t21_list.append(T21)

            # Save the reflection matrices to the layers
            # (for easier introspection/testing)
            st1.R12, st1.T12, st1.R21, st1.T21 = R12, T12, R21, T21

    # initiate (r)tnet as substrate top interface
        tnet_list = []
        rnet_list = []
        tnet      = t12_list[0]
        rnet      = r12_list[0]
        tnet_list.append(tnet)
        rnet_list.append(rnet)

        inv_t21_list   = []
        inv_t12_list   = []
        for i in range(1, len(self.layers) - 1):
            lay = self.layers[i]
    # through air layer at bottom of TF
            to_invert      = (I_air - r12_list[i]*rnet)
            inverted_t21   = np.linalg.solve(to_invert,t21_list[i])
            tnet           = tnet*inverted_t21
            rnet           = r21_list[i] + t12_list[i]*rnet*inverted_t21
            inv_t21_list.append(inverted_t21)
            tnet_list.append(tnet)
            rnet_list.append(rnet)
    # through TF layer
            P = lay.prop_fwd(self.heights_nm()[i-1]/self.period)
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

        self.R_net, self.T_net = rnet, tnet


        """ Calculate field expansions for all layers (including air) starting at top
            Ordering is now top to bottom (inverse of above)! ie f1 is superstrate (top)
            Calculate net downward energy flux in each infintesimal air layer & super/substrates
            (see appendix C in Dossou et al. JOSA 2012)
        """

        self.t_list = []
        self.r_list = []
        self.a_list = []
        num_prop_air    = self.layers[-1].air_ref().num_prop_pw_per_pol
        num_prop_in     = self.layers[-1].num_prop_pw_per_pol
        num_prop_out    = self.layers[0].num_prop_pw_per_pol
        out             = self.layers[0].specular_order

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

        if incoming_amplitudes is None:
            # Set the incident field to be a 0th order plane wave
            # in a given polarisation, from the semi-inf top layer
            d_minus = self.layers[-1].specular_incidence(pol)
        else:
            d_minus = incoming_amplitudes

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

        for i in range(len(self.layers) - 2):
            f1_plus = rnet_list[-2*i-2]*f1_minus
    # net downward flux in infintesimal air layer
            f_mat   = np.matrix(np.concatenate((f1_minus,f1_plus)))
            flux    = f_mat.H*U_mat*f_mat
            down_fluxes.append(flux)

            f2_minus = inv_t12_list[-i-1]*f1_minus
            f2_plus  = rnet_list[-2*i-3]*P_list[-i-1]*f2_minus

            f1_minus = inv_t21_list[-i-2]*P_list[-i-1]*f2_minus

    # bottom air to semi-inf substrate
        f1_plus  = rnet_list[0]*f1_minus

        f2_minus = tnet_list[0]*f1_minus
        self.trans_vector = f2_minus
        flux_TE  = np.linalg.norm(f2_minus[0:num_prop_out])**2
        flux_TM  = np.linalg.norm(f2_minus[neq_PW:neq_PW+num_prop_out])**2
        down_fluxes.append(flux_TE + flux_TM)

    # calculate absorptance in each layer
        for i in range(1,len(down_fluxes)-1):
            a_layer = abs(abs(down_fluxes[i])-abs(down_fluxes[i+1]))
            self.a_list.append(a_layer)
        a_layer = abs(down_fluxes[0]-down_fluxes[-1]-up_flux[0])
        self.a_list.append(a_layer)
        # sum_abs = sum(a_list[-nu_TFs:])
        # a_list.append(a_layer - sum_abs)

    # calculate reflectance in each layer
        for i in range(1,len(up_flux)-1):
            r_layer = abs(abs(up_flux[i])/abs(down_flux[i]))
            self.r_list.append(r_layer)
        r_layer = abs(up_flux[0]/down_fluxes[0])
        self.r_list.append(r_layer)

    # calculate transmittance in each layer
        for i in range(0,len(down_fluxes)-2):
            t_layer = abs(abs(down_fluxes[i+2])/abs(down_fluxes[i]))
            self.t_list.append(t_layer)
        t_layer = abs(down_fluxes[-1]/down_fluxes[0])
        self.t_list.append(t_layer)


    # plot scattering matrices as grayscale images
        if self.layers[0].other_para.plot_scat_mats == True:
            im = np.real(abs(rnet_list[-1]))
            plt.matshow(im,cmap=plt.cm.gray)
            cbar = plt.colorbar(extend='neither')
            plt.xlabel('Incoming Orders')
            plt.ylabel('Outgoing Orders')
            plt.suptitle('Net Reflection Scattering Matrix')
            plt.savefig('Rnet_wl %f' % self.layers[0].light.wl_nm)
            # for i in range(len(rnet_list)):
            #     im = np.real(rnet_list[i])
            #     plt.matshow(im,cmap=plt.cm.gray)
            #     plt.savefig('0testmat%i' % i)

    def _check_periods_are_consistent(self):
        """ Raise an error if layers have different periods."""
        for lay in self.layers:
            assert lay.structure.period == self.period, \
                "All layers in a multilayer stack must have the same period."


def t_r_a_plots(stack_list):
    # plot t,r,a plots each containing results for each layer and total. 
    # Also save text t,r,a to files
    wavelengths = np.array([s.layers[0].light.wl_nm for s in stack_list])
    a_list = []
    t_list = []
    r_list = []
    for stack in stack_list:
        a_list.extend(stack.a_list)
        t_list.extend(stack.t_list)
        r_list.extend(stack.r_list)

    total_h = sum(stack_list[0].heights_nm())
    layers_plot('Lay_Absorb', a_list, wavelengths, total_h)
    layers_plot('Lay_Trans',  t_list, wavelengths, total_h)
    layers_plot('Lay_Reflec', r_list, wavelengths, total_h)


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