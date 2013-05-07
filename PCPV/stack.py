import numpy as np
from plotting import layers_plot
from objects import Anallo, Simmo

class Stack(object):
    """ Represents a stack of layers evaluated at one frequency.

        This includes the semi-infinite input and output layers.

        INPUTS:

          - `layers` : a tuple of :ThinFilm:s and :NanoStruct:s 
            ordered from top to bottom layer.

          - `heights` : a tuple of the heights of the inside layers,
            i.e., all layers except for the top and bottom. This 
            overrides any heights specified in the :ThinFilm: or
            :NanoStruct: objects.

    """
    def __init__(self, layers, heights = None):
        self.layers = tuple(layers)
        self.heights = heights

    def heights(self):
        if None != self.heights:
            return self.heights
        else:
            return (lay.structure.height for lay in self.layers[1:-1])

    def structures(self):
        return (lay.structure for lay in self.layers)

    def calc_scat(self):
        """ Calculate the transmission and reflection matrices of the stack"""
        #TODO: rewrite this using R and T for interfaces not layers

        # Check that all structures have the same period
        for cell in self.layers:
            if cell.structure.period == self.layers[0].structure.period:
                pass
            else:
                raise ValueError, "All layers in a multilayer stack must have the same period!"

        nu_intfaces     = 2*(len(self.layers)-1)
        neq_PW          = self.layers[0].structure.nu_tot_ords # assumes incident from homogeneous film
        PW_pols         = 2*neq_PW
        num_prop_air    = self.layers[-1].num_prop_air
        num_prop_in     = self.layers[-1].num_prop_TF
        num_prop_out    = self.layers[0].num_prop_TF
        zero_vec        = np.matrix(np.zeros((PW_pols, 1),complex))
        I_air           = np.matrix(np.eye(PW_pols),dtype='D')
        inc             = self.layers[-1].structure.set_ord_in 
        out             = self.layers[0].structure.set_ord_out

        self.t_list = []
        self.r_list = []
        self.a_list = []

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
            r12_list.append(st1.R12)
            r21_list.append(st1.R21)
            # potential to save on one transpose
            t12_list.append(st1.T12)
            t21_list.append(st1.T21)
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
            P = np.mat(np.diag(np.exp(1j*lay.beta * float(lay.structure.height_1)/lay.structure.period)))
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

        for i in range(len(self.layers) - 2):
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
            plt.savefig('Rnet_wl %f' % self.layers[0].light.Lambda)
            # for i in range(len(rnet_list)):
            #     im = np.real(rnet_list[i])
            #     plt.matshow(im,cmap=plt.cm.gray)
            #     plt.savefig('0testmat%i' % i)

    def total_height(self):
        return sum([l.structure.height_1 for l in self.layers[1:-1]])


def r_t_mat(lay1, lay2):
    """ Return R12, T12, R21, T21 at an interface between lay1
        and lay2.
    """
    if isinstance(lay1, Anallo) and isinstance(lay2, Anallo):
        return r_t_mat_anallo(lay1, lay2)
    elif isinstance(lay1, Anallo) and isinstance(lay2, Simmo):
        pass
    elif isinstance(lay1, Simmo) and isinstance(lay2, Simmo):
        pass
    elif isinstance(lay1, Simmo) and isinstance(lay2, Simmo):
        raise NotImplementedError















def t_r_a_plots(stack_list):
    # plot t,r,a plots each containing results for each layer and total. 
    #Also save text t,r,a to files
    wavelengths = np.array([s.layers[0].light.Lambda for s in stack_list])
    a_list = []
    t_list = []
    r_list = []
    for stack in stack_list:
        a_list.extend(stack.a_list)
        t_list.extend(stack.t_list)
        r_list.extend(stack.r_list)

    total_h = stack_list[0].total_height()
    layers_plot('Lay_Absorb', a_list, wavelengths, total_h)
    layers_plot('Lay_Trans',  t_list, wavelengths, total_h)
    layers_plot('Lay_Reflec', r_list, wavelengths, total_h)
    

def god_function(raw_layers, lights, simmoargs = {}):
    for lay in layers:
        # Make the mesh (if appropriate)
        lay.make_mesh()

    stack_list = []
    for li in lights:
        # Find modes of the layer
        layers = (lay.eval(li, simmoargs) for lay in layers)
        stack = Stack(layers)
        stack.calc_scat()
        stack.delete_working()
