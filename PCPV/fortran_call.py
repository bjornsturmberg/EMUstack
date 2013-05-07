# Description

import numpy as np
import subprocess
import sys
import multiprocessing   as mp
from scipy import sqrt
sys.path.append("../PCPV/")

import materials
import objects

from Fortran_pcpv import pcpv

pi = np.pi

class Modes(object):
    """ Super-class from which Simmo and Anallo inherit common functionality"""
    def k_pll_norm(self):
        return self.light.k_pll * self.structure.period

    def wl_norm(self):
        return float(self.light.Lambda) / self.structure.period

    def air_ref(self):
        """ Return an :Anallo: for air for the same :Light: as this."""
        return self.light._air_ref(self.structure.period, self.other_para)


class Anallo(Modes):
    """ Like a :Simmo:, but for a thin film, and calculated analytically."""
    def __init__(self, thin_film, light, other_para):
        self.structure = thin_film
        self.light = light
        self.other_para = other_para
        self.mode_pol       = None
        self.ordre_ls = other_para.max_order_PWs

    def calc_modes(self):
        #TODO: remove this in favour of calc_kz()?

        # set up equivalent plane waves to FEM calc        
        ref_an = self.air_ref()
        sort_order = ref_an.sort_order

        k0 = 2*pi / self.wl_norm()

        kzs_2 = self.calc_kz(self.k(), 
            self.other_para.x_order_in, self.other_para.x_order_out,
            self.other_para.y_order_in, self.other_para.y_order_out,
            sort_order)

        self.beta = np.append(kzs_2, kzs_2) # add 2nd polarisation
        self.structure.nu_tot_ords = len(kzs_2)

        # TODO: these can be de-separated
        # Calculate number of propagating plane waves in air
        kzs_0 = ref_an.beta[:len(ref_an.beta)/2]
        self.num_prop_air = (kzs_0.imag == 0).sum()
        # Calculate number of propagating plane waves in thin film
        self.num_prop_TF = (kzs_2.imag == 0).sum()

        # FIXME: Yes, this is ludicrous, but historically
        # 2 refers to the layer we're in; 1 is ref_an!!
        # calc_scat demands that R12 is from air to the thin-film
        r21, t21, r12, t12 = r_t_mat_anallo(self, ref_an)
        self.R12, self.T12, self.R21, self.T21 = r12, t12, r21, t21


    def calc_kz(self, bs_k,
            x_order_in, x_order_out, y_order_in, y_order_out,
            sort_order = None):
        """ Return a 1D array of grating orders' kz, ordered by 
            `sort_order`.
        """
        n = self.n()
        #TODO: get k from self.k() once calc_modes and calc_scat are fixed
        #k = self.k()
        k = bs_k
        d = 1 #TODO: lx, ly??

        ordre_ls = self.ordre_ls

        blochx, blochy = self.k_pll_norm()

        pxs = pys = np.arange(-ordre_ls, ordre_ls+1)
        # Prepare for an inner loop over y, not x
        pys_mesh, pxs_mesh = np.meshgrid(pys, pxs)
        low_p2 = (pxs_mesh**2 + pys_mesh**2 <= ordre_ls**2)

        # Calculate k_x and k_y components of scattered PWs
        # (using the grating equation)
        alphas = blochx + pxs * 2 * pi / d
        betas  = blochy + pys * 2 * pi / d

        # Calculate all wave vector components k_z
        alpha2_mesh, beta2_mesh = np.meshgrid(alphas**2, betas**2)
        k_zs_unsrt = sqrt((k**2 - alpha2_mesh - beta2_mesh)[low_p2])
        
        #TODO: move this into the reference anallo
        if None == sort_order:
            # Sort the modes from propagating to fastest decaying
            # k_z is real for propagating waves
            # This is consistent with FEM because TODO: we use it there too
            sort_order = np.argsort(-1*np.real(k_zs_unsrt) + np.imag(k_zs_unsrt))

            # TODO: And now for some stuff that should be (re)moved
            # Find element of k_zs_unsrt that corresponds to 
            # px = x_order_in*, py = y_order_*
            select_order_in = np.nonzero((pxs_mesh[low_p2] == x_order_in) *
                (pys_mesh[low_p2] == y_order_in))[0][0]
            select_order_out = np.nonzero((pxs_mesh[low_p2] == x_order_out) *
                (pys_mesh[low_p2] == y_order_out))[0][0]
            self.structure.set_ord_in  = np.nonzero(sort_order==select_order_in)[0][0]
            self.structure.set_ord_out = np.nonzero(sort_order==select_order_out)[0][0]

            return k_zs_unsrt[sort_order], sort_order
        else:
            return k_zs_unsrt[sort_order]

    def n(self):
        if self.structure.loss:
            return self.structure.material.n(self.light.Lambda)
        else:
            return self.structure.material.n(self.light.Lambda).real


    def k(self):
        """ Return the normalised wavenumber in the background material"""
        return 2 * pi * self.n() / self.wl_norm()

    def Z(self):
        """ Return the wave impedance as a 1D array."""
        # Zcr is relative characteristic impedance Zc / Z0
        # Zcr = 1/n assumes that relative permeability is 1
        # Otherwise, use Zcr = \sqrt(epsilon_r / mu_r)
        Zcr = 1./self.n()

        # self.beta repeats itself halfway through
        # First half is for TE pol, second is for TM
        num_pw2 = len(self.beta) / 2

        # Calculate the (relative) wave impedances Z
        # TE (E in interface plane): Z = Zcr * k/k_z
        # TM (H in interface plane): Z = Zcr / (k/k_z)
        k_on_kz = self.k() / self.beta
        return np.concatenate((Zcr * k_on_kz[:num_pw2], Zcr / k_on_kz[num_pw2:]))



class Simmo(Modes):
    """docstring for Simmo"""
    def __init__(self, structure, light, other_para):
        self.structure      = structure
        self.light          = light
        self.other_para     = other_para
        self.max_order_PWs  = other_para.max_order_PWs
        self.prop_consts    = None
        self.mode_pol       = None
        self.T12            = None
        self.R12            = None
        self.T21            = None
        self.R21            = None

    def run(self, num_BM):
        """Return a string that runs the simmo, to execute in a shell"""
        struc = self.structure
        wl = self.light.Lambda
        n_effs = np.array([struc.superstrate.n(wl), struc.substrate.n(wl), 
            struc.background.n(wl), struc.inclusion_a.n(wl), 
            struc.inclusion_b.n(wl)])
        n_effs = n_effs[:struc.nb_typ_el]

        if self.structure.loss == False:
            n_effs = n_effs.real

        self.num_BM = num_BM
        ordre_ls  = self.max_order_PWs

        #TODO: break this stuff off into subfunctions
        pw_ords_x_1d = np.arange(-ordre_ls, ordre_ls + 1)
        pw_ords_y_1d = np.arange(-ordre_ls, ordre_ls + 1)
        # Y is inner loop in fortran
        #FIXME: make X the inner loop?
        pw_ords_y, pw_ords_x = np.meshgrid(pw_ords_y_1d, pw_ords_x_1d)
        sum_sq_ords = pw_ords_x**2 + pw_ords_y**2
        neq_pw = (sum_sq_ords <= ordre_ls**2).sum()
        # Fortran counter starts at 1
        zeroth_order = sum_sq_ords.reshape(-1).argmin() + 1

        # Avoid hitting Wood anomalies
        d = self.structure.period
        norm_wl = self.light.Lambda / d
        if self.light.Lambda % d == 0:
            norm_wl += 1e-15

        # Prepare for the mesh
        with open("../PCPV/Data/"+self.structure.mesh_file) as f:
            n_msh_pts, n_msh_el = [int(i) for i in f.readline().split()]

        resm = pcpv.calc_modes(
            norm_wl, self.num_BM, 
            ordre_ls, d, self.other_para.debug, 
            self.structure.mesh_file, self.structure.mesh_format, 
            n_msh_pts, n_msh_el,
            n_effs, self.light.k_pll, 
            self.structure.lx, self.structure.ly, self.other_para.tol, 
            self.other_para.E_H_field, self.other_para.i_cond, 
            self.other_para.itermax, #self.light.pol_for_fortran(), 
            self.other_para.PropModes, 
            self.other_para.PrintSolution, 
            self.other_para.PrintOmega, self.other_para.PrintAll, 
            self.other_para.Checks, self.other_para.q_average, 
            self.other_para.plot_real, self.other_para.plot_imag, 
            self.other_para.plot_abs, 
            self.structure.loss,
            neq_pw,
        )

        (   self.prop_consts, self.sol1, self.sol2, self.mode_pol, 
            type_el, table_nod, x_arr, pp, qq
            ) = resm

        self.beta = self.prop_consts

        # TODO: the following should be calculated later?
        ress = pcpv.calc_scat(
            norm_wl, 
            ordre_ls, self.other_para.debug, 
            n_effs, self.light.k_pll, 
            self.structure.lx, self.structure.ly,
            self.other_para.PrintAll, 
            self.other_para.Checks, 
            neq_pw, zeroth_order,
            self.sol1, self.sol2, 
            type_el, table_nod, x_arr, pp, qq
        )

        self.T12, self.R12, self.T21, self.R21 = [np.mat(x) for x in ress]


def r_t_mat_anallo(an1, an2):
    """ Return reflection and trans matrices at an1-an2 interface.

        Returns R12, T12, R21, T21.

        The sign of elements in T12 and T21 is fixed to be positive,
        in the eyes of `numpy.sign`
    """
    if len(an1.beta) != len(an2.beta):
        raise ValueError, "Need the same number of plane waves in \
        Anallos %(an1)s and %(an2)s" % {'an1' : an1, 'an2' : an2}

    Z1 = an1.Z()
    Z2 = an2.Z()

    R12 = np.mat(np.diag((Z2 - Z1)/(Z2 + Z1)))
    # Alas, we have a branch choice problem.
    # This stems from the desire for unit flux normalisation.
    # If we do not normalise field amplitudes by
    # $chi^\pm 1 = sqrt(k_z/k)$, then the numerator of T12 is
    # instead 2*Z_2, as per most expressions for impedance mismatch

    sqrt_Z2_Z1 = sqrt(Z2*Z1)
    # The correct branch is the one of the same sign as Z2 and Z1
    # (if they are the same sign)
    # scipy.sqrt is supposed to pick it, but sometimes doesn't.

    # Here we choose so that the real parts of T12 are positive
    sqrt_Z2_Z1 *= np.sign(sqrt_Z2_Z1) * np.sign(Z2+Z1)
    T12 = np.mat(np.diag(2.*sqrt_Z2_Z1/(Z2 + Z1)))
    R21 = -R12
    T21 = T12

    return R12, T12, R21, T21
