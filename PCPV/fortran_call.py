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
        self.is_air_ref = False

    def calc_modes(self):
        #TODO: switch to just using calc_kz()?

        kzs = self.calc_kz()

        self.beta = np.append(kzs, kzs) # add 2nd polarisation
        self.structure.nu_tot_ords = len(kzs)

        # FIXME: Yes, this is ludicrous, but historically 2 refers
        # to the layer we're in; 1 is ref_an!!
        # calc_scat demands that R12 is from air to the thin-film
        r21, t21, r12, t12 = r_t_mat_anallo(self, self.air_ref())
        self.R12, self.T12, self.R21, self.T21 = r12, t12, r21, t21


    def calc_kz(self):
        """ Return a sorted 1D array of grating orders' kz."""
        d = 1 #TODO: are lx, ly relevant here??
        ordre_ls = self.ordre_ls
        # Create arrays of grating order indexes p
        pxs = pys = np.arange(-ordre_ls, ordre_ls+1)
        # The inner loop in the fortran is over y, not x
        # So we call meshgrid with y first
        pys_mesh, pxs_mesh = np.meshgrid(pys, pxs)
        # Which elements of pys_mesh and pxs_mesh correspond to
        # orders low enough that we're interested in?
        low_ord = (pxs_mesh**2 + pys_mesh**2 <= ordre_ls**2)

        # Calculate k_x and k_y components of scattered PWs
        # (using the grating equation)
        alpha0, beta0 = self.k_pll_norm()
        alphas = alpha0 + pxs * 2 * pi / d
        betas  = beta0 + pys * 2 * pi / d

        # Calculate all wave vector components k_z
        alpha2_mesh, beta2_mesh = np.meshgrid(alphas**2, betas**2)
        k_zs_unsrt = sqrt((self.k()**2 - alpha2_mesh - beta2_mesh)[low_ord])
        
        if self.is_air_ref:
            assert not hasattr(self, 'sort_order'), \
                "Are you sure you want to reset the sort_order?"
            # Sort the modes from propagating to fastest decaying
            # k_z is real for propagating waves
            # This is consistent with FEM because TODO: we use it there too
            s = np.argsort(-1*k_zs_unsrt.real + k_zs_unsrt.imag)
            self.sort_order = s

            #TODO: make the following redundant?
            # Find element of k_zs_unsrt corresponding to designated pws
            bs = self.other_para
            px_in, px_out = bs.x_order_in, bs.x_order_out
            py_in, py_out = bs.y_order_in, bs.y_order_out

            pxm, pym = pxs_mesh[low_ord][s], pys_mesh[low_ord][s]
            self.structure.set_ord_in  = np.nonzero((pxm == px_in) *
                (pym == py_in))[0][0]
            self.structure.set_ord_out = np.nonzero((pxm == px_out) *
                (pym == py_out))[0][0]

        else:
            s = self.air_ref().sort_order
            assert s.shape == k_zs_unsrt.shape, (s.shape, 
                k_zs_unsrt.shape)

        # Calculate number of propagating plane waves in thin film
        # FIXME: surely this should include both pols?
        self.num_prop_pw = (k_zs_unsrt.imag == 0).sum()

        return k_zs_unsrt[s]

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
        #TODO: make X the inner loop?
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

        self.J, self.J_dag, self.T12, self.R12, self.T21, self.R21 = [
                        np.mat(x) for x in ress]


def r_t_mat_anallo(an1, an2):
    """ Returns R12, T12, R21, T21 at an interface between thin films.

        R12 is the reflection matrix from Anallo 1 off Anallo 2

        The sign of elements in T12 and T21 is fixed to be positive,
        in the eyes of `numpy.sign`
    """
    if len(an1.beta) != len(an2.beta):
        raise ValueError, "Need the same number of plane waves in \
        Anallos %(an1)s and %(an2)s" % {'an1' : an1, 'an2' : an2}

    Z1 = an1.Z()
    Z2 = an2.Z()

    R12 = np.mat(np.diag((Z2 - Z1)/(Z2 + Z1)))
    # N.B. there is potentially a branch choice problem here, stemming
    # from the normalisation to unit flux.
    # We normalise each field amplitude by
    # $chi^{\pm 1/2} = sqrt(k_z/k)^{\pm 1} = sqrt(Z/Zc)^{\pm 1}$
    # The choice of branch in those square roots must be the same as the
    # choice in the related square roots that we are about to take:
    T12 = np.mat(np.diag(2.*sqrt(Z2)*sqrt(Z1)/(Z2+Z1)))
    R21 = -R12
    T21 = T12

    return R12, T12, R21, T21

def r_t_mat_tf_ns(an1, sim2):
    """ Returns R12, T12, R21, T21 at an1-sim2 interface"""
    pass
