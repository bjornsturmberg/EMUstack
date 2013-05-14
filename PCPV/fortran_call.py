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
        wl = float(self.light.wl_nm) / self.structure.period
        if self.light.wl_nm % self.structure.period == 0:
            wl += 1e-15
        return wl

    def air_ref(self):
        """ Return an :Anallo: for air for the same :Light: as this."""
        return self.light._air_ref(self.structure.period, self.other_para)

    def calc_grating_orders(self, max_order):
        """ Return the grating order indices px and py, unsorted."""
        # Create arrays of grating order indexes (-p, ..., p)
        pxs = pys = np.arange(-max_order, max_order + 1)
        # The inner loop in the fortran is over y, not x
        # So we call meshgrid with y first
        pys_mesh, pxs_mesh = np.meshgrid(pys, pxs)
        # Which elements of pys_mesh and pxs_mesh correspond to
        # orders low enough that we're interested in?
        low_ord = (pxs_mesh**2 + pys_mesh**2 <= max_order**2)

        return pxs_mesh[low_ord], pys_mesh[low_ord]


class Anallo(Modes):
    """ Like a :Simmo:, but for a thin film, and calculated analytically."""
    def __init__(self, thin_film, light, other_para):
        self.structure  = thin_film
        self.light      = light
        self.other_para = other_para
        self.max_order_PWs   = other_para.max_order_PWs
        self.is_air_ref = False

    def calc_modes(self):
        #TODO: switch to just using calc_kz()?

        kzs = self.calc_kz()

        self.k_z = np.append(kzs, kzs) # add 2nd polarisation
        self.structure.num_pw_per_pol = len(kzs)

    def calc_kz(self):
        """ Return a sorted 1D array of grating orders' kz."""
        d = 1 #TODO: are lx, ly relevant here??

        # Calculate vectors of pxs and pys of all orders
        # with px^2 + py^2 <= self.max_order_PWs
        pxs, pys = self.calc_grating_orders(self.max_order_PWs)

        # Calculate k_x and k_y components of scattered PWs
        # (using the grating equation)
        alpha0, beta0 = self.k_pll_norm()
        alphas = alpha0 + pxs * 2 * pi / d
        betas  = beta0 + pys * 2 * pi / d

        k_z_unsrt = sqrt(self.k()**2 - alphas**2 - betas**2)

        if self.is_air_ref:
            assert not hasattr(self, 'sort_order'), \
                "Are you sure you want to reset the sort_order?"
            # Sort the modes from propagating to fastest decaying
            # k_z is real for propagating waves
            # This must be done consistently
            s = np.argsort(-1*k_z_unsrt.real + k_z_unsrt.imag)
            self.sort_order = s
        else:
            s = self.air_ref().sort_order
            assert s.shape == k_z_unsrt.shape, (s.shape, 
                k_z_unsrt.shape)

        # Find element of k_z_unsrt corresponding to zeroth order
        self.specular_order = np.nonzero((pxs[s] == 0) * (pys[s] == 0))[0][0]

        # Calculate number of propagating plane waves in thin film
        self.num_prop_pw_per_pol = (k_z_unsrt.imag == 0).sum()

        return k_z_unsrt[s]

    def n(self):
        if self.structure.loss:
            return self.structure.material.n(self.light.wl_nm)
        else:
            return self.structure.material.n(self.light.wl_nm).real

    def k(self):
        """ Return the normalised wavenumber in the background material"""
        return 2 * pi * self.n() / self.wl_norm()

    def Z(self):
        """ Return the wave impedance as a 1D array."""
        # Zcr is relative characteristic impedance Zc / Z0
        # Zcr = 1/n assumes that relative permeability is 1
        # Otherwise, use Zcr = \sqrt(epsilon_r / mu_r)
        Zcr = 1./self.n()

        # self.k_z repeats itself halfway through
        # First half is for TE pol, second is for TM
        num_pw2 = len(self.k_z) / 2
        k_z = self.k_z[:num_pw2]
        assert (k_z == self.k_z[num_pw2:]).all()

        # Calculate the (relative) wave impedances Z
        # TE (E in interface plane): Z = Zcr * k/k_z
        # TM (H in interface plane): Z = Zcr / (k/k_z)
        k_on_kz = self.k() / k_z

        # TE is always represented first
        return np.concatenate((Zcr * k_on_kz, Zcr / k_on_kz))


class Simmo(Modes):
    """docstring for Simmo"""
    def __init__(self, structure, light, other_para):
        self.structure      = structure
        self.light          = light
        self.other_para     = other_para
        self.max_order_PWs  = other_para.max_order_PWs
        self.prop_consts    = None
        self.mode_pol       = None

    def run(self, num_BM):
        """ Run the FEM in Fortran"""
        st = self.structure
        wl = self.light.wl_nm
        # 1st and 2nd elements of n_eff are deprecated
        # and _hopefully_ do nothing now (no guarantees)
        # previously, 1st element was superstrate index,
        # 2nd was substrate.
        n_effs = np.array([1., 1., st.background.n(wl), st.inclusion_a.n(wl), 
            st.inclusion_b.n(wl)])
        n_effs = n_effs[:st.nb_typ_el]

        if self.structure.loss == False:
            n_effs = n_effs.real

        pxs, pys = self.calc_grating_orders(self.max_order_PWs)
        num_pw_per_pol = pxs.size
        self.num_BM = num_BM

        d = self.structure.period

        # Prepare for the mesh
        with open("../PCPV/Data/"+self.structure.mesh_file) as f:
            n_msh_pts, n_msh_el = [int(i) for i in f.readline().split()]

        # Size of Fortran's complex superarray (scales with mesh)
        # In theory could do some python-based preprocessing
        # on the mesh file to work out RAM requirements
        cmplx_max = 2**25

        resm = pcpv.calc_modes(
            self.wl_norm(), self.num_BM, self.max_order_PWs, 
            self.structure.period, self.other_para.debug, 
            self.structure.mesh_file, self.structure.mesh_format, 
            n_msh_pts, n_msh_el,
            n_effs, self.k_pll_norm(), 
            self.structure.lx, self.structure.ly, self.other_para.tol, 
            self.other_para.E_H_field, self.other_para.i_cond, 
            self.other_para.itermax, 
            self.other_para.PropModes, 
            self.other_para.PrintSolution, 
            self.other_para.PrintOmega, self.other_para.PrintAll, 
            self.other_para.Checks, self.other_para.q_average, 
            self.other_para.plot_real, self.other_para.plot_imag, 
            self.other_para.plot_abs, 
            self.structure.loss,
            num_pw_per_pol, cmplx_max
        )

        self.k_z, J, J_dag, self.sol1, self.sol2, self.mode_pol = resm

        self.J, self.J_dag = np.mat(J), np.mat(J_dag)

def simulate_and_stack_and_calc(light, layers, pol, options = None):
    """ Simulate each layer, create a stack and run calc_scat.

        INPUTS:

        - `light`   : :Light: object for the stack.

        - `layers`  : A list of the stack's layer objects 
                        (:ThinFilm: or :NanoStruct:), ordered from
                        bottom to top.

        - `pol`     : Polarisation, to pass to :Stack.calc_scat:

        - `options` : A list of dictionaries of options to hand to
                        :Layer.calc_modes:, e.g.:
                        `[{}, {num_BM = 30}, {}]`
    """
    if None == options:
        slayers = [lay.calc_modes(light, simo_para) for lay in layers]
    else:
        slayers = [lay.calc_modes(light, simo_para, **opt) 
                        for lay, opt in zip(layers, options)]
