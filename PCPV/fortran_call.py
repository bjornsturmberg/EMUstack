# Description

import numpy as np
import subprocess
import sys
import multiprocessing   as mp
from scipy import sqrt
sys.path.append("../PCPV/")

import materials
import objects
import temporary_bullshit as bs

from Fortran_pcpv import pcpv

pi = np.pi

class Modes(object):
    """ Super-class from which Simmo and Anallo inherit common functionality"""
    def bloch_vec(self):
        # FIXME: Felix believes the following line is wrong but it's 
        # how it's always been done
        n_eff = self.structure.superstrate.n(self.light.Lambda)
        return self.light.bloch_vec(self.structure.period, n_eff)


class Anallo(Modes):
    """ Like a :Simmo:, but for a thin film, and calculated analytically."""
    def __init__(self, thin_film, light, other_para):
        self.structure = thin_film
        self.light = light
        self.other_para = other_para
        self.mode_pol       = None
        # self.T12            = None
        # self.R12            = None
        # self.T21            = None
        # self.R21            = None

    def calc_modes(self):
        self.structure.num_prop_air = []
        self.structure.num_prop_TF  = []

        wl = self.light.Lambda
        # refractive indeces listed so that thin film is 2nd,
        # that way can order plane waves in material as in infintesimal air layers,
        # which is ordered consistent with FEM.
        n_1 = self.structure.substrate.n(wl)
        n_2 = self.structure.material.n(wl)
        n_3 = self.structure.superstrate.n(wl)

        n = np.array([n_1, n_2, n_3])
        if self.structure.loss == False:
            n = np.real(n)
        # print n

        #set up equivalent plane waves to FEM calc
        # normalise to lattice constant equal 1 as in FEM
        d_in_nm  = self.structure.period
        d_norm   = 1
        wl_in_nm = wl
        wl_norm  = float(wl_in_nm/d_in_nm)
        light_angles = self.light
        k_list   = 2*pi*n/wl_norm
        k_perp   = self.calc_k_perp(n, k_list, d_norm,
            light_angles.theta, light_angles.phi, self.other_para.max_order_PWs,
            self.other_para.x_order_in, self.other_para.x_order_out,
            self.other_para.y_order_in, self.other_para.y_order_out)
        k_film   = k_perp[1]
        #FIXME: is this really correct??
        self.beta = np.append(k_film, k_film) # add 2nd polarisation
        # print 'k_perp[0]', k_perp[0]
        # print 'k_film', k_film

        # Impedance method only holds when pereability = 1 (good approx for Ag Al etc)
        num_ks, num_pw  = np.shape(k_perp)
        self.structure.nu_tot_ords = num_pw
        matrix_size  = 2*num_pw
        wave_imp_mat = np.zeros((num_ks,matrix_size),complex)
        for i in range(num_ks):
            impedance = np.ones(num_pw)/n[i]
            k         = np.ones(num_pw)*k_list[i]
            wave_imp  = np.zeros(matrix_size, complex)
            # for TE and TM polarisations NEED TO SORT OUT WHICH IS WHICH
            wave_imp[:num_pw] = impedance*k/k_perp[i]
            wave_imp[num_pw:] = impedance*k_perp[i]/k 
            wave_imp_mat[i,:] = wave_imp #numpy array rather than matrix!

        # Scattering matrices from wave impedances
        Z_1 = wave_imp_mat[0,:] # for substrate
        Z_2 = wave_imp_mat[1,:] # for thin film
        Z_3 = wave_imp_mat[2,:] # for superstrate
        r12 = np.mat(np.diag((Z_2-Z_1)/(Z_2 + Z_1)))
        r23 = np.mat(np.diag((Z_3-Z_2)/(Z_3 + Z_2)))
        t12 = np.mat(np.diag((2.*sqrt(Z_2*Z_1))/(Z_2 + Z_1)))
        t23 = np.mat(np.diag((2.*sqrt(Z_3*Z_2))/(Z_3 + Z_2)))


        # Store the matrices
        self.T12 = t12
        self.R12 = r12
        self.T21 = t23
        self.R21 = r23

        if self.structure.height_1 != 'semi_inf':
            # layer thickness in units of period d in nanometers
            h_normed = float(self.structure.height_1)/d_in_nm
            P_array  = np.exp(1j*np.array(k_film, dtype='complex128')*h_normed)
            P_array  = np.append(P_array, P_array) # add 2nd polarisation
            P        = np.matrix(np.diag(P_array),dtype='D')


    def calc_k_perp(self, n, k_list, d, theta, phi, ordre_ls, 
            x_order_in, x_order_out, y_order_in, y_order_out):
        k_perp = []
        # zero_ord = 0
        k_el = 0
        self.num_prop_air = 0
        self.num_prop_TF  = 0
        for k in k_list:
            common_factor = np.sin(theta*pi/180.0)*k
            bloch1 = common_factor*np.cos(phi*pi/180.0)
            bloch2 = common_factor*np.sin(phi*pi/180.0)
            vec_kx = 2.0*pi/d
            vec_ky = 2.0*pi/d

            raw_beta_z_pw = np.array([])
            for px in np.linspace(-ordre_ls, ordre_ls, 2*ordre_ls +1):
                for py in np.linspace(-ordre_ls, ordre_ls, 2*ordre_ls +1):
                    if (px**2 + py**2) <= ordre_ls**2:
                        alpha = bloch1 + vec_kx*px  # Bloch vector along x
                        beta  = bloch2 + vec_ky*py  # Bloch vector along y
                        z_tmp = k**2 - alpha**2 - beta**2
                        sqtmp = sqrt(z_tmp)
                        raw_beta_z_pw  = np.append(raw_beta_z_pw,sqtmp)
                        # number of propagating plane waves in thin film layer
                        if k_el == 0:
                            if z_tmp > 0.0e-5:
                                self.num_prop_air += 1
                            # if px == py == 0:
                            #     zero_ord = len(raw_beta_z_pw)-1
                            if px == x_order_in and py == y_order_in:
                                select_order_in  = len(raw_beta_z_pw)-1
                            if px == x_order_out and py == y_order_out:
                                select_order_out = len(raw_beta_z_pw)-1
                        if k_el == 1 and z_tmp > 0.0e-5:
                            self.num_prop_TF += 1


            # sort plane waves as [propagating big -> small, evanescent small -> big]
            # which is consistent with FEM
            # to be consistent in impedances Z, sub/superstrate sorted to order of thin film
            if k_el == 0: # air superstrates as medium to consistently sort by
                # if np.imag(n[0]) != 0.0:
                #     rev_ind = np.argsort(1*np.imag(raw_beta_z_pw))
                # else:
                rev_ind = np.argsort(-1*np.real(raw_beta_z_pw) + np.imag(raw_beta_z_pw))
                # layer.zero_ord    = int(np.where(rev_ind==zero_ord)[0])
                self.structure.set_ord_in  = int(np.where(rev_ind==select_order_in)[0])
                self.structure.set_ord_out = int(np.where(rev_ind==select_order_out)[0])
            sorted_beta_z_pw = raw_beta_z_pw[rev_ind]
            k_perp.append(sorted_beta_z_pw)
            k_el += 1

        # print 'k_perp[0]', k_perp[0]
        # print 'k_perp[1]', k_perp[1]
        # print 'k_perp[2]', k_perp[2]
        return k_perp


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
            n_effs, self.bloch_vec(), 
            self.structure.lx, self.structure.ly, self.other_para.tol, 
            self.other_para.E_H_field, self.other_para.i_cond, 
            self.other_para.itermax, #self.light.pol_for_fortran(), 
            self.other_para.PropModes, 
            self.other_para.PrintSolution, self.other_para.PrintSupModes, 
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

        # TODO: the following should be calculated later
        ress = pcpv.calc_scat(
            norm_wl, 
            ordre_ls, self.other_para.debug, 
            n_effs, self.bloch_vec(), 
            self.structure.lx, self.structure.ly,
            self.other_para.PrintAll, 
            self.other_para.Checks, 
            neq_pw, zeroth_order,
            self.sol1, self.sol2, 
            type_el, table_nod, x_arr, pp, qq
        )

        self.T12, self.R12, self.T21, self.R21 = [np.mat(x) for x in ress]
        if 1 == self.other_para.PrintOmega:
            P = np.diag(np.exp(1j*self.prop_consts*self.structure.height_1/d))

