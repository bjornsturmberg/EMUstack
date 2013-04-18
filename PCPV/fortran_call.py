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

pcpv_location = '../PCPV/'

pi = np.pi

class Simmo(object):
    """docstring for Simmo"""
    def __init__(self, structure, light, other_para, p):
        self.structure     = structure
        self.light         = light
        self.other_para    = other_para
        self.max_num_BMs   = structure.max_num_BMs
        self.var_BM_min    = other_para.var_BM_min
        self.max_order_PWs = other_para.max_order_PWs
        self.p             = p


    def inclusion_a_n(self):
        if isinstance(self.structure.inclusion_a, complex):
            return self.structure.inclusion_a
        elif isinstance(self.structure.inclusion_a,materials.Material):
            return self.structure.inclusion_a.n(self.light.Lambda)
        else:
            raise  NotImplementedError, "inclusion_a must either be complex number or material instance from materials.py"

    def inclusion_b_n(self):
        if isinstance(self.structure.inclusion_b, complex):
            return self.structure.inclusion_b
        elif isinstance(self.structure.inclusion_b,materials.Material):
            return self.structure.inclusion_b.n(self.light.Lambda)
        else:
            raise  NotImplementedError, "inclusion_b must either be complex number or material instance from materials.py"

    def background_n(self):
        if isinstance(self.structure.background, complex):
            return self.structure.background
        elif isinstance(self.structure.background,materials.Material):
            return self.structure.background.n(self.light.Lambda)
        else:
            raise  NotImplementedError, "background must either be complex number or material instance from materials.py"

    def superstrate_n(self):
        if isinstance(self.structure.superstrate, complex):
            return self.structure.superstrate
        elif isinstance(self.structure.superstrate,materials.Material):
            return self.structure.superstrate.n(self.light.Lambda)
        else:
            raise  NotImplementedError, "superstrate must either be complex number or material instance from materials.py"

    def substrate_n(self):
        if isinstance(self.structure.substrate, complex):
            return self.structure.substrate
        elif isinstance(self.structure.substrate,materials.Material):
            return self.structure.substrate.n(self.light.Lambda)
        else:
            raise  NotImplementedError, "substrate must either be complex number or material instance from materials.py"
        
    def run_fortran(self):
        """Return a string that runs the simmo, to execute in a shell"""
        inclusion_a_n = self.inclusion_a_n()
        inclusion_b_n = self.inclusion_b_n()
        background_n  = self.background_n()
        superstrate_n = self.superstrate_n()
        substrate_n   = self.substrate_n()

        n_effs = np.array([self.superstrate_n(), self.substrate_n(), self.background_n(), 
            self.inclusion_a_n(), self.inclusion_b_n()])
        n_effs = n_effs[:self.structure.nb_typ_el]

        if self.structure.loss == False:
            n_effs = n_effs.real
        if self.light.Lambda > self.var_BM_min:
            if isinstance(self.structure.inclusion_a, complex):
                max_n = self.structure.inclusion_a.real
            else:
                max_n = self.structure.inclusion_a.max_n()
            nval_tmp      = self.max_num_BMs*inclusion_a_n.real/max_n
            nval          = round(nval_tmp)
            ordre_ls      = self.max_order_PWs 
            # # in single layer calc scale # PWs by # of BMs to be used, in multilayer need all equal
            # ordre_ls      = round(self.max_order_PWs*nval/self.max_num_BMs)
        else:
            nval          = self.max_num_BMs
            ordre_ls      = self.max_order_PWs          

        # Run the fortran!
        pcpv.main(
            np.array([self.p, self.p]), self.light.Lambda, nval, 
            ordre_ls, self.structure.period, self.other_para.debug, 
            self.structure.mesh_file, self.structure.mesh_format, 
            n_effs, self.structure.has_substrate, self.light.theta, 
            self.light.phi, self.structure.height_1, 
            self.structure.height_2, self.structure.num_h, 
            self.structure.lx, self.structure.ly, self.other_para.tol, 
            self.other_para.E_H_field, self.other_para.i_cond, 
            self.other_para.itermax, self.light.pol, 
            self.other_para.traLambda, self.other_para.PropModes, 
            self.other_para.PrintSolution, self.other_para.PrintSupModes, 
            self.other_para.PrintOmega, self.other_para.PrintAll, 
            self.other_para.Checks, self.other_para.q_average, 
            self.other_para.plot_real, self.other_para.plot_imag, 
            self.other_para.plot_abs, self.other_para.incident, 
            self.other_para.what4incident, self.other_para.out4incident, 
            self.structure.loss, self.structure.label_nu
        )


# def run_command_in_shell(command):
#     """Submit full string of fortran call to command line"""
#     print command
#     return subprocess.call(command, shell = True)

def calc_k_perp(layer, n, k_list, d, theta, phi, ordre_ls, 
        x_order_in, x_order_out, y_order_in, y_order_out):
    k_perp = []
    # zero_ord = 0
    k_el = 0
    num_prop_air = 0
    num_prop_TF  = 0
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
                            num_prop_air += 1
                        # if px == py == 0:
                        #     zero_ord = len(raw_beta_z_pw)-1
                        if px == x_order_in and py == y_order_in:
                            select_order_in  = len(raw_beta_z_pw)-1
                        if px == x_order_out and py == y_order_out:
                            select_order_out = len(raw_beta_z_pw)-1
                    if k_el == 1 and z_tmp > 0.0e-5:
                        num_prop_TF += 1


        # sort plane waves as [propagating big -> small, evanescent small -> big]
        # which is consistent with FEM
        # to be consistent in impedances Z, sub/superstrate sorted to order of thin film
        if k_el == 0: # air superstrates as medium to consistently sort by
            # if np.imag(n[0]) != 0.0:
            #     rev_ind = np.argsort(1*np.imag(raw_beta_z_pw))
            # else:
            rev_ind = np.argsort(-1*np.real(raw_beta_z_pw) + np.imag(raw_beta_z_pw))
            # layer.zero_ord    = int(np.where(rev_ind==zero_ord)[0])
            layer.set_ord_in  = int(np.where(rev_ind==select_order_in)[0])
            layer.set_ord_out = int(np.where(rev_ind==select_order_out)[0])
        sorted_beta_z_pw = raw_beta_z_pw[rev_ind]
        k_perp.append(sorted_beta_z_pw)
        k_el += 1

    layer.num_prop_air.append(num_prop_air)
    layer.num_prop_TF.append(num_prop_TF)

    # print 'k_perp[0]', k_perp[0]
    # print 'k_perp[1]', k_perp[1]
    # print 'k_perp[2]', k_perp[2]
    return k_perp


def save_scat_mat(matrix, name, st, p, num_pw):
            # reshape matrices to be consistent with pcpv.exe output
            format_label_nu = '%04d' % st
            format_p        = '%04d' % p

            file_name = "st%(st)s_wl%(wl)s_%(mat_name)s.txt" % {
                'st' : format_label_nu, 'wl' : format_p, 'mat_name' : name }
            with open(file_name, 'w') as outfile:
                for k in range(num_pw):
                    for i in range(num_pw):
                        data = [i+1,  k+1, np.real(matrix[i,k]), np.imag(matrix[i,k]),
                            np.abs(matrix[i,k])**2]
                        data = np.reshape(data, (1,5))
                        np.savetxt(outfile, data, fmt=['%4i','%4i','%25.17G','%25.17G','%25.17G'], delimiter='')

def save_k_perp(matrix, name, st, p, num_pw, wl):
            format_label_nu = '%04d' % st
            format_p        = '%04d' % p

            file_name = "st%(st)s_wl%(wl)s_%(mat_name)s.txt" % {
                'st' : format_label_nu, 'wl' : format_p, 'mat_name' : name }
            data_list = [num_pw, wl]
            with open(file_name, 'w') as outfile:
                for k in range(num_pw):
                    # data = [[float(np.real(matrix[k])), float(np.imag(matrix[k]))]]
                    # np.savetxt(outfile, data, fmt=['%25.17G', '%25.17G'], delimiter='')
                    data_list.append(float(np.real(matrix[k])))
                    data_list.append(float(np.imag(matrix[k])))
                np.savetxt(outfile, np.atleast_2d(data_list), fmt='%25.17G', delimiter='')




def scat_mats(layer, light_list, simo_para):
    """Run the calculation of this structures/films scattering matrices"""

    print "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv"
    print "LAYER %s" % layer.label_nu
    print "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"

    wavelengths = np.array([l.Lambda for l in light_list], float)

# For nanostructured layers run pcpv.exe to find scattering matrices
    if isinstance(layer,objects.NanoStruct):
        # Interpolate refractive indecies over wavelength array
        # materials.interp_all(wavelengths)
        materials.interp_needed(wavelengths, layer.inclusion_a)
        materials.interp_needed(wavelengths, layer.inclusion_b)
        materials.interp_needed(wavelengths, layer.background)
        materials.interp_needed(wavelengths, layer.superstrate)
        materials.interp_needed(wavelengths, layer.substrate)

        # List of simulations to calculate, with full arguments
        simmo_list = []
        p          = 1
        for light in light_list:
            simmo_list += [Simmo(layer, light, simo_para, p)]
            p += 1

        # # Launch simos using pool to limit simultaneous instances 
        # command_list = [s.fortran_command_str() for s in simmo_list]
        # pool = mp.Pool(simo_para.num_cores)
        # pool.map(run_command_in_shell, command_list)

        # Run simmos sequentially for now
        for s in simmo_list:
            s.run_fortran()


# For homogeneous layers calculate scattering matrices analytically
    elif isinstance(layer,objects.ThinFilm):
        materials.interp_needed(wavelengths, layer.film_material)
        materials.interp_needed(wavelengths, layer.superstrate)
        materials.interp_needed(wavelengths, layer.substrate)

        layer.num_prop_air = []
        layer.num_prop_TF  = []

        p          = 1
        for wl in wavelengths:
            # refractive indeces listed so that thin film is 2nd,
            # that way can order plane waves in material as in infintesimal air layers,
            # which is ordered consistent with FEM.
            if isinstance(layer.substrate, complex):
                n_1 = layer.substrate
            elif isinstance(layer.substrate,materials.Material):
                n_1 = layer.substrate.n(wl)
            else:
                raise  NotImplementedError, "substrate must either be complex number or material instance from materials.py"

            if isinstance(layer.film_material, complex):
                n_2 = layer.film_material
            elif isinstance(layer.film_material,materials.Material):
                n_2 = layer.film_material.n(wl)
            else:
                raise  NotImplementedError, "film_material must either be complex number or material instance from materials.py"

            if isinstance(layer.superstrate, complex):
                n_3 = layer.superstrate
            elif isinstance(layer.superstrate,materials.Material):
                n_3 = layer.superstrate.n(wl)
            else:
                raise  NotImplementedError, "superstrate must either be complex number or material instance from materials.py"

            n = np.array([n_1, n_2, n_3])
            if layer.loss == False:
                n = np.real(n)
            # print n

            #set up equivalent plane waves to FEM calc
            # normalise to lattice constant equal 1 as in FEM
            d_in_nm  = layer.period
            d_norm   = 1
            wl_in_nm = wl
            wl_norm  = float(wl_in_nm/d_in_nm)
            light_angles = light_list[0]
            k_list   = 2*pi*n/wl_norm
            k_perp   = calc_k_perp(layer, n, k_list, d_norm,
                light_angles.theta, light_angles.phi, simo_para.max_order_PWs,
                simo_para.x_order_in, simo_para.x_order_out,
                simo_para.y_order_in, simo_para.y_order_out)
            k_film   = k_perp[1]
            # print 'k_perp[0]', k_perp[0]
            # print 'k_film', k_film

            # Impedance method only holds when pereability = 1 (good approx for Ag Al etc)
            shape_k_perp = np.shape(k_perp)
            num_ks       = shape_k_perp[0]
            num_pw       = shape_k_perp[1]
            layer.nu_tot_ords = num_pw
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


            # saving matrices to file
            if layer.substrate == layer.superstrate:
                save_scat_mat(r12, 'R12', layer.label_nu, p, matrix_size)
                save_scat_mat(t12, 'T12', layer.label_nu, p, matrix_size)
                save_scat_mat(r23, 'R21', layer.label_nu, p, matrix_size)
                save_scat_mat(t23, 'T21', layer.label_nu, p, matrix_size)
            else:
                save_scat_mat(r12, 'R12', layer.label_nu, p, matrix_size)
                save_scat_mat(t12, 'T12', layer.label_nu, p, matrix_size)
                save_scat_mat(r23, 'R23', layer.label_nu, p, matrix_size)
                save_scat_mat(t23, 'T23', layer.label_nu, p, matrix_size)
            # save_scat_mat(np.mat(np.diag(Z_1)), 'Z1', layer.label_nu, p, matrix_size)
            # save_scat_mat(np.mat(np.diag(Z_2)), 'Z2', layer.label_nu, p, matrix_size)
            # save_scat_mat(np.mat(np.diag(Z_3)), 'Z3', layer.label_nu, p, matrix_size)

            if layer.height_1 != 'semi_inf':
                # layer thickness in units of period d in nanometers
                h_normed = float(layer.height_1)/d_in_nm
                P_array  = np.exp(1j*np.array(k_film, complex)*h_normed)
                P_array  = np.append(P_array, P_array) # add 2nd polarisation
                P        = np.matrix(np.diag(P_array),dtype='D')
                save_scat_mat(P, 'P', layer.label_nu, p, matrix_size)

            save_k_perp(k_film, 'beta', layer.label_nu, p, num_pw, wl)

            p += 1

    label_nu = layer.label_nu + 1
    return label_nu
