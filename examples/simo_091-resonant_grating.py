"""
    simo_101-resonant_grating.py is a simulation script template for EMUstack.

    Copyright (C) 2013  Bjorn Sturmberg

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

"""

"""



import time
import datetime
import numpy as np
import sys
from multiprocessing import Pool
sys.path.append("../backend/")

import objects
import materials
import plotting
from stack import *

start = time.time()
################ Simulation parameters ################

# Number of CPUs to use in simulation
num_cores = 5

# Remove results of previous simulations
plotting.clear_previous('.npz')
plotting.clear_previous('.txt')
plotting.clear_previous('.pdf')
plotting.clear_previous('.log')

################ Light parameters #####################
wl_1     = 900
wl_2     = 1200
no_wl_1  = 3
# Set up light objects
wavelengths = np.linspace(wl_1, wl_2, no_wl_1)
light_list  = [objects.Light(wl, max_order_PWs = 3) for wl in wavelengths]


# period must be consistent throughout simulation!!!
period = 120
num_BM = 90

superstrate = objects.ThinFilm(period, height_nm = 'semi_inf',
    material = materials.Material(3.5 + 0.0j), loss = True)

homo_film  = objects.ThinFilm(period, height_nm = 5,
    material = materials.Material(3.6 + 0.27j), loss = True)

substrate = objects.ThinFilm(period, height_nm = 'semi_inf',
    material = materials.Air, loss = False)

grating_diameter = 100
grating = objects.NanoStruct('1D_array', period, grating_diameter, height_nm = 25, 
    inclusion_a = materials.Ag, background = materials.Material(1.5 + 0.0j), loss = True,
    make_mesh_now = True, force_mesh = True, lc_bkg = 0.05, lc2= 4.0)

mirror = objects.ThinFilm(period, height_nm = 100,
    material = materials.Ag, loss = True)


def simulate_stack(light):
    ################ Evaluate each layer individually ##############
    sim_superstrate = superstrate.calc_modes(light)
    sim_homo_film   = homo_film.calc_modes(light)
    sim_substrate   = substrate.calc_modes(light)
    sim_grating     = grating.calc_modes(light, num_BM = num_BM)
    sim_mirror      = mirror.calc_modes(light)

    ################ Evaluate full solar cell structure ##############
    """ Now when defining full structure order is critical and
    solar_cell list MUST be ordered from bottom to top!
    """
    stack = Stack((sim_substrate, sim_mirror, sim_grating, sim_homo_film, sim_superstrate))
    stack.calc_scat(pol = 'TE')

    return stack



# Run in parallel across wavelengths.
pool = Pool(num_cores)
stacks_list = pool.map(simulate_stack, light_list)
# Save full simo data to timestamped .npz file for safe keeping!
simotime = str(time.strftime("%Y%m%d%H%M%S", time.localtime()))
np.savez('Simo_results'+simotime, stacks_list=stacks_list)
    

######################## Plotting ########################
last_light_object = light_list.pop()
param_layer = grating_1 # Specify the layer for which the parameters should be printed on figures.
params_string = plotting.gen_params_string(param_layer, last_light_object, max_num_BMs=num_BM)


active_layer_nu = 3 # Specify which layer is the active one (where absorption generates charge carriers).
# Plot total transmission, reflection, absorption & that of each layer. 
# Also calculate efficiency of active layer.
Efficiency = plotting.t_r_a_plots(stacks_list, wavelengths, params_string, active_layer_nu=active_layer_nu) 
# Dispersion
# plotting.omega_plot(stacks_list, wavelengths, params_string) 
# # Energy Concentration
# which_layer = 2
# which_modes = [1,2] # can be a single mode or multiple
# plotting.E_conc_plot(stacks_list, which_layer, which_modes, wavelengths, 
#     params_string)


######################## Single Wavelength Plotting ########################
# Plot transmission as a function of k vector.
# plotting.t_func_k_plot(stacks_list)


# # Visualise the Scattering Matrices
# for i in range(len(wavelengths)):
#     extra_title = 'R_net'
#     plotting.vis_scat_mats(stacks_list[i].R_net, i, extra_title)
#     extra_title = 'R_12'
#     plotting.vis_scat_mats(stacks_list[i].layers[2].T21, i, extra_title)







# betas = stacks_wl_list[0][0][0].layers[1].k_z
# print betas
# betas = stacks_wl_list[0][0][0].layers[0].k_z
# print betas

# Rnet = stacks_wl_list[0][0][0].R_net
# J_mat = stacks_wl_list[0][0][0].layers[1].J
# T_c = np.sum((np.abs(stacks_wl_list[0][0][0].layers[1].T12)), axis=1)
# print T_c
# print Rnet
# print J_mat
# print_fmt = zip(np.real(betas),np.imag(betas),T_c)
# np.savetxt('Coupling_beta.txt', print_fmt, fmt = '%7.4f')
# print_fmt = zip(np.real(betas),np.imag(betas))
# np.savetxt('Coupling_beta.txt', print_fmt)







# Wrapping up simulation by printing to screen and log file
print '\n*******************************************'
print 'The ultimate efficiency is %12.8f' % Efficiency
print '-------------------------------------------'

# Calculate and record the (real) time taken for simulation
elapsed = (time.time() - start)
hms     = str(datetime.timedelta(seconds=elapsed))
hms_string = 'Total time for simulation was \n \
    %(hms)s (%(elapsed)12.3f seconds)'% {
            'hms'       : hms,
            'elapsed'   : elapsed, }

python_log = open("python_log.log", "w")
python_log.write(hms_string)
python_log.close()

print hms_string
print '*******************************************'
print ''
