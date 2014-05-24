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

# Number of CPUs to use im simulation
num_cores = 8

# Remove results of previous simulations
plotting.clear_previous('.log')
plotting.clear_previous('.pdf')

################ Light parameters #####################
# wavelengths = np.linspace(800,1600,100)
# light_list  = [objects.Light(wl, max_order_PWs = 6, theta = 0.0, phi = 0.0) for wl in wavelengths]
wl = 1600
azi_angles = np.linspace(0,89,90)
light_list  = [objects.Light(wl, max_order_PWs = 6, theta = p, phi = 0.0) for p in azi_angles]

################ Grating parameters #####################
period = 760

superstrate = objects.ThinFilm(period, height_nm = 'semi_inf',
    material = materials.Air, loss = False)

substrate  = objects.ThinFilm(period, height_nm = 'semi_inf',
    material = materials.Air, loss = False)

grating_1 = objects.NanoStruct('1D_array', period, small_d=period/2,
    diameter1=int(round(0.25*period)), diameter2=int(round(0.25*period)), height_nm = 150, 
    inclusion_a = materials.Material(3.61 + 0.0j), inclusion_b = materials.Material(3.61 + 0.0j),
    background = materials.Material(1.46 + 0.0j), 
    loss = True, make_mesh_now = True, force_mesh = True, lc_bkg = 0.1, lc2= 3.0)

grating_2 = objects.NanoStruct('1D_array', period, int(round(0.75*period)), height_nm = 2900, 
    background = materials.Material(1.46 + 0.0j), inclusion_a = materials.Material(3.61 + 0.0j), 
    loss = True, make_mesh_now = True, force_mesh = True, lc_bkg = 0.1, lc2= 3.0)

num_BM = 250


def simulate_stack(light):
   
    ################ Evaluate each layer individually ##############
    sim_superstrate = superstrate.calc_modes(light)
    sim_substrate   = substrate.calc_modes(light)
    sim_grating_1   = grating_1.calc_modes(light, num_BM = num_BM)
    sim_grating_2   = grating_2.calc_modes(light, num_BM = num_BM)

    ################ Evaluate full solar cell structure ##############
    """ Now when defining full structure order is critical and
    solar_cell list MUST be ordered from bottom to top!
    """

    stack = Stack((sim_substrate, sim_grating_1, sim_grating_2, sim_superstrate))
    # stack = Stack((sim_substrate, sim_grating_2, sim_superstrate))
    stack.calc_scat(pol = 'TE')

    return stack


# Run in parallel across wavelengths.
pool = Pool(num_cores)
stacks_list = pool.map(simulate_stack, light_list)


# ### Plot as in Handmer Fig2
# # require phi == 0.0
last_light_object = light_list.pop()
n_H = 3.61 # high refractive index
min_k_label = 15
plotting.t_func_k_plot_1D(stacks_list, last_light_object, n_H, min_k_label)


### Plot as in Handmer Fig1
chosen_PW_order = [-1,0,1,2]
plotting.single_order_T(stacks_list, azi_angles, chosen_PW_order)


param_layer = grating_1 # Specify the layer for which the parameters should be printed on figures.
params_string = plotting.gen_params_string(param_layer, last_light_object, max_num_BMs=num_BM)

stack_wl_list = []
for i in range(len(azi_angles)):
# for i in range(len(wavelengths)):
    stack_wl_list.append(stacks_list[i])
active_layer_nu = 1

Efficiency = plotting.t_r_a_plots(stack_wl_list, azi_angles, params_string, 
    active_layer_nu=active_layer_nu) 
# Efficiency = plotting.t_r_a_plots(stack_wl_list, wavelengths, params_string, 
#     active_layer_nu=active_layer_nu) 


# select_stack = stacks_list[-1]
# plot_mat = select_stack.T_net
# num_prop_PWs = select_stack.layers[0].num_prop_pw_per_pol
# plotting.vis_scat_mats(plot_mat,num_prop_PWs,extra_title='Transmission')

######################## Wrapping up ########################
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
