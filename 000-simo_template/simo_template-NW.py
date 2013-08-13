"""
The Ueber python script, the only one that needs to be edited to set up all 
simulation parameters.
Uses other python scripts to prime the simulation (interpolate raw data over chosen 
wavelengths etc.), then calls the fortran routine pcpv.exe for each wavelength giving 
it all the required details. It does this by spanning a new process for each wavelength,
keeping the total running instances to a maximum number (num_cores_to_use). Finally all 
results are collected in text files and the spectra are plotted. A log file is found in
python_log.txt
"""

import time
import datetime
import numpy as np
import sys
from multiprocessing import Pool
# import multiprocessing   as mp
sys.path.append("../PCPV/")

import clear_previous
import objects
import materials
import plotting
from stack import *


# The following should be in the function "setup_module()",
# but unfortunately simulate_stack is defined in a lazy-but-easy
# way: the structures are inherited rather than passed in.

start = time.time()
################ Simulation parameters ################

# Number of CPUs to use im simulation
num_cores = 5
# # Alternatively specify the number of CPUs to leave free on machine
# leave_cpus = 4 
# num_cores = mp.cpu_count() - leave_cpus

# Remove results of previous simulations
clear_previous.clean('.txt')
clear_previous.clean('.pdf')
clear_previous.clean('.gif')
clear_previous.clean('.log')

################ Light parameters #####################
wl_1     = 310
wl_2     = 1127
no_wl_1  = 3
# Set up light objects
wavelengths = np.linspace(wl_1, wl_2, no_wl_1)
light_list  = [objects.Light(wl, max_order_PWs = 3) for wl in wavelengths]
# Single wavelength run
# wl_super = 450
# wavelengths = np.array([wl_super])
# light_list  = [objects.Light(wl) for wl in wavelengths]


# period must be consistent throughout simulation!!!
period = 600.0

cover  = objects.ThinFilm(period = period, height_nm = 'semi_inf',
    material = materials.Air, loss = True)

bottom1  = objects.ThinFilm(period = period, height_nm = 'semi_inf',
    material = materials.Si_c, loss = False)
TH2  = objects.ThinFilm(period = period, height_nm = 100,
    material = materials.SiO2_a, loss = True)
TF4  = objects.ThinFilm(period = period, height_nm = 200,
    material = materials.Si_c, loss = True)
bottom3  = objects.ThinFilm(period = period, height_nm = 'semi_inf',
    material = materials.TiO2, loss = False)

NWs = objects.NanoStruct('NW_array', period, 60, height_nm = 2330, 
    inclusion_a = materials.Si_c, background = materials.Air, 
    loss = True, nb_typ_el = 2, 
    make_mesh_now = True, force_mesh = True, lc_bkg = 0.2, lc2= 1.0)

# Find num_BM for each simulation in a somewhat arbitrary way
# Maybe roll this out into a Bjorn-specific function
max_n = max([NWs.inclusion_a.n(wl).real for wl in wavelengths])
max_num_BMs = 200

def simulate_stack(light):
    num_BM = round(max_num_BMs * NWs.inclusion_a.n(light.wl_nm).real/max_n)
    # num_BM = max_num_BMs
    
    ################ Evaluate each layer individually ##############
    sim_cover = cover.calc_modes(light)
    sim_bot1  = bottom1.calc_modes(light)
    sim_TH2   = TH2.calc_modes(light)
    sim_bot3  = bottom3.calc_modes(light)
    sim_TF4   = TF4.calc_modes(light)
    sim_NWs   = NWs.calc_modes(light, num_BM = num_BM)

    ################ Evaluate full solar cell structure ##############
    """ Now when defining full structure order is critical and
    solar_cell list MUST be ordered from bottom to top!
    """

    stack0 = Stack((sim_bot1, sim_cover))
    stack1 = Stack((sim_bot1, sim_TH2, sim_NWs, sim_cover))
    # stack1 = Stack((sim_bot1, sim_TH2, sim_TF4, sim_cover))
    stack0.calc_scat(pol = 'TE')
    stack1.calc_scat(pol = 'TE')

# multiple heights for sim_TF4
    stack2_indiv_hs = []
    average_t = 0
    average_r = 0
    average_a = 0

    num_h = 10
    for h in np.linspace(100,2000,num_h):
        stack2 = Stack((sim_bot1,sim_TH2,sim_TF4, sim_cover), heights_nm = ([TH2.height_nm,h]))
        stack2.calc_scat(pol = 'TE')

        stack2_indiv_hs.append(stack2)

        average_t += stack2.t_list[-1]/num_h
        average_r += stack2.r_list[-1]/num_h
        average_a += stack2.a_list[-1]/num_h
    stack2.t_list[-1] = average_t
    stack2.r_list[-1] = average_r
    stack2.a_list[-1] = average_a

# stack2 contains info on the last height and the average spectra
    return [stack0,stack1,stack2,stack2_indiv_hs] 
    # return stack1


# Run in parallel across wavelengths.
pool = Pool(num_cores)
stacks_wl_list = pool.map(simulate_stack, light_list)
# Run one at a time
# stacks_wl_list = map(simulate_stack, light_list)



# Pull apart simultaneously simulated stakes into single stack, all wls arrays.
# Unnecissary if just returning a single stack
np.array(stacks_wl_list)


######################## Plotting ########################
last_light_object = light_list.pop()




# #### Example 0: simple interface.
# param_layer = bottom1 # Specify the layer for which the parameters should be printed on figures.
# params_string = plotting.gen_params_string(param_layer, last_light_object)
# stack_label = 0 # Specify which stack you are dealing with.
# stack0_wl_list = []
# for i in range(len(wavelengths)):
#     stack0_wl_list.append(stacks_wl_list[i][stack_label])
# # Plot total transmission, reflection, absorption & that of each layer.
# Efficiency = plotting.t_r_a_plots(stack0_wl_list, wavelengths, params_string, 
#     stack_label=stack_label) 



param_layer = TF4 # Specify the layer for which the parameters should be printed on figures.
params_string = plotting.gen_params_string(param_layer, last_light_object, max_num_BMs=max_num_BMs)

#### Example 1: simple multilayered stack.
stack_label = 1 # Specify which stack you are dealing with.
stack1_wl_list = []
for i in range(len(wavelengths)):
    stack1_wl_list.append(stacks_wl_list[i][stack_label])
active_layer_nu = 2 # Specify which layer is the active one (where absorption generates charge carriers).
# Plot total transmission, reflection, absorption & that of each layer. 
# Also calculate efficiency of active layer.
Efficiency = plotting.t_r_a_plots(stack1_wl_list, wavelengths, params_string, 
    active_layer_nu=active_layer_nu, stack_label=stack_label) 
# Dispersion
plotting.omega_plot(stack1_wl_list, wavelengths, params_string, stack_label=stack_label) 
# # Energy Concentration
# which_layer = 2
# which_modes = [1,2] # can be a single mode or multiple modes
# plotting.E_conc_plot(stack1_wl_list, which_layer, which_modes, wavelengths, 
#     params_string, stack_label=stack_label)





# #### Example 2: averaged multilayered stack where one layer has many heights.
# stack_label = 2
# active_layer_nu = 2
# stack2_wl_list = []
# for i in range(len(wavelengths)):
#     stack2_wl_list.append(stacks_wl_list[i][stack_label])
# Efficiency = plotting.t_r_a_plots(stack2_wl_list, wavelengths, params_string, 
#     active_layer_nu=active_layer_nu, stack_label=stack_label)



# #### Example 3: individual spectra of multilayered stack where one layer has many heights.
# stack_label = 3
# active_layer_nu = 2
# number_of_hs = len(stacks_wl_list[0][stack_label])
# for h in range(number_of_hs):
#     gen_name = '_h-'
#     h_name = str(h)
#     additional_name = gen_name+h_name # You can add an arbitry string onto the end of the spectra filenames.
#     stack3_hs_wl_list = []
#     for i in range(len(wavelengths)):
#         stack3_hs_wl_list.append(stacks_wl_list[i][stack_label][h])
#     Efficiency = plotting.t_r_a_plots(stack3_hs_wl_list, wavelengths, params_string, 
#         active_layer_nu=active_layer_nu, stack_label=stack_label, add_name = additional_name)
# # Animate spectra as a function of heights.
# from os import system as ossys
# delay = 5 # delay between images in gif in seconds
# names = 'Lay_Absorb_stack'+str(stack_label)+gen_name
# gif_cmd = 'convert -delay %(d)i +dither -layers Optimize -colors 16 %(n)s*.pdf %(n)s.gif'% {
# 'd' : delay, 'n' : names}
# ossys(gif_cmd)
# opt_cmd = 'gifsicle -O2 %(n)s.gif -o %(n)s-opt.gif'% {'n' : names}
# ossys(opt_cmd)
# rm_cmd = 'rm %(n)s.gif'% {'n' : names}
# ossys(rm_cmd)



######################## Single Wavelength Plotting ########################
# Plot transmission as a function of k vector.
# plotting.t_func_k_plot(stack3_wl_list)

# # Visualise the Scattering Matrices
# for i in range(len(wavelengths)):
#     extra_title = 'R_net'
#     plotting.vis_scat_mats(stack1_wl_list[i].R_net, i, extra_title)
#     extra_title = 'R_12'
#     plotting.vis_scat_mats(stack1_wl_list[i].layers[2].T21, i, extra_title)




######################## Wrapping up ########################
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
