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

simo_para  = objects.Controls(debug = False,max_order_PWs = 5, num_cores = 3,
    PrintAll = 0, Checks = 0, PrintSolution = 0, PrintSupModes = 0)
# Remove results of previous simulations
clear_previous.clean('.txt')
clear_previous.clean('.pdf')
clear_previous.clean('.gif')
clear_previous.clean('.log')

################ Light parameters #####################
wl_1     = 900
wl_2     = 1200
no_wl_1  = 3
# Set up light objects
wavelengths = np.linspace(wl_1, wl_2, no_wl_1)
light_list  = [objects.Light(wl) for wl in wavelengths]
# # Single wavelength run
# wl_super = 1000
# wavelengths = np.array([wl_super])
# light_list  = [objects.Light(wl) for wl in wavelengths]


# period must be consistent throughout simulation!!!
period = 120

cover  = objects.ThinFilm(period = period, height_nm = 'semi_inf',
    material = materials.Material(3.5 + 0.0j), loss = True)

homo_film  = objects.ThinFilm(period = period, height_nm = 5,
    material = materials.Material(3.6 + 0.27j), loss = True)

bottom = objects.ThinFilm(period = period, height_nm = 'semi_inf',
    material = materials.Air, loss = False)

grating_1 = objects.NanoStruct('1D_grating', period, 100, height_nm = 25, 
    inclusion_a = materials.Ag, background = materials.Material(1.5 + 0.0j), loss = True, nb_typ_el = 4, 
    make_mesh_now = True, force_mesh = True, lc_bkg = 0.05, lc2= 4.0)

mirror = objects.ThinFilm(period = period, height_nm = 100,
    material = materials.Ag, loss = True)


# Find num_BM for each simulation in a somewhat arbitrary way
# Maybe roll this out into a Bjorn-specific function
# max_n = max([grating_1.inclusion_a.n(wl).real for wl in wavelengths])
max_num_BMs = 163

def simulate_stack(light):
    num_BM = max_num_BMs#round(max_num_BMs * grating_1.inclusion_a.n(light.wl_nm).real/max_n)
    
    ################ Evaluate each layer individually ##############
    sim_cover = cover.calc_modes(light, simo_para)
    sim_homo_film = homo_film.calc_modes(light, simo_para)
    sim_bot = bottom.calc_modes(light, simo_para)
    sim_grat1 = grating_1.calc_modes(light, simo_para, num_BM = num_BM)
    sim_mirror = mirror.calc_modes(light, simo_para)

    ################ Evaluate full solar cell structure ##############
    """ Now when defining full structure order is critical and
    solar_cell list MUST be ordered from bottom to top!
    """
    stack = Stack((sim_bot, sim_mirror, sim_grat1, sim_homo_film, sim_cover))
    stack.calc_scat(pol = 'TE')

    return stack



# Run in parallel across wavelengths.
pool = Pool(simo_para.num_cores)
stack_list = pool.map(simulate_stack, light_list)
# # Run one at a time
# module.stack_list = map(simulate_stack, light_list)
    

######################## Plotting ########################
last_light_object = light_list.pop()
param_layer = grating_1 # Specify the layer for which the parameters should be printed on figures.
params_string = plotting.gen_params_string(param_layer, last_light_object, simo_para.max_order_PWs,
    max_num_BMs=max_num_BMs)


active_layer_nu = 2 # Specify which layer is the active one (where absorption generates charge carriers).
# Plot total transmission, reflection, absorption & that of each layer. 
# Also calculate efficiency of active layer.
Efficiency = plotting.t_r_a_plots(stack_list, wavelengths, params_string, 
    active_layer_nu=active_layer_nu) 
# Dispersion
plotting.omega_plot(stack_list, wavelengths, params_string) 
# # Energy Concentration
# which_layer = 2
# which_modes = [1,2] # can be a single mode or multiple
# plotting.E_conc_plot(stack_list, which_layer, which_modes, wavelengths, 
#     params_string)




######################## Single Wavelength Plotting ########################
# Plot transmission as a function of k vector.
plotting.t_func_k_plot(stack_list)



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
