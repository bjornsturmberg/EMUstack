"""
    simmo_template-many_in_one.py is a simulation script template for EMUstack.

    Copyright (C) 2013  Bjorn Sturmberg, Kokou Dossou, Felix Lawrence

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
Template python script file to execute a simulation. To start, open a terminal and change
directory to the directory containing this file (which must be in the same directory as 
the EMUstack directory). Run this script file by executing the following in the command line

$ python simo_template-many_in_one.py

This will use num_cores worth of your CPUs, and by default return you in the command
line, having printed results and saved plots to file as specified towards the end of 
this file. If instead you wish to have direct access to the simulation results (for 
further manipulation, debugging etc.) run this script with

$ python -i simo_template-many_in_one.py

which, after the calculations are complete, will return you into an interactive session 
of python, in which all simulation objects are accessible. In this session you can access
the docstrings of objects/classes/methods by typing

>>> from pydoc import help
>>> help(objects.Light)

where we have accessed the docstring of the Light class from objects.py


In real simulation scripts replace this docstring with a brief description of the 
simulation, eg.
`Simulating NW array with period 600 nm and NW diameter 120 nm, placed ontop of 
different substrates.'
"""

import time
import datetime
import numpy as np
import sys
from multiprocessing import Pool
sys.path.append("../EMUstack_backend/")

import objects
import materials
import plotting
from stack import *

start = time.time()
################ Simulation parameters ################

# Number of CPUs to use im simulation
num_cores = 7
# # Alternatively specify the number of CPUs to leave free on machine
# leave_cpus = 4 
# num_cores = mp.cpu_count() - leave_cpus

# Remove results of previous simulations
plotting.clear_previous('.txt')
plotting.clear_previous('.pdf')
plotting.clear_previous('.gif')
plotting.clear_previous('.log')

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
period = 600
max_num_BMs = 200

cover  = objects.ThinFilm(period, height_nm = 'semi_inf',
    material = materials.Air, loss = True)

bottom1  = objects.ThinFilm(period, height_nm = 'semi_inf',
    material = materials.Si_c, loss = False)
TH2  = objects.ThinFilm(period, height_nm = 100,
    material = materials.SiO2_a, loss = True)
TF4  = objects.ThinFilm(period, height_nm = 200,
    material = materials.Si_c, loss = True)
bottom3  = objects.ThinFilm(period, height_nm = 'semi_inf',
    material = materials.TiO2, loss = False)

NW_diameter = 120
NWs = objects.NanoStruct('NW_array', period, NW_diameter, height_nm = 2330, 
    inclusion_a = materials.Si_c, background = materials.Air, loss = True,    
    make_mesh_now = True, force_mesh = True, lc_bkg = 0.2, lc2= 1.0)

# Find num_BM for each simulation (wl) as num decreases w decreasing index contrast.
max_n = max([NWs.inclusion_a.n(wl).real for wl in wavelengths])

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
    stack1 = Stack((sim_bot1, sim_NWs, sim_cover))
    # stack1 = Stack((sim_bot1, sim_TH2, sim_NWs, sim_cover))
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





param_layer = NWs # Specify the layer for which the parameters should be printed on figures.
params_string = plotting.gen_params_string(param_layer, last_light_object, max_num_BMs=max_num_BMs)

#### Example 1: simple multilayered stack.
stack_label = 1 # Specify which stack you are dealing with.
stack1_wl_list = []
for i in range(len(wavelengths)):
    stack1_wl_list.append(stacks_wl_list[i][stack_label])
active_layer_nu = 1

Efficiency = plotting.t_r_a_plots(stack1_wl_list, wavelengths, params_string, 
    active_layer_nu=active_layer_nu, stack_label=stack_label) 
# Dispersion
plotting.omega_plot(stack1_wl_list, wavelengths, params_string, stack_label=stack_label) 














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




# param_layer = TF4 # Specify the layer for which the parameters should be printed on figures.
# params_string = plotting.gen_params_string(param_layer, last_light_object, max_num_BMs=max_num_BMs)

# #### Example 1: simple multilayered stack.
# stack_label = 1 # Specify which stack you are dealing with.
# stack1_wl_list = []
# for i in range(len(wavelengths)):
#     stack1_wl_list.append(stacks_wl_list[i][stack_label])
# active_layer_nu = 2 # Specify which layer is the active one (where absorption generates charge carriers).
# # Plot total transmission, reflection, absorption & that of each layer. 
# # Also calculate efficiency of active layer.
# Efficiency = plotting.t_r_a_plots(stack1_wl_list, wavelengths, params_string, 
#     active_layer_nu=active_layer_nu, stack_label=stack_label) 
# # Dispersion
# plotting.omega_plot(stack1_wl_list, wavelengths, params_string, stack_label=stack_label) 
# # # Energy Concentration
# # which_layer = 2
# # which_modes = [1,2] # can be a single mode or multiple modes
# # plotting.E_conc_plot(stack1_wl_list, which_layer, which_modes, wavelengths, 
# #     params_string, stack_label=stack_label)





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
# delay = 5 # delay between images in gif in hundredths of a second
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
