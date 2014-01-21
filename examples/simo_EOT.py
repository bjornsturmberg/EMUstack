"""
    simmo_EOT.py is a simulation script template for EMUstack.

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
Simulating Extraordinary Optical Transmission 
as in H. Liu, P. Lalanne, Nature 452 2008 doi:10.1038/nature06762

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
num_cores = 16
# # Alternatively specify the number of CPUs to leave free on machine
# leave_cpus = 4 
# num_cores = mp.cpu_count() - leave_cpus

# Remove results of previous simulations
plotting.clear_previous('.txt')
plotting.clear_previous('.pdf')
plotting.clear_previous('.gif')
plotting.clear_previous('.log')

################ Light parameters #####################
wl_1     = 0.85*940
wl_2     = 1.15*940
no_wl_1  = 600
# Set up light objects
wavelengths = np.linspace(wl_1, wl_2, no_wl_1)
# wavelengths = np.array([785,788,790,792,795])
light_list  = [objects.Light(wl, max_order_PWs = 4) for wl in wavelengths]
# Single wavelength run
# wl_super = 750
# wavelengths = np.array([wl_super])
# light_list  = [objects.Light(wl, max_order_PWs = 1) for wl in wavelengths]


#period must be consistent throughout simulation!!!
period = 940
diam1 = 266
NHs = objects.NanoStruct('NW_array', period, diam1, height_nm = 200, 
    inclusion_a = materials.Air, background = materials.Au, loss = True,
    square = True,    
    make_mesh_now = True, force_mesh = True, lc_bkg = 0.12, lc2= 5.0, lc3= 3.0)#lc_bkg = 0.08, lc2= 5.0)

cover  = objects.ThinFilm(period = period, height_nm = 'semi_inf',
    material = materials.Air, loss = False)
sub = objects.ThinFilm(period = period, height_nm = 'semi_inf',
    material = materials.Air, loss = False)

max_num_BMs = 50

NH_heights = [200]
# num_h = 21
# NH_heights = np.linspace(50,3000,num_h)

def simulate_stack(light):  
    num_BM = max_num_BMs

    
    ################ Evaluate each layer individually ##############
    sim_NHs = NHs.calc_modes(light, num_BM = num_BM)
    sim_cover  = cover.calc_modes(light)
    sim_sub    = sub.calc_modes(light)

# Loop over heights
    height_list = []
    for h in NH_heights:
        stackSub = Stack((sim_sub, sim_NHs, sim_cover), heights_nm = ([h]))
        stackSub.calc_scat(pol = 'TE')
        height_list.append(stackSub)

    return [height_list]


# Run in parallel across wavelengths.
pool = Pool(num_cores)
stacks_wl_list = pool.map(simulate_stack, light_list)
# Run one at a time
# stacks_wl_list = map(simulate_stack, light_list)

######################## Plotting ########################
last_light_object = light_list.pop()


param_layer = NHs # Specify the layer for which the parameters should be printed on figures.
params_string = plotting.gen_params_string(param_layer, last_light_object, max_num_BMs=max_num_BMs)

wls_normed = wavelengths/period

for h in range(len(NH_heights)):
    height = NH_heights[h]
    wl_list = []
    stack_label = 0
    for wl in range(len(wavelengths)):
        wl_list.append(stacks_wl_list[wl][stack_label][h])
    mess_name = '_h%(h)i'% {'h'   : h, }
    plotting.EOT_plot(wl_list, wls_normed, params_string, add_name = mess_name) 
# Dispersion
plotting.omega_plot(wl_list, wavelengths, params_string) 


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

print '*******************************************'
print hms_string
print '*******************************************'
