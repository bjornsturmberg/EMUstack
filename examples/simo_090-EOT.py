"""
    simo_090-EOT.py is a simulation script template for EMUstack.

    Copyright (C) 2015  Bjorn Sturmberg

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
sys.path.append("../backend/")

import objects
import materials
import plotting
from stack import *

start = time.time()
################ Simulation parameters ################

# Number of CPUs to use in simulation
num_cores = 16

# Remove results of previous simulations
plotting.clear_previous()

################ Light parameters #####################
wl_1     = 0.85*940
wl_2     = 1.15*940
no_wl_1  = 600
# Set up light objects
wavelengths = np.linspace(wl_1, wl_2, no_wl_1)
# wavelengths = np.array([785,788,790,792,795])
light_list  = [objects.Light(wl, max_order_PWs = 5, theta = 0.0, phi = 0.0) for wl in wavelengths]


#period must be consistent throughout simulation!!!
period = 940
diam1 = 266
NHs = objects.NanoStruct('2D_array', period, diam1, height_nm = 200,
    inclusion_a = materials.Air, background = materials.Au, loss = True,
    inc_shape = 'square',
    make_mesh_now = True, force_mesh = True, lc_bkg = 0.12, lc2= 5.0, lc3= 3.0)#lc_bkg = 0.08, lc2= 5.0)

strate  = objects.ThinFilm(period = period, height_nm = 'semi_inf',
    material = materials.Air, loss = False)

NH_heights = [200]
# num_h = 21
# NH_heights = np.linspace(50,3000,num_h)

def simulate_stack(light):
    ################ Evaluate each layer individually ##############
    sim_NHs    = NHs.calc_modes(light)
    sim_strate = strate.calc_modes(light)

# Loop over heights
    height_list = []
    for h in NH_heights:
        stackSub = Stack((sim_strate, sim_NHs, sim_strate), heights_nm = ([h]))
        stackSub.calc_scat(pol = 'TE')
        height_list.append(stackSub)

    return [height_list]


# Run in parallel across wavelengths.
pool = Pool(num_cores)
stacks_list = pool.map(simulate_stack, light_list)
# Save full simo data to .npz file for safe keeping!
np.savez('Simo_results', stacks_list=stacks_list)



######################## Plotting ########################
last_light_object = light_list.pop()

wls_normed = wavelengths/period

for h in range(len(NH_heights)):
    height = NH_heights[h]
    wl_list = []
    stack_label = 0
    for wl in range(len(wavelengths)):
        wl_list.append(stacks_list[wl][stack_label][h])
    mess_name = '_h%(h)i'% {'h'   : h, }
    plotting.EOT_plot(wl_list, wls_normed, add_name = mess_name, savetxt = True)
# Dispersion
plotting.omega_plot(wl_list, wavelengths)


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
print ''
