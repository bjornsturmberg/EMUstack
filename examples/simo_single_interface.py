"""
    simmo_single_interface.py is a simulation example for EMUstack.

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
Calculating transmission, and reflection of a simple interface between 2 homogeneous media. 
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
num_cores = 7

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


# period must be consistent throughout simulation!!!
period = 600

superstrate = objects.ThinFilm(period, height_nm = 'semi_inf',
    material = materials.Air, loss = True)
substrate   = objects.ThinFilm(period, height_nm = 'semi_inf',
    material = materials.Si_c, loss = False)

def simulate_stack(light):    
    ################ Evaluate each layer individually ##############
    sim_superstrate = superstrate.calc_modes(light)
    sim_substrate   = substrate.calc_modes(light)
    ################ Evaluate full solar cell structure ##############
    """ Now when defining full structure order is critical and
    solar_cell list MUST be ordered from bottom to top!
    """

    stack0 = Stack((sim_substrate, sim_superstrate))
    stack0.calc_scat(pol = 'TE')

    return [stack0]


# Run in parallel across wavelengths.
pool = Pool(num_cores)
stacks_wl_list = pool.map(simulate_stack, light_list)


# Pull apart simultaneously simulated stakes into single stack, all wls arrays.
# Unnecissary if just returning a single stack
np.array(stacks_wl_list)


######################## Plotting ########################
last_light_object = light_list.pop()


param_layer = substrate # Specify the layer for which the parameters should be printed on figures.
params_string = plotting.gen_params_string(param_layer, last_light_object)
stack_label = 0 # Specify which stack you are dealing with.
stack_wl_list = []
for i in range(len(wavelengths)):
    stack_wl_list.append(stacks_wl_list[i][stack_label])
# Plot total transmission, reflection, absorption & that of each layer.
Efficiency = plotting.t_r_a_plots(stack_wl_list, wavelengths, params_string, 
    stack_label=stack_label) 



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
print ''