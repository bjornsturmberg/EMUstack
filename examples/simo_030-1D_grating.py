"""
    simmo_NW_array.py is a simulation example for EMUstack.

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
Template python script file to execute a simulation. To start, open a terminal and change
directory to the directory containing this file (which must be in the same directory as 
the EMUstack directory). Run this script file by executing the following in the command line

$ python simmo_NW_array.py

This will use num_cores worth of your CPUs, and by default return you in the command
line, having printed results and saved plots to file as specified towards the end of 
this file. If instead you wish to have direct access to the simulation results (for 
further manipulation, debugging etc.) run this script with

$ python -i simmo_NW_array.py

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
max_num_BMs = 200

superstrate = objects.ThinFilm(period, height_nm = 'semi_inf',
    material = materials.Air, loss = False)

substrate  = objects.ThinFilm(period, height_nm = 'semi_inf',
    material = materials.SiO2_a, loss = False)

NW_diameter = 120
NW_array = objects.NanoStruct('2D_array', period, NW_diameter, height_nm = 2330, 
    inclusion_a = materials.Si_c, background = materials.Air, loss = True,    
    make_mesh_now = True, force_mesh = True, lc_bkg = 0.1, lc2= 2.0)

# Find num_BM for each simulation (wl) as num decreases w decreasing index contrast.
max_n = max([NW_array.inclusion_a.n(wl).real for wl in wavelengths])

def simulate_stack(light):
    num_BM = round(max_num_BMs * NW_array.inclusion_a.n(light.wl_nm).real/max_n)
    # num_BM = max_num_BMs
    
    ################ Evaluate each layer individually ##############
    sim_superstrate = superstrate.calc_modes(light)
    sim_substrate   = substrate.calc_modes(light)
    sim_NWs         = NW_array.calc_modes(light, num_BM = num_BM)

    ################ Evaluate full solar cell structure ##############
    """ Now when defining full structure order is critical and
    solar_cell list MUST be ordered from bottom to top!
    """

    stack = Stack((sim_substrate, sim_NWs, sim_superstrate))
    stack.calc_scat(pol = 'TE')

    return [stack] 


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



param_layer = NW_array # Specify the layer for which the parameters should be printed on figures.
params_string = plotting.gen_params_string(param_layer, last_light_object, max_num_BMs=max_num_BMs)

#### Example 1: simple multilayered stack.
stack_label = 0 # Specify which stack you are dealing with.
stack_wl_list = []
for i in range(len(wavelengths)):
    stack_wl_list.append(stacks_wl_list[i][stack_label])
active_layer_nu = 1

Efficiency = plotting.t_r_a_plots(stack_wl_list, wavelengths, params_string, 
    active_layer_nu=active_layer_nu, stack_label=stack_label) 
# Dispersion
plotting.omega_plot(stack_wl_list, wavelengths, params_string, stack_label=stack_label) 


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
