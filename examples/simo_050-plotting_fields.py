"""
    simo_050-plotting_fields.py is a simulation example for EMUstack.

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
plotting.clear_previous('.npz')
plotting.clear_previous('.txt')
plotting.clear_previous('.pdf')
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

    return stack


# Run in parallel across wavelengths.
pool = Pool(num_cores)
stacks_list = pool.map(simulate_stack, light_list)
# Save full simo data to .npz file for safe keeping!
np.savez('Simo_results', stacks_list=stacks_list)


######################## Plotting ########################
#### Example 1: simple multilayered stack.

plotting.t_r_a_plots(stack_wl_list) 

# Dispersion
plotting.omega_plot(stack_wl_list, wavelengths) 


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
