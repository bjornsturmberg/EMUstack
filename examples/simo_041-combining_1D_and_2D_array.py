"""
    simo_041-combining_1D_and_2D_array.py is a simulation example for EMUstack.

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
Combining 1D gratings and 2D arrays in the same stack.
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
num_cores = 7

# Remove results of previous simulations
plotting.clear_previous()

################ Light parameters #####################
wl_1     = 310
wl_2     = 1127
no_wl_1  = 3
# Set up light objects
wavelengths = np.linspace(wl_1, wl_2, no_wl_1)
light_list  = [objects.Light(wl, max_order_PWs = 4, theta = 0.0, phi = 0.0) \
    for wl in wavelengths]


# Period must be consistent throughout simulation!!!
period = 600

# In this example we set the number of Bloch modes to use in the simulation
# Be default it is set to be slightly greater than the number of PWs.
num_BMs = 200

superstrate = objects.ThinFilm(period, height_nm = 'semi_inf',
    material = materials.Air, loss = False)

substrate  = objects.ThinFilm(period, height_nm = 'semi_inf',
    material = materials.SiO2_a, loss = False)

NW_diameter = 120
NW_array = objects.NanoStruct('2D_array', period, NW_diameter, height_nm = 2330,
    inclusion_a = materials.Si_c, background = materials.Air, loss = True,
    make_mesh_now = True, force_mesh = True, lc_bkg = 0.1, lc2= 2.0)

# We now create a 1D grating that is periodic in x only, but whose scattering
# matrices are constructed with of the 2D plane wave basis. This allows this layer
# to be combined with 2D_arrays.
grating = objects.NanoStruct('1D_array', period, int(round(0.75*period)), height_nm = 2900,
    background = materials.Material(1.46 + 0.0j), inclusion_a = materials.Material(5.0 + 0.0j),
    loss = True, lc_bkg = 0.01, world_1d = False)


def simulate_stack(light):

    ################ Evaluate each layer individually ##############
    sim_superstrate = superstrate.calc_modes(light)
    sim_substrate   = substrate.calc_modes(light)
    sim_NWs         = NW_array.calc_modes(light, num_BMs=num_BMs)

    ###################### Evaluate structure ######################
    """ Now define full structure. Here order is critical and
        stack list MUST be ordered from bottom to top!
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

# Plot the transmission, reflection and absorption.
plotting.t_r_a_plots(stacks_list, active_layer_nu=1)

# We also plot the dispersion relation for the NW layer.
plotting.omega_plot(stacks_list, wavelengths)

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
