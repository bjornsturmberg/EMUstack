"""
    simo_053-plotting_k_space.py is a simulation example for EMUstack.

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
Here we investigate how efficiently a stack of 1D gratings excite diffraction orders.
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
num_cores = 5

# Remove results of previous simulations
plotting.clear_previous()

################ Light parameters #####################
azi_angles = np.linspace(-89,89,91)
wl = 1600
light_list  = [objects.Light(wl, max_order_PWs = 10, theta = p, phi = 0.0) \
    for p in azi_angles]


################ Grating parameters #####################
# The period must be consistent throughout a simulation!
period = 700

superstrate = objects.ThinFilm(period, height_nm = 'semi_inf',
    material = materials.Air, loss = False, world_1d =True)

substrate  = objects.ThinFilm(period, height_nm = 'semi_inf',
    material = materials.Air, loss = False, world_1d =True)

absorber    = objects.ThinFilm(period, height_nm = 10,
    material = materials.Material(1.0 + 0.05j), loss = True, world_1d =True)

grating_1 = objects.NanoStruct('1D_array', period, int(round(0.75*period)),
    height_nm = 2900, background = materials.Material(1.46 + 0.0j),
    inclusion_a = materials.Material(3.61 + 0.0j), loss = True,
    lc_bkg = 0.005)


def simulate_stack(light):

    ################ Evaluate each layer individually ##############
    sim_superstrate = superstrate.calc_modes(light)
    sim_substrate   = substrate.calc_modes(light)
    sim_absorber    = absorber.calc_modes(light)
    sim_grating_1   = grating_1.calc_modes(light)

    ###################### Evaluate structure ######################
    """ Now define full structure. Here order is critical and
        stack list MUST be ordered from bottom to top!
    """

    stack = Stack((sim_substrate, sim_absorber, sim_grating_1, sim_superstrate))
    stack.calc_scat(pol = 'TE')

    return stack


# Run in parallel across wavelengths.
pool = Pool(num_cores)
stacks_list = pool.map(simulate_stack, light_list)
# Save full simo data to .npz file for safe keeping!
np.savez('Simo_results', stacks_list=stacks_list)

######################## Post Processing ########################
# We can represent the strength with which different orders are excited
# in k-space.
plotting.t_func_k_plot_1D(stacks_list)
# This corresponds to Fig 2 of Handmer et al. Optics Lett. 35, 2010.
# (The PW_amplitudes plots correspond to Fig 1 of this paper).

# Lastly we also plot the transmission, reflection and absorption of each
# layer and the stack.
plotting.t_r_a_plots(stacks_list, xvalues=azi_angles)


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
