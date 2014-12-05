"""
    simo_020-thin_film_mulilatered_stack.py is a simulation example for EMUstack.

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
EUMstack loves metal \m/
However, as we saw in the previous example the substrate layer must be lossless,
so that we can distinguish propagating waves from evanescent ones.
To terminate the stack with a metalic mirror we must make it finite, but very thick.
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

# Remove results of previous simulations.
plotting.clear_previous()

################ Simulation parameters ################
# Select the number of CPUs to use in simulation.
num_cores = 2

################ Light parameters #####################
wl_1     = 400
wl_2     = 800
no_wl_1  = 4
wavelengths = np.linspace(wl_1, wl_2, no_wl_1)
light_list  = [objects.Light(wl, max_order_PWs = 1, theta = 0.0, phi = 0.0)\
    for wl in wavelengths]

# The period must be consistent throughout a simulation!
period = 300

# Define each layer of the structure, as in last example.
superstrate = objects.ThinFilm(period, height_nm = 'semi_inf',
    material = materials.Air)
TF_2 = objects.ThinFilm(period, height_nm = 5e6,
    material = materials.InP, loss=False)
TF_3 = objects.ThinFilm(period, height_nm = 52,
    material = materials.Si_a)
# Realistically a few micron thick mirror would do the trick,
# but EMUstack is height agnostic.... so what the hell.
mirror = objects.ThinFilm(period, height_nm = 1e5,
    material = materials.Ag)
substrate   = objects.ThinFilm(period, height_nm = 'semi_inf',
    material = materials.Air)

def simulate_stack(light):
    ################ Evaluate each layer individually ##############
    sim_superstrate = superstrate.calc_modes(light)
    sim_mirror = mirror.calc_modes(light)
    sim_TF_2 = TF_2.calc_modes(light)
    sim_TF_3 = TF_3.calc_modes(light)
    sim_substrate   = substrate.calc_modes(light)
    ###################### Evaluate structure ######################
    """ Now define full structure. Here order is critical and
        stack list MUST be ordered from bottom to top!
    """
# Put semi-inf substrate below thick mirror so that propagating energy is defined.
    stack = Stack((sim_substrate, sim_mirror, sim_TF_3, sim_TF_2, sim_superstrate))
    stack.calc_scat(pol = 'TM')

    return stack

pool = Pool(num_cores)
stacks_list = pool.map(simulate_stack, light_list)
# Save full simo data to .npz file for safe keeping!
np.savez('Simo_results', stacks_list=stacks_list)

######################## Post Processing ########################
# The total transmission should be zero.
plotting.t_r_a_plots(stacks_list)

######################## Wrapping up ########################
print '\n*******************************************'
# Calculate and record the (real) time taken for simulation,
elapsed = (time.time() - start)
hms     = str(datetime.timedelta(seconds=elapsed))
hms_string = 'Total time for simulation was \n \
    %(hms)s (%(elapsed)12.3f seconds)'% {
            'hms'       : hms,
            'elapsed'   : elapsed, }
print hms_string
print '*******************************************'
print ''

# and store this info.
python_log = open("python_log.log", "w")
python_log.write(hms_string)
python_log.close()
