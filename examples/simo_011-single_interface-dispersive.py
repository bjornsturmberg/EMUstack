"""
    simo_011-single_interface-dispersive.py is a simulation example for EMUstack.

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
Simulating an interface between 2 homogeneous, dispersive media.
We use multiple CPUs.
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

# We begin by remove all results of previous simulations.
plotting.clear_previous()

################ Simulation parameters ################
# Select the number of CPUs to use in simulation.
num_cores = 2

################ Light parameters #####################
wl_1     = 400
wl_2     = 800
no_wl_1  = 4
# Set up light objects (no need to specifiy n_inc as light incident from
# Air with n_inc = 1.0).
wavelengths = np.linspace(wl_1, wl_2, no_wl_1)
light_list  = [objects.Light(wl, max_order_PWs = 1, theta = 0.0, phi = 0.0) \
    for wl in wavelengths]

# The period must be consistent throughout a simulation!
period = 300

# Define each layer of the structure, now with dispersive media.
# The refractive indices are interpolated from tabulated data.
superstrate = objects.ThinFilm(period, height_nm = 'semi_inf',
    material = materials.Air)
substrate   = objects.ThinFilm(period, height_nm = 'semi_inf',
    material = materials.SiO2_a) # Amorphous silica

def simulate_stack(light):
    ################ Evaluate each layer individually ##############
    sim_superstrate = superstrate.calc_modes(light)
    sim_substrate   = substrate.calc_modes(light)
    ###################### Evaluate structure ######################
    """ Now define full structure. Here order is critical and
        stack list MUST be ordered from bottom to top!
    """

    stack = Stack((sim_substrate, sim_superstrate))
    stack.calc_scat(pol = 'TM') # This time TM polarised light is incident.

    return stack

# Run wavelengths in parallel across num_cores CPUs using multiprocessing package.
pool = Pool(num_cores)
stacks_list = pool.map(simulate_stack, light_list)
# Save full simo data to .npz file for safe keeping!
np.savez('Simo_results', stacks_list=stacks_list)

######################## Post Processing ########################
# This time let's visualise the net Transmission scattering matrix,
# which describes the propagation of light all the way from the superstrate into
# the substrate. When studying diffractive layers it is useful to know how many
# of theplane waves of the substrate are propagating, so lets include this.
wl_num = -1
T_net = stacks_list[wl_num].T_net
nu_prop = stacks_list[wl_num].layers[0].num_prop_pw_per_pol
plotting.vis_scat_mats(T_net, nu_prop_PWs=nu_prop)

# Let's just plot the spectra and see the effect of changing refractive indices.
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
