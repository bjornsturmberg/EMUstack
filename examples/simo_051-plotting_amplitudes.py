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
num_cores = 5

# Remove results of previous simulations
plotting.clear_previous('.txt')
plotting.clear_previous('.pdf')
plotting.clear_previous('.log')

################ Light parameters #####################
azi_angles = np.linspace(0,89,5)
wl = 1600
light_list  = [objects.Light(wl, max_order_PWs = 4, theta = p, phi = 0.0) for p in azi_angles]


################ Grating parameters #####################
# The period must be consistent throughout a simulation!
period = 700

superstrate = objects.ThinFilm(period, height_nm = 'semi_inf',
    material = materials.Air, loss = False)

substrate  = objects.ThinFilm(period, height_nm = 'semi_inf',
    material = materials.Air, loss = False)

absorber    = objects.ThinFilm(period, height_nm = 10,
    material = materials.Material(2.0 + 0.05j), loss = True)

grating_1 = objects.NanoStruct('1D_array', period, int(round(0.75*period)), height_nm = 2900, 
    background = materials.Material(1.46 + 0.0j), inclusion_a = materials.Material(3.61 + 0.0j), 
    loss = True, make_mesh_now = True, force_mesh = False, lc_bkg = 0.1, lc2= 3.0)

num_BM = 100


def simulate_stack(light):
   
    ################ Evaluate each layer individually ##############
    sim_superstrate = superstrate.calc_modes(light)
    sim_substrate   = substrate.calc_modes(light)
    sim_absorber    = absorber.calc_modes(light)
    sim_grating_1   = grating_1.calc_modes(light, num_BM = num_BM)

    ################ Evaluate full solar cell structure ##############
    """ Now when defining full structure order is critical and
    stack MUST be ordered from bottom to top!
    """

    stack = Stack((sim_substrate, sim_absorber, sim_grating_1, sim_superstrate))
    stack.calc_scat(pol = 'TE')

    return stack


# Run in parallel across wavelengths.
pool = Pool(num_cores)
stacks_list = pool.map(simulate_stack, light_list)

######################## Post Processing ########################
# We wish to know how evanescent orders are so we note the high refractive 
# index of grating.
n_H = 3.61

# We can plot the amplitudes of each transmitted plane wave order as a function of angle.
plotting.amps_of_orders(stacks_list, azi_angles, add_title='-default_substrate')
# By default this will plot the amplitudes in the substrate, however we can also give
# the index in the stack of a different homogeneous layer and calculate them here.
plotting.amps_of_orders(stacks_list, azi_angles, lay_interest=1)

# When many plane wave orders are included these last plots can become confusing,
# so instead one may wish to sum together the amplitudes of all propagating orders,
# of all evanescent orders, and all far-evanescent orders (which have in plane k>n_H * k0).
plotting.evanescent_merit(stacks_list, azi_angles, n_H, lay_interest=0)

# We can represent the strength with which different orders are excited in k-space.
plotting.t_func_k_plot_1D(stacks_list, n_H)
# This corresponds to Fig 2 of Handmer et al. Optics Lett. 35, 2010.
# (The amps_of_orders plots correspond to Fig 1 of this paper).

# Lastly we also plot the transmission, reflection and absorption of each layer and the stack.
plotting.t_r_a_plots(stacks_list, azi_angles)


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
