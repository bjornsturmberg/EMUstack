"""
    simo_030-1D_grating.py is a simulation example for EMUstack.

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
Simulating a lamellar grating that is periodic in x only.
For this simulation EMUstack uses the 1D diffraction orders for the basis
of the plane waves and carries out a 1D FEM calculation for the modes of
the grating.
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
no_wl_1  = 2
wavelengths = np.linspace(wl_1, wl_2, no_wl_1)
light_list  = [objects.Light(wl, max_order_PWs = 5, theta = 0.0, phi = 0.0) for wl in wavelengths]

# The period must be consistent throughout a simulation!
period = 300

# Define each layer of the structure
# We need to inform EMUstack at this point that all layers in the stack will
# be at most be periodic in one dimension (i.e. there are no '2D_arrays's).
superstrate = objects.ThinFilm(period, height_nm = 'semi_inf', world_1d=True,
    material = materials.Air)

substrate   = objects.ThinFilm(period, height_nm = 'semi_inf', world_1d=True,
    material = materials.Air)
# Define 1D grating that is periodic in x and contains 3 interleaved inclusions.
# Inclusion_a is in the center of the unit cell. Inclusions 2 and 3 have
# diameters diameter2, diameter3, and are of material inclusion_b.
# Inclusion 1 is still centered in the center and by default all inclusions are
# seperated by period/(# inclusions) so in this case perid/3.
# See Fortran Backends section of tutorial for more details.
grating = objects.NanoStruct('1D_array', period, int(round(0.05*period)),
    diameter2 = int(round(0.17*period)), diameter3 = int(round(0.03*period)),
    diameter4 = int(round(0.07*period)), height_nm = 2900,
    background = materials.Material(1.46 + 0.0j), inclusion_a = materials.Material(5.0 + 0.0j),
    inclusion_b = materials.Material(3.0 + 0.0j),
    loss = True, lc_bkg = 0.0071)
# To instead seperate the inclusions with an equal distance between their edges use
# the Keyword Arg edge_spacing = True.
grating = objects.NanoStruct('1D_array', period, int(round(0.15*period)),
    diameter2 = int(round(0.27*period)), diameter3 = int(round(0.03*period)),
    edge_spacing = True, height_nm = 2900,
    background = materials.Material(1.46 + 0.0j), inclusion_a = materials.Material(5.0 + 0.0j),
    inclusion_b = materials.Material(3.0 + 0.0j),
    loss = True, lc_bkg = 0.0071)

def simulate_stack(light):
    ################ Evaluate each layer individually ##############
    sim_superstrate = superstrate.calc_modes(light)
    sim_grating     = grating.calc_modes(light)
    sim_substrate   = substrate.calc_modes(light)
    ###################### Evaluate structure ######################
    """ Now define full structure. Here order is critical and
        stack list MUST be ordered from bottom to top!
    """

    stack = Stack((sim_substrate, sim_grating, sim_superstrate))
    stack.calc_scat(pol = 'TE')

    return stack

pool = Pool(num_cores)
# stacks_list = pool.map(simulate_stack, light_list)
stacks_list = map(simulate_stack, light_list)
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
