"""
    simo_051-plotting_fields_2d.py is a simulation example for EMUstack.

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
Show how to plot electric fields.
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
wl     = 615
light_list  = [objects.Light(wl, max_order_PWs = 15, theta = 0.0, phi = 0.0)]

# Period must be consistent throughout simulation!!!
period = 600

superstrate = objects.ThinFilm(period, height_nm = 'semi_inf',
    material = materials.Air, loss = False)

substrate  = objects.ThinFilm(period, height_nm = 'semi_inf',
    material = materials.Air, loss = False)

spacer  = objects.ThinFilm(period, height_nm = 200,
    material = materials.SiO2_a, loss = True)

NW_diameter = 120
NW_array = objects.NanoStruct('2D_array', period, NW_diameter,
    height_nm = 2330, inclusion_a = materials.Si_c, background = materials.Air,
    loss = True, make_mesh_now = True, force_mesh = True, lc_bkg = 0.1,
    lc2= 2.0, plotting_fields = True)


def simulate_stack(light):

    ################ Evaluate each layer individually ##############
    sim_superstrate = superstrate.calc_modes(light)
    sim_substrate   = substrate.calc_modes(light)
    sim_NWs         = NW_array.calc_modes(light)
    sim_spacer      = spacer.calc_modes(light)

    ###################### Evaluate structure ######################
    """ Now define full structure. Here order is critical and
        stack list MUST be ordered from bottom to top!
    """

    stack = Stack((sim_substrate, sim_spacer, sim_NWs, sim_superstrate))
    stack.calc_scat(pol = 'TE')

    return stack


# Run in parallel across wavelengths.
pool = Pool(num_cores)
stacks_list = pool.map(simulate_stack, light_list)
# Save full simo data to .npz file for safe keeping!
np.savez('Simo_results', stacks_list = stacks_list)


######################## Plotting ########################

# Plot fields on slices through stack along the x & y axis,
# and along the diagonals.
# This is done through all layers of the stack and saved as png files.
#
# Note that all field plots of previous simulations are deleted! Move any
# results that you wish to keep into a different folder, ideally copying the
# whole simo directory to future reference to simo parameters.
#
plotting.fields_vertically(stacks_list)

# Plot fields in the x-y plane at a list of specified heights.
plotting.fields_in_plane(stacks_list, lay_interest = 2, z_values = [0.0, 2.0])
plotting.fields_in_plane(stacks_list, lay_interest = 1, z_values = [1.0, 3.2])

# Plot fields inside nanostructures in 3D which are viewed using gmsh.
plotting.fields_3d(stacks_list, lay_interest = 2)

# Save electric field values (all components) at a list of selected point.
plotting.field_values(stacks_list, lay_interest = 0, xyz_values = [(4.0, 2.5, 7.0), (1.0, 1.5, 3.0)])

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
