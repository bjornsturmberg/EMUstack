"""
    simo_011-single_interface.py is a simulation example for EMUstack.

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
Simulating an interface between 2 homogeneous, non-dispersive media.
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

################ Light parameters #####################
wl_1     = 500
wl_2     = 600
no_wl_1  = 4
# Set up light objects, starting with the wavelengths,
wavelengths = np.linspace(wl_1, wl_2, no_wl_1)
# and also specifying angles of incidence and refractive medium of semi-infinite
# layer that the light is incident upon (default value is n_inc = 1.0).
# Fields in homogeneous layers are expressed in a Fourier series of diffraction
# orders,where all orders within a radius of max_order_PWs in k-space are included.
light_list  = [objects.Light(wl, max_order_PWs = 1, theta = 0.0, phi = 0.0, \
    n_inc=1.5) for wl in wavelengths]

# Our structure must have a period, even if this is artificially imposed
# on a homogeneous thin film. What's more,
# it is critical that the period be consistent throughout a simulation!
period = 300

# Define each layer of the structure.
superstrate = objects.ThinFilm(period, height_nm = 'semi_inf',
    material = materials.Material(1.5 + 0.0j))
substrate   = objects.ThinFilm(period, height_nm = 'semi_inf',
    material = materials.Material(3.0 + 0.0j))

def simulate_stack(light):
    ################ Evaluate each layer individually ##############
    sim_superstrate = superstrate.calc_modes(light)
    sim_substrate   = substrate.calc_modes(light)
    ###################### Evaluate structure ######################
    """ Now define full structure. Here order is critical and
        stack list MUST be ordered from bottom to top!
    """

    stack = Stack((sim_substrate, sim_superstrate))
    # Calculate scattering matrices of the stack (for all polarisations).
    stack.calc_scat(pol = 'TE') # Incident light has TE polarisation,
    # which only effects the net transmission etc, not the matrices.

    return stack

stacks_list = map(simulate_stack, light_list)
# Save full simo data to .npz file for safe keeping!
np.savez('Simo_results', stacks_list=stacks_list)

# Calculation of the modes and scattering matrices of each layer
# as well as the scattering matrices of the interfaces of the stack
# is complete.
# From here on we can print, plot or manipulate the results.

# Alternatively, you may wish to finish the simo file here,
# and be output into an interactive python instance were you
# have access to all simulation objects and results for further
# manipulation. In this case you run this file as
# $ python -i simo_010-single_interface.py
# In this session the docstrings of objects/classes/methods
# can be accessed by typing

# >>> from pydoc import help
# >>> help(objects.Light)

# where we have accessed the docstring of the Light class from objects.py


######################## Post Processing ########################
# We can retrieve the propagation constants (k_z) of each layer.
# Let's print the values at the short wavelength in the superstrate,
wl_num = 0
lay = 1
betas = stacks_list[wl_num].layers[lay].k_z
print 'k_z of superstrate \n', betas
# and save the values for the longest wavelength for the substrate.
wl_num = -1
lay = 0
betas = stacks_list[wl_num].layers[lay].k_z
np.savetxt('Substrate_k_zs.txt', betas.view(float).reshape(-1, 2))
# Note that saving to txt files is slower than saving data as .npz
# However txt files may be easily read by other programs...


# We can also access the scattering matrices of individual layers,
# and of interfaces of the stack.
# For instance the reflection scattering matrix off the top
# of the substrate when considered as an isolated layer.
wl_num = -1
lay = 0
R12_sub = stacks_list[wl_num].layers[lay].R12
print 'R12 of substrate \n', R12_sub

# The reflection matrix for the reflection off the top of the
# superstrate-substrate interface meanwhile is a property of the stack.
R_interface = stacks_list[wl_num].R_net
# Let us plot this matrix in greyscale.
plotting.vis_scat_mats(R_interface)
# Since all layers are homogeneous this matrix should only have non-zero
# entries on the diagonal.

# Lastly, we can also plot the transmission, reflection, absorption
# of each layer and of the stack as a whole.
plotting.t_r_a_plots(stacks_list)

# p.s. we'll keep an eye on the time...
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
