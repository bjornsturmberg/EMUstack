"""
    simo_070-ultrathin_limit.py is a simulation example for EMUstack.

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
Simulating an ultrathin film with a range of real and imaginary refractive
indices. Can we reach the theoretical limit of 0.5 absorption?
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
num_cores = 8

################ Light parameters #####################
wl     = 700
light  = objects.Light(wl, max_order_PWs = 0, theta = 0.0, phi = 0.0)

# The period must be consistent throughout a simulation!
period = 660

# Define each layer of the structure.
superstrate = objects.ThinFilm(period, height_nm = 'semi_inf',
    material = materials.Air, world_1d=True)
substrate   = objects.ThinFilm(period, height_nm = 'semi_inf',
    material = materials.Air, loss=False, world_1d=True)

n_min = 1
n_max = 10
num_n_re = 51
num_n_im = num_n_re
Re_n = np.linspace(n_min,n_max,num_n_re)
Im_n = np.linspace(n_max,n_min,num_n_im)
# Having lists run this way will ease plotting, as matshow plots from top left

def simulate_stack(Re):
    ################ Evaluate each layer individually ##############
    sim_superstrate = superstrate.calc_modes(light)
    sim_substrate   = substrate.calc_modes(light)

    # Re_stack = []
    # for Re in Re_n:
    Im_stack = []
    for Im in Im_n:
        TF_1 = objects.ThinFilm(period, height_nm = 10,
            material = materials.Material(Re + Im*1j))
        sim_TF_1 = TF_1.calc_modes(light)

        stack = Stack((sim_substrate, sim_TF_1, sim_superstrate))
        stack.calc_scat(pol = 'TM')

        Im_stack.append(stack)
        # Re_stack.append(Im_stack)

    return Im_stack

# Run wavelengths in parallel across num_cores CPUs using multiprocessing package.
pool = Pool(num_cores)
stacks_list = pool.map(simulate_stack, Re_n)
# # Save full simo data to .npz file for safe keeping!
# np.savez('Simo_results', stacks_list=stacks_list)

######################## Post Processing ########################

abs_mat = np.zeros((num_n_im,num_n_re))
for i in range(num_n_re):
    for j in range(num_n_im):
        abs_mat[j,i] = stacks_list[i][j].a_list[-1]

# Now plot as a function of Real and Imaginary refractive index.
# Requires a bit of manipulation of axis...
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
fig = plt.figure()
linesstrength = 3
font = 18
ax1 = fig.add_subplot(1,1,1)
mat = ax1.matshow(abs_mat,cmap=plt.cm.hot)
cbar1 = plt.colorbar(mat, extend='neither',alpha=1)
ax1.xaxis.set_ticks_position('bottom')
ax1.set_xticks(np.linspace(n_min,(num_n-1),n_max))
ax1.set_yticks(np.linspace(n_min,(num_n-1),n_max))
ax1.set_xticklabels([str(i) for i in np.linspace(n_min,n_max,n_max-n_min+1)])
ax1.set_yticklabels([str(i) for i in np.linspace(n_max,n_min,n_max-n_min+1)])
ax1.set_xlabel('Re(n)',fontsize=font)
ax1.set_ylabel('Im(n)',fontsize=font)
plt.title('Absorption of %(h)5.1f nm thick film @ wl = %(wl)5.1f'% \
    {'h' : stacks_list[0][0].heights_nm()[0],'wl' : wl})
plt.savefig('ultrathin_limit')



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
