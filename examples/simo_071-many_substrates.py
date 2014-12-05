"""
    simo_070-many_substrates.py is a simulation script template for EMUstack.

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
Simulating solar cell efficiency of nanohole array as a function of
substrate refractive indices (keeping geometry fixed).
We also average over a range of thicknesses to remove sharp Fabry-Perot resonances.
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
num_cores = 4

# Remove results of previous simulations
plotting.clear_previous()

################ Light parameters #####################
wl_1     = 310
wl_2     = 1127
no_wl_1  = 3
# Set up light objects
wavelengths = np.linspace(wl_1, wl_2, no_wl_1)
light_list  = [objects.Light(wl, max_order_PWs = 2, theta = 0.0, phi = 0.0) \
    for wl in wavelengths]


# Period must be consistent throughout simulation!!!
period = 550

cover  = objects.ThinFilm(period = period, height_nm = 'semi_inf',
    material = materials.Air, loss = True)

sub_ns = np.linspace(1.0,4.0,100)

NW_diameter = 480
NWs = objects.NanoStruct('1D_array', period, NW_diameter, height_nm = 2330,
    inclusion_a = materials.Si_c, background = materials.Air, loss = True,
    make_mesh_now = True, force_mesh = True, lc_bkg = 0.17, lc2= 2.5)

def simulate_stack(light):
    ################ Evaluate each layer individually ##############
    sim_cover = cover.calc_modes(light)
    sim_NWs   = NWs.calc_modes(light)

    # Loop over substrates
    stack_list = []
    for s in sub_ns:
        sub = objects.ThinFilm(period = period, height_nm = 'semi_inf',
        material = materials.Material(s + 0.0j), loss = False)
        sim_sub = sub.calc_modes(light)

        # Loop over heights to wash out sharp FP resonances
        average_t = 0
        average_r = 0
        average_a = 0

        num_h = 21
        for h in np.linspace(2180,2480,num_h):
            stackSub = Stack((sim_sub, sim_NWs, sim_cover), heights_nm = ([h]))
            stackSub.calc_scat(pol = 'TE')
            average_t += stackSub.t_list[-1]/num_h
            average_r += stackSub.r_list[-1]/num_h
            average_a += stackSub.a_list[-1]/num_h
        stackSub.t_list[-1] = average_t
        stackSub.r_list[-1] = average_r
        stackSub.a_list[-1] = average_a
        stack_list.append(stackSub)

    return stack_list


# Run in parallel across wavelengths.
pool = Pool(num_cores)
stacks_list = pool.map(simulate_stack, light_list)
# Save full simo data to .npz file for safe keeping!
np.savez('Simo_results', stacks_list=stacks_list)



######################## Plotting ########################
eta = []
for s in range(len(sub_ns)):
    stack_label = s # Specify which stack you are dealing with.
    stack1_wl_list = []
    for i in range(len(wavelengths)):
        stack1_wl_list.append(stacks_list[i][stack_label])
    sub_n = sub_ns[s]
    Efficiency = plotting.t_r_a_plots(stack1_wl_list, ult_eta=True,
        stack_label=stack_label, add_name = str(s))
    eta.append(100.0*Efficiency[0])
    # Dispersion of structured layer is the same for all cases.
    if s == 0:
        plotting.omega_plot(stack1_wl_list, wavelengths, stack_label=stack_label)

# Now plot as a function of substrate refractive index.
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
fig = plt.figure()
linesstrength = 3
font = 18
ax1 = fig.add_subplot(1,1,1)
ax1.plot(sub_ns,eta, 'k-o', linewidth=linesstrength)
ax1.set_xlabel('Substrate refractive index',fontsize=font)
ax1.set_ylabel(r'$\eta$ (%)',fontsize=font)
plt.savefig('eta_substrates')

# Animate spectra as a function of substrates.
from os import system as ossys
delay = 30 # delay between images in gif in hundredths of a second
names = 'Total_Spectra_stack'
gif_cmd = 'convert -delay %(d)i +dither -layers Optimize -colors 16 \
%(n)s*.pdf %(n)s.gif'% {
'd' : delay, 'n' : names}
ossys(gif_cmd)
opt_cmd = 'gifsicle -O2 %(n)s.gif -o %(n)s-opt.gif'% {'n' : names}
ossys(opt_cmd)
rm_cmd = 'rm %(n)s.gif'% {'n' : names}
ossys(rm_cmd)



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

print '*******************************************'
print hms_string
print '*******************************************'
