"""
    simmo_circular_dichroism.py is a simulation script template for EMUstack.

    Copyright (C) 2013  Bjorn Sturmberg, Kokou Dossou, Felix Lawrence

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
Simulating circular dichroism effect in elliptic nano hole arrays
as in T Cao1 and Martin J Cryan doi:10.1088/2040-8978/14/8/085101.
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
num_cores = 15

# Remove results of previous simulations
plotting.clear_previous('.txt')
plotting.clear_previous('.pdf')
plotting.clear_previous('.gif')
plotting.clear_previous('.log')

################ Light parameters #####################
wl_1     = 300
wl_2     = 1000 
no_wl_1  = 21
# Set up light objects
wavelengths = np.linspace(wl_1, wl_2, no_wl_1)
light_list  = [objects.Light(wl,theta = 45, phi = 45, max_order_PWs = 2) for wl in wavelengths]



#period must be consistent throughout simulation!!!
period = 165
diam1  = 140
diam2  = 60
ellipticity = (float(diam1-diam2))/float(diam1)

Au_NHs = objects.NanoStruct('2D_array', period, diam1,ellipticity = ellipticity,height_nm = 60, 
    inclusion_a = materials.Air, background = materials.Au_drude, loss = True,    
    make_mesh_now = True, force_mesh = True, lc_bkg = 0.2, lc2= 5.0)

superstrate = objects.ThinFilm(period = period, height_nm = 'semi_inf',
    material = materials.Air, loss = True)
substrate = objects.ThinFilm(period = period, height_nm = 'semi_inf',
    material = materials.Air, loss = False)

num_BM = 50

def simulate_stack(light):  
    num_h = 21
    NH_heights = [60]#np.linspace(10,100,num_h)

    ################ Evaluate each layer individually ##############
    sim_superstrate = superstrate.calc_modes(light)
    sim_Au   = Au_NHs.calc_modes(light, num_BM = num_BM)
    sim_substrate = substrate.calc_modes(light)

# Loop over heights
    height_list = []
    # for h in NH_heights:
    stackSub = Stack((sim_substrate, sim_Au, sim_superstrate))#, heights_nm = ([h]))
    stackSub.calc_scat(pol = 'R Circ')
    stackSub2 = Stack((sim_substrate, sim_Au, sim_superstrate))#, heights_nm = ([h]))
    stackSub2.calc_scat(pol = 'L Circ')
    saveStack = Stack((sim_substrate, sim_Au, sim_superstrate))#, heights_nm = ([h]))
    
    a_CD = []
    t_CD = []
    r_CD = []
    for i in range(len(stackSub.a_list)):
        a_CD.append(stackSub.a_list.pop() - stackSub2.a_list.pop())
    for i in range(len(stackSub.t_list)):
        t_CD.append(stackSub.t_list.pop() - stackSub2.t_list.pop())
    for i in range(len(stackSub.r_list)):
        r_CD.append(stackSub.r_list.pop() - stackSub2.r_list.pop())
    saveStack.a_list = a_CD
    saveStack.t_list = t_CD
    saveStack.r_list = r_CD
    height_list.append(saveStack)
    # height_list.append(stackSub)

    return height_list


# Run in parallel across wavelengths.
pool = Pool(num_cores)
stacks_list = pool.map(simulate_stack, light_list)
# Save full simo data to .npz file for safe keeping!
simotime = str(time.strftime("%Y%m%d%H%M%S", time.localtime()))
np.savez('Simo_results'+simotime, stacks_list=stacks_list)

######################## Plotting ########################
last_light_object = light_list.pop()


#### Individual spectra of multilayered stack where one layer has many heights.
stack_label = 0
active_layer_nu = 1
# for h in range(num_h):
h = 0
gen_name = '_h-'
h_name = str(h)
additional_name = gen_name+h_name # You can add an arbitry string onto the end of the spectra filenames.
stack3_hs_wl_list = []
for i in range(len(wavelengths)):
    stack3_hs_wl_list.append(stacks_list[i][h])
plotting.t_r_a_plots(stack3_hs_wl_list, stack_label=stack_label, add_name = additional_name)


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
print ''