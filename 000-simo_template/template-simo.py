# Des

import time
import datetime
import numpy as np
import multiprocessing 	 as mp
import subprocess
import sys
sys.path.append("../PCPV/")

import clear_previous
import objects
import materials
from   fortran_call     import Simmo
import cat_n_clean
import plotting

start = time.time()

################ Simo specific parameters ################
# solar cell parameters
radius1 = 150
radius2 = 150
period  = 1190
ff      = 0.20
# light parameters
wl_1     = 310
wl_2     = 410   # < 400 Im(n(cSi)) large
wl_3     = 1010  # > 1010 Im(n(cSi)) = 0
wl_4     = 1110
wl_5     = 1127
no_wl_1  = 46    # 2 nm steps
no_wl_2  = 610   # 1 nm steps
no_wl_3  = 10    # 10 nm steps - total 667 wavelengths
# simulation parameters
max_num_BMs   = 200
var_BM_min    = 370
max_order_PWs = 3
# number of cpus to use
num_cores_to_use = 10
# number of cpus to leave free
# num_cores_sahand_gets = 8
################ Simo specific parameters ################


# Remove results of previous simulations
clear_previous.clean('.txt')
clear_previous.clean('.pdf')
clear_previous.clean('.log')

# Set up solar cell
mesh = '../PCPV/Data/2by2_ff20_tau_01.mail'
solar_cell  = objects.SolarCell(radius1, radius2, period, ff, mesh)#, inclusion=materials.air, substrate=materials.SiO2)

# Set up light objects
wl_array_1  = np.linspace(wl_1, wl_2, no_wl_1)
wl_array_2  = np.linspace(wl_2+1, wl_3, no_wl_2)
wl_array_3  = np.linspace(wl_3+1, wl_4, no_wl_3)
wavelengths = np.concatenate((wl_array_1,wl_array_2, wl_array_3, [wl_5]))
light_list  = [objects.Light(wl) for wl in wavelengths]

# Simulation controls
other_para  = objects.Controls()

# Interpolate refractive indecies over wavelength array
materials.interp_all(wavelengths)

# List of simulations to calculate, with full arguments
simmo_list = []
p 		   = 1
for light in light_list:
	simmo_list += [Simmo(solar_cell, light, other_para, 
		max_num_BMs, var_BM_min, max_order_PWs, p)]
	p += 1

# # Start all simulations running
# [s.run()  for s in simmo_list]
# # Wait until every simulation is finished
# [s.wait() for s in simmo_list]

def run_command_in_shell(command):
	return subprocess.call(command, shell = True)

command_list = [s.fortran_command_str() for s in simmo_list]
# pool = mp.Pool(mp.cpu_count() - num_cores_sahand_gets)
pool = mp.Pool(num_cores_to_use)
pool.map(run_command_in_shell, command_list)


# Concatinate results
if other_para.traLambda  == 1:
	cat_n_clean.c_c_tra()
if other_para.PrintOmega == 1:
	cat_n_clean.c_c_omega()
if other_para.PropModes  == 1:
	cat_n_clean.c_c_detA()
	if solar_cell.loss   == False:
		cat_n_clean.c_c_prop_modes()

# Interpolate solar spectrum and calculate efficiency
Efficiency = plotting.irradiance('Absorptance','../PCPV/Data/ASTM_1_5_spectrum','Weighted_Abs')
# Plot sprectra
spec_list = ['Weighted_Abs', 'Absorptance', 'Transmittance', 'Reflectance']
last_light_object = light_list.pop()
plotting.tra_plot(spec_list, solar_cell, last_light_object, max_num_BMs, max_order_PWs, Efficiency)



# Wraping up simulation by printing to screen and log file
print '\n*******************************************'
print 'The ultimate efficiency is %12.8f' % Efficiency
print '-------------------------------------------'

# Calculate and record the time taken for simulation
elapsed = (time.time() - start)
hms     = str(datetime.timedelta(seconds=elapsed))
hms_string = 'Total time for simulation was \n \
	%(hms)s (%(elapsed)12.3f seconds)'% {
			'hms' 	    : hms,
			'elapsed'	: elapsed, }

python_log = open("python_log.log", "w")
python_log.write(hms_string)
python_log.close()

print hms_string
print '*******************************************'
