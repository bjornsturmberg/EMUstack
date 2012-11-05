"""
The Ueber python script, the only one that needs to be edited to set up all 
simulation parameters.
Uses other python scripts to prime the simulation (interpolate raw data over chosen 
wavelengths etc.), then calls the fortran routine pcpv.exe for each wavelength giving 
it all the required details. It does this by spanning a new process for each wavelength,
keeping the total running instances to a maximum number (num_cores_to_use). Finally all 
results are collected in text files and the spectra are plotted. A log file is found in
python_log.txt
"""

import time
import datetime
import numpy as np
import subprocess
import sys
# import multiprocessing   as mp
sys.path.append("../PCPV/")

from calculate_ff      import calculate_ff
import clear_previous
import objects
import materials
from fortran_call      import scat_mats
from combine_layers    import net_scat_mats
import cat_n_clean
import plotting

start = time.time()
label_nu = 0
################ Simulation parameters ################

simo_para  = objects.Controls(debug = 0,max_order_PWs = 4, num_cores = 6)

# Remove results of previous simulations
clear_previous.clean('.txt')
clear_previous.clean('.pdf')
clear_previous.clean('.log')

################ Light parameters #####################
# wl_1     = 310
# wl_2     = 410   # < 400 Im(n(cSi)) large
# wl_3     = 600#1010  # > 1010 Im(n(cSi)) = 0
# wl_4     = 800#1110
# wl_5     = 990#1127
# no_wl_1  = 46    # 2 nm steps
# no_wl_2  = 610   # 1 nm steps
# no_wl_3  = 10    # 10 nm steps - total 667 wavelengths
# # Set up light objects
# wl_array_1  = np.linspace(wl_1, wl_2, no_wl_1)
# wl_array_2  = np.linspace(wl_2+1, wl_3, no_wl_2)
# wl_array_3  = np.linspace(wl_3+1, wl_4, no_wl_3)
# wavelengths = np.concatenate((wl_array_1,wl_array_2, wl_array_3, [wl_5]))
# light_list  = [objects.Light(wl) for wl in wavelengths]

# Single wavelength run
wl_super = 700
wavelengths = np.array([wl_super])
# wavelengths = np.linspace(310, 900, 10)
light_list  = [objects.Light(wl) for wl in wavelengths]


################ Scattering matrices (for distinct layers) ##############
""" Calculate scattering matrices for each distinct layer.
Calculated in the order listed below, however this does not influence final 
structure which is defined later
"""

# period must be consistent throughout simulation!!!
period      = 600

cover  = objects.ThinFilm(simo_period = period, height_1 = 'semi_inf',
    film_material = materials.Air, superstrate = materials.Air, 
    substrate = materials.Air,loss = True, label_nu = label_nu)
label_nu = scat_mats(cover, light_list, simo_para)


radius1     = 60
ff          = calculate_ff(period,radius1)
max_num_BMs = 50
# mesh = 'bj_can_a60_d600.mail'
grating_1 = objects.NanoStruct(period, ff, radius1,
    # mesh_file = mesh,
    set_ff = False, height_1 = 2200, height_2 = 2400, num_h = 1,
    inclusion_a = materials.SiO2_a, background = materials.Si_c,
    superstrate = materials.Air, substrate = materials.Air, nb_typ_el = 4, 
    make_mesh_now = True, force_mesh = False,
    lc_bkg = 0.1, lc2_surf= 1.9, lc3_inc1 = 1.9, lc4_inc2 = 1.9, lc5_inc3 = 1.9, lc6_inc4 = 1.9,
    posx = 0, posy = 0, max_num_BMs = max_num_BMs, label_nu = label_nu)
label_nu = scat_mats(grating_1, light_list, simo_para)



homo_film  = objects.ThinFilm(simo_period = period, height_1 = 500, height_2 = 1000, num_h = 1,
    film_material = materials.Si_c, superstrate = materials.Air, 
    substrate = materials.Air,loss = True, label_nu = label_nu)
# label_nu = scat_mats(homo_film, light_list, simo_para)


# will only ever use top scattering matrices for the bottom layer
bottom = objects.ThinFilm(simo_period = period, height_1 = 'semi_inf',
    film_material = materials.SiO2_a, superstrate = materials.Air, 
    loss = True, label_nu = label_nu)
label_nu = scat_mats(bottom, light_list, simo_para)



################ Construct full solar cell structure ##############
""" Now when defining full structure order is critical and
solar_cell list MUST be ordered from bottom to top!
"""
# starting from bottom, building upwards
solar_cell = [bottom, grating_1, cover]
# test that all structures have the same period!

net_scat_mats(solar_cell,wavelengths)


# ################ Concatinate results ################ 
# last_light_object = light_list.pop()
# if simo_para.traLambda  == 1:
#     cat_n_clean.c_c_tra(pol = last_light_object.pol)
# if simo_para.PrintOmega == 1:
#     cat_n_clean.c_c_omega()
# if simo_para.PropModes  == 1:
#     cat_n_clean.c_c_detA()
#     if solar_cell.loss   == False:
#         cat_n_clean.c_c_prop_modes()

# # # Plotting
# # plotting.average_spec('Absorptance','Av_Absorb', len(wavelengths), solar_cell.num_h)
# # plotting.average_spec('Transmittance','Av_Trans', len(wavelengths), solar_cell.num_h)
# # plotting.average_spec('Reflectance','Av_Reflec', len(wavelengths), solar_cell.num_h)
# # # Interpolate solar spectrum and calculate efficiency
# # Efficiency = plotting.irradiance('Av_Absorb', 'Weighted_Absorb', 'Av_Trans', 'Weighted_Trans',
# #  'Av_Reflec', 'Weighted_Reflec', radius1, radius2, period, ff)
# # # Plot averaged sprectra
# # spec_list = ['Av_Absorb', 'Av_Trans', 'Av_Reflec']
# # plotting.tra_plot('Spectra', spec_list, solar_cell, last_light_object,
# #     max_num_BMs, max_order_PWs, Efficiency)
# # # Plot weighted averaged sprectra
# # spec_list = ['Weighted_Absorb', 'Weighted_Trans', 'Weighted_Reflec']
# # plotting.tra_plot('Spectra_weighted', spec_list, solar_cell, last_light_object, 
# #     max_num_BMs, max_order_PWs, Efficiency)
# # # Plot dispersion diagrams
# # plotting.omega_plot(solar_cell, last_light_object, 
# #     max_num_BMs, max_order_PWs, Efficiency)

# # if solar_cell.num_h != 1:
# # # Calculate and plot efficiency as a function of height
# #     plotting.efficiency_h('Absorptance','Efficiency_h',wavelengths,len(wavelengths),
# #     solar_cell.num_h, simo_para.Animate)
# # # Plot absorptance as funtion of height
# #     plotting.height_plot('Spectra_height', 'Absorptance', solar_cell, last_light_object,
# #     max_num_BMs, max_order_PWs, 'Efficiency_h', Efficiency, solar_cell.num_h)


# # # Wraping up simulation by printing to screen and log file
# # print '\n*******************************************'
# # print 'The ultimate efficiency is %12.8f' % Efficiency
# # print '-------------------------------------------'

# # # Calculate and record the (real) time taken for simulation
# # elapsed = (time.time() - start)
# # hms     = str(datetime.timedelta(seconds=elapsed))
# # hms_string = 'Total time for simulation was \n \
# #     %(hms)s (%(elapsed)12.3f seconds)'% {
# #             'hms'       : hms,
# #             'elapsed'   : elapsed, }

# # python_log = open("python_log.log", "w")
# # python_log.write(hms_string)
# # python_log.close()

# # print hms_string
# # print '*******************************************'