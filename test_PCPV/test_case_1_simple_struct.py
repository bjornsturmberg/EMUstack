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

simo_para  = objects.Controls(debug = 0, max_order_PWs = 0, num_cores = 1)

# Remove results of previous simulations
clear_previous.clean('.txt')
clear_previous.clean('.pdf')
clear_previous.clean('.log')

################ Light parameters #####################

# # Set up light objects
# wavelengths = [310, 410, 531.2007, 707.495,  881.786, 987.9632]
# light_list  = [objects.Light(wl) for wl in wavelengths]
# Single wavelength run
wl_super =  500.0
wavelengths = np.array([wl_super])
light_list  = [objects.Light(wl) for wl in wavelengths]


################ Scattering matrices (for distinct layers) ##############
""" Calculate scattering matrices for each distinct layer.
Calculated in the order listed below, however this does not influence final 
structure which is defined later
"""

# period must be consistent throughout simulation!!!
period  = 600

cover  = objects.ThinFilm(period = period, height_1 = 'semi_inf',
    film_material = materials.Air, superstrate = materials.Air, 
    substrate = materials.Air,loss = False, label_nu = 0)
scat_mats(cover, light_list, simo_para)

radius1 = 60
max_num_BMs = 30
grating_1 = objects.NanoStruct('NW_array', period, radius1, square = False,
    set_ff = False, height_1 = 1000, height_2 = 1000, num_h = 1,
    inclusion_a = materials.Si_c, background = materials.Air,
    superstrate = materials.Air, substrate = materials.Air, loss = True, nb_typ_el = 4, 
    make_mesh_now = True, force_mesh = True,
    lc_bkg = 0.15, lc2= 1.5, lc3= 1.5,
    max_num_BMs = max_num_BMs, label_nu = 1)
simmo_list = scat_mats(grating_1, light_list, simo_para)


# will only ever use top scattering matrices for the bottom layer
bottom = objects.ThinFilm(period = period, height_1 = 'semi_inf',
    film_material = materials.SiO2_a, superstrate = materials.Air, 
    loss = False, label_nu = 2)
scat_mats(bottom, light_list, simo_para)



################ Construct & solve for full solar cell structure ##############
""" Now when defining full structure order is critical and
solar_cell list MUST be ordered from bottom to top!
"""
solar_cell = [bottom, grating_1, cover]

net_scat_mats(solar_cell, wavelengths, simo_para)

def results_match_reference(filename):
    reference = np.loadtxt("ref/case_1/" + filename)
    result    = np.loadtxt(filename)
    np.testing.assert_allclose(result, reference, 1e-7, 1e-6, filename)

def test_results():
    result_files = ("Absorptance.txt", "Lay_Absorb_0.txt",
                    "Lay_Trans_0.txt", "Reflectance.txt",
                    "Transmittance.txt", "beta_st0000.txt",
                    "beta_st0002.txt", "omega_Ft_st0001.txt",
                    "omega_Fz_st0001.txt", "omega_pol_st0001.txt",
                    "omega_st0001.txt",)
    for f in result_files:
        yield results_match_reference, f
