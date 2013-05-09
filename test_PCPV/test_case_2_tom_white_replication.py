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
import plotting
from stack import *
from numpy.testing import assert_allclose as assert_ac
from numpy.testing import assert_equal
import testing

start = time.time()
label_nu = 0
################ Simulation parameters ################

simo_para  = objects.Controls(debug = 0,max_order_PWs = 5, num_cores = 7,
    PrintAll = 0, Checks = 0, PrintSolution = 0, PrintSupModes = 0)
# Remove results of previous simulations
clear_previous.clean('.txt')
clear_previous.clean('.pdf')
clear_previous.clean('.log')

################ Light parameters #####################
wl_1     = 900
wl_2     = 1200
no_wl_1  = 2#8
# Set up light objects
wavelengths = np.linspace(wl_1, wl_2, no_wl_1)
light_list  = [objects.Light(wl) for wl in wavelengths]
# # Single wavelength run
# wl_super = 1000
# wavelengths = np.array([wl_super])
# light_list  = [objects.Light(wl) for wl in wavelengths]


# period must be consistent throughout simulation!!!
period = 120

cover  = objects.ThinFilm(period = period, height_1 = 'semi_inf',
    material = materials.Material(3.5 + 0.0j), loss = True)

homo_film  = objects.ThinFilm(period = period, height_1 = 5,
    material = materials.Material(3.6 + 0.27j), loss = True)

bottom = objects.ThinFilm(period = period, height_1 = 'semi_inf',
    material = materials.Air, loss = False)

grating_1 = objects.NanoStruct('1D_grating', period, 100, height_1 = 25, 
    inclusion_a = materials.Ag, background = materials.Material(1.5 + 0.0j), loss = True, nb_typ_el = 4, 
    make_mesh_now = True, force_mesh = True, lc_bkg = 0.1, lc2= 4.0)

mirror = objects.ThinFilm(period = period, height_1 = 100,
    material = materials.Ag, loss = True)


stack_list = []
# Find num_BM for each simulation in a somewhat arbitrary way
# Maybe roll this out into a Bjorn-specific function
max_n = max([grating_1.inclusion_a.n(wl).real for wl in wavelengths])
max_num_BMs = 120

for light in light_list:

    num_BM = round(max_num_BMs * grating_1.inclusion_a.n(light.Lambda).real/max_n)

    ################ Scattering matrices (for distinct layers) ##############
    """ Calculate scattering matrices for each distinct layer.
    Calculated in the order listed below, however this does not influence final 
    structure which is defined later
    """


    sim_cover = cover.calc_modes(light, simo_para)
    sim_homo_film = homo_film.calc_modes(light, simo_para)
    sim_bot = bottom.calc_modes(light, simo_para)
    sim_grat1 = grating_1.calc_modes(light, simo_para, num_BM = num_BM)
    sim_mirror = mirror.calc_modes(light, simo_para)


    ################ Construct & solve for full solar cell structure ##############
    """ Now when defining full structure order is critical and
    solar_cell list MUST be ordered from bottom to top!
    """
    stack = Stack((sim_bot, sim_mirror, sim_grat1, sim_homo_film, sim_cover))
    stack.calc_scat(pol = 'TE')
    stack_list.append(stack)

t_r_a_plots(stack_list)


# solar_cell = [bottom, mirror, grating_1, homo_film, cover]
# specify which layer is the active one (where absorption generates charge carriers)
active_layer = homo_film
act_lay_nu   = 0
# specify which layer for which the parameters should be printed on figures
lay_print_params = grating_1


# net_scat_mats(solar_cell, wavelengths, simo_para)

################# Efficiency & weighted spectra for active layer ################
plotting.average_spec('Lay_Absorb_%d' % act_lay_nu, 'Av_Absorb',  
    len(wavelengths), 1)
plotting.average_spec('Lay_Trans_%d'  % act_lay_nu, 'Av_Trans',   
    len(wavelengths), 1)
plotting.average_spec('Reflectance', 'Av_Reflec',  
    len(wavelengths), 1)
# Interpolate solar spectrum and calculate efficiency
Efficiency = plotting.irradiance('Av_Absorb', 'Weighted_Absorb', 'Av_Trans', 'Weighted_Trans',
 'Av_Reflec', 'Weighted_Reflec', lay_print_params.radius1, period, ff = lay_print_params.ff)
# Plot averaged sprectra
last_light_object = light_list.pop()
spec_list = ['Av_Absorb', 'Av_Trans', 'Av_Reflec']
plotting.tra_plot('Spectra', spec_list, lay_print_params, last_light_object,
    max_num_BMs, simo_para.max_order_PWs, Efficiency)
# Plot weighted averaged sprectra
spec_list = ['Weighted_Absorb', 'Weighted_Trans', 'Weighted_Reflec']
plotting.tra_plot('Spectra_weighted', spec_list, lay_print_params, last_light_object, 
    max_num_BMs, simo_para.max_order_PWs, Efficiency)
# # Plot dispersion diagrams for each layer
# plotting.omega_plot(solar_cell, lay_print_params, last_light_object, 
#     max_num_BMs, simo_para.max_order_PWs, Efficiency)




# Wrapping up simulation by printing to screen and log file
print '\n*******************************************'
print 'The ultimate efficiency is %12.8f' % Efficiency
print '-------------------------------------------'

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



# # SAVE DATA AS REFERENCE
# # Only run this after changing what is simulated - this
# # generates a new set of reference answers to check against
# # in the future
# testing.save_reference_data("case_2", stack_list)



def results_match_reference(filename):
    reference = np.loadtxt("ref/case_2/" + filename)
    result    = np.loadtxt(filename)
    np.testing.assert_allclose(result, reference, 1e-6, 1e-6, filename)

def test_txt_results():
    result_files = (
        "Lay_Trans_0.txt",
        "Av_Absorb.txt",
        "Lay_Trans_1.txt",      "Reflectance.txt",
        "Av_Reflec.txt",        "Efficiency.txt",
        "Lay_Trans_2.txt",      "Transmittance.txt",
        "Av_Trans.txt",         "Lay_Absorb_0.txt",
        "Weighted_Absorb.txt",
        "Lay_Absorb_1.txt",
        "Weighted_Reflec.txt",
        "Lay_Absorb_2.txt",
        "Weighted_Trans.txt",
    )
    for f in result_files:
        yield results_match_reference, f

def test_stack_list_matches_saved(casefile_name = 'case_2', stack_list = stack_list, rtol = 1e-6, atol = 4e-6):
    ref = np.load("ref/%s.npz" % casefile_name)
    yield assert_equal, len(stack_list), len(ref['stack_list'])
    for stack, rstack in zip(stack_list, ref['stack_list']):
        yield assert_equal, len(stack.layers), len(rstack['layers'])
        lbl_s = "wl = %f, " % stack.layers[0].light.Lambda
        for i, (lay, rlay) in enumerate(zip(stack.layers, rstack['layers'])):
            lbl_l = lbl_s + "lay %i, " % i
            yield assert_ac, lay.R12,  rlay['R12'],  rtol, atol, lbl_l + 'R12'
            yield assert_ac, lay.T12,  rlay['T12'],  rtol, atol, lbl_l + 'T12'
            yield assert_ac, lay.R21,  rlay['R21'],  rtol, atol, lbl_l + 'R21'
            yield assert_ac, lay.T21,  rlay['T21'],  rtol, atol, lbl_l + 'T21'
            yield assert_ac, lay.beta, rlay['beta'], rtol, atol, lbl_l + 'beta'
            #TODO: assert_ac(lay.sol1, rlay['sol1'])
        #TODO: assert_ac(stack.R_net, rstack['R_net'])
        #TODO: assert_ac(stack.T_net, rstack['T_net'])

def test_final_absorptance_for_bjorn():
    results_match_reference("Absorptance.txt")
