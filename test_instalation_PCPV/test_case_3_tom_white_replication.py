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
import sys
from multiprocessing import Pool
sys.path.append("../PCPV/")

import clear_previous
import objects
import materials
import plotting
from stack import *
from numpy.testing import assert_allclose as assert_ac
from numpy.testing import assert_equal
import testing


# The following should be in the function "setup_module()",
# but unfortunately simulate_stack is defined in a lazy-but-easy
# way: the structures are inherited rather than passed in.

################ Simulation parameters ################

# Number of CPUs to use im simulation
num_cores = 7
# # Alternatively specify the number of CPUs to leave free on machine
# leave_cpus = 4 
# num_cores = mp.cpu_count() - leave_cpus

# Remove results of previous simulations
clear_previous.clean('.txt')
clear_previous.clean('.pdf')
# clear_previous.clean('.log')

################ Light parameters #####################
# wl_1     = 900
# wl_2     = 1050
# no_wl_1  = 2
# # Set up light objects
# wavelengths = np.linspace(wl_1, wl_2, no_wl_1)
# light_list  = [objects.Light(wl, max_order_PWs = 3) for wl in wavelengths]
# # Single wavelength run
wl_super = 1000
wavelengths = np.array([wl_super])
light_list  = [objects.Light(wl) for wl in wavelengths]


# period must be consistent throughout simulation!!!
period = 120

cover  = objects.ThinFilm(period = period, height_nm = 'semi_inf',
    material = materials.Material(3.5 + 0.0j), loss = True)

homo_film  = objects.ThinFilm(period = period, height_nm = 5,
    material = materials.Material(3.6 + 0.27j), loss = True)

bottom = objects.ThinFilm(period = period, height_nm = 'semi_inf',
    material = materials.Air, loss = False)

grating_diameter = 100
grating_1 = objects.NanoStruct('1D_grating', period, grating_diameter, height_nm = 25, 
    inclusion_a = materials.Ag, background = materials.Material(1.5 + 0.0j), loss = True,
    make_mesh_now = True, force_mesh = True, lc_bkg = 0.05, lc2= 4.0)

mirror = objects.ThinFilm(period = period, height_nm = 100,
    material = materials.Ag, loss = True)


stack_list = []
num_BM = 90

def simulate_stack(light):
    
    ################ Evaluate each layer individually ##############
    sim_cover = cover.calc_modes(light)
    sim_homo_film = homo_film.calc_modes(light)
    sim_bot = bottom.calc_modes(light)
    sim_grat1 = grating_1.calc_modes(light, num_BM = num_BM)
    sim_mirror = mirror.calc_modes(light)

    ################ Evaluate full solar cell structure ##############
    """ Now when defining full structure order is critical and
    solar_cell list MUST be ordered from bottom to top!
    """
    stack = Stack((sim_bot, sim_mirror, sim_grat1, sim_homo_film, sim_cover))
    stack.calc_scat(pol = 'TE')

    return stack

def setup_module(module):
    start = time.time()

    # Run in parallel across wavelengths.
    # This has to be in a setup_module otherwise nosetests will crash :(
    pool = Pool(2)
    module.stack_list = pool.map(simulate_stack, light_list)
    # # Run one at a time
    # module.stack_list = map(simulate_stack, light_list)

    last_light_object = light_list.pop()
    param_layer = grating_1 # Specify the layer for which the parameters should be printed on figures.
    params_string = plotting.gen_params_string(param_layer, last_light_object, max_num_BMs=num_BM)

    active_layer_nu = 3 # Specify which layer is the active one (where absorption generates charge carriers).
    # Plot total transmission, reflection, absorption & that of each layer. 
    # Also calculate efficiency of active layer.
    Efficiency = plotting.t_r_a_plots(stack_list, wavelengths, params_string, 
        active_layer_nu=active_layer_nu) 



    # # SAVE DATA AS REFERENCE
    # # Only run this after changing what is simulated - this
    # # generates a new set of reference answers to check against
    # # in the future
    # testing.save_reference_data("case_3", stack_list)



def results_match_reference(filename):
    rtol = 1e-6
    atol = 1e-2
    reference = np.loadtxt("ref/case_3/" + filename)
    result    = np.loadtxt(filename)
    np.testing.assert_allclose(result, reference, rtol, atol, filename)

def test_txt_results():
    result_files = (
        "Absorptance_stack1.txt",
        "Lay_Absorb_0_stack1.txt",
        "Lay_Absorb_1_stack1.txt",
        "Lay_Absorb_2_stack1.txt",
        "Lay_Trans_0_stack1.txt",
        "Lay_Trans_1_stack1.txt",
        "Lay_Trans_2_stack1.txt",
        "Reflectance_stack1.txt",
        "Transmittance_stack1.txt",
    )
    for f in result_files:
        yield results_match_reference, f

def test_stack_list_matches_saved(casefile_name = 'case_3'):
    rtol = 1e-0
    atol = 1e-0
    ref = np.load("ref/%s.npz" % casefile_name)
    yield assert_equal, len(stack_list), len(ref['stack_list'])
    for stack, rstack in zip(stack_list, ref['stack_list']):
        yield assert_equal, len(stack.layers), len(rstack['layers'])
        lbl_s = "wl = %f, " % stack.layers[0].light.wl_nm
        for i, (lay, rlay) in enumerate(zip(stack.layers, rstack['layers'])):
            lbl_l = lbl_s + "lay %i, " % i
            yield assert_ac, lay.R12, rlay['R12'], rtol, atol, lbl_l + 'R12'
            yield assert_ac, lay.T12, rlay['T12'], rtol, atol, lbl_l + 'T12'
            yield assert_ac, lay.R21, rlay['R21'], rtol, atol, lbl_l + 'R21'
            yield assert_ac, lay.T21, rlay['T21'], rtol, atol, lbl_l + 'T21'
            yield assert_ac, lay.k_z, rlay['k_z'], rtol, atol, lbl_l + 'k_z'
            # TODO: yield assert_ac, lay.sol1, rlay['sol1'], rtol, atol, lbl_l + 'k_z'
        yield assert_ac, stack.R_net, rstack['R_net'], rtol, atol, lbl_s + 'R_net'
        yield assert_ac, stack.T_net, rstack['T_net'], rtol, atol, lbl_s + 'T_net'

# def test_final_absorptance_last():
#     results_match_reference("Absorptance.txt")