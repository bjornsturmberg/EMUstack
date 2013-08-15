"""
Test simulation of a very simple structure;
a stack of lossy homogeneous dielectric films.
NOTE: This calculation is entirely analytical & should produce excellent 
agreement on all machines.
"""

import time
import datetime
import numpy as np
import sys
from multiprocessing import Pool
sys.path.append("../PCPV/")

import objects
import materials
import plotting
from stack import *
import testing
from numpy.testing import assert_allclose as assert_ac
from numpy.testing import assert_equal

start = time.time()
################ Simulation parameters ################

# Number of CPUs to use im simulation
num_cores = 3
# # Alternatively specify the number of CPUs to leave free on machine
# leave_cpus = 4 
# num_cores = mp.cpu_count() - leave_cpus

# Remove results of previous simulations
plotting.clear_previous('.txt')
plotting.clear_previous('.pdf')
# plotting.clear_previous('.log')

################ Light parameters #####################

# Set up light objects
wl_1     = 400
wl_2     = 1000
no_wl_1  = 3
wavelengths = np.linspace(wl_1, wl_2, no_wl_1)
light_list  = [objects.Light(wl, max_order_PWs = 1) for wl in wavelengths]
# Single wavelength run
# wl_super =  500.0
# wavelengths = np.array([wl_super])
# light_list  = [objects.Light(wl, max_order_PWs = 1) for wl in wavelengths]
# light = light_list[0]


################ Scattering matrices (for distinct layers) ##############
""" Calculate scattering matrices for each distinct layer.
Calculated in the order listed below, however this does not influence final 
structure which is defined later
"""

# period must be consistent throughout simulation!!!
period = 1

cover  = objects.ThinFilm(period = period, height_nm = 'semi_inf',
    material = materials.Material(3.5 + 0.0j), loss = False)

homo_film1  = objects.ThinFilm(period = period, height_nm = 50,
    material = materials.Material(3.6 + 0.27j), loss = True)

homo_film2  = objects.ThinFilm(period = period, height_nm = 200,
    material = materials.Si_c, loss = True)

mirror = objects.ThinFilm(period = period, height_nm = 100,
    material = materials.Ag, loss = True)

bottom = objects.ThinFilm(period = period, height_nm = 'semi_inf',
    material = materials.Air, loss = False)


stack_list = []

def simulate_stack(light):
    
    ################ Evaluate each layer individually ##############
    sim_cover = cover.calc_modes(light)
    sim_homo_film1 = homo_film1.calc_modes(light)
    sim_homo_film2 = homo_film2.calc_modes(light)
    sim_mirror = mirror.calc_modes(light)
    sim_bot = bottom.calc_modes(light)

    ################ Evaluate full solar cell structure ##############
    """ Now when defining full structure order is critical and
    solar_cell list MUST be ordered from bottom to top!
    """
    stack = Stack((sim_bot,sim_mirror,sim_homo_film1,sim_homo_film2,sim_homo_film1,sim_homo_film2, sim_cover))
    stack.calc_scat(pol = 'TM')

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
    param_layer = homo_film2 # Specify the layer for which the parameters should be printed on figures.
    params_string = plotting.gen_params_string(param_layer, last_light_object)

    active_layer_nu = 3 # Specify which layer is the active one (where absorption generates charge carriers).
    # Plot total transmission, reflection, absorption & that of each layer. 
    # Also calculate efficiency of active layer.
    Efficiency = plotting.t_r_a_plots(stack_list, wavelengths, params_string, 
        active_layer_nu=active_layer_nu) 



    # # SAVE DATA AS REFERENCE
    # # Only run this after changing what is simulated - this
    # # generates a new set of reference answers to check against
    # # in the future
    # testing.save_reference_data("case_0", stack_list)


def results_match_reference(filename):
    rtol = 1e-6
    atol = 1e-6
    reference = np.loadtxt("ref/case_0/" + filename)
    result    = np.loadtxt(filename)
    np.testing.assert_allclose(result, reference, rtol, atol, filename)

def test_txt_results():
    result_files = (
        "Absorptance_stack1.txt",
        "Lay_Absorb_0_stack1.txt",
        "Lay_Absorb_1_stack1.txt",
        "Lay_Absorb_2_stack1.txt",
        "Lay_Absorb_3_stack1.txt",
        "Lay_Absorb_4_stack1.txt",
        "Lay_Trans_0_stack1.txt",
        "Lay_Trans_1_stack1.txt",
        "Lay_Trans_2_stack1.txt",
        "Lay_Trans_3_stack1.txt",
        "Lay_Trans_4_stack1.txt",
        "Reflectance_stack1.txt",
        "Transmittance_stack1.txt",
        )
    for f in result_files:
        yield results_match_reference, f

def test_stack_list_matches_saved(casefile_name = 'case_0'):
    rtol = 1e-6
    atol = 1e-6
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
            #TODO: yield assert_ac, lay.sol1, rlay['sol1']
        yield assert_ac, stack.R_net, rstack['R_net'], rtol, atol, lbl_s + 'R_net'
        yield assert_ac, stack.T_net, rstack['T_net'], rtol, atol, lbl_s + 'T_net'