"""
Test simulation of a relatively simple structure;
a dilute silicon nanowire array, at non-normal incidence.
"""

import time
import datetime
import numpy as np
import sys
sys.path.append("../PCPV/")

import objects
import materials
import plotting
from stack import *
import testing
from numpy.testing import assert_allclose as assert_ac
from numpy.testing import assert_equal

def setup_module(module):
    ################ Simulation parameters ################


    # Number of CPUs to use im simulation
    num_cores = 1
    # # Alternatively specify the number of CPUs to leave free on machine
    # leave_cpus = 4 
    # num_cores = mp.cpu_count() - leave_cpus

    # Remove results of previous simulations
    plotting.clear_previous('.txt')
    plotting.clear_previous('.pdf')
    # plotting.clear_previous('.log')
    ################ Light parameters #####################

    # # Set up light objects
    # wavelengths = [310, 410, 531.2007, 707.495,  881.786, 987.9632]
    # light_list  = [objects.Light(wl) for wl in wavelengths]
    # Single wavelength run
    wl_super =  500.0
    wavelengths = np.array([wl_super])
    light_list  = [objects.Light(wl, theta = 20, phi = 40, max_order_PWs = 1) for wl in wavelengths]
    light = light_list[0]


    ################ Scattering matrices (for distinct layers) ##############
    """ Calculate scattering matrices for each distinct layer.
    Calculated in the order listed below, however this does not influence final 
    structure which is defined later
    """

    # period must be consistent throughout simulation!!!
    period  = 600

    cover  = objects.ThinFilm(period = period, height_nm = 'semi_inf',
        material = materials.Air, loss = False)
    sim_cover = cover.calc_modes(light)

    NW_diameter = 120
    num_BM = 20
    grating_1 = objects.NanoStruct('NW_array', period, NW_diameter, height_nm = 2330,
        inclusion_a = materials.Si_c, background = materials.Air,
        loss = True, make_mesh_now = True, force_mesh = True,
        lc_bkg = 0.07, lc2= 1.5, lc3= 2.0)
    sim_grat1 = grating_1.calc_modes(light, num_BM = num_BM)

    # will only ever use top scattering matrices for the bottom layer
    bottom = objects.ThinFilm(period = period, height_nm = 'semi_inf',
        material = materials.SiO2_a, loss = False)
    sim_bot = bottom.calc_modes(light)



    ################ Construct & solve for full solar cell structure ##############
    """ Now when defining full structure order is critical and
    solar_cell list MUST be ordered from bottom to top!
    """
    stack = Stack((sim_bot, sim_grat1, sim_cover))
    stack.calc_scat()
    module.stack_list = [stack]

    last_light_object = light_list.pop()
    param_layer = grating_1 # Specify the layer for which the parameters should be printed on figures.
    params_string = plotting.gen_params_string(param_layer, last_light_object, max_num_BMs=num_BM)
    active_layer_nu = 1
    Efficiency = plotting.t_r_a_plots(stack_list, wavelengths, params_string, 
        active_layer_nu=active_layer_nu)


    # # SAVE DATA AS REFERENCE
    # # Only run this after changing what is simulated - this
    # # generates a new set of reference answers to check against
    # # in the future
    # testing.save_reference_data("case_2", stack_list)

def results_match_reference(filename):
    rtol = 1e-6
    atol = 1e-2
    reference = np.loadtxt("ref/case_2/" + filename)
    result    = np.loadtxt(filename)
    np.testing.assert_allclose(result, reference, rtol, atol, filename)

def test_txt_results():
    result_files = (
        "Absorptance_stack0001.txt",
        "Lay_Absorb_0_stack0001.txt",
        "Lay_Trans_0_stack0001.txt",
        "Reflectance_stack0001.txt",
        "Transmittance_stack0001.txt",
        )
    for f in result_files:
        yield results_match_reference, f

def test_stack_list_matches_saved(casefile_name = 'case_2'):
    rtol = 1e-1
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
            #TODO: yield assert_ac, lay.sol1, rlay['sol1']
        yield assert_ac, stack.R_net, rstack['R_net'], rtol, atol, lbl_s + 'R_net'
        yield assert_ac, stack.T_net, rstack['T_net'], rtol, atol, lbl_s + 'T_net'
