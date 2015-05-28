"""
    test_case_2_simple_non_norm_inc.py is a simulation example for EMUstack.

    Copyright (C) 2015  Bjorn Sturmberg, Kokou Dossou, Felix Lawrence

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
Test simulation of a relatively simple structure;
a dilute silicon nanowire array, at non-normal incidence.
Uses .mail file from repository (to avoid meshing discrepancies).
"""

import time
import datetime
import numpy as np
import sys
sys.path.append("../backend/")

import objects
import materials
import plotting
from stack import *
import testing
from numpy.testing import assert_allclose as assert_ac
from numpy.testing import assert_equal

def setup_module(module):
    ################ Light parameters #####################
    # Set up light objects
    wl_super =  500.0
    wavelengths = np.array([wl_super])
    light_list  = [objects.Light(wl, max_order_PWs = 1, theta = 20, phi = 40) for wl in wavelengths]
    light = light_list[0]


    ################ Scattering matrices (for distinct layers) ##############
    """ Calculate scattering matrices for each distinct layer.
    Calculated in the order listed below, however this does not influence final
    structure which is defined later
    """

    # period must be consistent throughout simulation!!!
    period  = 600

    NW_diameter = 120
    num_BMs = 20
    NW_array = objects.NanoStruct('2D_array', period, NW_diameter, height_nm = 2330,
        inclusion_a = materials.Si_c, background = materials.Air,
        loss = True, make_mesh_now = False, mesh_file='4testing-600_120.mail')
    sim_NW_array = NW_array.calc_modes(light, num_BMs = num_BMs)

    superstrate  = objects.ThinFilm(period = period, height_nm = 'semi_inf',
        material = materials.Air, loss = False)
    sim_superstrate = superstrate.calc_modes(light)

    substrate = objects.ThinFilm(period = period, height_nm = 'semi_inf',
        material = materials.SiO2, loss = False)
    sim_substrate = substrate.calc_modes(light)



    ################ Construct & solve for full solar cell structure ##############
    """ Now when defining full structure order is critical and
    solar_cell list MUST be ordered from bottom to top!
    """
    stack = Stack((sim_substrate, sim_NW_array, sim_superstrate))
    stack.calc_scat()
    module.stack_list = [stack]

    plotting.t_r_a_write_files(stack_list, wavelengths)


    # # SAVE DATA AS REFERENCE
    # # Only run this after changing what is simulated - this
    # # generates a new set of reference answers to check against
    # # in the future
    # testing.save_reference_data("case_2", stack_list)

def results_match_reference(filename):
    rtol = 1e-6
    atol = 1e-6
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
    rtol = 1e-4
    atol = 1e-4
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
