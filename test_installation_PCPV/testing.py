import numpy as np
from numpy.testing import assert_allclose as assert_ac

from objects import Simmo
from os import system as ossys

def save_reference_data(casefile_name, stack_list):
    ref_stack_list = []
    for stack in stack_list:
        rstack = {'layers' : []}
        for lay in stack.layers:
            rlay = {}
            rlay['R12']  = lay.R12
            rlay['T12']  = lay.T12
            rlay['R21']  = lay.R21
            rlay['T21']  = lay.T21
            rlay['k_z']  = lay.k_z
            if isinstance(rlay, Simmo):
                rlay['sol1'] = lay.sol1
            rstack['layers'].append(rlay)
        rstack['R_net'] = stack.R_net
        rstack['T_net'] = stack.T_net
        ref_stack_list.append(rstack)
    np.savez_compressed("ref/%s.npz" % casefile_name, 
        stack_list = ref_stack_list)

    cp_cmd = 'cp *.txt ref/%s/' % casefile_name
    ossys(cp_cmd)

    assert False, "Reference results saved successfully, \
but tests will now pass trivially so let's not run them now."
