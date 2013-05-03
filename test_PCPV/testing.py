import numpy as np
from numpy.testing import assert_allclose as assert_ac


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
            rlay['beta'] = lay.beta
            rstack['layers'].append(rlay)
        #TODO: rstack['R_net'] = stack.R_net
        ref_stack_list.append(rstack)
    np.savez_compressed("ref/%s.npz" % casefile_name, 
        stack_list = ref_stack_list)
    raise ValueError, "Reference results saved successfully, \
    but tests will now pass trivially so let's not run them now."
