import unittest
import numpy
import numpy.testing as nt
#from numpy.ma.testutils import assert_array_approx_equal as arr_aeq
arr_eq = nt.assert_array_almost_equal
import os
import sys
import glob

class BjornTest(unittest.TestCase):
    """docstring for BjornTest"""

    def test_file_1(self):
        arr_eq(a, b, decimal=6, err_msg='Failed in %s' % test_folder) 
        #arrayalmostequal(a, b, atol = 1.e-5).all()
        #assert_arr_almost_eq(a, b, atol = 1.e-5, msg = "Help me! file=%s" % filename)
	
if __name__ == '__main__':
    unittest.main()
