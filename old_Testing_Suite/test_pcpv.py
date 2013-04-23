# Testing script for pcpv.exe
import os
import glob
import subprocess
import numpy
import unittest
import numpy.testing as nt
#from numpy.ma.testutils import assert_array_approx_equal as arr_aeq
arr_eq = nt.assert_array_almost_equal
from check_arrays_equal import BjornTest

test_cases = ['Air', 'Fano', 'Povinelli', 'Povinelli_Substrate']

result_file_list = ['VFEM_2D/Output/*.txt',
                   'VFEM_2D/Output/Dispersion/*.txt',
                   'VFEM_2D/Normed/*.txt',
                   'VFEM_2D/Normed/*.log',
                   'VFEM_2D/Matrices/*.txt',
                   'VFEM_2D/Matrices/*/*.txt',
                   ]

# delete previous test results
for test_folder in test_list:
    for result_location in result_file_list:
        for f in glob.glob(os.path.join(test_folder, result_location)):
            print f
#                os.remove(f)

# call pcpv for all different directories (geometries)
ps = []
#for test_folder in test_list:
os.chdir('Air/VFEM_2D/')
# p=subprocess.Popen('PWD=Air/VFEM_2D/ pcpv.exe', shell=True)
# ps.append(p)

# for p in ps:
    # p.wait()

print 'MOVING-ON-NOW-----------'


# test_list = ['Air']#, 'Fano', 'Povinelli', 'Povinelli_Substrate']
# result_file_list = ['Output/*.txt',
#                    'Output/Dispersion/*.txt',
#                   #'Normed/*.txt',
#                    'Matrices/*.txt',
#                    ]

os.chdir('../../')
# for test_folder in test_list:
#         for result_location in result_file_list:
#             ref_list = glob.glob(os.path.join(test_folder, 'Reference', result_location))
#             res_list = glob.glob(os.path.join(test_folder, 'VFEM_2D', result_location))
#             for ref_file in ref_list:
#                 a = numpy.loadtxt(ref_file)
#             for res_file in res_list:
#                 b = numpy.loadtxt(res_file)

# a = numpy.loadtxt('Air/Reference/Output/d_0600p_0000_A_Lambda.txt')
# b = numpy.loadtxt('Air/VFEM_2D/Output/d_0600p_0000_A_Lambda.txt')

arr_eq(a, b, decimal=6, err_msg='Failed in %s' % test_folder, verbose='False') 

# suite = unittest.TestLoader().loadTestsFromTestCase(BjornTest)
# unittest.TextTestRunner(verbosity=2).run(suite)
# test results by comparing difference to established results
#BjornTest.run()
