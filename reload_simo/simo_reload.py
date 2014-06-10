"""
Reload simulation results for further manipulation.
"""

import numpy as np
import sys
sys.path.append("../backend/")

import objects
import materials
import plotting
from stack import *

directory = '/home/bjorn/Results/14_06-Ev_FoMs/fig_of_merit-900-1100/'
npz_file  = 'Simo_results20140610060241.npz'
data = np.load(directory+npz_file)
stacks_list = data['stacks_list']

""" stacks_list is now just like in the original simo.py and we can access/plot as before. """

### Plot transmission, reflection and absorption spectra
# plotting.t_r_a_plots(stacks_list) 

# plotting.evanescent_merit(stacks_list, lay_interest=0)
# plotting.evanescent_merit(stacks_list, lay_interest=1)
