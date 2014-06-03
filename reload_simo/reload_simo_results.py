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


data = np.load('../fig_of_merit/Simo_results20140603173129.npz')
stacks_list = data['stacks_list']

""" stacks_list is now just like in the original simo.py and we can access/plot as before. """

### Plot transmission, reflection and absorption spectra
plotting.t_r_a_plots(stacks_list) 
