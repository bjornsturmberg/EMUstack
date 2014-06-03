"""
Reload all simulation results from specified simulation
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
