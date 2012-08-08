import numpy as np
import random
from calculate_ff     import calculate_ff

period = 600
supercell = 4
radius1 = 60
radius2 = 60
radius3 = 60
radius4 = 60

ff = 0.3
ff_tol = 0.001
min_a = 50
max_a = (period/2)/np.sqrt(supercell)
rad_array = [radius1,radius2,radius3,radius4]
for i in range(0,supercell):
    rad_array[i] = random.randint(min_a, max_a)
    print rad_array[i]
test_ff = calculate_ff(period,radius1,radius2,radius3,radius4,0,0,0,0,0,0,0,0,0,0,0,0)
while test_ff > ff+ff_tol or test_ff < ff-ff_tol:
    if test_ff > ff+ff_tol:
        rad_2_change = rad_array[random.randint(0,supercell-1)]
        rad_2_change = rad_2_change - random.randint(0,10)
    elif test_ff < ff-ff_tol:
        rad_2_change = rad_array[random.randint(0,supercell-1)]
        rad_2_change = rad_2_change + random.randint(0,10)
    test_ff = calculate_ff(period,radius1,radius2,radius3,radius4,0,0,0,0,0,0,0,0,0,0,0,0)
    print test_ff
ff = test_ff

print radius1
print radius2
print radius3
print radius4
print ff