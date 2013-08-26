Template python script file to execute a simulation. To start, open a terminal and change
directory to the directory containing this file (which must be in the same directory as 
the PCPV directory). Run this script file by executing the following in the command line

$ python simo_template-grating.py

This will use num_cores worth of your CPUs, and by default return you in the command
line, having printed results and saved plots to file as specified towards the end of 
this file. If instead you wish to have direct access to the simulation results (for 
further manipulation, debugging etc.) run this script with

$ python -i simo_template-grating.py

which, after the calculations are complete, will return you into an interactive session 
of python, in which all simulation objects are accessible. In this session you can access
the docstrings of objects/classes/methods by typing

>>> from pydoc import help
>>> help(objects.Light)

where we have accessed the docstring of the Light class from objects.py


In real simulation scripts replace this docstring with a brief description of the 
simulation, eg.
`Simulating the coupling of normally incident light into evanescent orders through a 
metallic grating of period 120 nm. Included 3 PW orders.'


-------------


-------------



-------------
|           |
|           |
|           |
|           |
|           |
-------------



      _
    (   )
  (       )
 (         )
  (       )
    ( _ )
      
'gmsh '+ data_location + msh_name + '.msh'
'gmsh '+ data_location + msh_name + '.geo'




##### Notes for running PCPV on NCI supercomputers #####

Make sure 

intel-fc
intel-cc
intel0mkl
python
python...-matplotlib

are all loaded. This is best done by adding

module load intel-fc/12.1.9.293
module load intel-cc/12.1.9.293
module load intel-mkl/12.1.9.293
module load python/2.7.3
module load python/2.7.3-matplotlib

to the bottom of your .login, .profile, and .cshrc/.bashrc files (or different versions numbers).