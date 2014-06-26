2D FEM Mode Solver
====================

2D Mesh
--------

2D FEM mesh are created using the open source program `gmsh <http://geuz.org/gmsh/>`_.
In general they are created automatically by EMUstack using the templates files for each inclusion shape. These are stored in backend/fortran/msh. For an up to date list of templates see the 'inc_shape' entry in the NanoStruct docstring. 

An advantage of using the FEM to calculate the modes of layers is that there is absolutely no constraints on the content of the unit cell. If you wish to create a different structure this can be done using gmsh, which is also used to view the mesh files (select files with the extension .msh).

Note that the area of the unit cell must always be unity! This has been assumed throughout the theoretical derivations.



FEM Errors
-----------

There are 2 errors that can be easily triggered within the Fortran FEM routines. These both cause them to simulation to abort and the terminal to be unresponsive (until you kill python or the screen session as described in :ref:`screen_sesh`).

The first of these is ::

	Error with _naupd, info_32 =           -3
	Check the documentation in _naupd.
	Aborting...

Long story short, this indicates that the FEM mesh is too coarse for solutions for higher order Bloch modes (Eigenvaules) to converge. To see this run the simulation with FEM_debug = 1 (in mode_calcs.py) and it will print the number of converged Eigenvalues nconv != nval.
This error is easily fixed by increasing the mesh resolution. Decrease 'lc_bkg' and/or increase 'lc2' etc.


The second error is :: 

	Error with _naupd, info_32 =           -8
	Check the documentation in _naupd.
	Aborting...

This is the opposite problem, when the mesh is so fine that the simulation is overloading the memory of the machine. More accurately the memory depends on the number of Eigenvalues being calculated as well as the number of FEM mesh points.
The best solution to this is to increase 'lc_bkg' and/or decrease 'lc2' etc.