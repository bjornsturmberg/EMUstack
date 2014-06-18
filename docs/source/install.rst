
The source code for EMUstack is hosted  `here on Github <https://github.com/bjornsturmberg/EMUstack>`_. The most convenient method of installing EMUstack is to download the installation script as follows (the alternative being to download direct from Github).

Open a terminal and navigate to the directory where you wish to install and run EMUstack from. Then run the following commands::

	$ wget -O EMUstack_setup.sh \ 
	http://www.physics.usyd.edu.au/emustack/setup_scripts/EMUstack_setup.sh

	$ chmod +x EMUstack_setup.sh
	
	$ bash EMUstack_setup.sh


This will download the latest release of EMUstack, install all dependencies (such as gfortran and linear algebra libraries) and compile EMUstack. Finally it will run some test simulations, after which you are all ready to go! 