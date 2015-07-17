Welcome to EMUstack!
--------------------

EMUstack is an open-source simulation package for calculating light propagation through multi-layered stacks of dispersive, lossy, nanostructured, optical media. It implements a generalised scattering matrix method, which extends the physical intuition of thin film optics to complex structures.


ORIGIN
------

EMUstack is the product of many years research, and is designed primarily as a research tool.
The underpinnings of EMUstack were majoritively developed within CUDOS (the ARC Centre of Excellence for Ultra-high bandwidth Devices for Optical Systems), at its University of Sydney (USyd) and University of Technology Sydney (UTS) nodes.

The scattering matrix method formalism was developed by Lindsay Botten, Ara Asatryan, and Kokou Dossou at UTS. The FEM routine was developed by Kokou Dossou while at the Université du Québec en Outaouais, Université Laval and UTS.

EMUstack was written by Björn Sturmberg during his Ph.D. at USyd, which was supported by the Australian Renewable Energy Agency. THE FEM routine was written by Kokou Dossou, and Felix Lawrence created the smooth f2py interface of Fortran and Python while at UTS.


LICENSE
-------

EMUstack is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License in the EMUstack/ directory. If not, see `here <http://www.gnu.org/copyleft/gpl.html>`_


OFFICIALLY SUPPORTED PLATFORMS
------------------------------

EMUstack has been developed for use on Linux and Unix-like operating systems. It may be easily ported to other operating systems, but there are no current plans for doing so. If you are willing and able to do so, please get in contact!


SETUP
-------

EMUstack has been developed on Ubuntu and is easiest to install on this platform. Simply sudo apt-get install the packages listed in the dependencies.txt file and then run setup.sh.
On other linux distributions you may have to install some packages manually, hints for which are given in the documentation.


TESTING
-------

EMUstack comes with a range of tests to ensure it is running correctly. These are found in the EMUstack/tests directories.


FILE MAP
--------

**backend** contains all the modules that make up EMUstack.

____ **data** contains refractive index data for materials and solar irradiance spectra.

____ **fortran** the Fortran code for the FEM solvers and the scattering matrix calculations.

________ **lib** libraries used by the Fortran routines.

________ **msh** location of Gmsh FEM mesh files, including template files.

**docs** contains the documentation.

____ **build** the documentation you want to read.

________ **doctrees** created by sphynix in compilation of docs.

________ **html** in web format.

________ **latex** in pdf format.

____ **source** the documentation you want to edit.

**examples** contains example simulations. Numbered files form a tutorial.

**reload_simo** contains python script demonstrating how to load a .npz file.

**tests** test files, best run with nosetests.

____ **ref** reference data for tests.

*dependencies.txt* list of linux packages needed to run EMUstack.

*LICENSE.txt* GNU GPL software license.

*README.rst* this file...

*setup.sh* bash script that compiles Fortran code and runs tests.


HELP!
-----

The growing human-readable documentation for EMUstack lives in EMUstack/docs/ and online at `readthedocs <http://emustack.readthedocs.org/en/latest/index.html>`_.org.

Limited community support is `available on the GitHub site <https://github.com/bjornsturmberg/EMUstack>`_, and via `the EMUstack mailing list <https://groups.google.com/forum/#!forum/emustack>`_. Feel free to ask questions.


CONTRIBUTING
------------

If you make an improvement to EMUstack, please share it with others by contributing back to the project via the GitHub site.
