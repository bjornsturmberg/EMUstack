Installation
================

The source code for EMUstack is hosted `here on Github <https://github.com/bjornsturmberg/EMUstack>`_. Please download the latest release from here.

EMUstack has been developed on Ubuntu and is easiest to install on this platform. Simply 'sudo apt-get install' the packages listed in the dependencies.txt file and then run setup.sh. ::

    $ sudo apt-get update
    $ sudo apt-get -y install <dependencies>
    $ /setup.sh

UPDATE: the current version of SuiteSparse is not fully compatible with 64 bit Linux... a solution to this is to backport `SuiteSparse 3.4 from Ubuntu 12.04 <http://packages.ubuntu.com/source/precise/suitesparse>`_ using the method described `here <https://help.ubuntu.com/community/PinningHowto#Example_.231:_Pinning_the_ubuntu-x-swat.2BAC8-q-lts-backport-precise_PPA>`_. Alternatively the pre-compiled libraries have been shown to work on Ubuntu 14.04

On other linux distributions either use the pre-compiled libraries of install them from the package manager or manually.

All that is required to use the pre-compiled libraries is to switch to a slightly modified Makefile and then run setup.sh. ::

    $ cd backend/fortran/
    $ mv Makefile Makefile_ubuntu
    $ mv Makefile-pre_compiled_libs Makefile
    $ cd ../../
    $ /setup.sh

The Fortran components (EMUstack source code and libraries) have been successfully compiled with intel's ifortran as well as open-source gfortran. In this documentation we use gfortran.

NOTE: different versions of gmsh can give errors in the final test. This is okay, provided the test simulation ran, i.e. the test gives E rather than F.

SuiteSparse
----------------

The FEM routine used in EMUstack makes use of the highly optimised `UMFPACK <https://www.cise.ufl.edu/research/sparse/umfpack/>`_ (Unsymmetric MultiFrontal Package) direct solver for sparse matrices developed by Prof. Timothy A. Davis. This is distributed as part of the  SuiteSparse libraries under a GPL license. It can be downloaded from `https://www.cise.ufl.edu/research/sparse/SuiteSparse/ <https://www.cise.ufl.edu/research/sparse/SuiteSparse/>`_

This is the process I followed in my installations. They are provided as little more than tips...

Unpack SuiteSparse into EMUstack/backend/fortran/, it should create a directory there; SuiteSparse/
Make a directory where you want SuiteSparse installed, in my case SS_installed ::

    $ mkdir SS_installed/

edit SuiteSparse/SuiteSparse\_config/SuiteSparse\_config.mk for consistency across the whole build; i.e. if using intel fortran compiler ::

    line 75 F77 = gfortran --> ifort

set path to install folder::

    line 85 INSTALL_LIB = /$Path_to_EMustack/EMUstack/backend/fortran/SS_install/lib
    line 86 INSTALL_INCLUDE = /$Path_to_EMustack/EMUstack/backend/fortran/SS_install/include

line 290ish commenting out all other references to these::

    F77 = ifort
    CC = icc
    BLAS   = -L/apps/intel-ct/12.1.9.293/mkl/lib/intel64 -lmkl_rt
    LAPACK = -L/apps/intel-ct/12.1.9.293/mkl/lib/intel64 -lmkl_rt

Now make new directories for the paths you gave 2 steps back::

    $ mkdir SS_installed/lib SS_installed/include

Download `metis-4.0 <http://glaros.dtc.umn.edu/gkhome/fsroot/sw/metis/OLD>`_ and unpack metis into SuiteSparse/ Now move to the metis directory::

    $ cd SuiteSparse/metis-4.0

Optionally edit metis-4.0/Makefile.in as per SuiteSparse/README.txt plus with -fPIC::

    CC = gcc
    or
    CC = icc
    OPTFLAGS = -O3 -fPIC

Now make metis (still in SuiteSparse/metis-4.0/)::

    $ make

Now move back to EMUstack/backend/fortran/ ::

    $ cp SuiteSparse/metis-4.0/libmetis.a SS_install/lib/

and then move to SuiteSparse/ and execute the following::

    $ make library
    $ make install
    $ cd SuiteSparse/UMFPACK/Demo
    $ make fortran64
    $ cp SuiteSparse/UMFPACK/Demo/umf4_f77zwrapper64.o into SS_install/lib/

Copy the libraries into EMUstack/backend/fortran/Lib/ so that EMUstack/ is a complete package that can be moved across machine without alteration. This will override the pre-compiled libraries from the release (you may wish to save these somewhere).::

    $ cp SS_install/lib/*.a EMUstack/backend/fortran/Lib/
    $ cp SS_install/lib/umf4_f77zwrapper64.o EMUstack/backend/fortran/Lib/





EMUstack Makefile
-------------------

Edit EMUstack/backend/fortran/Makefile to reflect what compiler you are using and how you installed the libraries. The Makefile has further details.

Then finally run the setup.sh script!
