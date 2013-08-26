
#################################################
PCPV is a cross platform, open source 
#################################################

Easy installation
-------------------------------------------------
* On Ubuntu the following may be installed from the 'Ubunutu software centre'
  libsuitesparse-metis-3.1.0
  libsuitesparse-metis-dev
  In which case you can skip down to the ARPACK step.




SuiteSparse (of which we particularly need UMFPACK)
-------------------------------------------------
* Download SuiteSparse https://www.cise.ufl.edu/research/sparse/SuiteSparse/
  Unpack SuiteSparse into PCPV/Fortran_pcpv/, it should create a directory there, SuiteSparse/
  $ mkdir where you want SuiteSparse installed, in my case PCPV/Fortran_pcpv/SS_installed
  $ mkdir  SS_installed/lib SS_installed/include
  edit SuiteSparse/SuiteSparse_config/SuiteSparse_config.mk
  for consitency across the whole build, use intel fortran compiler
  line 75 F77 = gfortran --> ifort
  set path to install folder
  line 85 INSTALL_LIB = /suphys/bjorn/Usyd_Running/f2pytest/PCPV/Fortran_pcpv/SS_install/lib
  line 86 INSTALL_INCLUDE = /suphys/bjorn/Usyd_Running/f2pytest/PCPV/Fortran_pcpv/SS_install/include
  
  line 290ish commenting out all other references to these
  F77 = ifort
  CC = icc
  BLAS   = -L/apps/intel-ct/12.1.9.293/mkl/lib/intel64 -lmkl_rt
  LAPACK = -L/apps/intel-ct/12.1.9.293/mkl/lib/intel64 -lmkl_rt
  
  $ mkdir the INSTALL_LIB and INSTALL_INCLUDE dirs
  

* Download metis-4.0 http://glaros.dtc.umn.edu/gkhome/fsroot/sw/metis/OLD
  Unpack metis into SuiteSparse/
  $ cd into SuiteSparse/metis-4.0
  optionally edit metis-4.0/Makefile.in as per SuiteSparse/ README.txt plus with -fPIC
  CC = gcc
  or 
  CC = icc
  OPTFLAGS = -O3 -fPIC
  $ make
  $ cp f2pytest/PCPV/Fortran_pcpv/SuiteSparse/metis-4.0/libmetis.a f2pytest/PCPV/Fortran_pcpv/SS_install/lib/


* in SuiteSparse/
  $ make library
  $ make install
  $ cd SuiteSparse/UMFPACK/Demo
  $ make fortran64
  $ cp SuiteSparse/UMFPACK/Demo/umf4_f77zwrapper64.o into SS_install/lib/
  Copy the libraries into PCPV/Fortran_pcpv/Lib/ so that PCPV/ is a complete package that can be moved across machine without alteration.
  $ cp SS_install/lib/*.a PCPV/Fortran_pcpv/Lib/
  $ cp SS_install/lib/umf4_f77zwrapper64.o PCPV/Fortran_pcpv/Lib/




ARPACK
-------------------------------------------------

* IF you have made SuiteSparse yourself, then edit PCPV/Fortran_pcpv/zarpack_util.f 
  selecting fortran compiler you are using on lines 768-771
  C     gfortran likes the following
  C      INTRINSIC          ETIME
  C     ifort likes the following
        EXTERNAL           ETIME

  ELSE IF, when running PCPV there are errors involving zneupd or znaupd then 
  $ mv zarpack.f zarpack.f.bak
  $ mv zarpack_util.f zarpack_util.f.bak
  and build ARPACK as described below.
 

* ELSE IF you are using Ubuntu package versions of SuiteSparse, then there is some issue with zarpack_util.f
  You will need to build ARPACK yourself. Download it here http://www.caam.rice.edu/software/ARPACK/download.html (choosing the stable arpack96 version).
  Unpack into same directory as SuiteSparse/ in my case Fortran_pcpv/
  edit ARmake.inc
  line 28 home = $(HOME)/ARPACK --> /suphys/bjorn/Usyd_Running/f2pytest/PCPV/Fortran_pcpv/ARPACK
  line 104 FC      = f77 --> ifort or gfort
  line 105 FFLAGS = -O -cg89 --> -O -fPIC
  line 115 MAKE    = /bin/make --> /usr/bin/make
  $ make lib




Make PCPV!!!
-------------------------------------------------

* Edit Fortran_pcpv/Makefile to reflect what compiler you are using and how you installed the libraries.
  The Makefile has further details.



#################################################
The installation of PCPV is complete.
However you almost certainly wish to have the freedom to create FEM mesh for your chosen parameters. To do this you will need to install the meshing program Gmsh and compile a simple Fortran routine.
#################################################



Install Gmsh
-------------------------------------------------
* IF using Ubuntu, then install gmsh from the 'Ubunutu software centre'. Done.

* ELSE, download a stable version of gmsh from http://geuz.org/gmsh/ 
  To date, PCPV has been successfully used with gmsh versions 2.5.1, 2.6.1

* Install 



Compile gmsh_conversion Fortran routine
-------------------------------------------------
* edit PCPV/Data/gmsh_conversion/Makefile selecting your Fortran compiler on line 10.
  
  # FC = gfortran
  FC = ifort

  then,

  $ cd PCPV/Data/gmsh_conversion/
  $ make




#################################################
The installation of PCPV is now fully complete. 
#################################################



Test the installation
-------------------------------------------------
* PCPV ships with 2 sets of tests; one to test an installation on a new machine, and another to test against when making modifications on the same machine. These tests are located in test_installation_PCPV/ and test_local_PCPV/ respectively. The actual test calculations performed in these sets are identical, and on downloading the reference data tested against is also identical.

This structure has been designed in response to installations on different machines giving slightly different results, particularly when different versions of gmsh have been used. To test a new installation please first run the tests in test_installation_PCPV/. To do this

$ cd test_installation_PCPV/
$ nosetests

During testing, individual test results are displayed with
. = pass
F = FAIL.

Once you are satisfied with you installation and the test_installation_PCPV/ tests have all passed, update the test reference data in test_local_PCPV/. To do this uncomment the line beginning with testing.save_reference_data in test_local_PCPV/test_case*.py and then

$ cd test_local_PCPV/
$ nosetests

and then re-comment those lines. Note that the test itself will fail giving the following message;
'AssertionError: Reference results saved successfully, but tests will now pass trivially so let's not run them now.'

When making updates and modifications on PCPV you should check that your results are still consistent with the test cases to the accuracy as set in test_local_PCPV/. To do this run

$ cd test_local_PCPV/
$ nosetests

with the testing.save_reference_data line commented out.





#################################################
Enjoy PCPV!!!
#################################################