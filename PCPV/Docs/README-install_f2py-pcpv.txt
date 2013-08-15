##### HOW TO install the pcpv with f2py #####

1) take the day off, go to the beach, enjoy yourself.

2) get your favourite work beverage, or 30...

3) bzr checkout ...

4) Download SuiteSparse https://www.cise.ufl.edu/research/sparse/SuiteSparse/
Unpack SuiteSparse into PCPV/Fortran_pcpv/, it should create a directory there, SuiteSparse/
mkdir where you want SuiteSparse installed, in my case PCPV/Fortran_pcpv/SS_installed
mkdir  SS_installed/lib SS_installed/include
edit SuiteSparse/SuiteSparse_config/SuiteSparse_config.mk
for consitency across the whole build, use intel fortran compiler
line 75 F77 = gfortran --> ifort
set path to install folder
line 85 INSTALL_LIB = /suphys/bjorn/Usyd_Running/f2pytest/PCPV/Fortran_pcpv/SS_install/lib
line 86 INSTALL_INCLUDE = /suphys/bjorn/Usyd_Running/f2pytest/PCPV/Fortran_pcpv/SS_install/include


5) Download metis-4.0 http://glaros.dtc.umn.edu/gkhome/fsroot/sw/metis/OLD
Unpack metis into SuiteSparse/
cd into SuiteSparse/metis-4.0
optionally edit metis-4.0/Makefile.in as per SuiteSparse/ README.txt plus with -fPIC
CC = gcc
OPTFLAGS = -O3 -fPIC
make
cp f2pytest/PCPV/Fortran_pcpv/SuiteSparse/metis-4.0/libmetis.a f2pytest/PCPV/Fortran_pcpv/SS_install/lib/


6) in SuiteSparse/
make library
make install
cd SuiteSparse/UMFPACK/Demo
make fortran64
cp SuiteSparse/UMFPACK/Demo/umf4_f77zwrapper64.o into SS_install/lib/


7) Download ARPACK http://www.caam.rice.edu/software/ARPACK/download.html
Unpack into same directory as SuiteSparse/ in my case Fortran_pcpv/
edit ARmake.inc
line 28 home = $(HOME)/ARPACK --> /suphys/bjorn/Usyd_Running/f2pytest/PCPV/Fortran_pcpv/ARPACK
line 35 PLAT = nickname of platform compiling for my case SILLIAC
line 104 FC      = f77 --> ifort
line 105 FFLAGS	= -O -cg89 --> -O -fPIC
line 115 MAKE    = /bin/make --> /usr/bin/make
make lib


zarpack_util.f 
line 768-771
C     gfortran likes the following
C      INTRINSIC          ETIME
C     ifort likes the following
      EXTERNAL           ETIME



8) Edit Fortran_pcpv/Makefile update paths and names to yours
SS_LIB_LOCATION = /suphys/bjorn/Usyd_Running/f2pytest/PCPV/Fortran_pcpv/SSparse_install/lib
ARPACK_LIB_LOCATION = /suphys/bjorn/Usyd_Running/f2pytest/PCPV/Fortran_pcpv/ARPACK
UMFPACK_NAME = umf4_f77zwrapper.o
#uncomment line 62
#comment line 63 (this removes timing routines)


9) Edit PCPV/Data/gmsh_conversion/Makefile
selecting your fortran compiler on line 10.

# FC = gfortran
FC = ifort

Make in this directory.


9) make (paying close care to cross all your digits)


11) 2 set of tests test_instalation installation test_local

Open test_PCPV/test_case_* and make adjustments near the bottom as described there.
Run tests once with loose constraints, if these all pass copy your results into ref/ and reset hard constraints for your future comparisons.

10) cd all the way back to base of pcpv package in this case f2pytest/
then cd test_PCPV

python nosetests
. = pass
F = FAIL

12) That's all there is, there isn't any more.








#######################################

Re Makefiles - vayu links to arpack automatically if you load python 2.7.3
module load intel-fc/12.1.6.233 intel-cc/12.1.6.233 intel-mkl/12.1.6.233 python/2.7.3 python/2.7.3-matplotlib


"check if suitesparse-metis-dev, arpack-dev (or whatever it's called), etc are installed on your machine already. If so, you're golden and skip to step X"


When I suggested installing to a different location, I didn't mean inside the PCPV folder - I meant in ~/lib/ or on the computer's hard drive or something else. I think putting the libraries in the PCPV folder is a bad idea, but I don't know enough to know why.


I haven't bothered with grammar or rephrasing, but I will pass on the advice once given to me by a Swiss mathematician:
If you write an "if", then you must always follow it by a "then".
I've found that rule is very useful for technical writing.