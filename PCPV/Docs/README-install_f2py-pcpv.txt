##### HOW TO install the pcpv with f2py #####

1) take the day off, go to the beach, enjoy yourself.

2) get your favorite work bevarage, or 30...

3) bzr checkout ...

4) Download SuiteSparse SuiteSparse https://www.cise.ufl.edu/research/sparse/SuiteSparse/
Unpack SuiteSparse into PCPV/Fortran_pcpv/, it should create a directory there SuiteSparse/
mkdir where you want installed SuiteSparse in my case PCPV/Fortran_pcpv/SS_installed
mkdirs  SS_installed/lib SS_installed/include
edit SuiteSparse/SuiteSparse_config/SuiteSparse_config.mk
for consistence across whole build use intel fortran compiler
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


8) Edit Fortran_pcpv/Makefile update paths and names to yours
SS_LIB_LOCATION = /suphys/bjorn/Usyd_Running/f2pytest/PCPV/Fortran_pcpv/SSparse_install/lib
ARPACK_LIB_LOCATION = /suphys/bjorn/Usyd_Running/f2pytest/PCPV/Fortran_pcpv/ARPACK
ARPACK_NAME = SILLIAC
UMFPACK_NAME = umf4_f77zwrapper.o
#uncomment line 62
#comment line 63 (this removes timing routines)


9) make (paying close care to cross all your digits)

10) cd all the way back to base of pcpv package in this case f2pytest/
the cd test_PCPV
(find more digits to cross, steal from others by amputation if necissary)
python nosetests
. = pass
F = FAIL

11) There may be many fails, these are probably related to using different versions of gmsh
If tests all complete but with fails you will need to update the test references for your compilation.

Open test_PCPV/test_case_* and make adjustments near the bottom as described there.
Rune tests once with loose constraints, if these all pass copy your results into ref/ and reset hard constraints for your future comparisions.

12) That's all there is, there isn't any more.


