Installation
================

The source code for EMUstack is hosted `here on Github <https://github.com/bjornsturmberg/EMUstack>`_. Please download the latest release from here.

EMUstack has been developed on Ubuntu and is easiest to install on this platform. Simply sudo apt-get install the packages listed in the dependencies.txt file and then run setup.sh.
On other linux distributions either install the packages from the distributions package manager or follow the process below to manually download and install the packages that may not be included in your platform.

The Fortran components (EMUstack source code and libraries) have been successfully compiled with intel's ifortran as well as open-source gfortran. In this documentation we use gfortran.

SuiteSparse
----------------

The FEM routine used in EMUstack makes use of the highly optimised `UMFPACK <https://www.cise.ufl.edu/research/sparse/umfpack/>`_ (Unsymmetric MultiFrontal Package) direct solver for sparse matrices developed by Prof. Timothy A. Davis. This is distributed as part of the  SuiteSparse libraries under a GPL license. It can be downloaded from `https://www.cise.ufl.edu/research/sparse/SuiteSparse/ <https://www.cise.ufl.edu/research/sparse/SuiteSparse/>`_ 

Unpack SuiteSparse into EMUstack/backend/fortran/, it should create a directory there, SuiteSparse/
mkdir where you want SuiteSparse installed, in my case SS_installed
\begin{lstlisting}
$ mkdir  SS_installed/lib SS_installed/include
\end{lstlisting}

edit SuiteSparse/SuiteSparse\_config/SuiteSparse\_config.mk for consistency across the whole build, if using intel fortran compiler

line 75 F77 = gfortran --> ifort

set path to install folder

line 85 INSTALL_LIB = /$Path_to_EMustack/EMUstack/backend/fortran/SS_install/lib
line 86 INSTALL_INCLUDE = /$Path_to_EMustack/EMUstack/backend/fortran/SS_install/include


line 290ish commenting out all other references to these
\begin{lstlisting}
F77 = ifort
CC = icc
BLAS   = -L/apps/intel-ct/12.1.9.293/mkl/lib/intel64 -lmkl_rt
LAPACK = -L/apps/intel-ct/12.1.9.293/mkl/lib/intel64 -lmkl_rt
\end{lstlisting}

\begin{lstlisting}
$ mkdir the INSTALL_LIB and INSTALL_INCLUDE dirs
\end{lstlisting}


\item Download \href{http://glaros.dtc.umn.edu/gkhome/fsroot/sw/metis/OLD}{metis-4.0}
Unpack metis into SuiteSparse/
\begin{lstlisting}
$ cd SuiteSparse/metis-4.0
\end{lstlisting}
optionally edit metis-4.0/Makefile.in as per SuiteSparse/ README.txt plus with -fPIC
\begin{lstlisting}
CC = gcc
or 
CC = icc
OPTFLAGS = -O3 -fPIC
\end{lstlisting}
\begin{lstlisting}
$ make
$ cp /$Path_to_EMustack/EMUstack/backend/fortran/SuiteSparse/metis-4.0/libmetis.a /$Path_to_EMustack/EMUstack/backend/fortran/SS_install/lib/
\end{lstlisting}


\item in SuiteSparse/
\begin{lstlisting}
$ make library
$ make install
$ cd SuiteSparse/UMFPACK/Demo
$ make fortran64
$ cp SuiteSparse/UMFPACK/Demo/umf4_f77zwrapper64.o into SS_install/lib/
\end{lstlisting}
Copy the libraries into EMUstack/backend/fortran/Lib/ so that EMUstack/ is a complete package that can be moved across machine without alteration.

\begin{lstlisting}
$ cp SS_install/lib/*.a EMUstack/backend/fortran/Lib/
$ cp SS_install/lib/umf4_f77zwrapper64.o EMUstack/backend/fortran/Lib/
\end{lstlisting}
\end{enumerate}





ARPACK
-----------


IF you have made SuiteSparse yourself, then edit EMUstack/backend/fortran/zarpack\_util.f 
specifying whether or not your fortran compiler has ETIME is a built in function of your compiler (INTRINSIC). On lines 768-771,
\begin{lstlisting}
C     gfortran likes the following
C      INTRINSIC          ETIME
C     ifort likes the following
      EXTERNAL           ETIME
\end{lstlisting}

ELSE IF, when running EMUstack there are errors involving zneupd or znaupd then
\begin{lstlisting}
$ mv zarpack.f zarpack.f.bak
$ mv zarpack_util.f zarpack_util.f.bak
\end{lstlisting}
and build ARPACK as described below.
 

\item ELSE IF you are using Ubuntu package versions of SuiteSparse, then there is some issue with zarpack\_util.f
You will need to build ARPACK yourself. Download it here \href{http://www.caam.rice.edu/software/ARPACK/download.html}{http://www.caam.rice.edu/software/ARPACK/download.html} (choosing the stable arpack96 version).\
\textbf{ARPACK plus ARPACK patch}
Unpack into same directory as SuiteSparse/ in my case EMUstack/backend/fortran/
edit ARmake.inc

line 28 

home = $(HOME)/ARPACK --> EMUstack/backend/fortran/ARPACK
line 104 FC      = f77 --> ifort or gfort
line 105 FFLAGS = -O -cg89 --> -O -fPIC
line 115 MAKE    = /bin/make --> /usr/bin/make

$ make lib




EMUstack Makefile
-------------------

Edit EMUstack/backend/fortran/Makefile to reflect what compiler you are using and how you installed the libraries. The Makefile has further details.

Then finally run the setup.sh script!