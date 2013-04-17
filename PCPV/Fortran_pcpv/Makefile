#
SHELL = /bin/sh

CODE =  $(wildcard *.f)
OBJS = $(CODE:.f=.o)

BIN = pcpv.exe

FC = gfortran

# Libraries
# LIBS = -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_sequential \
# 	-lmkl_core
#LIBRARY_PATH = /usr/lib
# LIBS = -l:libamd.so.2.2.0 -l:libumfpack.so.5.4.0 -l:libcholmod.so.1.7.1 \
# 	-l:libcolamd.so.2.7.1 -l:liblapack.so.3gf -l:libblas.so.3gf
# LIBS = -lamd -lumfpack -lcholmod -lcolamd -lcamd \
# 	-llapack -lblas #-lmetis #-lcerbla
# LIBS = -lblas

#CFLAGS = -c -w -O2 -I./ -traceback -check
CFLAGS = -c -O3
EXTRA = ./Lib/umf4_f77zwrapper64.o
#CFLAGS = -c -w -traceback -check all
LIBRARY_PATH = ./Lib
LIBS = -lumfpack -lcholmod -lcolamd -lcamd -lmetis -lcerbla \
	-llapack -lblas -lamd -lccolamd -lsuitesparseconfig

# LFLAGS = ./Lib/umf4_f77zwrapper64.o ./Lib/libumfpack.a ./Lib/libamd.a\
# 	 ./Lib/libcerbla.a ./Lib/libcholmod.a\
# 	 ./Lib/libcolamd.a ./Lib/libcamd.a ./Lib/libmetis.a\
# 	 ./Lib/libccolamd.a ./Lib/libumfpack.a #./Lib/blas.lib
# Lib/umf4_f77zwrapper.o -LLib -lumfpack -lamd #./Lib/blas_LINUX.a

$(BIN):	$(OBJS)
	$(FC) $(OBJS) $(EXTRA) -L$(LIBRARY_PATH) $(LIBS) $(LFLAGS) -o $(BIN)
	 # $(FC) $(OBJS) $(LFLAGS) -o $(BIN)
.f.o:	$(CODE)	
	$(FC) $(CFLAGS) $?

clean:
	rm *.o





