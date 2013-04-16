#
SHELL = /bin/sh

CODE =  $(wildcard *.f)
OBJS = $(CODE:.f=.o)

BIN = pcpv.exe

FC = gfortran

# Libraries
# LIBS = -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_sequential \
# 	-lmkl_core
LIBRARY_PATH = /usr/lib/
LIBS = -lumfpack -lcholmod -lamd -lcolamd -lcamd \
	-llapack -lblas #-lcerbla -lmetis
# LIBS = -lblas

#CFLAGS = -c -w -O2 -I./ -traceback -check
CFLAGS = -c -O3
EXTRA = ./Lib/umf4_f77zwrapper64.o
#CFLAGS = -c -w -traceback -check all
# LFLAGS = ./Lib/umf4_f77zwrapper64.o ./Lib/libumfpack.a ./Lib/libamd.a\
#	 ./Lib/libcerbla.a ./Lib/libcholmod.a\
#	 ./Lib/libcolamd.a ./Lib/libcamd.a ./Lib/libmetis.a\
#	 ./Lib/libccolamd.a ./Lib/blas.lib ./Lib/libumfpack.a
# Lib/umf4_f77zwrapper.o -LLib -lumfpack -lamd #./Lib/blas_LINUX.a

$(BIN):	$(OBJS)
	$(FC) $(OBJS) $(EXTRA) -L$(LIBRARY_PATH) $(LIBS) $(LFLAGS) -o $(BIN)
	 # $(FC) $(OBJS) $(LFLAGS) -o $(BIN)
.f.o:	$(CODE)	
	$(FC) $(CFLAGS) $?

clean:
	rm *.o





