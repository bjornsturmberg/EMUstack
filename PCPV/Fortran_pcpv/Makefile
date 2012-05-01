#
SHELL = /bin/sh

CODE =  $(wildcard *.f)
OBJS = $(CODE:.f=.o)

BIN = pcpv.exe

FC = ifort

# Libraries
LIBS = -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_sequential \
	-lmkl_core

CFLAGS = -c -w -O2 -I./ -traceback -check
#CFLAGS = -c -O3
#CFLAGS = -c -w -traceback -check all
LFLAGS = ./Lib/umf4_f77zwrapper64.o ./Lib/libumfpack.a ./Lib/libamd.a\
	 ./Lib/libcerbla.a ./Lib/libcholmod.a\
	 ./Lib/libcolamd.a ./Lib/libcamd.a ./Lib/libmetis.a\
	 ./Lib/libccolamd.a
# Lib/umf4_f77zwrapper.o -LLib -lumfpack -lamd #./Lib/blas_LINUX.a

$(BIN):	$(OBJS)
	$(FC) $(OBJS) -L$(LIBRARY_PATH) -I$(FPATH) $(LIBS) $(LFLAGS) -o $(BIN)
.f.o:	$(CODE)	
	$(FC) $(CFLAGS) $?

clean:
	rm *.o





