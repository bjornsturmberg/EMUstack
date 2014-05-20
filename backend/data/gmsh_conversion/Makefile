#
SHELL = /bin/sh

#CODE =  $(wildcard *.f)
CODE =  $(wildcard conv_gmsh_py.f)
OBJS = $(CODE:.f=.o)

BIN = conv_gmsh_py.exe

FC = gfortran
# FC = ifort

# Libraries
LIBS = 

#CFLAGS = -c -w -O2 -I./ -traceback -check
CFLAGS = -c -w -O2 -I./ 
LFLAGS = 

$(BIN):	$(OBJS)
	$(FC) $(OBJS) -L$(LIBRARY_PATH) -I$(FPATH) -o $(BIN)
.f.o:	$(CODE)	
	$(FC) $(CFLAGS) $?

clean:
	rm *.o





