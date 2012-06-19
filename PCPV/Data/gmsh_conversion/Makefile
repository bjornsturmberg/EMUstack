#
SHELL = /bin/sh

#CODE =  $(wildcard *.f)
CODE =  $(wildcard conv_gmsh_py.f)
OBJS = $(CODE:.f=.o)

BIN = conv_gmsh_py.exe

FC = ifort

# Libraries
LIBS = 

CFLAGS = -c -w -O2 -I./ -traceback -check
LFLAGS = 

$(BIN):	$(OBJS)
	$(FC) $(OBJS) -L$(LIBRARY_PATH) -I$(FPATH) -o $(BIN)
.f.o:	$(CODE)	
	$(FC) $(CFLAGS) $?

clean:
	rm *.o





