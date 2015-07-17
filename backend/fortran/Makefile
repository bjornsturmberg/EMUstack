# Makefile for EMUstack. Requires some adjustment to reflect your installation.

COMPILER_VENDOR = gnu95
#COMPILER_VENDOR = intelem

SUBROUTINES_FOR_PYTHON = py_calc_modes.f py_calc_modes_1d.f conv_gmsh.f \
	gmsh_plot_field_3d.f gmsh_plot_field.f gmsh_plot_PW.f field_value_plane.f \
	gmsh_plot_slice.f

LIB_LOCATION = lib
UMFPACK_NAME = umf4_f77zwrapper64.o



### --- IF using Ubuntu SuiteSparse --------------------------------------
EMUstack.so: EMUstack.pyf lib/umf4_f77zwrapper64.o *.f
	f2py -c EMUstack.pyf *.f --fcompiler=$(COMPILER_VENDOR) \
	--link-lapack_opt --link-blas_opt  \
	$(LIB_LOCATION)/$(UMFPACK_NAME) --link-umfpack

### With installed libsuitesparse v. < 4.0.0
lib/umf4_f77zwrapper64.o: lib/umf4_f77zwrapper_SS_3.7.1.c
	gcc -O5 -DZLONG -I/usr/include/suitesparse -c lib/umf4_f77zwrapper_SS_3.7.1.c \
		-fPIC -o lib/umf4_f77zwrapper64.o
### With installed libsuitesparse v. >= 4.0.0
### 12/2014 libsuitesparse >= 4.0.0 does not support 64 bit linux umf4_f77zwrapper64
#lib/umf4_f77zwrapper64.o: lib/umf4_f77zwrapper_SS_4.2.0.c
#	gcc -O5 -DZLONG -I/usr/include/suitesparse -c lib/umf4_f77zwrapper_SS_4.2.0.c \
#       -fPIC -o lib/umf4_f77zwrapper64.o

EMUstack.pyf: $(SUBROUTINES_FOR_PYTHON)
	f2py -m EMUstack -h EMUstack.pyf $(SUBROUTINES_FOR_PYTHON) --overwrite-signature

purge:
	-rm EMUstack.pyf
	-rm EMUstack.so
	-rm lib/umf4_f77zwrapper64.o
clean:
	-rm EMUstack.so
###------------------------------------------------------------------


# ### --- ELSE IF using self built SuiteSparse ------------------------
# EMUstack.so: EMUstack.pyf *.f
# 	f2py -c EMUstack.pyf *.f --fcompiler=$(COMPILER_VENDOR) \
#         --link-lapack_opt --link-blas_opt \
#         $(LIB_LOCATION)/$(UMFPACK_NAME) \
#         -L$(LIB_LOCATION) -lumfpack -lamd -lcholmod \
#         -lcolamd -lcamd -lmetis -lccolamd -lamd -lsuitesparseconfig -lrt

# EMUstack.pyf: $(SUBROUTINES_FOR_PYTHON)
# 	f2py -m EMUstack -h EMUstack.pyf $(SUBROUTINES_FOR_PYTHON) --overwrite-signature

# purge:
# 	-rm EMUstack.pyf
# 	-rm EMUstack.so
# clean:
# 	-rm EMUstack.so
# ###------------------------------------------------------------------
