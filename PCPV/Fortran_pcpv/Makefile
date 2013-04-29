COMPILER_VENDOR = g95

SUBROUTINES_FOR_PYTHON = main_64.f calc_modes.f calc_scat.f


pcpv.so: pcpv.pyf Lib/umf4_f77zwrapper64.o
	f2py -c pcpv.pyf *.f Lib/umf4_f77zwrapper64.o --link-lapack_opt \
		--link-umfpack --link-blas_opt

pcpv.pyf: $(SUBROUTINES_FOR_PYTHON)
	f2py -m pcpv -h pcpv.pyf $(SUBROUTINES_FOR_PYTHON) --overwrite-signature

# With installed libsuitesparse v. < 4.0.0
Lib/umf4_f77zwrapper64.o: Lib/umf4_f77zwrapper_SS_3.7.1.c
	gcc -DZLONG -I/usr/include/suitesparse -c Lib/umf4_f77zwrapper_SS_3.7.1.c \
		-fPIC -o Lib/umf4_f77zwrapper64.o

# With installed libsuitesparse v. >= 4.0.0
# Lib/umf4_f77zwrapper64.o: Lib/umf4_f77zwrapper_SS_4.2.0.c
# 	gcc -DZLONG -I/usr/include/suitesparse -c Lib/umf4_f77zwrapper_SS_4.2.0.c \
# 		-fPIC -o umf4_f77zwrapper64.o

purge:
	-rm pcpv.pyf
	-rm pcpv.so
	-rm Lib/umf4_f77zwrapper64.o

clean:
	-rm pcpv.so
