COMPILER_VENDOR = g95


pcpv.so: pcpv.pyf
	f2py -c pcpv.pyf *.f Lib/umf4_f77zwrapper64.o --link-lapack_opt \
		--link-umfpack --link-blas_opt

pcpv.pyf: main_64.f
	f2py -m pcpv -h pcpv.pyf main_64.f --overwrite-signature

purge:
	rm pcpv.pyf
	rm pcpv.so

clean:
	rm pcpv.so
