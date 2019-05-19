
fortran_files = hwm14.f90
mod = hwm14py
only = hwm14

clean:
	rm *~ *.x *.mod *.pyf out* *.so

sig:
	f2py -m $(mod) -h sig_file.pyf $(fortran_files) only: $(only) : | tee out1

mod:
	f2py -c sig_file.pyf $(fortran_files) --f77flags="-std=legacy" | tee out2

source:
	make sig;

compile: source
	make mod;


