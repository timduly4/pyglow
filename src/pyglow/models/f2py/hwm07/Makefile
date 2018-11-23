fortran_files = hwm07e_modified.f90 apexcord.f90
only = hwmqt
mod = hwm07py


clean:
	rm *~ *.x *.mod *.pyf out* *.so hwm07e_modified.f90

hwm07_test:
	gfortran-4.5 -ffree-form hwm07_driver.f hwm07e.f90 dwm07b.f90 apexcord.f90 -o hwm07_driver.x

checkhwm07:
	gfortran checkhwm07.f90 hwm07e_modified.f90 dwm07b.f90 apexcord.f90 -o checkhwm.x

checkhwm07_g95:
	g95 -cpp checkhwm07.f90 hwm07e.f90 dwm07b.f90 apexcord.f90 -o checkhwm.x

patch_hwm07e:
	perl -p -e 's/\r//' < hwm07e.f90 > hwm07e_unix.f90
	patch hwm07e_unix.f90 -i hwm07e.patch -o hwm07e_modified.f90

sig:
	f2py -m $(mod) -h sig_file.pyf $(fortran_files) only: $(only) : | tee out1

mod:
	f2py -c sig_file.pyf $(fortran_files)  --fcompiler=gfortran | tee out2
remove_decode_errors:
	python remove_decode_errors.py
source:
	make remove_decode_errors;
	make patch_hwm07e;
	make sig;

compile: source
	make mod;

install:
	cp $(mod).so ../model_atmosphere/modules/

# info about converting from dos file to unix file, for perl command in patch_hwm07e: http://kb.iu.edu/data/acux.html
