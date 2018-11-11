fortran_files = nrlmsise00_sub_patched.for
mod = msis00py
only = gtd7

clean:
	rm -f *~ *.x *.mod *.pyf out* *.so

msis_test:
	gfortran -std=legacy -fno-automatic -O2 msis00_driver.f nrlmsise00.f -o msis_test.x

sig:
	f2py -m $(mod) -h sig_file.pyf $(fortran_files) only: $(only) : | tee out1

patch_sig:
	patch sig_file.pyf -i sig.patch -o sig_file_patched.pyf

patch_msis:
	patch nrlmsise00_sub.for -i nrlmsise00_sub.patch -o nrlmsise00_sub_patched.for

mod:
	f2py -c sig_file_patched.pyf $(fortran_files) --f77flags="-std=legacy" | tee out2


source:
	make patch_msis;
	make sig;
	make patch_sig;

compile: source
	make mod;

install:
	cp $(mod).so ../model_atmosphere/modules/
