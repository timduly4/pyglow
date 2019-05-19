fortran_files = igrf11_modified.f
mod = igrf11py
only = igrf11syn

clean:
	rm *~ *.x *.mod *.pyf out* *.so igrf11_modified.f

igrf:
	gfortran igrf11.f -o igrf11.x

igrf_test:
	gfortran igrf11_driver.f igrf11_cut.f -o igrf_test.x



patch_igrf11:
	patch igrf11.f -i igrf11.patch -o igrf11_modified.f

sig:
	f2py -m $(mod) -h sig_file.pyf $(fortran_files) only: $(only) : | tee out1

patch_sig:
	patch sig_file.pyf -i sig.patch -o sig_file_patched.pyf

mod:
	f2py -c sig_file_patched.pyf $(fortran_files) | tee out2

source:
	make patch_igrf11;	
	make sig;
	make patch_sig;

compile: source
	make mod;

install:
	cp $(mod).so ../model_atmosphere/modules
