fortran_files = igrf12_modified.f
mod = igrf12py
only = igrf12syn

clean:
	rm *~ *.x *.mod *.pyf out* *.so igrf12_modified.f

igrf:
	gfortran igrf12.f -o igrf12.x

igrf_test:
	gfortran igrf12_driver.f igrf12_cut.f -o igrf_test.x



patch_igrf12:
	patch igrf12.f -i igrf12.patch -o igrf12_modified.f

sig:
	f2py -m $(mod) -h sig_file.pyf $(fortran_files) only: $(only) : | tee out1

patch_sig:
	patch sig_file.pyf -i sig.patch -o sig_file_patched.pyf

mod:
	f2py -c sig_file_patched.pyf $(fortran_files) | tee out2

source:
	make patch_igrf12;	
	make sig;
	make patch_sig;

compile: source
	make mod;

install:
	cp $(mod).so ../model_atmosphere/modules
