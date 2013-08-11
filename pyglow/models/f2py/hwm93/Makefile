fortran_files = hwm93_modified.f
mod = hwm93py
only = gws5

clean:
	rm *~ *.x *.mod *.pyf out* *.so hwm93_modified.f

patch_hwm93:
	patch hwm93.f -i hwm93.patch -o hwm93_modified.f

sig:
	f2py -m $(mod) -h sig_file.pyf $(fortran_files) only: $(only) : | tee out1

patch_sig:
	patch sig_file.pyf -i sig.patch -o sig_file_patched.pyf

mod:
	f2py -c sig_file_patched.pyf $(fortran_files) --f77flags="-std=legacy" | tee out2



compile:
	make patch_hwm93;
	make sig;
	make patch_sig;
	make mod;

install:
	cp $(mod).so ../model_atmosphere/modules/




