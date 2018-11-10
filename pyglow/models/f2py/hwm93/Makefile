fortran_files = hwm93_modified.f
mod = hwm93py
only = gws5

clean:
	rm -f *~ *.x *.mod *.pyf out* *.so hwm93_modified.f

hwm93_modified.f:
	patch hwm93.f -i hwm93.patch -o hwm93_modified.f

sig_file.pyf: $(fortran_files)
	f2py -m $(mod) -h sig_file.pyf $^ only: $(only) : | tee out1

sig_file_patched.pyf: sig_file.pyf
	patch sig_file.pyf -i sig.patch -o sig_file_patched.pyf

mod: sig_file_patched.pyf $(fortran_files) 
	f2py -c $^ --f77flags="-std=legacy" | tee out2



compile: mod

source: sig_file_patched.pyf $(fortran_files)

install:
	cp $(mod).so ../model_atmosphere/modules/




