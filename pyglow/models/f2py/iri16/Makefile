fortran_files = cira.for igrf.for iridreg.for irifun.for irisub.for iritec.for iriflip_modified.for cosd_sind.for
only = iri_sub read_ig_rz readapf107
mod = iri16py

f77flags = -std=legacy -w -O2 -fbacktrace -fno-automatic -fPIC

clean:
	rm -f *~
	rm -f *.x
	rm -f *.mod
	rm -f *.pyf
	rm -f *.so
	rm -f *.rej
	rm -f iriflip_modified.for
	rm -f irisub_modified.for

test_gfortran:
	make patch_iriflip;
	make patch_irisub;
	gfortran $(f77flags) -o test_gfortran.x simple_iri_driver.for $(fortran_files)

patch_iriflip:
	python delete_iriflip_comments.py

sig:
	f2py -m $(mod) -h sig_file.pyf $(fortran_files) only: $(only) : | tee out1

patch_sig:
	patch sig_file.pyf -i sig.patch -o sig_file_patched.pyf

mod:
	# This works well [1]
	f2py -c sig_file_patched.pyf $(fortran_files) --fcompiler=gfortran --f77flags="$(f77flags)" | tee out2

compile:
	make patch_iriflip;
	make sig;
	make patch_sig;
	make mod;

install:
	cp $(mod).so ../model_atmosphere/modules/

#   [1]  I had a huge bug where iri was not repeating itself on python calls. I found
# this github https://github.com/sdelarquier/iri and the compiler options
#   	f2py --overwrite-signature irisub.for -m ${MODULE} -h iri.pyf
#   	gfortran -w -O2 -fbacktrace -fno-automatic -fPIC -c *.for
# 	    f2py --f77flags="-w" -c iri.pyf *.for
# makes it so it works! I need to look into these compiler options and
# how they relate to IRI

#remove_comments:
#	perl -pi -e 's/![\.|-].*$//g' iriflip.for
