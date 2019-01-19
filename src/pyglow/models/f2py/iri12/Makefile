fortran_files = cira.for igrf.for iridreg_modified.for irifun.for irisub.for iritec.for iriflip.for
only = iri_sub
mod = iri12py

clean:
	rm *~ *.x *.mod *.pyf *.so
	#iriflip_modified.for

test_gfortran:
	gfortran -std=legacy -w -O2 --backtrace -fno-automatic -fPIC -o test_gfortran.x tims_iri11_driver.for irisub.for irifun.for iritec.for iridreg.for igrf.for cira.for iriflip.for

test_g95:
	g95 -fPIC -fexceptions -o test_g95.x tims_iri11_driver.for irisub.for irifun.for iritec.for iridreg.for igrf.for cira.for iriflip.for

test_ifort:
	ifort -o test_ifort.x tims_iri11_driver.for irisub.for irifun.for iritec.for iridreg.for igrf.for cira.for iriflip.for

#test_f77:
#	f77 -o test_gfortran.x tims_iri11_driver.for irisub.for irifun.for iritec.for iridreg.for igrf.for cira.for iriflip.for

patch_iridreg:
	patch iridreg.for -i iridreg.patch -o iridreg_modified.for

patch_iriflip:
	python3 delete_iriflip_comments.py
	#patch iriflip.for -i iriflip.patch -o iriflip_modified.for

sig:
	f2py -m $(mod) -h sig_file.pyf $(fortran_files) only: $(only) : | tee out1

patch_sig:
	patch sig_file.pyf -i sig.patch -o sig_file_patched.pyf

mod:
#	f2py -c sig_file_patched.pyf $(fortran_files) --fcompiler=f77 --f77flags="-fPIC -fexceptions" | tee out2
	# I KNOW THIS ONE WORKS:
	#f2py -c sig_file_patched.pyf $(fortran_files) --fcompiler=gfortran --f77flags="-std=legacy" | tee out2
	# BUT THIS ONE WORKS BETTER [1]
	f2py -c sig_file_patched.pyf $(fortran_files) --fcompiler=gfortran --f77flags="-std=legacy -w -O2 -fbacktrace -fno-automatic -fPIC" | tee out2


source:
	make patch_iridreg;
	make patch_iriflip;
	sed -i'' -e 's/cm-3//g' iriflip.for
	make sig;
	make patch_sig;


compile: source
	make mod;

install:
	cp $(mod).so ../model_atmosphere/modules/



#   [1]  I had a huge bug where iri was not repeating itself on Python calls. I found
# this github https://github.com/sdelarquier/iri and the compiler options
#   	f2py --overwrite-signature irisub.for -m ${MODULE} -h iri.pyf
#   	gfortran -w -O2 -fbacktrace -fno-automatic -fPIC -c *.for
# 	    f2py --f77flags="-w" -c iri.pyf *.for
# makes it so it works! I need to look into these compiler options and
# how they relate to IRI

#remove_comments:
#	perl -pi -e 's/![\.|-].*$//g' iriflip.for
