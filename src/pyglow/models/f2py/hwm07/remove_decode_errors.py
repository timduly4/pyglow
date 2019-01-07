import os

import glob

forfiles = glob.glob('*.f90')

for f in forfiles:
	content = open(f,'r',errors='replace').read()
	open(f,'w').write(content)


# this perl command deletes all the comments in iriflip...
# my compiler was having trouble with them, for some reason


