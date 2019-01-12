import glob
import sys

if sys.version_info[0] >= 3:
    forfiles = glob.glob('*.f90')
    for f in forfiles:
        content = open(f, 'r', errors='replace').read()
        open(f, 'w').write(content)
