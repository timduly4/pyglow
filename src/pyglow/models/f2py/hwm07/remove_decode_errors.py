import glob

forfiles = glob.glob('*.f90')

for f in forfiles:
    content = open(f, 'r').read()
    open(f, 'w').write(unicode(content, errors='ignore'))
