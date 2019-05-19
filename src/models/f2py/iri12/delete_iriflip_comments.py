from __future__ import print_function
import os

# this perl command deletes all the comments in iriflip...
# my compiler was having trouble with them, for some reason
cmd = "perl -pi -e 's/![\.|-].*$//g' iriflip.for"
print(cmd)
os.system(cmd)

