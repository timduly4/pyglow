#!/bin/bash


#git clone git://github.com/timduly4/pyglow.git; 
#cd pyglow/; 

cd ./pyglow/models/; 
make all; 
cd ../../; 

# On a Mac, we require --prefix=
# (so spaces afterwards)
# Reference: https://goo.gl/7AG01z
# (tmd, 9/17/16)
python ./setup.py install --user --prefix=;

# Run update indices:
cd ~/; 
#python -c "from pyglow import pyglow; pyglow.update_indices()"

