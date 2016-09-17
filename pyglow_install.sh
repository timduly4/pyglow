#!/bin/bash


#git clone git://github.com/timduly4/pyglow.git; 
#cd pyglow/; 

cd ./pyglow/models/; 
make all; 
cd ../../; 

# (tmd, 9/17/16)
# TODO: on a mac, we require:
# --prefix=
# to be appended for this installation command:
python ./setup.py install --user;
# How to check to see if we're on a mac?
# Or, does this command work OK with Linux and we can
# just add it by default?
# Reference: https://goo.gl/7AG01z

# Run update indices:
cd ~/; 
python -c "from pyglow import pyglow; pyglow.update_indices()"

