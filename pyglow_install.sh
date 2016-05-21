#!/bin/bash


#git clone git://github.com/timduly4/pyglow.git; 
#cd pyglow/; 

cd ./pyglow/models/; 
make all; 
cd ../../; 
python ./setup.py install --user;

# Run update indices:
cd ~/; 
python -c "from pyglow import pyglow; pyglow.update_indices()"

