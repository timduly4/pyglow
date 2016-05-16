#!/bin/bash


#git clone git://github.com/timduly4/pyglow.git; 

cd pyglow/; 
cd ./pyglow/models/; 
make all; 
cd ../../; 
python ./setup.py install --user;

