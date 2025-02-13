#!/bin/bash

PROCESSOR_COUNTER=20

echo "Starting flux generation..."
date 
cd FLUX
./ffs_flux.py -c $PROCESSOR_COUNTER
cd ..


echo "Starting interface..."
date 
cd IF1
./ffs_shoot.py -c $PROCESSOR_COUNTER
cd ..;


echo "Starting interface..."
date 
cd IF2
./ffs_shoot.py -c $PROCESSOR_COUNTER
cd ..;


echo "Starting interface..."
date 
cd IF3
./ffs_shoot.py -c $PROCESSOR_COUNTER
cd ..


echo "Starting interface..."
date 
cd IF4
./ffs_shoot.py -c $PROCESSOR_COUNTER
cd ..


echo "Starting interface..."
date 
cd IF5
./ffs_shoot.py -c $PROCESSOR_COUNTER
cd ..



echo "Finished ..."
date 
