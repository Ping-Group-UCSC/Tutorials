#!/usr/bin/env bash

############################################
# wrapper script for converting qe charge
# density to jdftx charge density and for
# converting qe atomic positions to jdftx
############################################

# python script for converting charge density
conv=/home/tjsmart/Programs/Ping-Group/qe-jdftx-convert/conv.py
# number of atoms
nat=97

for d in Q*/ ; do
    cd $d
    $conv temp/BN.save/ bnjdftx         # convert charge
    pw2j.sh $nat relax.out > atom.pos   # generate atom.pos file
    cd ..
done
