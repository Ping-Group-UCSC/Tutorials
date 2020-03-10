#!/bin/sh
####################################################
# This is a sample script to run scf total-energy
# calculations on a unit cell of Al using four 
# different values for the Monkhorst-Pack grid divisions
# and two different values of smearing parameter degauss
#
# Written by Shobhana Narasimhan
#
# You should copy this file and modify it as 
# appropriate for the tutorial.
####################################################
# Notes:
#
# 1. You can loop over a variable by using the 
#    'for...do...done' construction. 
# 2. Variables can be referred to within the script 
#    by typing the variable name preceded by the '$' 
#    sign. 
#
####################################################
# Important initial variables for the code
# (change these as necessary)
####################################################

NAME1='degauss'
NAME2='nk'

####################################################

for DEG in 0.02 0.04 0.06 0.08 0.10
do
for NKDIV in 6 8 12 16
do
cat > al_$NAME1.${DEG}_$NAME2.${NKDIV}_scf.in << EOF
 &control
    calculation = 'scf',
    verbosity = 'high'
    prefix = 'Al_exc3'
    pseudo_dir = './'
    outdir = 'temp'
 /
 &system
    ibrav =  2, 
    celldm(1) = 7.5, 
    nat =  1, 
    ntyp = 1,
    ecutwfc = 12.0,
    occupations = 'smearing',
    smearing = 'marzari-vanderbilt',
    degauss = $DEG
 /
 &electrons
    mixing_beta = 0.7
 /

ATOMIC_SPECIES
 
Al 26.98  Al_ONCV_PBE-1.0.upf

ATOMIC_POSITIONS (alat)
 Al 0.0 0.0 0.0

K_POINTS (automatic)
  $NKDIV $NKDIV $NKDIV 1 1 1
EOF

pw.x < al_$NAME1.${DEG}_$NAME2.${NKDIV}_scf.in >  al_$NAME1.${DEG}_$NAME2.${NKDIV}_scf.out

done

done
