#!/bin/bash
#SBATCH -p RM-small
#SBATCH -N 1
#SBATCH -t 01:00:00
#SBATCH --ntasks-per-node=1
#SBATCH -J jdftx%j
#SBATCH -o job%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jxu153@ucsc.edu

DIRY="/home/jxuucsc/yambo_codes/yambo-4.1.4/bin"
DIRPW="/home/jxuucsc/qe_codes/qe-6.1/bin"
MPICMD="mpirun -np  $SLURM_NTASKS"
MPICMDS="mpirun -np 1"

cd ./out/*.save
$MPICMDS $DIRY/p2y
cd ../..
mv ./out/*.save/SAVE ./
RL=`grep "Max WF components" ./out/*.save/l_stderr  | awk '{print $7}'`
RL=`echo "scale=2;$RL/(sqrt(45./15.))^3" | bc -l`
RL=`echo "$RL/1+1" | bc`
echo "RL = $RL"
sed "s/XXX/$RL/g" init.in0 > init.in
