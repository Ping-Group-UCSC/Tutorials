#!/bin/bash
#SBATCH -p RM-small
#SBATCH -N 2
#SBATCH -t 06:00:00
#SBATCH --ntasks-per-node=14
#SBATCH -J jdftx%j
#SBATCH -o job%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jxu153@ucsc.edu

DIRY="/home/jxuucsc/yambo_codes/yambo-4.1.4/bin"
DIRPW="/home/jxuucsc/qe_codes/qe-6.1/bin"
MPICMD="mpi -np  $SLURM_NTASKS"
MPICMDS="mpirun -np 1"

$MPICMD $DIRY/yambo -F rpa.in
