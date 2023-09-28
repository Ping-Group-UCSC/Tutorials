#!/bin/bash
#SBATCH -p RM-small
#SBATCH -N 2
#SBATCH -t 08:00:00
#SBATCH --ntasks-per-node=28
#SBATCH -J jdftx%j
#SBATCH -o job%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jxu153@ucsc.edu

DIRY="/home/jxuucsc/yambo_codes/yambo-4.1.4/bin"
DIRPW="/home/jxuucsc/qe_codes/qe-6.1/bin"
MPICMD="mpirun -np  $SLURM_NTASKS"
MPICMDS="mpirun -np 1"

$MPICMD $DIRPW/pw.x < scf.in > scf.out
$MPICMD $DIRPW/pw.x < nscf.in > nscf.out
