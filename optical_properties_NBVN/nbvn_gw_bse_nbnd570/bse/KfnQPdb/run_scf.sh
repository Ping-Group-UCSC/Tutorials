#!/bin/bash
#SBATCH -p skx-dev
#SBATCH -N 3
#SBATCH -t 01:00:00
#SBATCH --ntasks-per-node=48
#SBATCH -J pw%j
#SBATCH -o job%j.out
#SBATCH -e job%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jxu153@ucsc.edu

export OMP_NUM_THREADS=1
MPICMD="ibrun"
MPICMDS="mpirun -n 1"
DIRY="/home1/06235/tg855346/work/yambo_codes/yambo-4.1.4/bin"
#DIRY="/home1/06235/tg855346/work/yambo_codes/yambo-4.3.2/bin"
DIRPW="/home1/06235/tg855346/qe_codes/qe-6.1.0/bin"

$MPICMD $DIRPW/pw.x -nk 9 < scf.in > scf.out
$MPICMD $DIRPW/pw.x -nk 9 < nscf.in > nscf.out
