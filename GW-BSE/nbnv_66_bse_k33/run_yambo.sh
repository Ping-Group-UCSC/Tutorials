#!/bin/bash
#SBATCH --partition=long
#SBATCH --nodes=12
#SBATCH --ntasks-per-node=36
#SBATCH --time=24:00:00
#SBATCH --job-name=pw%j
#SBATCH --account=cfn37607
#SBATCH --mail-user=jxu153@ucsc.edu
#SBATCH --mail-type=FAIL

module load intel

export OMP_NUM_THREADS=1
MPICMD="srun -n $SLURM_NTASKS"
MPICMDS="mpirun -n 1"
DIRY="/sdcc/u/jxuucsc/yambo_codes/yambo-4.1.4/bin"
#DIRY="/sdcc/u/jxuucsc/yambo_codes/yambo-4.3.2/bin"
DIRPW="/sdcc/u/jxuucsc/qe_codes/qe-6.1/bin"

$MPICMD $DIRY/yambo -F eps0.in
#$MPICMD $DIRY/yambo -F bse.in
