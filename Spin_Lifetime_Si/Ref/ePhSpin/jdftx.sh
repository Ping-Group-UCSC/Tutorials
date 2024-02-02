#!/bin/bash
#SBATCH -p skx-dev
#SBATCH -N 4
#SBATCH -t 02:00:00
#SBATCH --ntasks-per-node=48
#SBATCH -J si_t1

module load gsl fftw3

MPICMD="ibrun"
DIRJ="/home1/06235/tg855346/jdftx_codes/jdftx-test/build"
DIRF="/home1/06235/tg855346/jdftx_codes/jdftx-test/build-FeynWann-mod"

${MPICMD} ${DIRF}/ePhSpinRelax -i ePhSpinRelax.in -G4 > ePhSpinRelax.out
