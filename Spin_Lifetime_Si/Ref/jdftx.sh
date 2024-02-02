#!/bin/bash
#SBATCH -p skx-dev
#SBATCH -N 1
#SBATCH -t 02:00:00
#SBATCH --ntasks-per-node=8
#SBATCH -J si_t1

module load gsl fftw3

MPICMD="ibrun"
DIRJ="/home1/06235/tg855346/jdftx_codes/jdftx-test/build"
DIRF="/home1/06235/tg855346/jdftx_codes/jdftx-test/build-FeynWann-mod"

${MPICMD} ${DIRJ}/jdftx -i totalE.in > totalE.out
${MPICMD} ${DIRJ}/jdftx -i bandstruct.in > bandstruct.out
${MPICMD} ${DIRJ}/phonon -ni phonon.in > split.out
