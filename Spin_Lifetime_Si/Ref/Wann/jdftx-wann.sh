#!/bin/bash
#SBATCH -p skx-dev
#SBATCH -N 1
#SBATCH -t 00:30:00
#SBATCH --ntasks-per-node=48
#SBATCH -J si_t1

module load gsl fftw3

MPICMD="ibrun"
DIRJ="/home1/06235/tg855346/jdftx_codes/jdftx-test/build"
DIRF="/home1/06235/tg855346/jdftx_codes/jdftx-test/build-FeynWann-mod"

prfx=wannier
#${MPICMD} ${DIRJ}/wannier -i ${prfx}.in > ${prfx}.out
#exit 0
python rand_wann-centers.py
cp wannier.in0 wannier.in
cat rand_wann-centers.dat >> wannier.in
rm rand_wann-centers.dat
${MPICMD} ${DIRJ}/wannier -i ${prfx}.in > ${prfx}.out
