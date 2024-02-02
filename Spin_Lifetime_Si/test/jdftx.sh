#!/bin/bash
##SBATCH -p debug
#SBATCH -N 18
#SBATCH -t 99:00:00
#SBATCH --ntasks-per-node=2
#SBATCH -J jdftx

source ~/.bashrc
module use /home/jxu153/modulefiles
module load myopenmpi-4.0.2_gcc-4.8.5
#module purge
#module load gnu openmpi mkl gsl
module list

MPICMD="mpirun -np $SLURM_NTASKS"
DIRJ="/export/data/share/jxu/jdftx_codes/jdftx-feynwann-stable-202109/build"
DIRF="/export/data/share/jxu/jdftx_codes/jdftx-feynwann-stable-202109/build-FeynWann"

${MPICMD} ${DIRJ}/jdftx -i scf.in > scf.out
${MPICMD} ${DIRJ}/jdftx -i totalE.in > totalE.out
#${MPICMD} ${DIRJ}/jdftx -i bandstruct.in > bandstruct.out
#${MPICMD} ${DIRJ}/phonon -ni phonon.in > split.out
