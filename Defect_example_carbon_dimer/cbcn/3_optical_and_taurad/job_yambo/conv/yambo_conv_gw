#!/bin/bash
#SBATCH --job-name=xxf
#SBATCH --output=qe.%j
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=32
#SBATCH --time=24:00:00
#SBATCH --partition=cpuq
# SBATCH --dependency=afterany:29392
#SBATCH --account=cfn307395
#SBATCH --mail-user=kli103@ucsc.edu
#SBATCH --mail-type=FAIL

####################### Kairay ###########################################
# echo "Start:"; date
# module add intel/17.0.5.239 impi/2017
# NORES=$(($SLURM_NTASKS_PER_NODE * $SLURM_JOB_NUM_NODES))
# MPICMD="mpirun -genv I_MPI_FABRICS shm:ofa -n $SLURM_NTASKS"
# PWDIR="/export/data/share/wufeng/programs-intel2017.5/qe-6.1-scal/bin"
# YAMDIR=/export/data/share/jxu/yambo-codes/yambo-4.1.4/bin
# YAMDIR=/export/data/share/jxu/yambo-codes/yambo-4.4.0/bin
##########################################################################

################### Stampede2 ###########################################
# echo "Start:"; date
# echo "Running program on $SLURM_JOB_NUM_NODES nodes with $SLURM_NTASKS total tasks, with each node getting $SLURM_NTASKS_PER_NODE running on cores."
# export OMP_NUM_THREADS=1
# MPICMD="ibrun"
# PWDIR=/home1/06931/kli1003/work/programs/qe-6.1.0/bin
# YAMDIR=/home1/06931/kli1003/work/yambo_install/yambo-4.1.4/bin
#########################################################################

####################################### Lux #########################################################################################################
 echo "Start:"; date;
 echo "Running program on $SLURM_JOB_NUM_NODES nodes with $SLURM_NTASKS total tasks, with each node getting $SLURM_NTASKS_PER_NODE running on cores."
 module load intel/impi
 export OMP_NUM_THREADS=1
 MPICMD="mpirun -n $SLURM_NTASKS --ppn 40"
 PWDIR="/data/users/jxu153/codes/qe/qe-6.1.0/bin"
 YAMDIR=/data/users/jxu153/codes/yambo/yambo-4.1.4/bin
# YAMDIR=/data/users/jxu153/codes/yambo/yambo-4.4.0/bin
# YAMDIR=/home/kli103/work/programs/yambo-4.4.0/bin
#####################################################################################################################################################

############################### BNL ###########################
# module load intel
# echo "Start:"; date
# export OMP_NUM_THREADS=1
# MPICMD="srun -n $SLURM_NTASKS"
# MPICMDS="mpirun -n 1"
# YAMDIR="/sdcc/u/kli/programs/yambo-4.1.4_feng/bin"
# PWDIR="/sdcc/u/kli/programs/qe-6.1.0/bin"
###############################################################

#######################################################################################################
bands='900 1200 1500 1900 2000'
blocks='1 2 3 4 5 6 7 8'
for i in ${bands}
do
  for j in ${blocks}
  do
    $MPICMD $YAMDIR/yambo -F gw_conv_$i'b_'$j'Ry.in' -J $i'b_'$j'Ry'
  done
done
echo "Done"
echo "End:"; date
#######################################################################################################
