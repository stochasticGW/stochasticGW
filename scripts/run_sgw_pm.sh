#!/bin/bash
#
# Slurm job script for running StochasticGW on Perlmutter-GPU
# Set '-n <#-MPI-proceses>' below; uses up to 4 MPI orocesses, 4 GPUs per node
#
#SBATCH -A m1759
#SBATCH -q debug
#SBATCH -n 1
#SBATCH -t 0:29:59

#SBATCH -C gpu
#SBATCH -c 32
#SBATCH --ntasks-per-node=4
#SBATCH --gpus-per-task=1
#SBATCH --gpu-bind=map_gpu:0,1,2,3
#SBATCH -J sGW
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err

module swap PrgEnv-gnu PrgEnv-nvidia
module load cray-fftw

export SLURM_CPU_BIND="cores"

srun sgw.x
