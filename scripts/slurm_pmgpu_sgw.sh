#!/bin/bash
#
# Slurm job script for running StochasticGW on Perlmutter-GPU
# Set '--ntasks=<#-MPI-proceses>' below; limit 4 MPI orocesses, 4 GPUs per node
#

#SBATCH --account=<your-account>
#SBATCH --constraint=gpu
#SBATCH --qos=debug

#SBATCH --time=00:29:59
#SBATCH --ntasks=4

#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=32
#SBATCH --gpus-per-task=1
#SBATCH --gpu-bind=map_gpu:0,1,2,3

#SBATCH -J sGW
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err

module swap PrgEnv-gnu PrgEnv-nvidia
module load cray-fftw

export SLURM_CPU_BIND="cores"

srun sgw_gpu.x

