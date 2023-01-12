#!/bin/bash
#
# This script runs StochasticGW on 1 GPU node, with 4 MPI ranks and 4 GPUs
#
#SBATCH -A <your_account_here>
#SBATCH -C gpu
#SBATCH -q debug
#SBATCH -t 00:29:59
#SBATCH -n 4
#SBATCH --ntasks-per-node=4
#SBATCH -c 32
#SBATCH --gpus-per-task=1
#SBATCH --gpu-bind=map_gpu:0,1,2,3

export SLURM_CPU_BIND="cores"

srun sgw.x > H2O_sgw.out

