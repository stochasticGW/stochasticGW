#!/bin/bash
#
# Slurm job script for running StochasticGW on Perlmutter-CPU
# Set '--ntasks=<#-MPI-proceses>' below; limit 128 MPI processes per node
# and set 'cpus-per-task' = 256 / 'ntasks-per-node'
#

#SBATCH --account=<your-account>
#SBATCH --constraint=cpu
#SBATCH --qos=debug

#SBATCH --time=00:29:59
#SBATCH --ntasks=128

#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=2

#SBATCH -J sGW
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err

module swap PrgEnv-nvidia PrgEnv-gnu
module load cray-fftw

export SLURM_CPU_BIND="cores"

srun sgw_cpu.x


