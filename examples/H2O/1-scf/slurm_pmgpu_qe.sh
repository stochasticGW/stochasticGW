#!/bin/bash
#SBATCH --job-name="H2O"
#SBATCH --time=00:09:59
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gpus-per-task=1
#SBATCH --cpus-per-task=32
#SBATCH --account=<your_account_here>
#SBATCH --constraint=gpu
#SBATCH --qos=debug

module load espresso

export OMP_NUM_THREADS=1

srun pw.x -in pwmt-h2o.in >& h2o.out
srun qe2sgw.x > qe2sgw.out
