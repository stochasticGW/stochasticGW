#!/bin/bash
#
# This script runs Quantum Espresso on 1 CPU + 1 GPU, followed by qe2sgw
# This is appropriate for small systems, since pp.x (called by qe2sgw.x)
# will generate the list of cube files for all molecular orbitals and the
# density in serial.
#
#SBATCH --account=<your-account>
#SBATCH --constraint=gpu
#SBATCH --qos=debug

#SBATCH --time=00:09:59
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

#SBATCH --cpus-per-task=32
#SBATCH --gpus-per-task=1
#SBATCH --gpu-bind=map_gpu:0,1,2,3

#SBATCH --job-name="H2O"

module load espresso

export OMP_NUM_THREADS=1

srun pw.x -in pwmt-h2o.in >& h2o.out
srun qe2sgw.x > qe2sgw.out
