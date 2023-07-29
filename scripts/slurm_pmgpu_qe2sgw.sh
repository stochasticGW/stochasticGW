#!/bin/bash
#
# This script runs Quantum Espresso's 'pp.x' utility to generate
# the cube files needed by 'qe2sgw_stepwise.x', for a list of orbitals
# and density given below. This is appropriate for large systems for
# which one wants to extract a small subset of the molecular orbitals
# to pass to sGW.
#
# One should run 'qe2sgw_stepwise.x' before calling this script to
# generate the input files listed (below) for 'pp.x'.
# One should then run 'qe2sgw_stepwise.x' again after this script 
# completes to process the cube files into 'sgwinp.txt'for stochasticGW.
#
#SBATCH --account=<your-account>
#SBATCH --constraint=gpu
#SBATCH --qos=debug

#SBATCH --time=00:29:59
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

#SBATCH --cpus-per-task=32
#SBATCH --gpus-per-task=1
#SBATCH --gpu-bind=map_gpu:0,1,2,3

#SBATCH --job-name="qe2sgw"

module load espresso
export OMP_NUM_THREADS=1

# List of orbitals + density to run
# (specific to H2O example)
srun pp.x -i pp_orb_1.in > pp_orb_1.out
srun pp.x -i pp_orb_2.in > pp_orb_2.out
srun pp.x -i pp_orb_3.in > pp_orb_3.out
srun pp.x -i pp_orb_4.in > pp_orb_4.out
srun pp.x -i pp_orb_5.in > pp_orb_5.out
srun pp.x -i pp_orb_6.in > pp_orb_6.out
srun pp.x -i pp_orb_7.in > pp_orb_7.out
srun pp.x -i pp_orb_8.in > pp_orb_8.out
srun pp.x -i pp_dens.in > pp_dens.out

