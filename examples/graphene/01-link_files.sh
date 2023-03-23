#!/bin/bash

# This serial script links files for running the graphene example
# Quantum Espresso is used for SCF, followed by sGW

## Link files in the scf run directory
cd 1-scf
# Pseudopotential file(s)
ln -sf ../../../PP/C.pw-mt_fhi.UPF .
# qe2sgw executable
ln -sf ../../../utils/qe2sgw/qe2sgw.x .
# Utility script for plotting geometry, orbitals, density
ln -sf ../../../utils/vis-qe2sgw/plotorbital.py .
cd ..

## Link files in the sGW run directory
cd 2-sgw
# QE processed output: coordinates file and wavefunction
ln -sf ../1-scf/cnt.ini .
ln -sf ../1-scf/wf.txt .
# sGW random seed file
ln -sf ../../common/random.inp .
# sGW counter input
ln -sf ../../common/counter.inp .
# Pseudopotential directory
ln -sf ../../../PP
# sGW executable
ln -sf ../../../src/sgw_gpu.x .
cd ..

