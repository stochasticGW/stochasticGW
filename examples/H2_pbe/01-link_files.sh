#!/bin/bash

# This serial script links files for running the H2 example
# Quantum Espresso is used for SCF, followed by sGW

## Link files in the scf run directory
cd 1-scf
# Pseudopotential file(s)
ln -sf ../../../PP/H.pw-mt_fhi.UPF .
# qe2sgw executable
ln -sf ../../../utils/qe2sgw/qe2sgw.x .
# Utility script for plotting geometry, orbitals, density
ln -sf ../../../utils/vis-qe2sgw/plotorbital.py .
cd ..

## Link files in the sGW-from-BGW
cd 2-sgw-bgw
# QE processed output: cnt.ini for atomic coordinates, wf.txt for orbita
ln -sf ../1-scf/WFN .
ln -sf ../1-scf/RHO .
# sGW random seed file
ln -sf ../../common/random_alt.inp random.inp
# sGW counter input
ln -sf ../../common/counter.inp .
# Pseudopotential directory
ln -sf ../../../PP
# sGW executable
ln -sf ../../../src/sgw_gpu.x .
cd ..

## Link files in the sGW-from-qe run directory
cd 2-sgw-qe
# QE processed output: cnt.ini for atomic coordinates, wf.txt for orbitals
ln -sf ../1-scf/cnt.ini .
ln -sf ../1-scf/wf.txt .
# sGW random seed file
ln -sf ../../common/random_alt.inp random.inp
# sGW counter input
ln -sf ../../common/counter.inp .
# Pseudopotential directory
ln -sf ../../../PP
# sGW executable
ln -sf ../../../src/sgw_gpu.x .
cd ..

