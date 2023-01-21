#!/bin/bash

# This serial script links files for running the C60 example
# sGW is run using full wf.bin file in project directory
# (should be modified if using local copy of wf.bin)

## Link files in the sGW run directory
cd 2-sgw
# Wavefunction path (modify if needed)
ln -sf /global/cfs/cdirs/nesap/BerkeleyGW/SGW/C60/wf.bin .
# sGW random seed file
ln -sf ../../common/random.inp .
# sGW counter input
ln -sf ../../common/counter.inp .
# Pseudopotential directory
ln -sf ../../../PP
# sGW executables
ln -sf ../../../src/sgw_cpu.x .
ln -sf ../../../src/sgw_gpu.x .
cd ..

