#!/bin/bash

# This serial script links files for running the naphthalene example

## Link files in the sGW run directory
cd 2-sgw
# sGW random seed file
ln -sf ../../common/random.inp .
# sGW counter input
ln -sf ../../common/counter.inp .
# Pseudopotential directory
ln -sf ../../../PP
# sGW executables
ln -sf ../../../src/sgw_gpu.x .
# Utility script for plotting filter coefficients, reconstructed filter
ln -sf ../../../utils/vis-filter/plotfilter.py .
cd ..

