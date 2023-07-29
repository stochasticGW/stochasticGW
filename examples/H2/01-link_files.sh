#!/bin/bash

# This serial script links files for running the standalone H2 example

## Link files in the sGW run directory
cd 2-sgw
# sGW random seed file
ln -sf ../../common/random_alt.inp random.inp
# sGW counter input
ln -sf ../../common/counter.inp .
# Pseudopotential directory
ln -sf ../../../PP
# sGW executable
ln -sf ../../../src/sgw_gpu.x .
cd ..

