# H2 PBE example, stochastic GW calculation direct from QE starting point

In this step we will compute the ionization potential of H2 using StochasticGW.
Here, StochasticGW uses the `cnt.ini` and `wf.txt` inputs from the `1-scf` 
directory which had been previously generated directly from their respective
Quantum Espresso outputs by the `qe2sgw` utility.

1. Examine the input file, `INPUT`. Note that, unlike in the standalone
   example, we specify the density functional via its LibXC code.

2. Run `./sgw.x` (or, on an HPC cluster, call `sgw_gpu.x` using an appropriate
   job script). Examine the output file and compare the result to the value
   in the Tutorial.

