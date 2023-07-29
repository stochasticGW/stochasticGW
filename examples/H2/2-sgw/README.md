# H2 standalone example, stochastic GW calculation

In this step we will compute the ionization potential of H2 using StochasticGW.
See "H2: Combined Stochastic and Deterministic treatment" in the StochasticGW
User manual and Tutorial for a detailed description of the calculation and
parameters.

1. Examine the input file, `INPUT`. Note also that the coordinates file,
   `cnt.ini`, the wavefunction file `wf.txt`, the random seed file,
   `random.inp`, and the counter file, `counter.inp` must be present.

2. Run `./sgw.x` (or, on an HPC cluster, call `sgw_gpu.x` using an appropriate
   job script). Examine the output file and compare the result to the value
   in the Tutorial.

