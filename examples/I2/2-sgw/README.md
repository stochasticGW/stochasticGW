# I2 example, stochastic GW calculation

In this step we will compute the orbital energies of each of the orbitals of I2
to test convergence of the starting energies with respect to the number of grid
points. Since only initial energies are desired and no stochastic GW calculation
is to be performed, several parameters ('nmctot', 'ntddft', 'gamma', 'nchb')
are set to be much smaller than in a production calculation to save numerical
effort. Note also that this job will throw an error during the stochastic GW
calculation when attempting to evaluate the Sigma operator for orbitals i and
j!=i due to an array having complex values; this is expected behavior.

1. Examine the input file, `INPUT`. Note that the parameter 'orbj_indx'
   contains the list of all orbitals to enforce reading them from `sgwinp.txt`.

2. Run `./sgw.x` (or, on an HPC system, using an appropriate job script calling
   `./sgw_gpu.x`).

3. Examine the energies of the orbitals generated from Quantum Espresso (as read
   from `sgwinp.txt`) by searching for the tag in the output file:
   'ORBJ(1) read: spin = 1; energy ='
   Compare these values with the energies computed with StochasticGW, which can
   be found by searching for the tag:
   'KS energy computed for orbital 1'
   Do the orbital energies computed with Quantum Espresso and with StochasticGW
   agree to better than 0.001 Hartree?

