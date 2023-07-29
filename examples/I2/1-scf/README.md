# I2 example, generating the orbitals

This step computes the orbitals using a series of SCF calculations. We will
vary the grid densities in the SCF runs and then compare the orbital energies
in StochasticGW with those in Quantum Espresso, increasing the density
of points until the energies converge in StochasticGW. We will use the
tool `qe2sgw.x` with the 'unified' option to convert the orbitals into
StochasticGW format after each SCF run.

1. Examine `pwmt-i2.in`. Note that we are computing 14 orbitals for the
   gamma k-point, which includes the 7 occupied orbitals plus 7 vacant
   orbitals. Note the values of the 'ecutwfc' and 'ecutrho' parameters (36.0
   and 36.00001 Ry, respectively). Run the SCF calculation using Quantum
   Espresso: `pw.x -input pwmt-i2 >& i2.out` (or, on an HPC cluster, run using
   an appropriate job script).

2. Examine `qe2sgw.in`. The 'output' setting, 'unified', means that a single
   file, `sgwinp.txt` will be generated. This file will contain all orbitals
   listed on the line beginning with 'orb', which in this case, includes all
   14 orbitals.

3. Run `./qe2sgw.x`. This utility generates and input for Quantum Espresso's
   `pp.x` utility and calls the utility to generate cube files. The cube files
   are then processed by `qe2sgw.x` to generate the required sGW files.

4. Go to the folder `../2-sgw` and follow the instructions there to run the
   StochasticGW calculation. Once you have finished, return to this directory.

5. Repeat steps 1-4, but, leaving 'ecutwfc' constant, increase 'ecutrho' to the
   following values in `pwmt-i2.in`:
   54.00001 Ry, 72.00001 Ry, 108.00001 Ry, 144.00001 Ry, 180.00001 Ry.
   Note that StochasticGW only supports even grid dimensions along each (x,y,z)
   dimension, so 'ecutrho' must be chosen in Quantum Espresso to give even 
   values of (nx,ny,nz) in `sgwinp.txt`.

