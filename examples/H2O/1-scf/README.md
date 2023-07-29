# H2O example, generating the orbitals

This step computes the orbitals using an SCF calculation.
The results are then converted into StochasticGW format using the tool
`qe2sgw.x`.

1. Examine `pwmt-h2o.in`. Note that we are computing 8 orbitals for the
   gamma k-point, which includes the 4 occupied orbitals plus 4 vacant
   orbitals. Run the SCF calculation using Quantum Espresso:
   `pw.x -input pwmt-h2o >& h2o.out` (or, on an HPC cluster, run using an
   appropriate job script).

2. Examine `qe2sgw.in`. The 'output' setting, 'full', means that two files
   will be generated for stochastic-GW. The file `cnt.ini` contains the 
   coordinates of the atoms while `wf.txt` contains all orbitals resulting
   from the SCF run. Note that the line beginning 'orb' should be left blank
   when 'output' is set to 'full'.

3. Run `./qe2sgw.x`. This utility generates and input for Quantum Espresso's
   `pp.x` utility and calls the utility to generate cube files. The cube files
   are then processed by `qe2sgw.x` to generate the required sGW files.

