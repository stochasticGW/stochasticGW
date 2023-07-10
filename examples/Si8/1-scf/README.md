# Si8 example, generating the orbitals

This step computes the orbitals using an SCF calculation.
The results are then converted into StochasticGW format using the tool
`qe2sgw.x`.

1. Examine `pwmt-si.in`. Note that we are computing 32 orbitals for the
   gamma k-point, which includes the 16 occupied orbitals plus 16 vacant
   orbitals. Run the SCF calculation using Quantum Espresso:
   `pw.x -input pwmt-si >& silicon.out` (or, on NERSC Perlmutter, run using an
   appropriate job script).

2. Examine `qe2sgw.in`. The 'output' setting, 'targeted', means that the HOMO
   and LUMO orbitals will be targeted for post-DFT treatment in stochastic-GW.
   Four files are generated: `cnt.ini` containing atomic coordinates, `dens.txt`
   containing the charge density, and `homo.txt` and `lumo.txt` containing
   the HOMO and LUMO orbitals, respectively. Note that the line beginning 'orb'
   should be left blank unless off-diagonal Sigma elements are to be computed 
   with StochasticGW.

3. Run `./qe2sgw.x`. This utility generates and input for Quantum Espresso's
   `pp.x` utility and calls the utility to generate cube files. The cube files
   are then processed by `qe2sgw.x` to generate the required sGW files.

