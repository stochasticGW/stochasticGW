# Graphene example, generating the orbitals

This step computes the orbitals using an SCF calculation.
The results are then converted into StochasticGW format using the tool
`qe2sgw.x`.

1. Examine `pwmt-graphene.in`. Note that we are computing 16 orbitals for the
   gamma k-point, which includes the 8 occupied orbitals plus 8 vacant
   orbitals. Run the SCF calculation using Quantum Espresso:
   `pw.x -input pwmt-graphene >& graphene.out` (or, on an HPC cluster,
   run using an appropriate job script).

2. Examine `qe2sgw.in`. The 'output' setting, 'full', will generate two files,
   `cnt.ini` containing atomic coordinates, and `wf.txt` containing the
   full set of orbitals.

3. Run `./qe2sgw.x`. This utility generates and input for Quantum Espresso's
   `pp.x` utility and calls the utility to generate cube files. The cube files
   are then processed by `qe2sgw.x` to generate the required sGW files.

