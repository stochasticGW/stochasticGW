# BN example, generating the orbitals

This step computes the orbitals using an SCF calculation.
The results are then converted into StochasticGW format using the tool
`qe2sgw.x`.

1. Examine `pwmt-bn.in`. Note that we are computing 32 orbitals for the
   gamma k-point, which includes the 16 occupied orbitals plus 16 vacant
   orbitals. Run the SCF calculation using Quantum Espresso:
   `pw.x -input pwmt-bn >& bn.out` (or, on an HPC cluster, run using an
   appropriate job script).

2. Copy the qe2sgw input file via `cp qe2sgw_full.in qe2sgw.in`. The 'output'
   setting, 'full', will generate two files, `cnt.ini` containing atomic 
   coordinates, and `wf.txt` containing the full set of orbitals.

3. Run `./qe2sgw.x`. This utility generates and input for Quantum Espresso's
   `pp.x` utility and calls the utility to generate cube files. The cube files
   are then processed by `qe2sgw.x` to generate the required sGW files.

4. The files generated in the previous steps are not in the correct format for
   the sGW filtering calculation, which requires either 'targeted' or 'unified'
   format. Copy the alternate qe2sgw input file by executing,
   `cp qe2sgw_unified.in qe2sgw.in`.

5. Rerun `./qe2sgw.x` to generate `sgwinp.txt` required by StochasticGW for 
   filtering.

