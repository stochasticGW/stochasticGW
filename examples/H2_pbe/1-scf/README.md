# H2 PBE example, generating the orbitals

This step computes the orbitals of H2 via an SCF calculation performed in
Quantum Espresso using the PBE functional. The results are then converted into
StochasticGW format using two different tool chains.

1. Examine `h2.in`. Note that we are computing 2 orbitals for the
   gamma k-point, one occupied plus one vacant. Run the SCF calculation 
   using Quantum Espresso: `pw.x -input h2 >& h2.out` (or, on an HPC cluster,
   run using an appropriate job script).

2. The first way to prepare data for StochasticGW is the direct method using
   the `qe2sgw` tool. Examine `qe2sgw.in`. The 'output' setting, 'full', means
   that two files will be generated for StochasticGW: `cnt.ini`
   contains the coordinates of the atoms, and `wf.txt` contains orbitals
   from the SCF run. Note that the line beginning with 'orb' must be left
   blank if 'output' is set to 'full'.

3. Run `./qe2sgw.x > qe2sgw.out`. This utility generates and input for 
   Quantum Espresso's `pp.x` utility and calls the utility to generate cube
   files. The cube files are then processed by `qe2sgw.x` to generate the 
   required sGW files. This example is continued in the `../2-sgw-qe`
   directory.

4. The second way to prepare data for StochasticGW is via BerkeleyGW. First,
   we use the `pw2bgw` tool to prepare the wavefunction and density in BGW
   format. Then we use the `bgw2sgw` tool to convert the data into sGW format.
   The first (`pw2bgw`) step is performed here. Examine `pw2bgw.in`.

5. Run `pw2bgw.x -i pw2bgw.in > pw2bgw.out`. Note that two files are generated
   for BerkeleyGW, `WFN` and `RHO`. This example is continued in the `../2-sgw-bgw`
   directory.

