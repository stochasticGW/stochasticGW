# Si214 example, generating the orbitals

This step computes the orbitals and density using an SCF calculation.
The results are then converted into StochasticGW format using the tool
`qe2sgw.x`.

1. Examine `Si214_fhi_scf.inp`. Note that we are computing 432 orbitals for the
   gamma k-point, which includes the 428 occupied orbitals plus a few vacant
   orbitals. Run the SCF calculation using Quantum Espresso:
   `pw.x -input Si214_fhi_scf.inp >& Si214_fhi_scf.out`

2. Examine `qe2sgw.in`. The 'output' setting, 'unified', means that a single
   file, `sgwinp.txt` will be generated. This file will contain all orbitals
   listed on the line beginning with 'orb'.

3. Run `./qe2sgw_stepwise.x`. This utility will attempt to construct the sGW
   input using cube files containing the orbitals listed in `qe2sgw.in`.
   Since the cube files do not yet exist, the utility will instead generate the
   inputs required to generate them from Quantum Espresso. After running this
   step there should be three new files: `pp_dens.in`, `pp_orb_428.in`, and
   `pp_orb_429.in`.

3. Use Quantum Espresso's post-processing utility to generate the appropriate
   cube files used to construct the sGW input by running:
   `pp.x -i pp_orb_428.in > pp_orb_428.out`
   `pp.x -i pp_orb_429.in > pp_orb_429.out`
   `pp.x -i pp_dens.in > pp_dens.out`

4. Run `./qe2sgw_stepwise.x` a second time. Now the cube files from the
   previous step have been processed, and the file `sgwinp.txt` should be
   generated.
