# Graphene example, stochastic GW calculation

In this step we will compute the orbital energies of each orbital of graphene
and compare the stochastic GW values with those from the preliminary DFT
calculation.

1. Examine the input file, `INPUT`. Note that the parameter 'orb_indx'
   indexes the lowest energy orbital (1).

2. Run `./sgw.x` (or, on an HPC system, using an appropriate job script calling
   `./sgw_gpu.x`).

3. Modify `INPUT` and repeat step 2 for 'orb_indx' values of 2,3...16 and
   examine the final quasiparticle energies. The corresponding DFT energies can
   be found by searching for the tag `vs. eorb read from file:` in the output.
   Which orbitals showed the largest change from the GW correction?

