# Si8 example, stochastic GW calculation

In this step we will run the stochastic GW calculation on a unit cell of bulk
silicon to compute the quasi-particle energies of the HOMO and LUMO orbitals.
This example uses gapped-filtering to prepare the initial stochastic states.
Since the output format of `qe2sgw` was set to 'targeted' in the previous step,
this example illustrates using symlinks to select the orbital of choice.

1. Examine the input file, `INPUT`. Note that the parameter 'orb_indx' is set
   to a value of '-1' indicating that the files `orb.txt` and `dens.txt` are to
   be read. In particular, `orb.txt` should be linked to `homo.txt` in the SCF
   directory as we will first run sGW to compute the QP energy of the HOMO
   orbital. Note also that the HOMO and LUMO energies must be provided in
   `INPUT` to use gapped-filtering in this example.

2. Run `./sgw.x` (or on NERSC-Perlmutter, `./sgw_gpu.x` using a job script
   which targets the Perlmutter GPU partition). Run the filter plotting script,
   `python plotfilter.py filtercheb.dat` to verify that the filter is converged
   with respect to the number of Chebyshev polynomials.

3. To compute the QP energy of the LUMO orbital, a separate calculation must be
   performed. Modify the existing link to the HOMO orbital file so that it
   points to the LUMO file by executing `ln -sfn ../1-scf/lumo.txt orb.txt`.

4. Run `./sgw.x` (or on NERSC-Perlmutter, `./sgw_gpu.x` using a job script
   which targets the Perlmutter GPU partition) to generate the QP energy
   of the LUMO.

