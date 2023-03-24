# BN example, stochastic GW calculation

In this step we will run the stochastic GW calculation on BN to compute the
HOMO-LUMO gap using three different levels of stochastic treatment.

1. Examine the input file, `INPUT`, noting the following parameters: 
   'orb_indx' is set to a value of '16' indicating that the HOMO orbital is
   selected; 'ntddft' is set to '-1' indicating that the action of W will be
   computed deterministically; 'projection' is set to 'T' indicating that 
   projection is used to generate the initial stochastic states.

2. First we will compute the HOMO-LUMO gap using deterministic evaluation of W.
   Run `./sgw.x` (or, on an HPC cluster, `./sgw_gpu.x` using a appropriate job
   script) to compute the energy of the HOMO. Then modify `INPUT` by setting
   'orb_indx' to '17' and rerun sGW to compute the LUMO energy.

3. Modify `INPUT` by setting 'ntddft' to '32' to evaluate the action of W
   stochastically instead of deterministically. Repeat the HOMO and LUMO
   calculations from step 2 using this setting.

4. Now we would like to repeat step 3, but using gapped-filtering instead of
   projection to generate the initial states. (Note that gapped filtering uses
   a different input format, so the atomic coordinates and orbitals will be
   read from `sgwinp.txt` instead of from `cnt.ini` and `wf.txt`, respectively).
   Modify `INPUT` by setting 'projection' to 'F' to activate filtering. Then
   run StochasticGW with this setting for 'orb_indx' values of '16' and '17'.

5. Compare the HOMO-LUMO gap computed with the three methods from steps 2-4.
   Which is the most accurate? Note that using deterministic action of W and
   projection are more efficient for small systems with few occupied orbitals;
   for large systems stochastic action of W and filtering show more favorable
   scaling.

