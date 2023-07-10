# Si214 example, stochastic GW calculation

In this step we will run the stochastic GW calculation, using gapped-filtering
to prepare the initial stochastic states.

1. Examine the input file, `INPUT`. In particular, note the values of the
   parameters 'nmctot' (number of Monte Carlo samples), 'gamma' (which
   inversely determines the length of the time propagation), and 'nchb'
   (maximum Chebyshev polynomial order used to fit the filter).

2. Run `./sgw.x`. This initial run generates a file, `filtercheb.dat`, which
   contains Chebyshev coefficients of the filter expansion. Run
   `python plotfilter.py filtercheb.dat` to generate plots showing the
   magnitudes of the filter coefficents and the quality of the fit. Increase
   the value of the 'nchb' parameter in `INPUT` and rerun sGW until the
   oscillations in the filter around the band gap become negligible.

3. Using a converged value for 'nchb', increase the value of 'nmctot' to 1024
   and decrease the value of 'gamma' to 0.06, and rerun `./sgw.x`. The
   resulting quasiparticle energy for orbital 428 (the HOMO) is the last
   instance of 'MC, E_qp' in the output file.

4. Change the value of 'orb_indx' to 429 in INPUT and rerun `./sgw.x` to obtain
   the quasiparticle energy of the LUMO orbital.
