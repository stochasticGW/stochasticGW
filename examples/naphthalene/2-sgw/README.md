# naphthalene example, stochastic GW calculation

In this step we will run the stochastic GW calculation, using Heaviside
filtering and gapped filtering to prepare the initial stochastic states.

1. Examine the input file, `INPUT`. In particular, note the values of the
   parameters 'nmctot' (number of Monte Carlo samples), 'gamma' (which
   inversely determines the length of the time propagation), and 'nchb'
   (maximum Chebyshev polynomial order used to fit the filter).

2. Run `./sgw.x` (or on NERSC-Perlmutter, `./sgw_gpu.x` using a job script
   which targets the Perlmutter GPU partition). This initial run uses
   Heaviside filtering and run generates a file, `filtercheb.dat`, which
   contains Chebyshev coefficients of the filter expansion. Run
   `python plotfilter.py filtercheb.dat` to generate plots showing the
   magnitudes of the filter coefficents and the quality of the fit.

3. Increase the value of the 'nchb' parameter in `INPUT` to 1800 and rerun sGW.
   Rerun `plotfilter.py` to regenerate the plot of filter coefficients. How did
   the plot of the Chebyshev fit of the filter change, relative to the previous
   run?

4. Using the larger value for 'nchb', increase the value of 'nmctot' to 1000
   and decrease the value of 'gamma' to 0.06, and rerun `./sgw.x`. The
   resulting quasiparticle energy for the HOMO orbital is the last
   instance of 'MC, E_qp' in the output file.

5. Comment out the label and value for the 'mu' parameter in INPUT and
   decrease the value of 'nchb' to 450 in INPUT and rerun `./sgw.x`.
   This switches the filtering algorithm from Heaviside filtering
   (default) to gapped filtering, which uses the HOMO and LUMO energies
   as filter parameters. Examine `filtercheb.dat`. Is the fit to the filter
   converged with respect to the number of Chebyshev polynomials? Is the
   final QP energy of the HOMO the same as in the Heaviside run?

