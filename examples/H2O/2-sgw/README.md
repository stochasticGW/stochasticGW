# H2O example, stochastic GW calculation

In this step we will run the stochastic GW calculation, using the 'orbj'
functionality to compute the matrix elements of the Sigma operator,
`<i|Sigma|i>` and `<i|Sigma|j>`, for i=4 (HOMO) and all j=1 through j=8.

1. Examine the input file, `INPUT`. The 'orb_indx' entry denotes the i-th
   orbital above, and the list of j-orbitals are provided in the entry
   'orbj_indx'.

2. Run `./sgw.x` (or, on an HPC cluster, call `sgw_gpu.x` using an appropriate
   job script). In the output, note that the
   results for `<i|Sigma|i>` are printed followed by those for `<i|Sigma|j>`.
   Do the values of the matrix elements for `<i|Sigma|j>` for j=4 agree with
   those printed for `<i|Sigma|i>`?

