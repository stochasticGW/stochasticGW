# C60 example, stochastic GW calculation

In this step we will compute the energy of the LUMO orbital using stochastic GW, 
first using the TDDFT effective potential, followed by an RPS.

1. Examine the input file, `INPUT`. The 'flgdyn' parameter should be set to 'T'
   to use the TDDFT potential, and 'useGPU' should be set to 'F' as the core 
   GW code makes MPI calls and is not yet implemented on GPUs.

2. Run `./sgw.x` (or on NERSC-Perlmutter, `./sgw_cpu.x` using a job script
   which targets the Perlmutter CPU partition). Note the final quasi-particle
   energy averaged over 840 samples.

3. Modify `INPUT` by setting 'flgdyn' to 'F' and 'useGPU' to 'T' and rerun
   `./sgw.x` (or on NERSC-Perlmutter, `./sgw_gpu.x` using a job script which
   targets the Perlmutter GPU partition), noting the final quasi-particle 
   energy.

