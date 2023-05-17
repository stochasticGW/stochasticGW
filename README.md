Open-source stochastic GW software. To install: clone the GitHub repository and follow the installation instructions in src/README.
The manual describes the installation and contains a tutorial.

Version Notes:

v3: GPU support (NVIDIA GPUs) has been added. Gapped-filtering method has been added for stochastic-GW on large systems. The wrapper to extract Quantum ESPRESSO wavefunctions has been improved to allow one to use QE to prepare orbitals for gapped filtering.

v2: Added the ability to extract wavefunctions from Quantum ESPRESSO output to be used as input for sGW. Version 2 also includes implementation of LDA and GGA functionals via the LIBXC library.
