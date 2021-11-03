# hons-thesis-code
Supplementary code for honours thesis.

This repository contains supplementary code for the Honours thesis 'Realising the GKP code using periodic driving'. The files included in this repository are:
* `Square code.ipynb` contains the code listed in Appendix I of the thesis, as well as some examples of using these functions. The functions include:
  * `decode()` decodes a quantum state into the square GKP codespace.
  * `squeezing()` computes the squeezing for an arbitrary state.
  * `SP()` defines the square pulse driving function.
  * `HD()` defines the harmonic driving function.
  * `Floquet_states()` gives the quasienergies, Floquet states, their corresponding squeezings and the decoded Floquet states for a given driving Hamiltonian.
  * `adiabatic_evolution()` evolves an initial state according to the adiabatic protocols defined in Chapter 3 of the thesis and returns a list of the evolved states at each period of evolution. It also has the option to compute the Floquet operator at each period of the evolution, in which case the function returns a list of Floquet operators as well as a list of the evolved states at each period of the evolution.
* `Hex code.ipynb` contains the analagous functions to those in `Square code.ipynb` for the hexagonal GKP code. This code was used to produce the figures in Appendix H of the thesis.
* `Optimised amplitudes.ipynb` contains the code used to find the optimised amplitudes for the harmonic driving function using SciPy's `optimize.minimize()` function.
* `Numerical integration.ipynb` contains the code used to numerically integrate the 1st and 2nd order terms in the Magnus expansion, which was used to produce the plot in Appendix D.
* `Thesis figures.ipynb` contains the code to produce all figures in the thesis using saved data from the `data` folder.
* `miscfuncs.py` contains the miscellaneous functions `displacement_fock()`, `gkp_measure_ops()`, and `fock_to_pos` used to define the GKP stabilisers, the Pauli-bin operators and convert from the Fock to the position basis to plot wave functions, respectively. The contents of this file were written by Dr. Grimsmo and Mackenzie Shaw.
* `pyplotsetup.py` contains code to set parameters for matplotlib figures. 
