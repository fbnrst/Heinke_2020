# Heinke_2020
Diploid hepatocytes drive physiological liver renewal in adult humans (https://doi.org/10.1101/2020.08.07.230086)

The model names for the raw files and in the Python code are slightly different:
| papaer | files and code |
| ----------- | ----------- |
|POP1 | POP1spline|
|POP2p | Rmspline|
|POP2p_stem | POP2p_stemspline|
|POP3p | POP3p_2x2n_NPspline|


The python package used for the simulation of the C14 concentration and the parameter estimation can be found here: 
https://github.com/rodjul42/pyC14

The cell based simulation is in the folder cellsim and is writen in C. It needs SWIG,Python and GSL to build.

The notebooks in this repo are for the corner plots and compression of the data and simulation.
The data for the corner plots are here: http://dx.doi.org/10.5281/zenodo.5938661
