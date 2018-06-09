# Readme
This repo conatins code necessary to reproduce the figures and analysis presented in the paper

Srinivasan, S., Cluett, W. R., and Mahadevan, R., Model-based Design of Bistable Cell Factories for Metabolic Engineering, Bioinformatics, 2017

All necessary scripts for generating the figures in the figureDrivers folder.
All other folders/files are necessary to generate most of these figures.
Feel free to mix and match functions used for generating different types of figures.

# Dependencies
The code was written to work in MATLAB 2014a

1. MATCONT_CL Command Line version of the Dynamical Systems Analysis Toolbox for MATLAB 
https://sourceforge.net/projects/matcont/

2. ADMAT Tool for algorithmic differentiation
http://www.cayugaresearch.com/admat.html

ADMAT can be obtained free of charge for upto one year for academic customers from Cayuga Research, Waterloo, ON, Canada

3. IBM ILOG CPLEX is also required to generate production boundary figures in the Supplementary Information

# Notes 
1. In case you want to use another LP solver for the FBA problem, make sure it does not conflict with ADMAT
2. OPTI suite of solvers (https://github.com/jonathancurrie/OPTI) conflict with ADMAT and cannot be used together
3. CasADi (https://github.com/casadi/casadi/wiki) is another alternative to ADMAT for generating gradients, Jacobians and Hessiansfor various MATLAB and Python functions. This has however not been implemented in this repository.
