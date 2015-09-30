clc

addpath(genpath('C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel'));
rxfname = 'C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel\ecoliN1.txt';

%create model structure
[FBAmodel,parameter,variable,nrxn,nmetab] = modelgen(rxfname);

%assign initial fluxes and calculate FBA fluxes for direction
FBAmodel = FBAfluxes(FBAmodel);

%sample initial metabolite concentrations for estimating kientic parameters
[mc,parameter] = parallel_sampling(FBAmodel,parameter);

%estimate kinetic parameters in an ensemble
ensb = parallel_ensemble(FBAmodel,parameter,mc);

%solve ODE of model to steady state
sol = IntegrateModel(FBAmodel,ensb,mc);
