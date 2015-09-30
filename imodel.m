function [model,solverP,saveData] = imodel(model,tmax)
if nargin < 2
    tmax = 500000;%s
end
solverP = struct();

%ODE Solver parameters
solverP.RabsTol = 1e-6;
solverP.PabsTol = 1e-6;
solverP.MabsTol = 1e-5;
solverP.RelTol = 1e-4;
solverP.MaxIter = 1000;    
solverP.MaxDataPoints = 200;
solverP.tmax = tmax;%s
solverP.tout = 0.01;

%data file save location/folder
saveData.filename = '';%sprintf('ExptCondition_%d',exptnum);
saveData.dirname =...
'C:\Users\shyam\Documents\Courses\CHE1125Project\Results\KModel';

%recheck growth rate
[gMax,~,flag] = estimateLPgrowth(model);
if flag
    fprintf('Maximum feasible growth rate = %2.3g h-1\n',-gMax);
else
    fprintf('Growth is Infeasible\n');
end
