function sol = IntegrateModel(model,ess_rxn,Vup_struct,ensb,mc,change_pos,change_neg)
%change in initial conditions
if nargin<7
    change_neg = ([]);
end
if nargin < 6
    change_pos = struct([]);
end
%check if concentrations are initialized
if nargin < 5
    %reinitialize concentrations
    imc = zeros(model.nt_metab,1);
%     imc(strcmpi(model.mets,'pep[c]')) = 0.002;
%     imc(strcmpi(model.mets,'pep[c]')) = 1e-5;
else
    imc = mc;
    imc(strcmpi(model.mets,'glc[e]')) = imc(strcmpi(model.mets,'glc[e]'));
%     imc(strcmpi(model.mets,'q8h2[c]')) = 100;
%     imc(strcmpi(model.mets,'pi[c]')) = 100;
end

if nargin<4
    error('getiest:NoA','No parameter vector');
else
    pvec = ensb{1,2};
end
    
if nargin<3
    Vup_struct([]);
end
if nargin<2
    ess_rxn = {};
end
%initialize solver properties
[model,solverP,saveData] = imodel(model,ess_rxn,Vup_struct,1.2);

% model.Vuptake = zeros(model.nt_rxn,1);
% h2o = find(strcmpi(model.rxns,'exH2O'));
% pi =  find(strcmpi(model.rxns,'exPI'));
h = find(strcmpi(model.rxns,'exH'));
% 
model.Vuptake([h]) = [1000];

%noramlize concentration vector to intial state
Nimc = imc./imc;
Nimc(imc==0) = 0;

%intorduce perturbation in initial conditions
% met.glc_e = 10;
% Nimc(strcmpi(model.mets,'glc[e]')) = 1.1;
if ~isempty(change_pos)
    Nimc = changeInitialCondition(model,Nimc,change_pos);
end
if ~isempty(change_neg)
    Nimc = changeInitialCondition(mdoel,Nimc,[],change_neg);
end

model.imc = imc;
model.imc(model.imc==0) = 1;
%calculate initial flux
flux = iflux(model,pvec,Nimc.*imc);
dXdt = ODEmodel(0,Nimc,[],model,pvec);

%integrate model
[sol,finalSS,status] = callODEsolver(model,pvec,Nimc,solverP);

%initialize solver properties
[model,solverP,saveData] = imodel(model,10.0);

%integrate model
[sol,finalSS,status] = callODEsolver(model,pvec,Nimc,solverP);

%initialize solver properties
[model,solverP,saveData] = imodel(model,50.0);

%integrate model
[sol,finalSS,status] = callODEsolver(model,pvec,Nimc,solverP);

%initialize solver properties
[model,solverP,saveData] = imodel(model,100.0);

%integrate model
[sol,finalSS,status] = callODEsolver(model,pvec,Nimc,solverP);

%initialize solver properties
[model,solverP,saveData] = imodel(model,200.0);

%integrate model
[sol,finalSS,status] = callODEsolver(model,pvec,Nimc,solverP);

%initialize solver properties
[model,solverP,saveData] = imodel(model,500.0);

%integrate model
[sol,finalSS,status] = callODEsolver(model,pvec,Nimc,solverP);

%initialize solver properties
[model,solverP,saveData] = imodel(model,1000.0);

%integrate model
[sol,finalSS,status] = callODEsolver(model,pvec,Nimc,solverP);

%initialize solver properties
[model,solverP,saveData] = imodel(model,2000.0);

%integrate model
[sol,finalSS,status] = callODEsolver(model,pvec,Nimc,solverP);

%initialize solver properties
[model,solverP,saveData] = imodel(model,5000.0);

%integrate model
[sol,finalSS,status] = callODEsolver(model,pvec,Nimc,solverP);

%initialize solver properties
[model,solverP,saveData] = imodel(model,10000.0);

%integrate model
[sol,finalSS,status] = callODEsolver(model,pvec,Nimc,solverP);

%initialize solver properties
[model,solverP,saveData] = imodel(model,15000.0);

%integrate model
[sol,finalSS,status] = callODEsolver(model,pvec,Nimc,solverP);

