function newmodel = setupMetLP_red(model,rxn_add,mc)
if nargin<3
    mc = [];
end
if nargin<2
    rxn_add = {};
end
%% %create new model structure
newmodel.S = model.S;
newmodel.Keq = model.Keq;
newmodel.Vss = model.Vss;
newmodel.rxns = model.rxns;
newmodel.CMPS = model.CMPS;

newmodel = addremoveRxns(newmodel);

%% %remove reactions with zero fluxes
newmodel.S(:,abs(newmodel.Vss)<1e-7) = [];
newmodel.Keq(abs(newmodel.Vss)<1e-7) = [];
newmodel.rxns(abs(newmodel.Vss)<1e-7) = [];
newmodel.CMPS(:,abs(newmodel.Vss)<1e-7) = [];
newmodel.Vss(abs(newmodel.Vss)<1e-7) = [];
newmodel.mets = model.mets;
if ~isempty(mc)
    newmodel.mc = mc;
end

%% %special reactions - reactions where h[c], h[e],pi[c] affect 
%thermodynamic equilibrium 
% vspl =  [];%[find(strcmpi(newmodel.rxns,'CYTBD'))];

%% %metabolites that do not affect thermodynamic equilibrium   
vmet= [];

%% %find all exchnage and transport reactions in newmodel
[Vind,VFex,Vex,bmrxn] = fluxIndex(newmodel);
%reactions to consider for kinetics other than Vind
Vind = addToVind(newmodel,Vind,rxn_add);
% Vind = [Vind vspl find(strcmpi(newmodel.rxns,'GLCpts'))]; 
Vex = setdiff(Vex,Vind);     

%% %remove h[c], h[e],pi[c]
%remove from all reactions except vspl
if ~isempty(vmet)
    newmodel.S(vmet,:) = [];
    newmodel.mets(vmet) = [];
    if ~isempty(mc)
        newmodel.mc(vmet) = [];
    end
end

%% %remove h2o[c]
vh2o = find(strcmpi(newmodel.mets,'h2o[c]'));

newmodel.S(vh2o,:) = [];
newmodel.mets(vh2o) = []; 
if ~isempty(mc)
    newmodel.mc(vh2o) = []; 
end

vhe = find(strcmpi(newmodel.mets,'h[e]'));
vhc = find(strcmpi(newmodel.mets,'h[c]'));   

vatps = find(strcmpi(newmodel.rxns,'ATPS4r'));

vother = setdiff(1:length(newmodel.rxns),vatps);
newmodel.S([vhe vhc],vother) = 0;

newmodel.S(:,[Vex VFex bmrxn]) = [];
newmodel.Keq([Vex VFex bmrxn]) = [];
newmodel.Vss([Vex VFex bmrxn]) = [];
newmodel.rxns([Vex VFex bmrxn]) = [];

%% %remove empty rows in S
newmodel.mets = newmodel.mets(logical(sum(logical(newmodel.S),2)));
newmodel.S = newmodel.S(logical(sum(logical(newmodel.S),2)),:);
if ~isempty(mc)
    newmodel.mc = newmodel.mc(logical(sum(logical(newmodel.S),2)));
end

vatpm = find(strcmpi(newmodel.rxns,'ATPM'));  

newmodel.Vss = newmodel.Vss(setdiff(1:size(newmodel.S,2),vatpm));
newmodel.rxns = newmodel.rxns(setdiff(1:size(newmodel.S,2),vatpm));
newmodel.Keq = newmodel.Keq(setdiff(1:size(newmodel.S,2),vatpm));
newmodel.S = newmodel.S(:,setdiff(1:size(newmodel.S,2),vatpm));

%% %nonzero reactions for thermodynamic analysis
nmet = size(newmodel.S,1);

%bounds
lb = zeros(nmet,1);
lb(lb==0) = log(1e-8);
ub = zeros(nmet,1);
ub(ub==0) = log(3e-2);

%concentrations in M = mole/L, Bennett et al., 2009
lb(strcmpi(newmodel.mets,'glc[e]')) = log(mc(strcmpi(model.mets,'glc[e]')));
lb(strcmpi(newmodel.mets,'pyr[c]')) = log(1e-6);
lb(strcmpi(newmodel.mets,'atp[c]')) = log(1e-6);
lb(strcmpi(newmodel.mets,'adp[c]')) = log(1e-6);
lb(strcmpi(newmodel.mets,'fdp[c]')) = log(1e-5);
lb(strcmpi(newmodel.mets,'dhap[c]')) = log(1e-4);
lb(strcmpi(newmodel.mets,'pi[c]')) = log(5e-4);
lb(strcmpi(newmodel.mets,'g3p[c]')) = log(5e-4);
lb(strcmpi(newmodel.mets,'h[c]')) = log(1.1e-7);
lb(strcmpi(newmodel.mets,'nadh[c]')) = log(1e-6);
% lb(strcmpi(newmodel.mets,'h[e]')) = log(1e-8);

ub(strcmpi(newmodel.mets,'glc[e]')) = log(mc(strcmpi(model.mets,'glc[e]')));
ub(strcmpi(newmodel.mets,'dhap[c]')) = log(3e-4);
ub(strcmpi(newmodel.mets,'h[c]')) = log(3e-7);
ub(strcmpi(newmodel.mets,'h[e]')) = log(3e-4);
ub(strcmpi(newmodel.mets,'pi[c]')) = log(5e-3);
ub(strcmpi(newmodel.mets,'q8h2[c]')) = log(1e-3);
ub(strcmpi(newmodel.mets,'nadh[c]')) = log(1e-3);
ub(strcmpi(newmodel.mets,'nad[c]')) = log(1e-2);

knwn_id = zeros(nmet,1);
for imet = 1:nmet
    if lb(imet) == ub(imet)
        knwn_id(imet)=1;
    end
end

%% %setup original problem
%Ax <=b 
A = newmodel.S';
A_ub = A(:,~logical(knwn_id));
newmodel.A = A_ub;
newmodel.A_kn = A(:,logical(knwn_id));

b_ub = log(newmodel.Keq)-A(:,logical(knwn_id))*lb(logical(knwn_id));
newmodel.b = b_ub;

if ~isempty(mc)
    newmodel.mc = mc(~logical(knwn_id));
end
newmodel.x_kn = lb(logical(knwn_id));%known concentrations
newmodel.lb = lb(~logical(knwn_id));
newmodel.ub = ub(~logical(knwn_id));
newmodel.mets_kn = newmodel.mets(logical(knwn_id));
newmodel.mets = newmodel.mets(~logical(knwn_id));





