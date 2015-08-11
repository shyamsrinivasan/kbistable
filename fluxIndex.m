function [Vind,Vuptake,VFup,VFex,Vex,bmrxn,Vup,Vdn] = fluxIndex(model,nt_rxn,newS)
%external mets
exind = ~cellfun('isempty',regexp(model.mets,'\w(?:\[e\])$'));

[~,allrxns] = find(newS(:,1:nt_rxn));
all_unqrxns = unique(allrxns);
nmetab_allrxns = histc(allrxns,all_unqrxns);
ex_rxn = all_unqrxns(nmetab_allrxns == 1);%one sided rxns

%FBA style Reactions A[e/c] <==> | <==> A[e/c]
%one sided rxns
[~,rxn] = find(newS(:,ex_rxn)<0);
VFex = ex_rxn(rxn);%Excrretion from the system
[~,rxn] = find(newS(:,ex_rxn)>0);
VFup = ex_rxn(rxn);%Uptake into the system
try
    VFext = [VFup VFex];
catch
    VFext = [VFup;VFex];
end

%Uptake Reactions: A[e] ---> A[c] | A[e] ---> B[c] | A[e] + B[c] ---> A[c]
%+ D[c]
noex_rxn = setdiff(1:nt_rxn,ex_rxn);
[~,rxn] = find(newS(exind,noex_rxn)<0);
Vuptake = noex_rxn(rxn);

%Matching VFup with Vuptake
Vup = [];
for ivf = 1:length(VFup)
    nmet = newS(:,VFup(ivf))>0;
    [~,rxn] = find(newS(nmet,noex_rxn)<0);
    Vup = union(Vup,noex_rxn(rxn));
end

%Excrertion Reaction: A[c] ---> A[e] | A[c] ---> B[e] | A[c] + B[c] --->
%A[e] + D[c]
[~,rxn] = find(newS(exind,noex_rxn)>0);
Vex = noex_rxn(rxn);
%Remove Vex from considering A[c] + B[e] ---> C[c] + D[c] type reactions


%matchinf VFex with Vex
Vdn = [];
for ive = 1:length(VFex)
    nmet = newS(:,VFex(ive))<0;
    [~,rxn] = find(newS(nmet,noex_rxn)>0);
    Vdn = d_uUnion(Vdn,noex_rxn(rxn));
end

%Identify biomass reaction
[~,bmrxn] = find(newS(strcmpi(model.mets,'Biomass'),:) > 0);
try
    Vind = setdiff(1:nt_rxn,[Vuptake;bmrxn;VFext;Vex']);%intracellular rxns
catch
    Vind = setdiff(1:nt_rxn,[Vuptake bmrxn VFext' Vex]);%intracellular rxns
end
end