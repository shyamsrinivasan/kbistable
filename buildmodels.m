function pvec = buildmodels(model,pvec,mc,rxn_add,rxn_excep)
if nargin < 5
    rxn_excep = {};
end
if nargin<4
    rxn_add = {};
end

%reactions to consider for kinetics other than Vind
Vind = addToVind(model,model.Vind,rxn_add,rxn_excep);

% Vind = [model.Vind find(strcmpi(model.rxns,'GLCpts'))];
% Vind = [Vind find(strcmpi(model.rxns,'NADTRHD'))];
% Vind = [Vind find(strcmpi(model.rxns,'THD2'))];
% Vind = [Vind find(strcmpi(model.rxns,'CYTBD'))];
% %         find(strcmpi(model.rxns,'NADH16'))...
% %         find(strcmpi(model.rxns,'ATPS4r'))];


%metabolites that do not affect thermodynamic equilibrium  
he = find(strcmpi(model.mets,'h[e]'));
hc = find(strcmpi(model.mets,'h[c]'));
h2o = find(strcmpi(model.mets,'h2o[c]'));
% pie = find(strcmpi(model.mets,'pi[e]'));
% pic = find(strcmpi(model.mets,'pi[c]'));
co2 = find(strcmpi(model.mets,'co2[c]'));
vmet = [he hc h2o co2];
        
% q8 = find(strcmpi(model.mets,'q8[c]'));
% q8h2 = find(strcmpi(model.mets,'q8h2[c]'));

nrxn = length(Vind);
ntrxn = model.nt_rxn;
ntmet = model.nt_metab;

%backup known parameters from pvec
newp = struct();
newp.K = pvec.K;
newp.Kind = sparse(ntmet,ntrxn);
newp.Vmax = pvec.Vmax;
newp.kfwd = pvec.kcat_fwd;
newp.kbkw = pvec.kcat_bkw;

pvec.Kin = sparse(ntmet,ntrxn);
S = model.S;
check = zeros(ntrxn,1);
flux = zeros(ntrxn,1);

for irxn = 1:nrxn
    %compensated species indices
%     sbcmp = zeros(length(model.mets),1);
%     prcmp = zeros(length(model.mets),1);
    nmet = size(S,1);

    sbid = S(:,Vind(irxn))<0;    
    prid = S(:,Vind(irxn))>0;
    
%     Kscol = zeros(length(find(sbid)),1);
%     Kpcol = zeros(length(find(prid)),1);
    Kscol = zeros(nmet,1);
    Kpcol = zeros(nmet,1);
    
    %no parameters for cofactors - assumed abundant 
    %cofactros are assumed as compensated species
    %hence
    if any(strcmpi(model.rxns{Vind(irxn)},'PPC')) ||...
        any(strcmpi(model.rxns{Vind(irxn)},'ATPS4r'))
        sbid(h2o) = 0;
        prid(h2o) = 0;       
    else
        if any(sbid)
            if ~any(model.CMPS(sbid,Vind(irxn))) 
                sbid(vmet) = 0;
                cmp_s = [];
            else
                sbid = find(sbid);
                cmp_s = sbid(logical(model.CMPS(sbid,Vind(irxn))));
                sbid = setdiff(sbid,cmp_s);
                sbid = setdiff(sbid,[he h2o]);
                sbid = logical(sparse(sbid,1,1,nmet,1));
            end
        end
        if any(prid)
            if ~any(model.CMPS(prid,Vind(irxn)))
                prid(vmet) = 0;
                cmp_p = [];
            else
                prid = find(prid);
                cmp_p = prid(logical(model.CMPS(prid,Vind(irxn))));
                prid = setdiff(prid,cmp_p);
                prid = setdiff(prid,[he h2o]);
                prid = logical(sparse(prid,1,1,nmet,1));
            end
        end
    end
    if ~any(strcmpi(model.rxns{Vind(irxn)},'PPC'))        
    else 
%         if ~any(model.CMPS(sbid,Vind(irxn)))  
%             sbid(vmet) = 0;
%             cmp_s = [];
%         else
%             sbid = find(sbid);
%             cmp_s = sbid(logical(model.CMPS(sbid,Vind(irxn))));
%             sbid = setdiff(sbid,cmp_s);
%             sbid = setdiff(sbid,[he h2o]);    
%             sbid = logical(sparse(sbid,1,1,nmet,1));
%         end
%         if ~any(model.CMPS(prid,Vind(irxn)))
%             prid(vmet) = 0;
%             cmp_p = [];
%         else
%             prid = find(prid);
%             cmp_p = prid(logical(model.CMPS(prid,Vind(irxn))));
%             prid = setdiff(prid,cmp_p);  
%             prid = setdiff(prid,[he h2o]);
%             prid = logical(sparse(prid,1,1,nmet,1));
%         end   
    end
%     if any(sbid)
%         if any(sbid(vmet))
%             sbcmp(vmet(logical(sbid(vmet)))) =...
%             sbid(vmet(logical(sbid(vmet))));    
%             sbid(vmet) = 0;
%         end
%     end
%     
%     if any(prid)
%         if any(prid(vmet))
%             prcmp(vmet(logical(prid(vmet)))) =...
%             prid(vmet(logical(prid(vmet))));
%             prid(vmet) = 0;
%         end
%     end
    pvec = estimateKm(pvec,sbid,prid,mc,Kscol,Kpcol,Vind(irxn));
    
    %forward and backward catalytic rates
    %kfwd and kbkw
    %kfwd or kbkw is sampled basedon reaction directionality from FBA for
    %thermodynamic consistency
    %sampling done only for unknown values
%     fprintf('%s \t delG = %3.6g \t Vflux = %3.6g\t',model.rxns{Vind(irxn)},...
%              pvec.delGr(Vind(irxn)),...
%              model.Vss(Vind(irxn)));    
    pvec = samplekcat(model,pvec,sbid,prid,Vind(irxn),mc);
    
    pvec.Vmax(Vind(irxn)) = 1;
    pvec.Vmax(model.Vss==0) = 0;
    
    %#check for vss and delGr direction    
    flux(Vind(irxn)) = CKinetics(model,pvec,mc,Vind(irxn));
    if pvec.delGr(Vind(irxn)) ~= 0
        if flux(Vind(irxn))*pvec.delGr(Vind(irxn))<0
            check(Vind(irxn)) = 1;
        else
            check(Vind(irxn)) = -1;
        end    
    else
        if flux(Vind(irxn))*pvec.delGr(Vind(irxn))<1e-6
            check(Vind(irxn)) = 1;
        else
            check(Vind(irxn)) = -1;
        end
    end    
end

%other reactions 
%transport reactions x[e] <==> x[c]
Vex = model.Vex;
Vex = setdiff(Vex,Vind); 
for irxn = 1:length(Vex)
    nmet = size(S,1);
    
    sbid = S(:,Vex(irxn))<0;    
    prid = S(:,Vex(irxn))>0;
    
    Kscol = zeros(nmet,1);
    Kpcol = zeros(nmet,1);
    
%     Kscol = zeros(length(find(sbid)),1);
%     Kpcol = zeros(length(find(prid)),1);
    if any(sbid)
        if ~any(model.CMPS(sbid,Vex(irxn))) 
            sbid([he hc]) = 0;                
        else
            sbid = find(sbid);
            cmp_s = sbid(logical(model.CMPS(sbid,Vex(irxn))));
            sbid = setdiff(sbid,cmp_s);
            sbid = logical(sparse(sbid,1,1,nmet,1));
        end
    end
    
    if any(prid)
        if ~any(model.CMPS(prid,Vex(irxn)))
            prid([he hc]) = 0;
        else
            prid = find(prid);
            cmp_p = prid(logical(model.CMPS(prid,Vex(irxn))));
            prid = setdiff(prid,cmp_p);
            prid = logical(sparse(prid,1,1,nmet,1));
%             prid = setdiff(prid,[he h2o]);
        end
    end
        
%     if any(sbid)
%         if ~any(strcmpi(model.rxns{Vex(irxn)},'h2ot'))
%             sbid(h2o) = 0;
%         end
%         if ~any(strcmpi(model.rxns{Vex(irxn)},'PIt2r'))
%             sbid([pie pic]) = 0;
%         end
%         if ~any(strcmpi(model.rxns{Vex(irxn)},'CO2t'))
%             sbid(co2) = 0;
%         end
%         
%         if ~any(model.CMPS(sbid,Vex(irxn))) 
%             sbid([he hc]) = 0;                
%         else
%             sbid = find(sbid);
%             cmp_s = sbid(logical(model.CMPS(sbid,Vex(irxn))));
%             sbid = setdiff(sbid,cmp_s);
%             sbid = logical(sparse(sbid,1,1,nmet,1));
% %             sbid = setdiff(sbid,[he h2o]);
%         end
%     end
%     if any(prid)
%         if ~any(strcmpi(model.rxns{Vex(irxn)},'h2ot'))
%             prid(h2o) = 0;
%         end
%         if ~any(strcmpi(model.rxns{Vex(irxn)},'PIt2r'))
%             prid([pie pic]) = 0;
%         end
%         if ~any(strcmpi(model.rxns{Vex(irxn)},'CO2t'))
%             prid(co2) = 0;
%         end
%         if ~any(model.CMPS(prid,Vex(irxn)))
%             prid([he hc]) = 0;
%         else
%             prid = find(prid);
%             cmp_p = prid(logical(model.CMPS(prid,Vex(irxn))));
%             prid = setdiff(prid,cmp_p);
%             prid = logical(sparse(prid,1,1,nmet,1));
% %             prid = setdiff(prid,[he h2o]);
%         end
%     end    
   pvec = estimateKm(pvec,sbid,prid,mc,Kscol,Kpcol,Vex(irxn));
end

%other reactions - redox balance
% pvec = samplekcatRedox(model,pvec,mc);
% 
% if any(isnan(pvec.kcat_fwd))
%     pvec.kcat_fwd(isnan(pvec.kcat_fwd)) = 1000;
% end
% if any(isnan(pvec.kcat_fwd))
%     pvec.kcat_fwd(isnan(pvec.kcat_fwd)) = 1000;
% end

%set irreversible kcats
for irxn = 1:length(model.rxns)
    if ~model.rev(irxn)
        pvec.kcat_bkw(irxn)=0;
    end
end

%exhcnage reactions
pvec.kcat_fwd(model.VFex) = 0;
pvec.kcat_bkw(model.VFex) = 0;

%restore Vmax from backup
pvec.Vmax = newp.Vmax;

%estimate Vmax
%for Vind
if all(check(Vind)>0)
    pvec.Vmax(pvec.delGr==0) = 0;
%     pvec = findVmax(model,pvec,mc);
    
    %simple vmax = vss/ck
    for irxn = 1:length(Vind)
        [~,ck] = CKinetics(model,pvec,mc,Vind(irxn));
        if ck
            pvec.Vmax(Vind(irxn)) = model.Vss(Vind(irxn))/(3600*ck);
        else
            pvec.Vmax(Vind(irxn)) = 1;
        end
    end
    
    %other reactions
    for irxn = 1:length(Vex)
        [~,tk] = TKinetics(model,pvec,mc,Vex(irxn));
        if tk
            pvec.Vmax(Vex(irxn)) = model.Vss(Vex(irxn))/(3600*tk);
        else
            pvec.Vmax(Vex(irxn)) = 1;
        end
    end  
    
    %atp maintanance
    atp = strcmpi(model.mets,'atp[c]');
    if any(atp) && any(strcmpi(model.rxns,'atpm'))
        pvec.Vmax(strcmpi(model.rxns,'atpm')) =...
        model.Vss(strcmpi(model.rxns,'atpm'))/(mc(atp)/1e-5/(1+mc(atp)/1e-5))/3600;
    end
    %for redox reactions    
%     [~,rk,vred] = RedoxKinetics(model,pvec,mc,flux);
%     pvec = getRKparameter(model,pvec,mc,vred);
%     for irxn = 1:length(vred)
%         if rk(irxn)
%             pvec.Vmax(vred(irxn)) = model.Vss(vred(irxn))/rk(irxn);
%         else
%             pvec.Vmax(vred(irxn)) = 1;
%         end
%     end   
    
    %for trasnport fluxes
%     Vex = model.Vex;
%     Vex = setdiff(Vex,Vind);    
%     pvec = getTKparameter(model,pvec,mc,Vex);
%     pvec.Vmax(Vex) = 1;
%     for irxn = 1:length(Vex)
%         [~,tk] = TKinetics(model,pvec,mc,Vex(irxn));
%         if tk
%             pvec.Vmax(Vex(irxn)) = model.Vss(Vex(irxn))/tk;
%         else
%             pvec.Vmax(Vex(irxn)) = 1;
%         end
%     end   
    
    pvec.Vmax(model.VFex) = 1;    
    pvec.Vmax(model.Vss==0) = 0;
    pvec.feasible = 1;
else
    fprintf('Thermodynamically infeasible parameters\n');
    fprintf('Discontinuing\n');
    pvec.feasible = 0;
    return
end



%check - calculate initial flux
flux = iflux(model,pvec,mc);

    
    
            



        