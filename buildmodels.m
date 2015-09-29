function ensb = buildmodels(model,pvec,mc)

%reactions to consider for kinetics other than Vind
Vind = [model.Vind...
        find(strcmpi(model.rxns,'GLCpts'))...
        find(strcmpi(model.rxns,'THD2'))];

Vind = setdiff(Vind,find(strcmpi(model.rxns,'ATPM')));

pic = find(strcmpi(model.mets,'pi[c]'));
pie = find(strcmpi(model.mets,'pi[e]'));
hc = find(strcmpi(model.mets,'h[c]'));
he = find(strcmpi(model.mets,'h[e]'));
q8 = find(strcmpi(model.mets,'q8[c]'));
q8h2 = find(strcmpi(model.mets,'q8h2[c]'));
h2o = find(strcmpi(model.mets,'h2o[c]'));

nrxn = length(Vind);
ntrxn = model.nt_rxn;
ntmet = model.nt_metab;

%backup known parameters from pvec
newp.K = pvec.K;
newp.Kind = sparse(ntmet,ntrxn);
newp.Vmax = pvec.Vmax;
newp.kfwd = pvec.kfwd;
newp.kbkw = pvec.bkw;

S = model.S;

for irxn = 1:nrxn
    sbid = S(:,Vind(irxn))<0;    
    prid = S(:,Vind(irxn))>0;
    
    Kscol = zeros(length(find(sbid)),1);
    Kpcol = zeros(length(find(prid)),1);
    
    %no parameters for cofactors - assumed abundant 
    sbid([pic pie hc he q8 q8h2 h2o]) = 0;
    prid([pic pie hc he q8 q8h2 h2o]) = 0;
    
    nsb = length(find(sbid));
    npr = length(find(prid));
    
    %enzyme saturation
    sigma = random(makedist('Uniform'),...
                   nsb+npr,...
                   1);
    
    %determine kientic parameters from concentration and sigma
    %substrates
    if ~ismepty(find(sbid,1))
        sb_rat = sigma(1:nsb)./(1-sigma(1:nsb));
        Ksb = mc(logical(sbid))./sb_rat;
        Kscol(logical(sbid),1) = pvec.K(logical(sbid),Vind(irxn));
        if any(Kscol==1)
            pvec.K(Kscol==1,Vind(irxn))=Ksb(pvec.K(logical(sbid),Vind(irxn))==1);
            pvec.Kind(Kscol==1,Vind(irxn))=1;
        end
    end
    
    %products
    if ~isempty(find(prid,1))
        pr_rat = sigma(nsb+1:nsb+npr)./(1-sigma(nsb+1:nsb+npr));
        Kpr = mc(logical(prid))./pr_rat;
        Kpcol(logical(prid),1) = pvec.K(logical(prid),Vind(irxn));
        if any(Kpcol==1)
            pvec.K(Kpcol==1,Vind(irxn))=Kpr(pvec.K(logical(prid),Vind(irxn))==1);
            pvec.Kind(Kpcol==1,Vind(irxn))=1;
        end
    end
    
    %forward and backward catalytic rates
    %kfwd and kbkw
    %kfwd = kcat
    %kbkw is sampled basedon reaction directionality from FBA for
    %thermodynamic consistency
    
    
            



        