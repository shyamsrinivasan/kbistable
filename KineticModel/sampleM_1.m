%function mconc = sampleM_1(model)
%Sample metabolites under thermodynamic constraints
%0 < qmA < Keq
%Sample qmA in [0, Keq] the caluclate [M] from qmA
function mconc = sampleM_1(model)
set(0,'RecursionLimit',750);
% n_rxn = model.n_rxn;
n_rxn = length(model.Vind);
Vind = model.Vind;
nmetab = model.nint_metab;
% accpt_mind = zeros(nmetab,1);
%Sample qmA
qmA = (model.Keq(Vind)).*betarnd(1.5,4.5,n_rxn,1);
%Determine metabolite concentration
mconc = zeros(nmetab,1);
mc_ind = zeros(nmetab,1);
%Set known metabolite concentrations
if any(model.MClow == model.MChigh)
    known_conc_indx = (model.MClow == model.MChigh);
    mconc(known_conc_indx) = model.MClow(known_conc_indx);
    mc_ind(known_conc_indx) = 1;
end
for irxn = 1:n_rxn
    subsind = model.S(:,Vind(irxn)) < 0;
    prodind = model.S(:,Vind(irxn)) > 0;    
    %start recursive call
    [mconc,mc_ind] =...
    r_c1(qmA(irxn),mconc,subsind,prodind,Vind(irxn),model,qmA,mc_ind);   
end
return

function [mconc,mc_ind] =...
         r_c1(qmA,mconc,subsind,prodind,irxn,model,allqmA,mc_ind)
if mc_ind(subsind)%if all substrates are available
    s_subs = -model.S(subsind,irxn);
    prod_prod = qmA*prod(mconc(subsind).^s_subs);    
    nprod = length(find(prodind));
    kn_prod = zeros(model.nint_metab,1);
    ukn_prod = zeros(model.nint_metab,1);
    kn_prod(prodind) = mc_ind(prodind);
    ukn_prod(prodind) = ~kn_prod(prodind); 
    kn_prod = logical(kn_prod);
    ukn_prod = logical(ukn_prod);
    skn_prod = model.S(kn_prod,irxn); 
    sukn_prod = model.S(ukn_prod,irxn);
    %if p1 to pn-1 is known
    if any(ukn_prod) && length(find(kn_prod)) == nprod-1
        %calculate pn from prod_prod
        mconc(ukn_prod) =...
        (prod_prod./prod(mconc(kn_prod).^skn_prod))^(1/sukn_prod);
        mc_ind(ukn_prod) = 1;
    elseif any(ukn_prod)
        %sample p1 to pn-1 
        ukn = find(ukn_prod);
        smpl_prod = zeros(model.nint_metab,1);
        smpl_prod(ukn(1:end-1)) = 1;        
        [mconc,mc_ind] = gen_newsample(mconc,smpl_prod,model,mc_ind);
    end
else%if one of the unknown substrates is the product of another reaction
%     nsubs = length(find(subsind));
    kn_sub = zeros(model.nint_metab,1);
    ukn_sub = zeros(model.nint_metab,1);
    kn_sub(subsind) = mc_ind(subsind);
    ukn_sub(subsind) = ~kn_sub(subsind); 
    %kn_sub = logical(kn_sub);
    ukn_sub = logical(ukn_sub);
%     s_sub_kn = -model.S(kn_sub,irxn);
%     s_sub_ukn = -model.S(ukn_sub,irxn);
%     s_prod = model.S(prodind,irxn);
%     [~,rxn_prod] = find(model.S(ukn_sub,1:model.n_rxn)>0);
    [~,rxn_prod] = find(model.S(ukn_sub,model.Vind)>0);
    rxn_prod = model.Vind(model.Vind(rxn_prod) ~= irxn);
%     rxn_prod = rxn_prod(rxn_prod ~= irxn);
    if any(rxn_prod)
        if length(rxn_prod) > 1
            rxn_prod = rxn_prod(1);
        end           
        subs_new = model.S(:,rxn_prod)<0;
        prod_new = model.S(:,rxn_prod)>0;
        [mconc,mc_ind] =...
        r_c1(allqmA(rxn_prod),mconc,subs_new,prod_new,rxn_prod,model,...
             allqmA,mc_ind);
    else
        %sample substrate
        [mconc,mc_ind] = gen_newsample(mconc,ukn_sub,model,mc_ind);
    end
    %calculate prod_prod
    [mconc,mc_ind] = r_c1(qmA,mconc,subsind,prodind,irxn,model,allqmA,mc_ind);
end
return

function [mconc,mc_ind] =...
         r_c(qmA,mconc,subsind,prodind,irxn,model,allqmA,mc_ind)
if mc_ind(subsind)    
    nprod = length(find(prodind));
    kn_prod = zeros(model.nint_metab,1);
    ukn_prod = zeros(model.nint_metab,1);
    kn_prod(prodind) = mc_ind(prodind);
    ukn_prod(prodind) = ~kn_prod(prodind); 
    kn_prod = logical(kn_prod);
    ukn_prod = logical(ukn_prod);
    s_subs = model.S(subsind,irxn);
    s_prod_kn = model.S(kn_prod,irxn);
    s_prod_ukn = model.S(ukn_prod,irxn);
    if ~any(ukn_prod) 
        %all concentrations have already been identified
        return
    end
    if length(find(kn_prod)) == nprod-1
        mconc(ukn_prod) = ((qmA*prod(mconc(subsind).^s_subs))/...
                          (prod(mconc(kn_prod).^s_prod_kn)))^(1/s_prod_ukn);
        mc_ind(ukn_prod) = 1;
    else
        [~,rxn_prod] = find(model.S(ukn_prod,:));
        rxn_prod = rxn_prod(rxn_prod ~= irxn);
        if any(rxn_prod) 
            for irxn1 = 1:length(rxn_prod)
                subsind_new = model.S(:,rxn_prod(irxn1)) < 0;
                prodind_new = model.S(:,rxn_prod(irxn1)) > 0;           
                [mconc,mc_ind] =...
                r_c(allqmA(rxn_prod(irxn1)),mconc,subsind_new,prodind_new,...
                    rxn_prod(irxn1),model,allqmA,mc_ind);
            end
        else
            %choose 1:n-1 products for sampling
            uknprods = find(ukn_prod);
            smpl_prod = zeros(nmetab,1);
            smpl_prod(uknprods(1:end-1)) = 1;
            %sample n-1 product metabolites            
            [mconc,mc_ind] = gen_newsample(mconc,smpl_prod,model,mc_ind);
        end
        %recalculate known & unknown indices
        kn_prod = zeros(model.nint_metab,1);
        ukn_prod = zeros(model.nint_metab,1);
        kn_prod(prodind) = mc_ind(prodind);
        ukn_prod(prodind) = ~kn_prod(prodind);  
        kn_prod = logical(kn_prod);
        ukn_prod = logical(ukn_prod);
        s_prod_kn = model.S(logical(kn_prod),irxn);
        s_prod_ukn = model.S(logical(ukn_prod),irxn);
        %calculate mconc(prod(1)) from 1:n-1 values
        mconc(ukn_prod) = ((qmA*prod(mconc(subsind).^s_subs))/...
                          (prod(mconc(kn_prod).^s_prod_kn)))^(1/s_prod_ukn);
        mc_ind(ukn_prod) = 1;
    end
elseif mc_ind(prodind)%if subs are unknown but prod known
    
    if ~any(ukn_sub)
        return
    end
    if length(find(kn_sub)) == nsubs-1
        mconc(ukn_sub) = (prod(mconc(prodind).^s_prod)/...
                         (qmA*prod(mconc(kn_sub).^s_sub_kn)))^(1/s_sub_ukn);
        mc_ind(ukn_sub) = 1;
    end    
else%if subs are unknown
    [~,rxn_prod] = find(model.S(prodind,:));
    rxn_prod = rxn_prod(rxn_prod ~= irxn);    
%     rxn_ind = zeros(model.n_rxn,1);
%     rxn_ind(irxn) = 1;
    if any(rxn_prod)
        rxn_prod = rxn_prod(1);
%         for irxn1 = 1:length(rxn_prod)
%             if ~rxn_ind(rxn_prod(irxn1))
                subsind_new = model.S(:,rxn_prod) < 0;
                prodind_new = model.S(:,rxn_prod) > 0;
                [mconc,mc_ind] =...
                r_c(allqmA(rxn_prod),mconc,subsind_new,prodind_new,...
                    rxn_prod,model,allqmA,mc_ind);
%             end
%         end
    else
        %sample substrate(s) metabolite
        [mconc,mc_ind] = gen_newsample(mconc,subsind,model,mc_ind);
    end
    %recalculate known & unknown indices
    nprod = length(find(prodind));
    kn_prod = zeros(model.nint_metab,1);
    ukn_prod = zeros(model.nint_metab,1);
    kn_prod(prodind) = mc_ind(prodind);
    ukn_prod(prodind) = ~kn_prod(prodind); 
    kn_prod = logical(kn_prod);
    ukn_prod = logical(ukn_prod);
    s_subs = model.S(subsind,irxn);
    s_prod_kn = model.S(kn_prod,irxn);
    s_prod_ukn = model.S(ukn_prod,irxn);
    if ~any(ukn_prod) 
        %all concentrations have already been identified
        return
    end
    if length(find(kn_prod)) == nprod-1 && all(mconc(subsind)~=0)        
    %calculate mconc(prod(1)) from 1:n-1 values
        mconc(ukn_prod) = ((qmA*prod(mconc(subsind).^s_subs))/...
                          (prod(mconc(kn_prod).^s_prod_kn)))^(1/s_prod_ukn);
        mc_ind(ukn_prod) = 1;
    elseif any(mconc(subsind)==0)
        [mconc,mc_ind] = gen_newsample(mconc,subsind,model,mc_ind);
        [mconc,mc_ind] =...
        r_c(qmA,mconc,subsind,prodind,irxn,model,allqmA,mc_ind);
    end
%     [mconc,mc_ind] = r_c(qmA,mconc,subsind,prodind,irxn,model,allqmA,mc_ind);   
end   
return

function [mconc,mc_ind] = gen_newsample(mconc,subsind,model,mc_ind)
subsind = logical(subsind);
  mconc(subsind) = model.MClow(subsind) +...
        (model.MChigh(subsind) - model.MClow(subsind)).*...
        betarnd(1.5,4.5,length(find(subsind)),1);
  mc_ind(subsind) = 1;
return

function [mconc,mc_ind] =...
          recursive_call(qmA,mconc,subsind,s_subs,prodind,s_prod_n,model,...
                         allqmA,mc_ind)
if mconc(subsind) ~= 0
    if ~mc_ind(prodind)
         mconc(prodind) = (qmA*prod(mconc(subsind).^s_subs)).^(1/s_prod_n);
         mc_ind(prodind) = 1;
    end
else
    rxn_prod = find(model.S(prodind,:));
    if any(rxn_prod)
        %recalculate subsind and prodind for reaction rxn_prod
        subsind_new = model.S(:,rxn_prod) < 0;
        prodind_new = model.S(:,rxn_prod) > 0;
        s_subs_n = model.S(subsind_new,irxn);
        s_prod_n = model.S(prodind_new,irxn);
        [mconc,mc_ind] =...
        recursive_call(allqmA(rxn_prod),mconc,subsind_new,s_subs_n,...
                       prodind_new,s_prod_n,model,allqmA,mc_ind); 
    else
        %sample the metabolite   
        mconc = gen_newsample(mconc,prodind,model);
    end
    mconc(prodind) = qmA*prod(mconc(subsind));
    mc_ind(prodind) = 1;
end
return


