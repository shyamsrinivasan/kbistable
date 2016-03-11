function [flux,vflux] = CKinetics(model,pvec,mc,Vind)
[~,nc] = size(mc);
allmc = mc;
S = model.S;
nrxn = model.nt_rxn;
K = pvec.K;
kfwd = pvec.kcat_fwd;
kbkw = pvec.kcat_bkw;
Vmax = pvec.Vmax;

vflux = zeros(nrxn,nc);   
flux = zeros(nrxn,nc);

%eliminate consideration for excess cofators
%pi[c],pi[e],h[c],h[e],h2o[c]
% find(strcmpi(model.mets,'pi[e]'))...
%         find(strcmpi(model.mets,'pi[c]'))...
he = find(strcmpi(model.mets,'h[e]'));
h2o = find(strcmpi(model.mets,'h2o[c]'));
vmet = [he...
        find(strcmpi(model.mets,'h[c]'))...        
        h2o...       
        find(strcmpi(model.mets,'co2[c]'))];


for ic = 1:nc   
    mc = allmc(:,ic);
    for irxn = 1:length(Vind)
        %compensated species indices
    %     sbcmp = zeros(length(model.mets),1);
    %     prcmp = zeros(length(model.mets),1);
        nr_flux = zeros(1,1);

        nmet = size(S,1);

        sbid = S(:,Vind(irxn))<0;
        prid = S(:,Vind(irxn))>0;    
        
        %if any(strcmpi(model.rxn{Vind(irxn)},'PPC')) ||...
             %any(strcmpi(model.rxn{Vind(irxn)},'ATPS4r'))
%              sbid(h2o) = 0;
%              prid(h2o) = 0;          
        if any(strcmpi(model.rxns{Vind(irxn)},'PPC')) ||...
             any(strcmpi(model.rxns{Vind(irxn)},'ATPS4r'))
         
            mc_alls = prod(logical(mc(sbid)));
            mc_allp = prod(logical(mc(prid)));
            sbid(h2o) = 0;
            prid(h2o) = 0;        
            if ~any(model.CMPS(sbid,Vind(irxn)))  
%                 sbid(vmet) = 0;
                cmp_s = [];
            else
                sbid = find(sbid);
                cmp_s = sbid(logical(model.CMPS(sbid,Vind(irxn))));
                sbid = setdiff(sbid,cmp_s);
                sbid = setdiff(sbid,[he h2o]);    
                sbid = logical(sparse(sbid,1,1,nmet,1));
            end
            if ~any(model.CMPS(prid,Vind(irxn)))
%                 prid(vmet) = 0;
                cmp_p = [];
            else
                prid = find(prid);
                cmp_p = prid(logical(model.CMPS(prid,Vind(irxn))));
                prid = setdiff(prid,cmp_p);  
                prid = setdiff(prid,[he h2o]);
                prid = logical(sparse(prid,1,1,nmet,1));
            end      
        else
            if any(sbid)
                mc_alls = prod(logical(mc(sbid)));
                if ~any(model.CMPS(sbid,Vind(irxn)))            
                    sbid(vmet) = 0;
                    cmp_s = [];
                else
                    sbid = find(sbid);
                    cmp_s = sbid(logical(model.CMPS(sbid,Vind(irxn))));
                    sbid = setdiff(sbid,cmp_s);
                    sbid = setdiff(sbid,[he h2o]);
                    sbid = logical(sparse(sbid,1,1,nmet,1));
        %             sbid(sbid==cmp_s) = [];                     
                end
            end

            if any(prid)
                mc_allp = prod(logical(mc(prid)));
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
%         if ~any(strcmpi(model.rxns{Vind(irxn)},'PPC'))||...
%                 ~any(strcmpi(model.rxns{Vind(irxn)},'ATPS4r'))                  
%         else                    
%         end
        if ~isempty(cmp_s)
            cmp_s = prod(mc(cmp_s).*(-model.S(cmp_s,Vind(irxn))));
            if cmp_s > 0
                cmp_s = 1;
            else
                cmp_s = 0;
            end
        else
            cmp_s = 1;
        end
        if ~isempty(cmp_p)
            cmp_p = prod(mc(cmp_p).*(model.S(cmp_p,Vind(irxn))));
            if cmp_p > 0
                cmp_p = 1;
            else
                cmp_p = 0;
            end
        else
            cmp_p = 1;
        end
        
        Sb = -S(sbid,Vind(irxn));
        Sp = S(prid,Vind(irxn));
        smet = model.mets(sbid);
        pmet = model.mets(prid);
        tfshc = strcmpi(smet,'h[c]');
        tfshe = strcmpi(smet,'h[e]');
        tfphc = strcmpi(pmet,'h[c]');
        tfphe = strcmpi(pmet,'h[e]');
        if any(tfshc)||any(tfshe)
            Sb(tfshc) = 1;
            Sb(tfshe) = 1;
        end            
        if any(tfphc) || any(tfphe)
            Sp(tfphc) = 1;
            Sp(tfphe) = 1;
        end
        if model.rev(Vind(irxn))          
            if all(mc(sbid)>0) && all(mc(prid)>0)
                nr_flux = mc_alls*kfwd(Vind(irxn))*cmp_s*prod((mc(sbid)./K(sbid,Vind(irxn))).^Sb) -...
                          mc_allp*kbkw(Vind(irxn))*cmp_p*prod((mc(prid)./K(prid,Vind(irxn))).^Sp);
            elseif all(mc(sbid)>0)
                nr_flux = mc_alls*kfwd(Vind(irxn))*cmp_s*prod((mc(sbid)./K(sbid,Vind(irxn))).^Sb);
            elseif all(mc(prid)>0)
                nr_flux = -mc_allp*kbkw(Vind(irxn))*cmp_p*prod((mc(prid)./K(prid,Vind(irxn))).^Sp);
            end    
            if any(sbid) && any(prid)
                %Denominator - 1.6
                dr_sb = 1+mc(sbid)./K(sbid,Vind(irxn));
                for j = 1:length(find(sbid))
                    for si = 2:Sb(j)
                        dr_sb(j) = dr_sb(j) + dr_sb(j)^si;                
                    end
                end
                %dr_pr
                dr_pr = 1+mc(prid)./K(prid,Vind(irxn));
                for j = 1:length(find(prid))
                    for si = 2:Sp(j)
                        dr_pr(j) = dr_pr(j)+dr_pr(j)^si;
                    end
                end
            else
                dr_sb = 0;
                dr_pr = 0;
            end
        elseif ~model.rev(Vind(irxn))                           
            if all(mc(sbid)>0)
                nr_flux =...
                mc_alls*kfwd(Vind(irxn))*cmp_s*prod((mc(sbid)./K(sbid,Vind(irxn))).^Sb);
            end
            if any(sbid) && any(prid)
                %Denominator - 1.6
                dr_sb = 1+mc(sbid)./K(sbid,Vind(irxn));
                for j = 1:length(find(sbid))
                    for si = 2:Sb(j)
                        dr_sb(j) = dr_sb(j) + dr_sb(j)^si;                
                    end
                end                  
            else
                dr_sb = 0;                
            end
            dr_pr = 0;
        end
        

        if any(sbid) && any(prid)

            %Denominator - 1.6
%             dr_sb = 1+mc(sbid)./K(sbid,Vind(irxn));
%             for j = 1:length(find(sbid))
%                 for si = 2:Sb(j)
%                     dr_sb(j) = dr_sb(j) + dr_sb(j)^si;                
%                 end
%             end
%             %dr_pr
%             dr_pr = 1+mc(prid)./K(prid,Vind(irxn));
%             for j = 1:length(find(prid))
%                 for si = 2:Sp(j)
%                     dr_pr(j) = dr_pr(j)+dr_pr(j)^si;
%                 end
%             end
            dr_flux = prod(dr_sb)+prod(dr_pr)-1;
            vflux(Vind(irxn),ic) = scale_flux(nr_flux/dr_flux);
        else
            vflux(Vind(irxn),ic) = 0;
        end
        flux(Vind(irxn),ic) = Vmax(Vind(irxn))*vflux(Vind(irxn),ic);
    end
end
flux = flux(Vind,:);
vflux = vflux(Vind,:);