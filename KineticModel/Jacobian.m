%Calculate Jacobian using a standalone function
%supercedes calcJacobian.m
function J = Jacobian(model,X)
nrxn = model.nrxn;
nmetab = model.nmetab;

jmetab = 0;
kreg = 0;

X = [flux;metabolite;enzyme];
flux = X(1:nrxn);
metab = X(nrxn+1:nrxn+nmetab);
enz = X(nrxn+nmetab+1:end);

[mindx,rindx] = find(model.S);
Jtherm = sparse(rindx,mindx,ones(length(rindx),1),nrxn,nmetab);
% relJtherm = sparse(rindx,mindx,ones(length(rindx),1),nrxn,nmetab);
Jsat = sparse(rindx,mindx,ones(length(rindx),1),nrxn,nmetab);
% relJsat = sparse(rindx,mindx,ones(length(rindx),1),nrxn,nmetab);

[regind,regrind] = find(model.SI);
Jreg = sparse(regrind,regind,ones(length(regrind),1),nrxn,nmetab);
% relJreg = sparse(regrind,regind,ones(length(regrind),1),nrxn,nmetab);

for irxn = 1:nrxn
% Initialization
    %indices for each reaction j
    subsind = model.S(:,irxn) < 0;%substrate
    prodind = model.S(:,irxn) > 0;%product
    actind = model.SI(:,irxn) > 0;%activator
    inhibind = model.SI(:,irxn) < 0;%inhibitor
    
    nsubs = length(find(subsind));
    nprod = length(find(prodind));
    nact = length(find(actind));
    ninb = length(find(inhibind));  

    %Parameters - K & KI are vectors 
    %Stoichiometric Coefficients
    %K = [Ksubs; Kprod];
    if nsubs ~= 0
        subs = metab(subsind);%Concentrations
        %elements in S are -ve for substrates 
        s_subs = -(model.S(subsind, irxn)); 
        Ksubs = model.K(jmetab + 1:jmetab + nsubs);              
    else
        Ksubs = [];
        s_subs = [];
        subs = [];
    end    
    if nprod ~= 0
        pruds = metab(prodind);
        s_prod = model.S(prodind, irxn);
        Kprod = model.K(jmetab + nsubs + 1:jmetab + nsubs + nprod);        
    else
        Kprod = [];
        s_prod = [];
        pruds = [];
    end
    jmetab = jmetab + nsubs + nprod;
    
    %KI = [KIact; KIinb]; 
    if nact ~= 0
        act = metab(actind);
        s_act = model.SI(actind, irxn);
        KIact = model.KI(kreg + 1:kreg + nact);        
    else
        KIact = [];
        s_act = [];
        act = [];
    end
    if ninb ~= 0
        inhib = metab(inhibind);
        %elements in SI are -ve for inhbition
        s_inhib = -(model.SI(inhibind, irxn));
        KIinb = model.KI(kreg + nact + 1:kreg + nact + ninb); 
    else
        KIinb = [];
        s_inhib = [];
        inhib = [];
    end
    kreg = kreg + nact + ninb;    


    if ~isempty(subs) && ~isempty(s_subs)
        gamma_subs = prod(subs.^s_subs);
        numrsat = prod(subs.^s_subs);
        drsatsubs = prod((1 + subs./Ksubs).^s_subs);
    else
        gamma_subs = [];
        numrsat = 1;
        drsatsubs = 1;
    end
    if ~isempty(pruds) && ~isempty(s_prod)
        gamma_prod = prod(pruds.^s_prod);
        drsatprod = prod((1 + pruds./Kprod).^s_prod);
    else
        gamma_prod = [];
        drsatprod = 1;
    end  
    
    % Thermodynamic Contribution
    gamma = gamma_subs/gamma_prod;
    if ~isempty(gamma(irxn))
        if model.Keq(irxn) ~= 0%avoid divide by zero error            
            gamma = gamma/model.Keq(irxn);
        end           
    end
    subsvec = (s_subs.*(subs.^(s_subs - 1)))./(subs.^s_subs);
    Jtherm_subs = gamma*subsvec;
%     relJtherm_subs = (gamma/(1-gamma))*(subs.*subsvec);
    prodvec = (s_prod.*(pruds.^(s_prod - 1)))./(pruds.^s_prod);
    Jtherm_prod = - gamma*prodvec;
%     relJtherm_prod = (gamma/(1-gamma))*(pruds.*prodvec);
    Jthermrx = [Jtherm_subs' Jtherm_prod'];
%     relJtherm = [relJtherm_subs' relJtherm_prod'];
      
    % Regulatory Contribution
    if ~isempty(act) && ~isempty(s_act)
        vact = prod(((act./KIact).^s_act)./(1 + (act./KIact).^s_act));
    else
        vact = 1;
    end
    
    if ~isempty(inhib) && ~isempty(s_inhib)
        vinhib = prod(1./(1 + (inhib./KIinb).^s_inhib));
    else
        vinhib = 1;
    end
    act_ratio = (act.^s_act)./(Kact.^s_act);
    Jreg_act = ((s_act./act).*(1-act_ratio./(1+act_ratio)));    
    Jreg_inhib = -(s_inhib.*(inhib.^s_inhib-1)./(Kinb.^s_inhib));    
    Jregrx = vact*vinhib*[Jreg_act' Jreg_inhib'];
%     relJreg = [(act.*Jreg_act)' (inhib.*Jreg_inhib)'];
  
    %Saturation Component
    subsvec1 = s_subs.*(subs.^(s_subs - 1))./(subs.^s_subs);
    subsvec2 = (s_subs./Ksubs).*(((1 + subs./Ksubs).^(s_subs-1))./((1 + subs./Ksubs).^s_subs));
    subsscalar = drsatsubs*drsatprod/(drsatsubs*drsatprod - 1);
    Jsat_subs = (numrsat/(drsatsubs*drsatprod - 1))*(subsvec1 - subsscalar*subsvec2);
%     relJsat_subs = subs.*(subsvec1 - subsscalar*subsvec2);
    prodscalar = drsatsubs*drsatprod;
    prodvector = (s_prod./Kprod).*(((1 + pruds./Kprod).^(s_prod-1))./((1 + pruds./Kprod).^s_prod));
    Jsat_prod = -(numrsat/(drsatsubs*drsatprod-1)^2)*prodscalar*prodvector;    
%     relJsat_prod = -(pruds/(drsatsubs*drsatprod - 1)).*(prodscalar*prodvector);
    Jsatrx = [Jsat_subs' Jsat_prod'];  
% relJsat = [relJsat_subs' relJsat_prod'];%relative jacobians - dlnv/dlnS

    Jtherm(irxn,subsind) = Jthermrx(1:nsubs);
    %relJtherm(irxn,subsind) = relJtherm(1:nsubs);
    Jtherm(irxn,prodind) = Jthermrx(nsubs+1:nsubs+nprod);
    %relJtherm(irxn,prodind) = relJthermrx(nsubs+1:nsubs+nprod);
    
    Jsat(irxn,subsind) = Jsatrx(1:nsubs);
    %relJsat(irxn,subsind) = relJsatrx(1:nsubs);
    Jsat(irxn,prodind) = Jsatrx(nsubs+1:nsubs+nprod);
    %relJsat(irxn,prodind) = relJsatrx(nsubs+1:nsubs+nprod);

    
     
    if ~isempty(act)
        Jreg(irxn,actind) = Jregrx(1:nact); 
%         relJreg(irxn,actind) = relJregrx(1:nact);
    end
    if ~isempty(inhib)
        Jreg(irxn,inhibind) = Jregrx(nact+1:nact+ninb);  
%         relJreg(irxn,inhibind) = Jregrx(nact+1:nact+ninb);  
    end   

end
J = Jsat + Jtherm + Jreg;
% relJ = relJsat + relJtherm + relJreg;
%relJ = relJsat + relJtherm + relJreg;
