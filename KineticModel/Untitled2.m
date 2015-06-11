%function [Vpts] = uptake_flux(Vuptake)

%in model formation
%check if any protein name is identical to metabolite name


%catabolite repression by PTS
%Metabolites external
Cglcx
Clacx
%Metabolites internal
Cglc
Clac
Cpep
Cpyr
%Enzymes
EI
EIIA
HPr
P
P_EI
P_EIIA
P_HPr
%Proteins/genes

%crp, cAMP system - gene+protein
crp %basal+activated by crp.cAMP + activated by (crp.cAMP)^2
%crp = alpha*RNAP*(1/b+freg)
%freg = crp.cAMP/(1+CRP.cAMP)
cya %basal 
%cya = alpha*RNAP*1/b
CrpA %protein from crp
cAMP %protein from cya
%crp.cAMP = k1*CRP*cAMP - k2*crp.cAMP

%PTS genes+proteins
ptsG %produces enzyme EIICB - will use existing TRN model
%ptsG = alpha*RNAP*(1/b+freg)
%freg = CRP.cAMP/(1+Mlc+CRP.cAMP)
ptsHI %operon produces enzyme EI, HPr, EIIA (minimal induction due to regulation)
%EI = beta*alpha*RNAP*(1/b+freg)
%freg = CRP.cAMP/(1+Mlc+CRP.cAMP)
%HPr = beta*alpha*RNAP*(1/b)
%above operon/gene regulated by Mlc and CRP.cAMP
%should include more CRP.cAMP mechanics

%PTS flux
Vpts %kinetic rate laws for PTS using enzyme & metabolite concentrations
%Vpts = f(CEIICB,CP_EIIA,Cglcx);

%uptake flux
dGlcdt = Vpts; %obtained from above

