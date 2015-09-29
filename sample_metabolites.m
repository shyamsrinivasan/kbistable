function mc = sample_metabolites(model,mc)
if nargin < 2
    mc = zeros(model.nt_metab,1);
end
delG = model.delSGr;
R = 0.008314;%kJ/mol.K
T = 298;%K
RT = R*T;
%indices - fluxes
vglc = strcmpi(model.rxns,'exGLC');
vpts = strcmpi(model.rxns,'glcpts');
vg6pd = strcmpi(model.rxns,'G6PDH2r');
vpgi = strcmpi(model.rxns,'pgi');
vpgl = strcmpi(model.rxns,'pgl');
vgnd = strcmpi(model.rxns,'gnd');
vrpe = strcmpi(model.rxns,'rpe');
vrpi = strcmpi(model.rxns,'rpi');
vtkt1 = strcmpi(model.rxns,'tkt1');
vtala = strcmpi(model.rxns,'tala');
vtkt2 = strcmpi(model.rxns,'tkt2');
vpfk = strcmpi(model.rxns,'pfk');
vfbp = strcmpi(model.rxns,'fbp');
vfba = strcmpi(model.rxns,'fba');
vtpi = strcmpi(model.rxns,'tpi');
vgapd = strcmpi(model.rxns,'gapd');
vpgk = strcmpi(model.rxns,'pgk');
vpgm = strcmpi(model.rxns,'pgm');
veno = strcmpi(model.rxns,'eno');
vacont = strcmpi(model.rxns,'aconta');
vcs = strcmpi(model.rxns,'cs');
vicl = strcmpi(model.rxns,'icl');
vicd = strcmpi(model.rxns,'icdhyr');
vakgd = strcmpi(model.rxns,'akgdh');
vsuca = strcmpi(model.rxns,'sucoas');
vsucd = strcmpi(model.rxns,'sucdi');
vfum = strcmpi(model.rxns,'fum');
vmals = strcmpi(model.rxns,'mals');
vme1 = strcmpi(model.rxns,'me1');
vmdh = strcmpi(model.rxns,'mdh');
vpdh = strcmpi(model.rxns,'pdh');
vpfl = strcmpi(model.rxns,'pfl');
vppc = strcmpi(model.rxns,'ppc');
vppck = strcmpi(model.rxns,'ppck');
vpyk = strcmpi(model.rxns,'pyk');
vldh = strcmpi(model.rxns,'ldh_d');
vpta = strcmpi(model.rxns,'ptar');
vack = strcmpi(model.rxns,'ackr');
vacald = strcmpi(model.rxns,'acald');
valcd2x = strcmpi(model.rxns,'alcd2x');
vatpm = strcmpi(model.rxns,'atpm');
vcytb = strcmpi(model.rxns,'cytbd');
vnadh16 = strcmpi(model.rxns,'NADH16');
vnadth = strcmpi(model.rxns,'nadtrhd');

%metabolites
glc = strcmpi(model.mets,'glc[e]');
g6p = strcmpi(model.mets,'g6p[c]');
f6p = strcmpi(model.mets,'f6p[c]');
pgl = strcmpi(model.mets,'6pgl[c]');
pgc = strcmpi(model.mets,'6pgc[c]');
r5p = strcmpi(model.mets,'r5p[c]');
ru5p = strcmpi(model.mets,'ru5p-D[c]');
fdp = strcmpi(model.mets,'fdp[c]');
dhap = strcmpi(model.mets,'dhap[c]');
g3p = strcmpi(model.mets,'g3p[c]');
dpg = strcmpi(model.mets,'13dpg[c]');
x5p = strcmpi(model.mets,'xu5p-D[c]');
e4p = strcmpi(model.mets,'e4p[c]');
s7p = strcmpi(model.mets,'s7p[c]');
pg3 = strcmpi(model.mets,'3pg[c]');
pg2 = strcmpi(model.mets,'2pg[c]');
pep = strcmpi(model.mets,'pep[c]');
pyr = strcmpi(model.mets,'pyr[c]');
accoa = strcmpi(model.mets,'accoa[c]');
coa = strcmpi(model.mets,'coa[c]');
cit = strcmpi(model.mets,'cit[c]');
icit = strcmpi(model.mets,'icit[c]');
glx = strcmpi(model.mets,'glx[c]');
akg = strcmpi(model.mets,'akg[c]');
sucoa = strcmpi(model.mets,'succoa[c]');
suc = strcmpi(model.mets,'succ[c]');
fum = strcmpi(model.mets,'fum[c]');
mal = strcmpi(model.mets,'mal[c]');
oaa = strcmpi(model.mets,'oaa[c]');
form = strcmpi(model.mets,'for[c]');
lac = strcmpi(model.mets,'lac[c]');
actp = strcmpi(model.mets,'actp[c]');
ac = strcmpi(model.mets,'ac[c]');
acald = strcmpi(model.mets,'acald[c]');
etoh = strcmpi(model.mets,'etoh[c]');
pic = strcmpi(model.mets,'pi[c]');
pie = strcmpi(model.mets,'pi[e]');
hc = strcmpi(model.mets,'h[c]');
he = strcmpi(model.mets,'h[e]');
o2 = strcmpi(model.mets,'o2[c]');


mclb = zeros(model.nt_metab,1);
mcub = zeros(model.nt_metab,1);
%assuymed concentrations in moles/L or M
mc(o2) = 3e-3;
mc(pic) = 1e-3;
mc(pie) = 1e-3;
mc(hc) = 1e-7;
mc(he) = 1e-7;
mc(coa) = 4.5e-3;
% mc(coa) = 1e-7;
nad_r = 1e3;%nad/nadh
nadp_r = 1e-3;%nad_r*exp(-delG(vnadth)/RT);%nadp/nadph



%atp/adp
% atp_adplb = max(mc(h)*mc(pic)*exp(delG(vatpm)/RT),...
%                 mc(fdp)*mc(h)/mc(f6p)*exp(delG(vpfk)/RT));
%adp/atp
% adp_atpub = min(1/(mc(h)*mc(pic))*exp(-delG(vatpm)/RT),...
%                 mc(f6p)/(mc(fdp)*mc(h))*exp(-delG(vpfk)/RT)); 
%1/(mc(hc)*mc(pic))*
adp_atplb = 0;
% adp_atpub  = 1/(mc(hc)*mc(pic))*min(exp((-delG(vppc)-delG(vppck))/RT),...
%                                    exp(-delG(vatpm)/RT)); 
adp_atpub  = min(exp((-delG(vppc)-delG(vppck))/RT),...
                                   exp(-delG(vatpm)/RT)); 
pd = makedist('Uniform','lower',adp_atplb,'upper',adp_atpub);
adp_atp = random(pd,1,1);
%q8h2/q8
% q8h2_q8lb = mc(he)/(mc(hc)*mc(o2))*exp(delG(vcytb)/RT);
% q8h2_q8ub = mc(hc)/(mc(he)*nad_r)*exp(-delG(vnadh16)/RT);
q8h2_q8lb = 1/mc(o2)*exp(delG(vcytb)/RT);
q8h2_q8ub = 1/1*nad_r*exp(-delG(vnadh16)/RT);
pd = makedist('Uniform','lower',q8h2_q8lb,'upper',q8h2_q8ub);
q8h2_q8 = random(pd,1,1);
%g6p
% mcub(g6p) = mc(glc)/(adp_atp*mc(hc))*exp((-delG(vpts)+delG(vpyk))/RT);
mcub(g6p) = mc(glc)/adp_atp*exp((-delG(vpts)+delG(vpyk))/RT);
pd = makedist('Uniform','lower',0,'upper',mcub(g6p));
mc(g6p) = random(pd,1,1);
%f6p
mcub(f6p) = mc(g6p)*exp(-delG(vpgi)/RT);
pd = makedist('Uniform','lower',0,'upper',mcub(f6p));
mc(f6p) = random(pd,1,1);
%fdp
% mcub(fdp) = mc(f6p)/(mc(hc)*adp_atp)*exp(-delG(vpfk)/RT);
% mclb(fdp) = mc(f6p)*mc(pic)*exp(delG(vfbp)/RT);
mcub(fdp) = mc(f6p)/adp_atp*exp(-delG(vpfk)/RT);
mclb(fdp) = mc(f6p)*exp(delG(vfbp)/RT);
pd = makedist('Uniform','lower',mclb(fdp),'upper',mcub(fdp));
mc(fdp) = random(pd,1,1);
% pd = makedist('Uniform');
% mc(fdp) = mclb(fdp) + (mcub(fdp)-mclb(fdp))*random(pd,1,1);           
%g3p
mcub(g3p) = sqrt(mc(fdp)*exp((-delG(vfba)-delG(vtpi))/RT));
pd = makedist('Uniform','lower',0,'upper',mcub(g3p));
mc(g3p) = random(pd,1,1);
%dhap
mcub(dhap) = mc(fdp)/mc(g3p)*exp(-delG(vfba)/RT);
mclb(dhap) = mc(g3p)*exp(delG(vtpi)/RT);
pd = makedist('Uniform','lower',mclb(dhap),'upper',mcub(dhap));
mc(dhap) = random(pd,1,1);
%6pgl
% mcub(pgl) = (mc(g6p)*nadp_r/mc(hc))*exp(-delG(vg6pd)/RT);
mcub(pgl) = mc(g6p)*nadp_r*exp(-delG(vg6pd)/RT);
pd = makedist('Uniform','lower',0,'upper',mcub(pgl));
mc(pgl) = random(pd,1,1);
%6pgc 
% mcub(pgc) = (mc(pgl)/mc(hc))*exp(-delG(vpgl)/RT);
mcub(pgc) = mc(pgl)*exp(-delG(vpgl)/RT);
pd = makedist('Uniform','lower',0,'upper',mcub(pgc));
mc(pgc) = random(pd,1,1);
%ru5p
mcub(ru5p) = mc(pgc)*nadp_r*exp(-delG(vgnd)/RT);
pd = makedist('Uniform','lower',0,'upper',mcub(ru5p));
mc(ru5p) = random(pd,1,1);
%r5p and x5p
if mc(ru5p) < 1e-6
    mc(r5p) = 1e-8;
    mc(x5p) = 1e-8;
else
    mcub(r5p) = mc(ru5p)*exp(delG(vrpi)/RT);
    pd = makedist('Uniform','lower',mcub(r5p)/10,'upper',mcub(r5p));
    mc(r5p) = random(pd,1,1);
    %xu5p
    mcub(x5p) = mc(ru5p)*exp(-delG(vrpe)/RT);
    pd = makedist('Uniform','lower',mcub(x5p)/10,'upper',mcub(x5p));
    mc(x5p) = random(pd,1,1);
end
%s7p
if mc(x5p)<1e-6 || mc(r5p)<1e-6
    mc(s7p) = 1e-8;
else
    mcub(s7p) = (mc(x5p)*mc(r5p)/mc(g3p))*exp(-delG(vtkt1)/RT);
    pd = makedist('Uniform','lower',mcub(s7p)/10,'upper',mcub(s7p));
    mc(s7p) = random(pd,1,1);
end
%e4p
if mc(s7p)<1e-6
    mc(e4p) = 1e-8;
else
    mclb(e4p) = (mc(f6p)*mc(g3p)/mc(x5p))*exp(delG(vtkt2)/RT);
    mcub(e4p) = mc(g3p)*mc(s7p)/mc(f6p)*exp(-delG(vtala)/RT);
    pd = makedist('Uniform','lower',mclb(e4p),'upper',mcub(e4p));
    mc(e4p) = random(pd,1,1);
end
%s7p
% mcub(s7p) = (mc(x5p)*mc(r5p)/mc(g3p))*exp(-delG(vtkt1)/RT);
% mclb(s7p) = (mc(e4p)*mc(f6p)/mc(g3p))*exp(delG(vtala)/RT);
% pd = makedist('Uniform','lower',mclb(s7p),'upper',mcub(s7p));
% mc(s7p) = random(pd,1,1);
%13dpg
% mcub(dpg) = (mc(g3p)*mc(pic)*nad_r/mc(hc))*exp(-delG(vgapd)/RT);
mcub(dpg) = (mc(g3p)*nad_r)*exp(-delG(vgapd)/RT);
pd = makedist('Uniform','lower',mcub(dpg)/10,'upper',mcub(dpg));
mc(dpg) = random(pd,1,1);
%3pg
mcub(pg3) = mc(dpg)*adp_atp*exp(-delG(vpgk)/RT);
pd = makedist('Uniform','lower',mcub(pg3)/10,'upper',mcub(pg3));
mc(pg3) = random(pd,1,1);
%2pg
mcub(pg2) = mc(pg3)*exp(-delG(vpgm)/RT);
pd = makedist('Uniform','lower',mcub(pg2)/10,'upper',mcub(pg2));
mc(pg2) = random(pd,1,1);
%pep
mcub(pep) = mc(pg2)*exp(-delG(veno)/RT);
pd = makedist('Uniform','lower',mcub(pep)/10,'upper',mcub(pep));
mc(pep) = random(pd,1,1);
%pyr
% mcub(pyr) = mc(pep)*mc(hc)*exp(-delG(vpyk)/RT)*adp_atp;
mcub(pyr) = mc(pep)*exp(-delG(vpyk)/RT)*adp_atp;
pd = makedist('Uniform','lower',mcub(pyr)/10,'upper',mcub(pyr));
mc(pyr) = random(pd,1,1);
%accoa
mcub(accoa) = mc(pyr)/nad_r*exp(-delG(vpdh)/RT);
pd = makedist('Uniform','lower',mcub(accoa)/10,'upper',mcub(accoa));
mc(accoa) = random(pd,1,1);
%coa
% mclb(coa) = mc(accoa)/(mc(pyr)*nad_r)*exp(delG(vpdh)/RT);
% pd = makedist('Uniform','lower',mclb(coa),'upper',mclb(coa)*10);
% mc(coa) = random(pd,1,1);

%TCA metabolites
%mal
mclb(mal) = mc(pyr)/nad_r*exp(delG(vme1)/RT);
mcub(mal) = mc(pyr)/nad_r*exp(delG(vme1)/RT);
pd = makedist('Uniform');
mc(mal) = mclb(mal) + (mcub(mal)-mclb(mal))*random(pd,1,1);
%oaa
% mcub(oaa) = mc(mal)/mc(hc)*nad_r*exp(delG(vmdh)/RT);
mcub(oaa) = mc(mal)*nad_r*exp(delG(vmdh)/RT);
pd = makedist('Uniform','lower',mcub(oaa)/10,'upper',mcub(oaa));
mc(oaa) = random(pd,1,1);
%cit
% mcub(cit) = mc(oaa)*mc(accoa)/(mc(coa)*mc(hc))*exp(-delG(vcs)/RT);
mcub(cit) = mc(oaa)*mc(accoa)*exp(-delG(vcs)/RT);
pd = makedist('Uniform','lower',mcub(cit)/10,'upper',mcub(cit));
mc(cit) = random(pd,1,1);
%icit
mcub(icit) = mc(cit)*exp(-delG(vacont)/RT);
pd = makedist('Uniform','lower',mcub(icit)/10,'upper',mcub(icit));
mc(icit) = random(pd,1,1);
%akg
mcub(akg) = mc(icit)*nadp_r*exp(-delG(vicd)/RT);
pd = makedist('Uniform','lower',mcub(akg)/10,'upper',mcub(akg));
mc(akg) = random(pd,1,1);
%glx
% mclb(glx) = mc(mal)*mc(coa)*mc(hc)/mc(accoa)*exp(delG(vmals)/RT);
mclb(glx) = mc(mal)/mc(accoa)*exp(delG(vmals)/RT);
mcub(glx) = mc(mal)/mc(accoa)*exp(delG(vmals)/RT);
pd = makedist('Uniform');
mc(glx) = mclb(glx)+(mclb(glx)-mclb(glx))*random(pd,1,1);
%succ
mclb(suc) = mc(icit)/mc(glx)*exp(-delG(vicl)/RT);
mcub(suc) = mc(icit)/mc(glx)*exp(-delG(vicl)/RT);
pd = makedist('Uniform');
mc(suc) = mclb(suc)+(mcub(suc)-mclb(suc))*random(pd,1,1);
%succoa
% mclb(sucoa) = mc(akg)*mc(coa)*nad_r*exp(-delG(vakgd)/RT);
mcub(sucoa) = mc(akg)*nad_r*exp(-delG(vakgd)/RT);
% mcub(sucoa) = mc(suc)*mc(coa)/(mc(pic)*adp_atp)*exp(-delG(vsuca)/RT);
mclb(sucoa) = mc(suc)/adp_atp*exp(-delG(vsuca)/RT);
pd = makedist('Uniform');
mc(sucoa) =  mclb(sucoa)+(mcub(sucoa)-mclb(sucoa))*random(pd,1,1);
%fum
% mcub(fum) = mc(suc)/q8h2_q8*exp(-delG(vsucd)/RT);
mcub(fum) = mc(suc)*exp(-delG(vsucd)/RT);
% mclb(fum) = mc(mal)*exp(delG(vfum)/RT);
pd = makedist('Uniform','lower',mcub(fum)/10,'upper',mcub(fum));
mc(fum) = random(pd,1,1);
%for
% mclb(form) = mc(pyr)*mc(coa)/mc(accoa)*exp(-delG(vpfl)/RT);
mclb(form) = mc(pyr)/mc(accoa)*exp(-delG(vpfl)/RT);
mcub(form) = mc(pyr)/mc(accoa)*exp(-delG(vpfl)/RT);
pd = makedist('Uniform');
mc(form) = mclb(form) + (mcub(form)-mclb(form))*random(pd,1,1);
%lac
% mcub(lac) = mc(pyr)*mc(hc)/nad_r*exp(-delG(vldh)/RT);
mcub(lac) = mc(pyr)/nad_r*exp(-delG(vldh)/RT);
pd = makedist('Uniform','lower',mcub(lac)/10,'upper',mcub(lac));
mc(lac) = random(pd,1,1);
%actp
% mcub(actp) = mc(accoa)*mc(pic)/mc(coa)*exp(-delG(vpta)/RT);
mcub(actp) = mc(accoa)*exp(-delG(vpta)/RT);
pd = makedist('Uniform','lower',mcub(actp)/10,'upper',mcub(actp));
mc(actp) = random(pd,1,1);
%ac
mclb(ac) = mc(actp)*adp_atp*exp(-delG(vack)/RT);
pd = makedist('Uniform','lower',mclb(ac),'upper',mclb(ac)*10);
mc(ac) = random(pd,1,1);
%acald
% mcub(acald) = mc(accoa)*mc(hc)/(mc(coa)*nad_r)*exp(delG(vacald)/RT);
mcub(acald) = mc(accoa)/nad_r*exp(delG(vacald)/RT);
pd = makedist('Uniform','lower',mcub(acald)/10,'upper',mcub(acald));
mc(acald) = random(pd,1,1);
%etoh
% mcub(etoh) = mc(acald)*mc(hc)/nad_r*exp(-delG(valcd2x)/RT);
mcub(etoh) = mc(acald)/nad_r*exp(-delG(valcd2x)/RT);
pd = makedist('Uniform','lower',mcub(etoh)/10,'upper',mcub(etoh));
mc(etoh) = random(pd,1,1);

%cit,icit,akg,succoa,succ,fum,mal,oaa,glx




