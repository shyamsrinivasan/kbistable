function newmodel = setupSlackVariables(model)

newmodel = model;

%setup slack problem
nconstr = length((newmodel.A(:,1)));
%add slack variables to all constraints
A_slack = sparse(1:nconstr,1:nconstr,1,nconstr,nconstr);

%adjust sign of slack variable coeffcient for variable to be >= 0
A_slack = repmat(sign(newmodel.Vss),1,nconstr).*A_slack;

newmodel.A = [newmodel.A A_slack];

lb_slack = zeros(size(newmodel.A,1),1);
lb_slack(lb_slack==0) = 1;%1.3201;%5e-1;
%-------
vpyk = strcmpi(newmodel.rxns,'pyk');
lb_slack(vpyk) = 14;
vmdh = strcmpi(newmodel.rxns,'mdh');
lb_slack(vmdh) = 7;
% vgpd = strcmpi(newmodel.rxns,'g6pdh2r');
% lb_slack(vgpd) = 1;
%-------
ub_slack = zeros(size(newmodel.A,1),1);
ub_slack(ub_slack==0) = Inf;

newmodel.lb = [newmodel.lb;lb_slack];
newmodel.ub = [newmodel.ub;ub_slack];