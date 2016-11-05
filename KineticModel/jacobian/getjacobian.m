function [Jxact,lambda,w] = getjacobian(x,pvec,model)
% get analytical jacobian based on modular rate law
nrxn = model.nt_rxn;
nmet = model.nt_metab;
DVX = zeros(nmet,nrxn);

Vind = model.Vind;
Vex = model.Vex;

allmet = [x;model.PM];


[~,~,DVX(:,Vind)] = CKinetics(model,pvec,allmet,Vind,1);
[~,~,DVX(:,Vex)] = TKinetics(model,pvec,allmet,Vex,1);

DVX(allmet==0,:) = 0;
DVX(allmet~=0,:) = DVX(allmet~=0,:)./repmat(allmet(allmet~=0),1,nrxn);
% remove mets in remid from jacobian - they do not control flux
% DVX(model.remid,:) = 0;
S = model.S;

% complete jacobian - build in column loops rowwise
% J = nmet x nmet
Jxact = sparse(S*DVX');
% for imet = 1:nmet
%     Jxact(:,imet) = S*DVX(imet,:)';
% end

% since extracellular concentrations are sometime zero the jacobian for
% internal metabolites alone is also given separately
% nint = model.nint_metab;
% Jint = Jxact(1:nint,1:nint);

% get jacobian using finite difference - this works and produces same
% result as analytical jacobians
% rhsfun = @(irow,mc)Svrow(irow,mc,model,pvec);
% Jnum = model_Jacobian(model,allmet,rhsfun);

% use ADMAT to calculate jacobians - see issue #18 on github
% admatfun = @(x)toyNLAE(x,model,pvec);
% xADMATobj = deriv(x,eye(size(x,1)));
% xADMATres = admatfun(xADMATobj);
% F = getval(xADMATres);
% J = getydot(xADMATres); 

% eigen value and eigen vector
[w,lambda] = eig(full(Jxact));
lambda = diag(lambda);
% 
% [w,lambda] = eig(full(Jint));
% lambda = diag(lambda);

% Jxact = Jint;
% lambda = [];
% w = [];