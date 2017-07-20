% get continuation diagrams on acetate for different fixed v4 enzyme parameters
% also get continuation diagrams on v4 parameters for
% different fixed acetate concentrations
% do a codimension-2 bifurcation analysis with acetate and these other enzyme
% parameters as the second parameter
runKotte
close all
figh = bifurcationPlot(data.x1,data.s1,data.f1,[4,1]);

% continue on acetate for different values of V4max
v4max = 0:.05:.35;
pvec(ap) = .01;
model.PM(ac-length(xeq)) = .01;
iguess = [.1 .1 .1];

npar = length(v4max);
nbifpts = zeros(npar,1);
nbifpts(1:4) = 2000;
nbifpts(4:6) = 10000;
nbifpts(6:end) = 30000;
allpvec = repmat(pvec,npar,1);
allpvec(:,11) = v4max;
tspan = 0:0.1:500;
[~,xeq1,~,feq1] = solveODEonly(npar,iguess',model,allpvec,opts,tspan);

for ip = 1:npar
    pvec = allpvec(ip,:);
    % run MATCONT
    [data,y,p] =...
    execMATCONT(@KotteMATCONT,@Kottecont_fluxcalc,...
                xeq1(:,ip),allpvec(ip,:),ap,fluxg,model,nbifpts(ip));
    s.(['pt' num2str(ip)]) = data;
%     % get the mss for y and p
%     if ~isempty(data)
%         [yss,iyval,fyval] = parseMATCONTresult(data.s1,y);
%         [pss,ipval,fpval] = parseMATCONTresult(data.s1,p);
%         [fss,ifval,ffval] = parseMATCONTresult(data.s1,data.flux);
%     end
end

%% check which solutions have mss
mssid = [];
nss = zeros(npts,1);
for ip = 1:npar
    if ~isempty(s.(['pt' num2str(ip)]))
        s1 = s.(['pt' num2str(ip)]).s1;
        nLP = size(s1,1)-2; % [initial lp1 lp2 final] 
        if nLP > 0
            fprintf('Vector %d has %d Limit Points\n',ip,nLP);
            mssid = union(mssid,ip);
            nss(ip) = nLP;
        end
    else
        fprintf('No convergence at %d\n',ip);
    end
end

%% print bifurcation diagrams for all mss pts
figh = figure;
for ip = 1:npar
    if ismember(ip,mssid)
        % draw bifurcation diagram
        str_data = s.(['pt' num2str(ip)]);
        figh = bifurcationPlot(str_data.x1,str_data.s1,str_data.f1,[4,1],[],[],figh);
    end
end


% find equilibrium solution










%% continue form liits points obtained on eq continuation
% apLP = [9 11];
% [LPdata_f,LPdata_b] =...
% execLPcont(@KotteMATCONT,xeq,pvec',apLP,ap,data,1000);
% 
% for id = 1:2
%     bifurcationPlot(LPdata_f{id}.x1,LPdata_f{id}.s1,LPdata_f{id}.f1,[4,1],[],1,figh);
%     bifurcationPlot(LPdata_b{id}.x1,LPdata_b{id}.s1,LPdata_b{id}.f1,[4,1],[],1,figh);
% end







