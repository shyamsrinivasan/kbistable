% continue on v4max for different acetate concentrations
runKotte
close all
figh = bifurcationPlot(data.x1,data.s1,data.f1,[4,1]);

% continue on acetate for different values of V4max
acetate = 0.01:.5:5;
% v4max = 0.8;
iguess = [.1 .1 .1];
tspan = 0:0.1:500;
npar = length(acetate);
ap = 11;

for ip = 1:npar
    pvec(ap) = 1;
    pvec(9) = acetate(ip);
    model.PM(ac-length(xeq)) = acetate(ip);
    [~,xeq1,~,feq1] = solveODEonly(1,iguess',model,pvec,opts,tspan);
    % run MATCONT
    [data,y,p] =...
    execMATCONT(@KotteMATCONT,@Kottecont_fluxcalc,...
                xeq1,pvec,ap,fluxg,model,800);
    s.(['pt' num2str(ip)]) = data;
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

%% 3-d bifurcation diagram
% acetate vs v4max vs pep/v4
recalcdata = struct();
figure
hold on
for ip = 1:npar    
    if ismember(ip,mssid)        
        s1 = s.(['pt' num2str(ip)]);
        % get only positive parameters
%         posid = s1.x1(end,:)>=0;
%         posend = find(posid,1,'last');
        % get all parameters
        npts = size(s1.x1,2);
        allpvec = repmat(pvec,npts,1);
        allpvec(:,9) = acetate(ip);        
        allpvec(:,ap) = s1.x1(end,:);
        alleqpts = s1.x1(1:3,:);
        recalcdata.x1 = alleqpts;
        recalcdata.flux = s1.flux;
        % recalculate stability info
        for ipos = 1:npts
            recalcdata.f1(:,ipos) =...
            stabilityInfo(@Kotte_givenNLAE,alleqpts(:,ipos)',model,allpvec(ipos,:));
        end
        % draw 3-d bifurcation plot
        bifurcationPlot([alleqpts;allpvec(:,9);allpvec(:,ap)],);
    end
end

