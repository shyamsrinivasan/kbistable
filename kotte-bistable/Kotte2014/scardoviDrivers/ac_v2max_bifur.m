% continue on v2max for different acetate concentrations
runKotte
close all
figh = bifurcationPlot(data.x1,data.s1,data.f1,[4,1]);

% continue on acetate for different values of V4max
acetate = 0.01:.5:5;
% v4max = 0.8;
iguess = [.1 .1 .1];
tspan = 0:0.1:500;
npar = length(acetate);
ap = 14;

for ip = 1:npar
    pvec(ap) = 2;
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
figh2 = figure;
for ip = 1:npar
    if ismember(ip,mssid)
        % draw bifurcation diagram
        str_data = s.(['pt' num2str(ip)]);
        figh = bifurcationPlot(str_data.x1,str_data.s1,str_data.f1,[4,1],[],[],figh);
        figh2 = bifurcationPlot([str_data.flux;str_data.x1(end,:)],...
                                str_data.s1,str_data.f1,[6,5],[],[],figh2);
    end
end
set(0,'CurrentFigure',figh);
line([0 0],[0 6],'LineStyle','--','Color',[0 0 1]);
line([0 1],[0 0],'LineStyle','--','Color',[0 0 1]);
set(0,'CurrentFigure',figh2);
line([0 0],[0 1],'LineStyle','--','Color',[0 0 1]);
line([0 1],[0 0],'LineStyle','--','Color',[0 0 1]);

%% 3-d bifurcation diagram - only positive parameter values
% acetate vs v4max vs pep/v4
recalcdata = struct();
hfig1 = figure;
hfig2 = figure;
hfig3 = figure;
hold on
for ip = 1:npar    
    if ismember(ip,mssid)        
        s1 = s.(['pt' num2str(ip)]);
        % get only positive parameters
        posid = s1.x1(end,:)>=0;
%         posend = find(posid,1,'last');
        % get all parameters
        npts = size(s1.x1,2);
        allpvec = repmat(pvec,npts,1);
        allpvec(:,9) = acetate(ip);   
        model.PM(ac-length(xeq)) = acetate(ip);
        allpvec(:,ap) = s1.x1(end,:);
        alleqpts = s1.x1(1:3,:);
        allflux = s1.flux(1:5,:);
        recalcdata.x1 = alleqpts;
        recalcdata.flux = s1.flux;
        % recalculate stability info
        for ipos = 1:npts
            recalcdata.f1(:,ipos) =...
            stabilityInfo(@Kotte_givenNLAE,alleqpts(:,ipos)',model,allpvec(ipos,:));
        end
        % draw 3-d bifurcation plot
%         bifurcationPlot([alleqpts;allpvec(:,9)';allpvec(:,ap)'],s1.s1,recalcdata.f1,[4 5 1],[],1,hfig);
        bifurcationPlot([alleqpts;allpvec(:,9)';allpvec(:,ap)'],s1.s1,recalcdata.f1,[4 5 1],[],1,hfig1);
        xlabel('acetate a.u.');
        ylabel('V2max a.u.');
        zlabel('pep a.u.');
        view([116 22]);
        grid on
        bifurcationPlot([allflux;allpvec(:,9)';allpvec(:,ap)'],s1.s1,recalcdata.f1,[6 7 5],[],1,hfig2);
        xlabel('acetate a.u.');
        ylabel('V2max a.u.');
        zlabel('v4 a.u.');
        view([116 22]);
        grid on
%         set(0,'CurrentFigure',hfig3);
%         set(gca,'NextPlot','add');
        bifurcationPlot([alleqpts;allflux(5,:);allpvec(:,ap)'],s1.s1,recalcdata.f1,[5 1 4],[],1,hfig3);
%         line(allpvec(posid,ap)',alleqpts(1,posid),allflux(5,posid));
        xlabel('V2max a.u.');
        ylabel('pep a.u.');
        zlabel('v4 a.u.');
        view([116 22]);
        grid on
    end
end
% axis([0 2.5 0 1.6 0 1.6]);

%% save figures
% gcf
% fname = 'C:\Users\shyam\Documents\Courses\CHE1125Project\Results\kottedesign\v2max_cont_curves\pep_v4_v2max_ac_july26';
% print('-depsc','-painters','-loose',fname)



