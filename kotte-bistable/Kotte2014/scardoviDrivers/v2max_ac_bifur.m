% get continuation diagrams on acetate for different fixed v4 enzyme parameters
% also get continuation diagrams on v4 parameters for
% different fixed acetate concentrations
% do a codimension-2 bifurcation analysis with acetate and these other enzyme
% parameters as the second parameter
runKotte
close all
figh = bifurcationPlot(data.x1,data.s1,data.f1,[4,1]);

% continue on acetate for different values of V4max
v2max = 4:-.05:0;
pvec(ap) = .01;
model.PM(ac-length(xeq)) = .01;
iguess = [.1 .1 .1];

npar = length(v2max);
nbifpts = zeros(npar,1);
nbifpts(1:4) = 800;
nbifpts(4:6) = 800;
nbifpts(7:end) = 1500;
allpvec = repmat(pvec,npar,1);
allpvec(:,14) = v2max;
tspan = 0:0.1:500;
[~,xeq1,~,feq1] = solveODEonly(npar,iguess',model,allpvec,opts,tspan);

for ip = 1:npar
%     pvec = allpvec(ip,:);
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
        allpvec(:,14) = v2max(ip);      
        allpvec(:,ap) = s1.x1(end,:);
        alleqpts = s1.x1(1:3,:);
        allflux = s1.flux(1:5,:);
        recalcdata.x1 = alleqpts;
        recalcdata.flux = s1.flux;
        % recalculate stability info
        for ipos = 1:npts
            model.PM(ac-length(xeq)) = allpvec(ipos,ap);
            recalcdata.f1(:,ipos) =...
            stabilityInfo(@Kotte_givenNLAE,alleqpts(:,ipos)',model,allpvec(ipos,:));
        end
        % draw 2-d bifurcation plot for different v4max values
%         bifurcationPlot([alleqpts;allpvec(:,ap)'],...
%                         s1.s1,recalcdata.f1,[4 1],@getKotteaxislabels);
%         bifurcationPlot([allflux;allpvec(:,ap)'],...
%                         s1.s1,recalcdata.f1,[6 5],@getKotteaxislabels);
        % draw 3-d bifurcation plot
%         bifurcationPlot([alleqpts;allpvec(:,9)';allpvec(:,ap)'],...
%                         s1.s1,recalcdata.f1,[4 5 1],[],1,hfig);
        bifurcationPlot([alleqpts;allpvec(:,14)';allpvec(:,ap)'],...
                        s1.s1,recalcdata.f1,[4 5 1],[],1,hfig1);
        xlabel('V2max a.u.');
        ylabel('acetate a.u.');
        zlabel('pep a.u.');
        view([116 22]);
        grid on
        bifurcationPlot([allflux;allpvec(:,14)';allpvec(:,ap)'],...
                        s1.s1,recalcdata.f1,[6 7 5],[],1,hfig2);
        xlabel('V2max a.u.');
        ylabel('acetate a.u.');
        zlabel('v4 a.u.');
        view([116 22]);
        grid on
%         set(0,'CurrentFigure',hfig3);
%         set(gca,'NextPlot','add');
        bifurcationPlot([alleqpts;allflux(5,:);allpvec(:,ap)'],...
                        s1.s1,recalcdata.f1,[5 1 4],[],1,hfig3);
%         line(allpvec(posid,ap)',alleqpts(1,posid),allflux(5,posid));
        xlabel('acetate a.u.');
        ylabel('pep a.u.');
        zlabel('v4 a.u.');
        view([116 22]);
        grid on
    end
end

%% save figures
% gcf
% fname = 'C:\Users\shyam\Documents\Courses\CHE1125Project\Results\kottedesign\cont_curves\pep_v4_ac_july25';
% print('-depsc','-painters','-loose',fname)











%% continue form liits points obtained on eq continuation
% apLP = [9 11];
% [LPdata_f,LPdata_b] =...
% execLPcont(@KotteMATCONT,xeq,pvec',apLP,ap,data,1000);
% 
% for id = 1:2
%     bifurcationPlot(LPdata_f{id}.x1,LPdata_f{id}.s1,LPdata_f{id}.f1,[4,1],[],1,figh);
%     bifurcationPlot(LPdata_b{id}.x1,LPdata_b{id}.s1,LPdata_b{id}.f1,[4,1],[],1,figh);
% end







