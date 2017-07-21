function hfig = posbifurplot(y,s1,f1,idx,axisfun,pid,hfig,addannot)
% idx - [x-axis y-axis z-axis] index values
if nargin<8
    addannot = [];
end
if nargin<7
    hfig = figure;
end
if nargin<6
    pid = 1;
end
if nargin<5
    axisfun = [];
end
global annot
annot = addannot;

% nvar = size(f1,1);
% npar = size(x1,1)-size(f1,1);

% y = x1(1:nvar,:);
% p = x1(nvar+1:nvar+npar,:);

if size(s1,1)>2
    xindex = cat(1,s1.index);
%     xindex = xindex(2:end-1);
end
evalRe = real(f1);
% get positive vals for y(end,:)
posid = find(y(end,:)<0,1,'first');
% posy = y(:,1:posid);

for kndex = 1:length(xindex)    
    % plot only if xindex is positive or within first positive index
    if xindex(kndex)<posid && xindex(kndex+1)<posid
        % plot data from initial value to first bifurcation
        if xindex(kndex)+1<=size(y,2)
            if any(evalRe(:,xindex(kndex))<0) && any(evalRe(:,xindex(kndex)+1)>=0)
                % entering unstable from stable region 
                LineP.LineStyle = '--';
                LineP.Color = 'k';
            elseif any(evalRe(:,xindex(kndex))>=0) && any(evalRe(:,xindex(kndex)+1)<0)
                % entering stable from unstable region
                LineP.LineStyle = '-';
                LineP.Color = 'k';
            elseif any(evalRe(:,xindex(kndex))<0) && any(evalRe(:,xindex(kndex)+1)<0)
                % staying in stable region
                LineP.LineStyle = '-';
                LineP.Color = 'k';
            elseif any(evalRe(:,xindex(kndex))>=0) && any(evalRe(:,xindex(kndex)+1)>=0)
                % staying in unstable region
                LineP.LineStyle = '--';
                LineP.Color = 'k';
            end
        end
        LineP.LineWidth = 1.5;
        LPval = y(:,xindex(2:end-1));
    
        if length(idx)<=2
            % plot2D
            if (kndex+1)<=length(xindex)
                xval = y(idx(1),xindex(kndex):xindex(kndex+1));
                yval = y(idx(2),xindex(kndex):xindex(kndex+1));
            else
                xval = y(idx(1),xindex(kndex):end);
                yval = y(idx(2),xindex(kndex):end);
            end
            % plot data
            hfig = plot2Dbifurcation(xval,yval,idx(1),idx(2),...
                                    LPval,pid,hfig,LineP,axisfun);
        elseif length(idx)<=3
            % plot3D 
            if (kndex+1)<=length(xindex)
                xval = y(idx(1),xindex(kndex):xindex(kndex+1));
                yval = y(idx(2),xindex(kndex):xindex(kndex+1));
                zval = y(idx(3),xindex(kndex):xindex(kndex+1));
            else
                xval = y(idx(1),xindex(kndex):end);
                yval = y(idx(2),xindex(kndex):end);
                zval = y(idx(3),xindex(kndex):end);
            end
            % get only positive data
%             posid = find(xval<0|yval<0|zval<0,1,'first');
%             posxval = xval(1:posid);
%             posyval = yval(1:posid);
%             poszval = zval(1:posid);
            % plot data
            hfig = plot3Dbifurcation(xval,yval,zval,idx(1),idx(2),idx(3),...
                                     LPval,pid,hfig,LineP,axisfun);
        end
    end
end
end

function hfig =...
        plot2Dbifurcation(xval,yval,xid,yid,LPval,pid,hfig,LineP,axisfun)
global annot
if ~isempty(hfig)
    set(0,'CurrentFigure',hfig);    
    set(gca,'NextPlot','add');
    set(hfig,'CurrentAxes',gca);
    hl = line(xval,yval);    
else
    hfig = figure;
    set(0,'CurrentFigure',hfig);  
    hl = line(xval,yval);
    hfig = gcf;
end

if ~isempty(LineP)
    set(hl,LineP);
else
    set(hl,'LineWidth',3);
end

if size(LPval,1)>4
    if xid > 5
    else
        if ~isempty(axisfun)
            [xlabel,ylabel] = axisfun(2,1,[xid yid]);
        else
            xlabel={};ylabel={};
        end
    end
else
    if xid > 3
        if ~isempty(axisfun)
            [xlabel,ylabel] = axisfun(2,3,[xid yid pid]);
        else
            xlabel = {};ylabel={};
        end
    else
        if ~isempty(axisfun)
            [xlabel,ylabel] = axisfun(2,2,[xid yid]);
        else
            xlabel={};ylabel={};
        end
    end
end
    

% plot bifurcation points
line(LPval(xid,:),LPval(yid,:),'LineStyle','none',...
                             'Marker','o','MarkerEdgeColor','r',...
                             'MarkerFaceColor','r','MarkerSize',6);
setproperties(2,gca,xlabel,ylabel);

if ~isempty(annot)
    text(xval(end),yval(end),annot.text,'FontSize',20);
end
end

function hfig = plot3Dbifurcation(xval,yval,zval,xid,yid,zid,...
                                LPval,pid,hfig,LineP,axisfun)
if ~isempty(hfig)
    figure(hfig);    
    hl = line(xval,yval,zval);
    set(get(hl,'Parent'),'NextPlot','add');
else
    hl = line(xval,yval,zval);
    hfig = gcf;
end

if ~isempty(LineP)
    set(hl,LineP);
else
    set(hl,'LineWidth',3);
end
if ~isempty(axisfun)
    [xlabel,ylabel,zlabel] = axisfun(3,2,[xid yid zid]);
else
    xlabel={};ylabel={};zlabel={};
end

line(LPval(xid,:),LPval(yid,:),LPval(zid,:),'LineStyle','none',...
                             'Marker','o','MarkerEdgeColor','r',...
                             'MarkerFaceColor','r','MarkerSize',6);
axlims = zeros(3,2);
% x-limits
if max(xval)>0
    axlims(1,:) = [0 max(xval)];
else
    axlims(1,:) = [0 1];
end
% y-limits
if max(yval)>0
    axlims(2,:) = [0 max(yval)];
else
    axlims(2,:) = [0 1];
end
% z-limit
if max(zval)>0
    axlims(3,:) = [0 max(zval)];
else
    axlims(3,:) = [0 1];
end
setproperties(3,gca,xlabel,ylabel,axlims,zlabel);  
end

function setproperties(plotype,hsfig,xlabel,ylabel,axlims,zlabel)
if nargin<5 
    axlims = [];
end
if nargin < 6
    zlabel = {};
end
% set axis properties
set(get(gca,'YLabel'),'String',ylabel);  
set(get(gca,'YLabel'),'FontName','Arial');   
set(get(gca,'YLabel'),'FontSize',22); 
set(get(gca,'XLabel'),'String',xlabel);  
set(get(gca,'XLabel'),'FontName','Arial');   
set(get(gca,'XLabel'),'FontSize',22);
if plotype == 3
    set(get(gca,'ZLabel'),'String',zlabel);  
    set(get(gca,'ZLabel'),'FontName','Arial');   
    set(get(gca,'ZLabel'),'FontSize',22);
    axesP.ZColor = [.1 .1 .1];
end
axesP.FontName  = 'Arial';
axesP.FontSize = 22;
axesP.LineWidth = 1.5;
axesP.TickLength = [0.01 0.01];
axesP.XColor = [.1 .1 .1];
axesP.YColor = [.1 .1 .1];
if ~isempty(axlims)
    axesP.XLim = axlims(1,:);
    axesP.YLim = axlims(2,:);
    if plotype==3
        axesP.ZLim = axlims(3,:);
    end
end
if plotype == 3
    set(gca,axesP);
else
    set(hsfig,axesP);
end
end


        


