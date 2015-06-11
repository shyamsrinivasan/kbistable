function [hsubfig,new_data] = steadystateplot(model,ng,varname,data,hfig,LineP,hsubfig)
if isfield(model,'nvar')
    nvar = model.nvar;
else
    nvar = size(data.y,1);
end
%data.x - parameter values (scalar)
nadd = 0;
Var = cell(nvar,1);
Var(1:ng(1)) = model.Gene;
Var(ng(1)+1:ng(1)+ng(2)) = model.Enzyme;
Var(ng(1)+ng(2)+1:ng(1)+ng(2)+ng(3)) = model.Metabolites(1:ng(3));
%Var(end) = model.Metabolites(end);
%Add protein name
for ivar = 1:length(varname)
    var_tf = strcmpi(varname{ivar},Var);
    if any(var_tf(1:ng(1)))        
        nadd = nadd+1;
        varname{length(varname)+1} = model.Regulators{logical(full(model.trate*var_tf(1:ng(1))))};        
    end    
end
nplots = length(varname);
if nargin < 7 || isempty(hsubfig)
    hsubfig = zeros(nplots,1);
end
if nargin < 6 || isempty(LineP)
    LineP = struct();
    LineP.Marker = 'o';
    LineP.MarkerSize = 10;
    LineP.MarkerEdgeColor = 'none';
    LineP.MarkerFaceColor = [0 .5 0];
end
if nargin < 5
    hfig = figure;    
else
    figure(hfig);    
end
if rem(nplots,2) == 0  
   n = nplots/2;
else
    n = (nplots+1)/2;
end
%Collecting data points
new_data = struct();
new_data.x = data.x;
new_data.y = [];
new_data.labels = {};
for ivar = 1:nplots    
    %Plotting data         
    var_tf = strcmpi(varname{ivar},Var);%TRN model/Integrated Model 
    if any(var_tf(1:ng(1)))||any(var_tf(ng(1)+1:ng(1)+ng(2)))%mRNA or protein        
        g_tf = strcmpi(varname{ivar},model.Gene);
        pg_tf = strcmpi(varname{ivar},model.Regulators);
        if any(g_tf)
            %mRNA
            y_label1 = sprintf('%s mRNA umole/gDCW \n Steady state',model.Gene{g_tf});
            new_data.labels = [new_data.labels;model.Gene{g_tf}];
            var_tf = logical([var_tf(1:ng(1));zeros(nvar-ng(1),1)]);     
%             LineP.Displayname = sprintf('%s',model.Gene{g_tf});          
        elseif any(pg_tf)
            %Protein            
            y_label1 = sprintf('%s umole/gDCW \n Steady state',model.Regulators{pg_tf});
            new_data.labels = [new_data.labels;model.Regulators{pg_tf}];
            var_tf = logical([zeros(ng(1),1);pg_tf;zeros(nvar-ng(1)-length(pg_tf),1)]);        
%             LineP.Displayname = sprintf('%s',model.Regulators{pg_tf});            
        end
    elseif any(var_tf(ng(1)+ng(2)+1:end))
        m_tf = strcmpi(varname{ivar},model.Metabolites);
%         if any(m_tf) && var_tf(end)%biomass
%             y_label1 = sprintf('%s (gDCW)',model.Metabolites{end});
%             LineP.Displayname = y_label1;
        if any(m_tf)%Metabolite
            y_label1 = sprintf('%s mmole/gDCW',model.Metabolites{m_tf});
            new_data.labels = [new_data.labels;model.Metabolites{m_tf}];
%             LineP.Displayname = sprintf('%s',model.Metabolites{m_tf});
        end
    else
        fprintf('Variable %s does not Exist\n',varname{ivar});
        continue
    end
    %Figure Properties, Data, etc.,
    figure(hfig);
    if hsubfig(ivar) ~= 0
        hca = findobj(hsubfig(ivar),'type','axes');
        set(hfig,'CurrentAxes',hca);  
    else%subplot is unassigned
        hsubfig(ivar) = subplot(2,n,ivar);       
        hca = gca;
        %Make sure more plots can be added at end of loop
        set(hca,'NextPlot','add');     
        set(hsubfig(ivar),'NextPlot','add');
    end
    hline1 = plot(hca,data.x,data.y(var_tf));    
    new_data.y = [new_data.y;data.y(var_tf)];    
%     if call %if called from trnperturbation
        setProperties(hfig,hsubfig(ivar),hca,hline1,data.x,data.y(var_tf),LineP);
%     else
%         set(hline1,LineP);
%     end
    x_label = sprintf('Parameter, %s',data.xlabel);
    % whitebg(hfig,[0 0 0]);
    set(get(hca,'XLabel'),'String',x_label); 
    set(get(hca,'XLabel'),'FontName','Lucida Sans');
    set(get(hca,'XLabel'),'FontSize',12);
    set(get(hca,'YLabel'),'String',y_label1); 
    set(get(hca,'YLabel'),'FontName','Lucida Sans');   
    set(get(hca,'YLabel'),'FontSize',12);
end
return
function setProperties(hfig,hsubfig,haxis,hgobj,XData,YData,LineP)
if nargin<7 
    LineP = struct();
    LineP.Marker = 'o';
    LineP.MarkerSize = 10;
    LineP.MarkerEdgeColor = 'none';
    LineP.MarkerFaceColor = [0 .5 0];
else
    if ~isfield(LineP,'Marker')
        LineP.Marker = 2.0;
    end
    if ~isfield(LineP,'MarkerSize')
        LineP.MarkerSize = 10;
    end
    if ~isfield(LineP,'MarkerEdgeColor')
        LineP.MarkerEdgeColor = 'none';
    end
    if ~isfield(LineP,'MarkerFaceColor')
        LineP.MarkerFaceColor = [0 .5 0];
    end
end
Xmax = max(XData);
Xdiv = sd_round(Xmax/10,1,1);
%Axis properties
AxisP = struct();
AxisP.XMinorTick = 'on';
AxisP.YMinorTick = 'on';
AxisP.TickLength = [0.04,0.04];
AxisP.FontName = 'Lucida Sans';
AxisP.FontSize = 12; 
AxisP.XColor = [.1 .1 .1];
AxisP.YColor = [.1 .1 .1];
% AxisP.XTick = 0:Xdiv:Xmax;
AxisP.XLim = [0 Xmax];

nlines = length(findobj(hsubfig,'type','line'));
if nlines > 1%Presence of 2 or more axes
  hgobj = findobj(haxis,'type','line'); 
  YDatam = zeros(nlines,1);
  YDatan = zeros(nlines,1);
  for iline = 1:nlines
    YDatam(iline) = max(get(hgobj(iline),'YData'));
    YDatan(iline) = min(get(hgobj(iline),'YData'));    
  end
  Ymax = sd_round(max(YDatam),3,5); 
  Ymin = sd_round(min(YDatan),3,5);
  Ydiv = sd_round((max(YDatam)-min(YDatan))/10,2,5);
else
  Ymin = sd_round(min(YData),3,5);
  Ymax = sd_round(max(YData),3,5);
  Ydiv = sd_round((max(YData)-min(YData))/10,2,5); 
end
if Ymin == Ymax 
%     AxisP.YLim = [Ymin-1 Ymax+1];
else
%     AxisP.YLim = [Ymin Ymax]; 
end
% AxisP.YTick = Ymin:Ydiv:Ymax;
set(haxis,AxisP);
set(get(haxis,'XLabel'),'FontName','Lucida Sans');
set(get(haxis,'XLabel'),'FontSize',12);
set(get(haxis,'YLabel'),'FontName','Lucida Sans');   
set(get(haxis,'YLabel'),'FontSize',12);
% legend(haxis,'show');
% legend(haxis,'boxoff');
%Common Graphic Object(line) Properties
if ~isempty(LineP)
    set(hgobj(1),LineP);
end
return