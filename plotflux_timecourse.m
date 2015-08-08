function plotflux_timecourse(flux,t,model,hsubfig,hfig)
if nargin < 5
    hfig = initialize();
end
if nargin < 4
    [~,hsubfig] = initialize(hfig,flux);
end
nflux = length(flux);
if rem(nflux,2)==0
    n = nflux/2;
else
    n = (nflux+1)/2;
end

if n > 4
% if n > 8
    if rem(n,2)==0
        nrows = n/2;
    else
        nrows=(n+1)/2;
    end

else
    nrows = 2;
%     nrows = 4;
end
    
for ivar = 1:nflux
    if ~isfield(model,'Regulators')
        y_label1 = sprintf('%s mmole/gDCW',model.rxns{ivar});
    end
    if hsubfig(ivar) ~= 0
        hca = findobj(hsubfig(ivar),'type','axes');
        set(hfig,'CurrentAxes',hca);  
    else%subplot is unassigned
        hsubfig(ivar) = subplot(nrows,n,ivar);       
        hca = gca;
        %Make sure more plots can be added at end of loop
        set(hca,'NextPlot','add');     
        set(hsubfig(ivar),'NextPlot','add');
    end
    hline = line(t,flux(ivar));
    %Object Properties
    set(hline,'LineStyle','-',...
              'Marker','o',...
              'MarkerSize',3,...
              'MarkerFaceColor',[1 0 0],...
              'MarkerEdgeColor',[1 0 0]);
    %Axis Properties
    set(get(gca,'YLabel'),'String',y_label1);  
    set(get(gca,'YLabel'),'FontName','Arabic Type Setting');  
    set(get(gca,'YLabel'),'FontSize',8); 
    
    xlabel = sprintf('Steady States');
    set(get(gca,'XLabel'),'String',xlabel);  
    set(get(gca,'XLabel'),'FontName','Arabic Type Setting');   
    set(get(gca,'XLabel'),'FontSize',12);
end

return
function [hfig,hsubfig] = initialize(hfig,flux,hsubfig)
if nargin < 1
    if isempty(findobj('type','figure','Name','Flux Course'))
                hfig = figure('Name','Flux Course'); 
    else
        hfig = findobj('type','figure','Name','Flux Course');
    end
    figure(hfig);
    return
end    

if nargin < 3
    hsubfig = findobj(hfig,'type','axes');
    if isempty(hsubfig)
        hsubfig = zeros(length(flux),1);
    end
end

return