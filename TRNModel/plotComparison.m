% function [allSolution,allfinalSS,trnmodel,plotData,h_subfig] =...
%          plotComparison(fname,ng,varname,varargin)
% Calculate and plot comparison for different simulation environments 
% steady state and time course profiles of the integrated network
function [allSolution,allfinalSS,trnmodel,plotData,h_subfig] =...
         plotComparison(fname,ng,varname,varargin)
     
allSolution = varargin{1};
trnmodel = varargin{2};
flag = 0;
if isempty(allSolution)
    flag = 1;
end
fileid = fopen(fname);
if fileid == -1
    fprintf('File %s cannot be opened.', fname);
    trnmodel = struct([]);
    return;
end

E = textscan(fileid, '%d%s%s%s%s%s%s', 'Delimiter', '\t',...
            'TreatAsEmpty', {'None'}, 'HeaderLines', 1);
fclose(fileid);
floc = 'C:\Users\shyam\Documents\Courses\CHE 1125 Project\TRN Model';
if flag
    if ~isempty(E)    
        allSolution = struct();
        nexpt = length(E{1});
        par_ = cell(nexpt,2);  
        initSolution = struct([]);
        for iexpt = 1:nexpt
            %Load initial model
            fprintf('Experiment %d\n',E{1}(iexpt));
            fname = [floc sprintf('\\%s.mat',E{2}{iexpt})];
            load(fname);
            fprintf('Model Loaded: %s\n',E{2}{iexpt});            
            %or build model
    %         [trnmodel,FBAmodel,defparval] = Tmodel(trnfname,regfname,FBAmodel,variable); 
            %Model Initialization
            sim_name = sprintf('sim_%d',iexpt);
            batch = struct();
            batch.init{1} = {'M1xt';'M2xt'};
            batch.init{2} = [1e2;10];%mmoles            
            batch.tpmax = 10;%h  
            %Assign plot variables
            E{3}{iexpt} = strtrim(strrep(E{3}{iexpt},'"',''));
            E{4}{iexpt} = strtrim(strrep(E{4}{iexpt},'"',''));
            E{5}{iexpt} = strtrim(strrep(E{5}{iexpt},'"',''));
            [varname] = ExtractData(2,E{3}{iexpt});
            %assign Parameters
            [trnmodel,parname,parval] = ExtractData(3,trnmodel,E{4}{iexpt},...
                                                               E{5}{iexpt},...
                                                               E{6}{iexpt});
            %Store Parameter Values used?
            par_{iexpt,1}{1} = parname{1};
            par_{iexpt,1}{2} = parval{1};
            %par_{iexpt,2}{1} = parname{2};
%            par_{iexpt,2}{2} = parval{2};           
            %ODE solver parameters
            solverP.RabsTol = 1e-6;
            solverP.PabsTol = 1e-6;
            solverP.MabsTol = 1e-6;
            solverP.RelTol = 1e-3;
            solverP.MaxIter = 1000;    
            solverP.MaxDataPoints = 500;         
            saveData.dirname = 'C:\Users\shyam\Documents\Courses\CHE 1125 Project\Results\TRN Model version 4';
            saveData.filename = sim_name;
            if isempty(initSolution)                 
                %calculate sol for model at gmax = 0.1 and set initSolution 
                if isfield(trnmodel,'gmax')
                    old_gmax = trnmodel.gmax;
                    trnmodel.gmax = 0.8;
                end
                batch.tmax = 432000;%s%batch time for initial solution
                [~,Solution,~,finalSS] =...
                trnperturbation(trnmodel,FBAmodel,batch,defparval,ng,...
                                solverP,initSolution,...
                                varname,saveData); 
                initSolution = Solution.initSS;
                allSolution.init = Solution.initSS;
                allfinalSS.init = finalSS;
                trnmodel.gmax = old_gmax;                    
            end
            %Solve Model ODE 
            batch.tmax = 36000;%s
            [~,Solution,status,finalSS] =...
            trnperturbation(trnmodel,FBAmodel,batch,defparval,ng,...
                            solverP,initSolution,...
                            varname,saveData);    
            if status < 0
                continue
            end
            allSolution.(sim_name) = Solution.initSS;
            allfinalSS.(sim_name) = finalSS;            
            savefile(Solution.initSS,sim_name,saveData);
            close all        
            fprintf('Completed Simulation #%d of %d\n',iexpt,nexpt);
        end
    end
end
%write data to excel file
selectData(allSolution,trnmodel,ng);
LineP = struct();
LineP.LineWidth = 2.0;
nexpt = length(fieldnames(allSolution));
ColorSpec = setLineP(E{7});
for iexpt = 1:nexpt
    %Plot initial solution curve
    if iexpt == 1
        sim_name = 'init'; 
        LineP.Color = 'Black';
        LineP.DisplayName = sprintf('Initial gmax = 0.8');
        LineP.LineStyle = '--';
    elseif iexpt > 1
        %Select model
        sim_name = sprintf('sim_%d',iexpt-1); 
        %Plot solution curve
        LineP.Color = ColorSpec{iexpt-1};
        LineP.LineStyle = '-';
        LineP.Displayname = sprintf('%s=%s',par_{iexpt-1,1}{1},par_{iexpt-1,1}{2});%,...
                                                   %par_{iexpt,2}{1},par_{iexpt,2}{2}); 
    end
    if exist('h_subfig','var')
        [h_fig,h_subfig] =...
        plotSims(sim_name,trnmodel,ng,varname,allSolution,LineP,h_subfig);
    else
        [h_fig,h_subfig] =...
        plotSims(sim_name,trnmodel,ng,varname,allSolution,LineP);
    end        
end
setProperties(h_fig,h_subfig,allSolution.(sim_name));
%Bifurcation Diagram
plotData = [];
% plotData = bifurcationPlot(ColorSpec,par_,nexpt,trnmodel,ng,varname,allfinalSS);
return
function [h_fig,h_subfig] =...
         plotSims(sim_name,trnmodel,ng,varname,allSolution,LineP,h_subfig)
%Plot solution curve                                                
    if isempty(findobj('type','figure','Name','Initial Solution'))
        h_fig = figure('Name','Initial Solution','Color',[1 1 1]);
    else
        h_fig = findobj('type','figure','Name','Initial Solution');
        figure(h_fig);
    end 
    if nargin > 6 || exist('h_subfig','var')
        [h_subfig] =...
        dynamicplot(trnmodel,ng,varname,allSolution.(sim_name),h_fig,LineP,h_subfig);
    else
        [h_subfig] =...
        dynamicplot(trnmodel,ng,varname,allSolution.(sim_name),h_fig,LineP);
    end
return
function [plotData,hb_subfig] =...
         bifurcationPlot(ColorSpec,par_,nexpt,trnmodel,ng,varname,allfinalSS)
plotData.val = zeros(nexpt,0);
for iexpt = 1:nexpt
    %Select model
    sim_name = sprintf('sim_%d',iexpt);
    %Select Data 
    data.x = str2double(par_{iexpt,1}{2});
    data.y = allfinalSS.(sim_name).y;    
    data.xlabel = par_{iexpt,1}{1};
    %Plot solution curve
    LineP.Marker = 'o';
    LineP.MarkerEdgeColor = 'none';
    LineP.MarkerFaceColor = ColorSpec{iexpt};
    LineP.Displayname = sprintf('%s=%s',par_{iexpt,1}{1},par_{iexpt,1}{2});
    if isempty(findobj('type','figure','Name','Bifurcation Diagram'))
        hb_fig = figure('Name','Bifurcation Diagram','Color',[1 1 1]);
    else
        hb_fig = findobj('type','figure','Name','Bifurcation Diagram');
        figure(hb_fig);
    end
    if exist('hb_subfig','var')        
        [hb_subfig,new_data] =...
        steadystateplot(trnmodel,ng,varname,data,hb_fig,LineP,hb_subfig);
    else
        [hb_subfig,new_data] =...
        steadystateplot(trnmodel,ng,varname,data,hb_fig,LineP);
    end  
    plotData.val(iexpt,1) = new_data.x;
    plotData.val(iexpt,2:length(new_data.y)+1) = new_data.y';    
end
plotData.label = [data.xlabel,new_data.labels'];
setProperties(hb_fig,hb_subfig,[],plotData);
return
function setProperties(hfig,hsubfig,Solution,finalSS)
if nargin < 4
    %Solution call
    Xmin = 0;
    Xmax = max(Solution.t(:)/3600);
else
    %bifurcation call
    Xmin = min(finalSS.val(:,1));
    Xmax = max(finalSS.val(:,1));
end
if ~isempty(hsubfig)
    for ivar = 1:length(hsubfig)
%         Xmax = max(Solution.t(:)/3600);
%         Xdiv = sd_round(max(Solution.t(:))/10,1,1);
        %Common Axis Properties
        AxisP = struct();
        AxisP.XMinorTick = 'on';        
%         AxisP.XTick = 0:Xdiv:Xmax;
        AxisP.YMinorTick = 'on';
        AxisP.TickLength = [0.04,0.04];
        AxisP.XColor = [.1 .1 .1];
        AxisP.YColor = [.1 .1 .1];
        AxisP.XLim = [Xmin Xmax];
        AxisP.FontName = 'Lucida Sans';
        AxisP.FontSize = 12; 
        AxisP.FontWeight = 'bold';
        
        %Variable Specific Properties
        hca = findobj(hsubfig(ivar),'type','axes');
        set(get(hca,'XLabel'),'FontName','Lucida Sans',...
                              'FontSize',12,...
                              'FontWeight','bold');
%         set(get(hca,'XLabel'),);
        set(get(hca,'YLabel'),'FontName','Lucida Sans',...
                              'FontSize',12,...
                              'FontWeight','bold');   
%         set(get(hca,'YLabel'),);
        set(hfig,'CurrentAxes',hca);
        hline = findobj(hca,'type','line');
        Y = zeros(length(hline),2);
        iline = 1;
        while iline <= length(hline)
            Y(iline,1) = min(get(hline(iline),'YData'));
            Y(iline,2) = max(get(hline(iline),'YData'));
            iline = iline+1;
        end  
%         Ymin = sd_round(min(Y(:,1)),3,5);
%         Ymax = sd_round(max(Y(:,2)),3,5);
%         Ydiv = sd_round((max(Y(:,2))-min(Y(:,1)))/10,2,5);
%         AxisP.YLim = [Ymin Ymax];           
%         AxisP.YTick = Ymin:Ydiv:Ymax;
        if ~isempty(AxisP)
            set(hca,AxisP);
        end
    end
    legend(findobj(hsubfig(1),'type','axes'),'show');
    legend(findobj(hsubfig(1),'type','axes'),'boxoff');
end
function ColorSpec = setLineP(cl_name)
n_color = length(cl_name);
ColorSpec = cell(n_color,1);
for icolor = 1:n_color
%     fac = icolor/(1+icolor);
%     ColorSpec{icolor} = [1/2*fac 2/3*fac 1/3*fac];  
    ColorSpec{icolor} = rgb(cl_name{icolor});
end
    