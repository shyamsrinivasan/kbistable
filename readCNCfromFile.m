% function [mc,model] = readCNCfromFile(fname,model)
% read concentrations from text file
function [mc,model] = readCNCfromFile(fname,model)
fileid = fopen(fname);
if fileid == -1
    fprintf('File %s cannot be opened.', fname);
    mc = [];
    return;
end

C = textscan(fileid, '%s%f',...
                     'Delimiter', '\t',...
                     'TreatAsEmpty', {'None'},...
                     'HeaderLines', 1);
mc = zeros(length(model.mets),1);
for ic = 1:length(C{1})
    tfm = strcmpi(C{1}{ic},model.mets);
    if any(tfm)
        mc(tfm) = C{2}(ic);
        if isfield(model,'MClow')
            model.MClow(tfm) = C{2}(ic);
        end
        if isfield(model,'MChigh')
            model.MChigh(tfm) = C{2}(ic);
        end
    end
end

%get metabolite structure suitable for input to iconcetration.m
h2oc = find(strcmpi(newmodel.mets,'h2o[c]'));
he = find(strcmpi(newmodel.mets,'h[e]'));
hc = find(strcmpi(newmodel.mets,'h[c]'));
o2e = find(strcmpi(newmodel.mets,'o2[e]'));

vmet = [he hc h2oc o2e];
newcell = cellfun(@(x)strcmpi(vmet,x),model.mets,'UniformOutput',false);
display(newcell);
% for im = 1:length(model.mets)
    

