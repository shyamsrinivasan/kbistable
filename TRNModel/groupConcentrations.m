function [newPIV,newColr] = groupConcentrations(petconc,conc,multiSS)
[nvar,nsamp] = size(petconc);
newPIV = cell(nvar,1);
newColr = cell(nvar,1);
load('C:\Users\shyam\Documents\Courses\CHE 1125 Project\Colors.mat');

for iv = 1:nvar
    pvi = petconc(iv,1:end);    
%     pviint = floor(pvi);
%     pviflt = roundsd(pvi-pviint,3,'ceil');
%     pvi = pviint+pviflt;
    
    ss = multiSS{iv}(:,1);
    nss = length(ss);
    newPIV{iv} = cell(nss,1);
    assigInd = zeros(nsamp,1);
    Clrnd = floor(1 + (76-1)*rand(nss,1));
    newColr{iv} = cell(length(Clrnd),1);
    iss = 1;
    while iss <= nss
        ssmat = repmat(ss(iss),1,nsamp);
        if ss(iss) > 1 && ss(iss) < 10
            ssr = abs(conc(iv,2:end)-ssmat)<1e-2;
            if any(ssr) 
                newPIV{iv}{iss} = pvi(setdiff(find(ssr),find(assigInd)));
                assigInd(ssr) = 1;
                try
                    newColr{iv}{iss} = rgb(Colors{Clrnd(iss)});
                catch
                    newColr{iv}{iss} = rgb('Black');
                end
            end
        elseif ss(iss) > 10
            ssr = abs(conc(iv,2:end)-ssmat)<1e-1;
            if any(ssr)
                newPIV{iv}{iss} = pvi(setdiff(find(ssr),find(assigInd)));
                assigInd(ssr) = 1;
                try
                    newColr{iv}{iss} = rgb(Colors{Clrnd(iss)});
                catch
                    newColr{iv}{iss} = rgb('Black');
                end
            end
        elseif ss(iss) > 1e-3
            ssr = abs(conc(iv,2:end)-ssmat)<1e-4;
            if any(ssr)               
                newPIV{iv}{iss} = pvi(setdiff(find(ssr),find(assigInd)));
                assigInd(ssr) = 1;
                try
                    newColr{iv}{iss} = rgb(Colors{Clrnd(iss)});
                catch
                    newColr{iv}{iss} = rgb('Black');
                end
            end 
        else
            ssr = abs(conc(iv,2:end)-ssmat)<1e-6;
            if any(ssr)
                newPIV{iv}{iss} = pvi(setdiff(find(ssr),find(assigInd))); 
                assigInd(ssr) = 1;
                try
                    newColr{iv}{iss} = rgb(Colors{Clrnd(iss)});
                catch
                    newColr{iv}{iss} = rgb('Black');
                end
            end
        end
%         if any(conc(ivar,:)-ssmat)
%             nagr = length(find(conc(ivar,:) == ssmat));
%             newPI(ivar,1:nagr) = petconc(conc(ivar,:) == ssmat);
%         end
        iss = iss + 1;
    end   
end