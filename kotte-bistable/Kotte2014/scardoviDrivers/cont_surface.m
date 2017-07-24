% contour plot of v4/pep/V4max
% choose 1 acetate
% x - v4max;
% y - pep
% z - v4

% 1. choose acetate and v4max values
% 2. solve NLAE to get pep and v4
% 3. plot them on separate surfaces
lambda = 0.01:.01:4;
iguess = [.01;.01;.01];
v4max = 0:.05:.34;
npar = length(v4max);
allcontpts = struct();
for ip = 1:npar
    pvec(11) = v4max(ip);
    [contpt,fluxpt] = continuation(model,pvec,lambda,9,iguess,opts);
    allcontpts.(['pt' num2str(ip)]).xc = contpt;
    allcontpts.(['pt' num2str(ip)]).flx = fluxpt;
end

%% figures
figure
hold on
plot3(repmat(v4max(7),1,400),lambda,allcontpts.pt7.flx(5,1:400),'r')
plot3(repmat(v4max(7),1,400),fliplr(lambda),allcontpts.pt7.flx(5,401:end))
figure
hold on
plot3(repmat(v4max(7),1,400),lambda,allcontpts.pt7.xc(1,1:400),'r')
plot3(repmat(v4max(7),1,400),fliplr(lambda),allcontpts.pt7.xc(1,401:end))


