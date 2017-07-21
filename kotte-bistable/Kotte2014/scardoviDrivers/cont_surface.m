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
    contpt = continuation(model,pvec,lambda,9,iguess,opts);
    allcontpts.(['pt' num2str(ip)]) = contpt;
end

