% run reginVolume dynamic simulations from sampled code
load('regionVolumeSamplevals.mat');
% integrate
options = [];
[ppival,npival,allxeq,ssid,allfeq] =...
getPequilibrium(rndivals(30001:end,:),model,pvec,options,opts,tspanf);
