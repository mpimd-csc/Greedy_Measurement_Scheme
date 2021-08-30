function [ Hdata, res, regCond ] = regProb( kminPercent, uInF, y, s, I )
%regProb infers frequency domain data from time domain data
% In
%   kminPercent     ... percentage of tie domain data to be used
%   uInF            ... truncated input vector
%   s               ... sampling points
%   I               ... index vector
% Out
%   Hdata           ... transfer functions evaluations
%   res             ... frobenius norm
%   regCond         ... condition numbers


%% compute Hdata
startIndx = length(s) - floor(length(s)*(kminPercent/100));
tSel = startIndx:length(s);

FF = (1/length(s))*exp((tSel-1)'*s(I));
D = diag(uInF(I));
condFF = cond(FF);
condD = cond(D);
FFxD = FF*D; %(1/length(t))*exp(t(tSel)'*s(I))*diag(uInF(I));
clear FF D;
condFFxD = cond(FFxD);
%FF = exp(t(tSel)'*s(I));
%disp(['cond ', num2str(cond(FF)), ' ', num2str(cond(D)), ' ', num2str(regCond)]);
Hdata = (FFxD)\(y(tSel).');
res = norm((FFxD)*Hdata - y(tSel).', 'fro');
regCond = [condFF condD condFFxD];

end

