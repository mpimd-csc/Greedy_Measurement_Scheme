function [ u ] = getInputSignal2( t, s, mNotTest, minFreq, maxFreq,SelectedIndices )
% getInputSignal Generates input signal
% In
%   t               ... time steps
%   s               ... sampling points
%   mNotTest        ... number of points that are NOT test points (=number
%                       of interpolation points)
%   minFreq         ... min of frequency range
%   maxFreq         ... max of frequency range
% Out
%   u               ... input signal

%% 
% select important modes
%I = selectImportantModes([], mNotTest, s, minFreq, maxFreq, 'equidistantConjugate');
 Isorted = sort(SelectedIndices, 'ascend');
 I = myReorderS(Isorted );
%I=SelectedIndices1;
assert(abs(sum(imag(s(I)))) < 1e-10);
coeff = 1/length(s)*ones(1, mNotTest) + 1/length(s)*sqrt(-1)*ones(1, mNotTest);
%u = sum(diag(coeff)*exp(s(I(1:mNotTest)).'*(0:length(s)-1)), 1);
nrChunks = ceil(length(s)/1e+6);
u = zeros(size(s));
for i=1:nrChunks
    curIndx = mygetchunk(1:length(s), nrChunks, i);
    u(curIndx) = sum(diag(coeff)*exp(s(I(1:mNotTest)).'*(curIndx-1)), 1);
end
% first input is zero
%u = u - u(1);


%generate input signal
myu = zeros(size(u));
for i=1:length(I)
    myu = myu + 1/length(s)*(1 + 1j)*(cos(imag(s(I(i)))*(0:length(t) - 1)) + 1j*sin(imag(s(I(i)))*(0:length(t)-1)));
%uInF2 = fftshift(fft(myu));
end
assert(norm(myu - u) < 1e-10);
u = myu;

% double check
uInF = fftshift(fft(u));
reorder=myReorderS(find(abs(uInF) > 1e-10)');
assert(norm(I' - myReorderS(find(abs(uInF) > 1e-10)')) == 0);

end

