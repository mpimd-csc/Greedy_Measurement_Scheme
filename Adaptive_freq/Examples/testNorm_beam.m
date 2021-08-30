% -------------------------------------------------------------------------
%   adaptive beam for noisy data
% -------------------------------------------------------------------------
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Affero General Public License for more details.
%
% You should have received a copy of the GNU Affero General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% -------------------------------------------------------------------------
%
% Copyright (C) 2019-2020 Karim Cherifi, Pawan Goyal, Peter Benner
% Contact: Karim Cherifi, cherifi@mpi-magdeburg.mpg.de

clearvars; clc; close all
addpath(genpath('../Functions/'))

%load the example
load('beam.mat'); 
A       = full(A);
E       = eye(size(A)); 
D       = 0;
[n,m]   = size(B);

% Define the original transfer function 
H_orig    = @(s) C*((s*E-A)\B) + D;

origsys.A = A; origsys.B = B; origsys.C = C; origsys.E = E; origsys.D = D;

%Array containing the number of interpolation points
ite_array = 10:4:58;    

%maximum and minimum frequency
Maxfreq    = log10(5) ;
Minfreq    = -1 ;

%bandwidth of the fitler
params.tol      = 1e-14;
params.bw       = 0.2;
params.NormFlag = 1;


for i=1:length(ite_array)
    % compute the Loewner model using only the chosen number of
    % interpolatin points
[IdenaModel,out,wt,outequi] = adaptive_freq(H_orig,ite_array(i),Maxfreq,Minfreq,params);

H2NormEqui(i)   = h2Norm(outequi,origsys,Maxfreq,Minfreq);

H2NormAdap(i)   = h2Norm(out,origsys,Maxfreq,Minfreq);


end

%% plot the H2 norm with respect to the number of chosen interpolation points
tikzflag = 0 ;

figure(1)
semilogy(ite_array(1:5),abs(H2NormAdap(1:5)),'b',ite_array,abs(H2NormEqui),'r')
legend('H2 error for Adaptive Loewner','H2 error for equidistant Loewner','Location','Best')


if tikzflag == 1
addpath('M2tikz/')
matlab2tikz('Exm_beam_H2norm.tikz','width','\fwidth','height','\fheight')
end
