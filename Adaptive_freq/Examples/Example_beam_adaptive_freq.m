% -------------------------------------------------------------------------
%   adaptive beam
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
% Copyright (C) 2019-2020 Karim Cherifi, Pawan Goyal, Peter Benner
% Contact: Karim Cherifi, cherifi@mpi-magdeburg.mpg.de
% -------------------------------------------------------------------------

clearvars; clc; close all
addpath(genpath('../Functions/'))

% Choosing the filter 'Fbell' as descirbed in the maunscript
Filter_type = 'Fbell'; %'notch';

params.bw = 0.2; %bandwidth of the bell filter

params.tol = 1e-12; % Tolerance in the truncation step of the Loewner framework

%% Define a system
load('beam.mat'); 
E       = eye(size(A)); 
[n,m]   = size(B);
D=0;
% Define the original transfer function 
H_orig  = @(s) C*((s*E-A)\B) + D;


%% Generate the frequency response data by taking Np in a defined frequency range
Maxfreq    = log10(5) ;
Minfreq    = -1 ;
Max_ite    = 350;

params. NormFlag   = 0;

[IdenaModel,out,wt,outequi] = adaptive_freq(H_orig,Max_ite,Maxfreq,Minfreq,params);


%% Compare the transfer functions of the original and identified systems
tikzflag = 0 ;

sH_orig = zeros(1,length(wt));
sH_Loew = zeros(1,length(wt));
sH_equi = zeros(1,length(wt));

for j = 1:length(wt)
    sH_orig(j) = H_orig(wt(j)*1i);
    sH_Loew(j) = out.LoewModel.TF(wt(j)*1i);
    sH_equi(j) = outequi.LoewModel.TF(wt(j)*1i);
    
end

w_equi = outequi.w_equi;
snew   = out.snew;

for j = 1:length(snew)
    orig_points(j) = H_orig(snew(j)*1i);
    int_points(j)  = IdenaModel.TF(snew(j)*1i);
    equi_points(j) = outequi.LoewModel.TF(w_equi(j)*1i);
end


figure(1)
subplot(1,2,1)
loglog(wt,abs(sH_orig),'b',wt,abs(sH_Loew),'--g',wt,abs(sH_equi),'-.r',snew,abs(int_points),'+k',w_equi,abs(equi_points),'om' )
legend('original','Adaptive Loewner','equidistant Loewner','adaptive int points','equidistant int points','Location','Best')
subplot(1,2,2)
loglog(wt,abs(sH_orig-sH_Loew),'--g',wt,abs(sH_orig-sH_equi),'-.r',snew,abs(int_points-orig_points),'+k')
legend('Adaptive Loewner','equidistant Loewner','Location','Best')

if tikzflag == 1
addpath('M2tikz/')
matlab2tikz('Exm_beam_TF.tikz','width','\fwidth','height','\fheight')
end
