%% Adaptive selection of interpolation points: Penzl example
% This example in section Example 5.1 in the paper
%
% Cherifi K., Goyal, P., and Benner, P., Adaptive selection of 
% interpolation points for data driven modelling, 2020.
%
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

clearvars; clc; close all
addpath(genpath('../Functions/'))

load('fom.mat'); 
E       = eye(size(A)); 
D       = 0;
[n,m]   = size(B);

bw1     = 0.1;
bw2     = 0.6;
bw3     = 3;


% Define the original transfer function 
H_orig  = @(s) C*((s*E-A)\B) + D;

%maximum and minimum frequency
Maxfreq    = 3 ;
Minfreq    = -1 ;
Max_ite    = 350;

params.tol      = 1e-14;
params.NormFlag = 0;

% compute the adaptive Loewner model using bw1
params.bw     = bw1;
[IdenaModel1,out1,wt1,outequi] = adaptive_freq(H_orig,Max_ite,Maxfreq,Minfreq,params);


% compute the adaptive Loewner model using bw2
params.bw     = bw2;
[IdenaModel2,out2,wt2,outequi] = adaptive_freq(H_orig,Max_ite,Maxfreq,Minfreq,params);


% compute the adaptive Loewner model using bw3
params.bw     = bw3;
[IdenaModel3,out3,wt3,outequi] = adaptive_freq(H_orig,Max_ite,Maxfreq,Minfreq,params);


%% compare and plot the bode plot of the different resulting systems
tikzflag = 0 ;

wt  =  sort(unique([wt1,wt2,wt3]));

sH_orig  = zeros(1,length(wt));
sH_Loew1 = zeros(1,length(wt));
sH_Loew2 = zeros(1,length(wt));
sH_Loew3 = zeros(1,length(wt));
sH_equi  = zeros(1,length(wt));

for j = 1:length(wt)
    sH_orig(j)  = H_orig(wt(j)*1i);
    sH_Loew1(j) = out1.LoewModel.TF(wt(j)*1i);
    sH_Loew2(j) = out2.LoewModel.TF(wt(j)*1i);
    sH_Loew3(j) = out3.LoewModel.TF(wt(j)*1i);
    sH_equi(j)  = outequi.LoewModel.TF(wt(j)*1i);
end


figure(1)
loglog(wt,abs(sH_orig-sH_Loew1),'--g',wt,abs(sH_orig-sH_Loew2),'-.r',wt,abs(sH_orig-sH_Loew3),'--m',wt,abs(sH_orig-sH_equi),'--c')
legend(['Adaptive Loewner bw=',num2str(bw1)],['Adaptive Loewner bw=',num2str(bw2)],['Adaptive Loewner bw=',num2str(bw3)],'equidistant Loewner','Location','Best')


if tikzflag == 1
addpath('M2tikz/')
matlab2tikz('Exm_fom_bw.tikz','width','\fwidth','height','\fheight')
end



