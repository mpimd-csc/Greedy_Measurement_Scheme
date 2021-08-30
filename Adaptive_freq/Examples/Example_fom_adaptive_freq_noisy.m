% -------------------------------------------------------------------------
% adaptive fom
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
Filter_type = 'Fbell'; 

Ntotal     = 1000;
Maxfreq    = 3 ;
Minfreq    = -1 ;

w = logspace(Minfreq ,Maxfreq ,Ntotal);

rng(0);

bw = 0.6 ; %bandwidth of the bell filter

Std0 =  0;        % noise free case 

Std1 =  10^-4;    % Standard deviation of the white noise
Std2 =  10^-5;    % Standard deviation of the white noise
Std3 =  10^-6;    % Standard deviation of the white noise

tol = 1e-14; % Tolerance in the truncation step of the Loewner framework

%% Define a system
load('fom.mat'); 
E       = eye(size(A)); 
[n,m]   = size(B);
D = 0;  A = full(A);
% Define the original transfer function 
H_orig  = @(s) C*((s*E-A)\B) + D;

origsys.A = A; origsys.B = B; origsys.C = C; origsys.E = E; origsys.D = D;

[sH_Idennew0,out0,snew0] = ModelwNoise(Std0,H_orig,D,tol,bw);
[sH_Idennew1,out1,snew1]  = ModelwNoise(Std1,H_orig,D,tol,bw);
[sH_Idennew2,out2,snew2]  = ModelwNoise(Std2,H_orig,D,tol,bw);
[sH_Idennew3,out3,snew3]  = ModelwNoise(Std3,H_orig,D,tol,bw);


%% construct model with logarithmically equidistant interpolation points
%w_equi = logspace(Minfreq ,Maxfreq ,length(snew)-2);
w_equi = logspace(Minfreq ,Maxfreq ,30);
for j = 1:length(w_equi)
   F_equi(j) = H_orig(1i*w_equi(j));  
end
Wn_real1 = Std1.*randn(1,length(F_equi)); % White noise with Standard deviation Std
Wn_imag1 = Std1.*randn(1,length(F_equi)); % White noise with Standard deviation Std

Wn_real2 = Std2.*randn(1,length(F_equi)); % White noise with Standard deviation Std
Wn_imag2 = Std2.*randn(1,length(F_equi)); % White noise with Standard deviation Std

Wn_real3 = Std3.*randn(1,length(F_equi)); % White noise with Standard deviation Std
Wn_imag3 = Std3.*randn(1,length(F_equi)); % White noise with Standard deviation Std

F_equi1 = (real(F_equi) + Wn_real1.*(abs(real(F_equi)))) + 1i *(imag(F_equi)+ Wn_imag1.*(abs(imag(F_equi)))); %Gaussian white noise W
F_equi2 = (real(F_equi) + Wn_real2.*(abs(real(F_equi)))) + 1i *(imag(F_equi)+ Wn_imag2.*(abs(imag(F_equi)))); %Gaussian white noise W
F_equi3 = (real(F_equi) + Wn_real3.*(abs(real(F_equi)))) + 1i *(imag(F_equi)+ Wn_imag3.*(abs(imag(F_equi)))); %Gaussian white noise W

[H_equi0, outequi0] = Construct_aModel(w_equi,F_equi,D,tol);
[H_equi1, outequi1] = Construct_aModel(w_equi,F_equi1,D,tol);
[H_equi2, outequi2] = Construct_aModel(w_equi,F_equi2,D,tol);
[H_equi3, outequi3] = Construct_aModel(w_equi,F_equi3,D,tol);


%[H_equi, outequi] = Construct_aModel(w_equi,F_equi,D,tol);


%% compute norms
sys_orig    = dss(A,B,C,0,E);

sysAdap_Std0   = dss(out0.LoewModel.A,out0.LoewModel.B,out0.LoewModel.C,0,out0.LoewModel.E);
sysAdap_Std1   = dss(out1.LoewModel.A,out1.LoewModel.B,out1.LoewModel.C,0,out1.LoewModel.E);
sysAdap_Std2   = dss(out2.LoewModel.A,out2.LoewModel.B,out2.LoewModel.C,0,out2.LoewModel.E);
sysAdap_Std3   = dss(out3.LoewModel.A,out3.LoewModel.B,out3.LoewModel.C,0,out3.LoewModel.E);

sysEqui_Std0   = dss(outequi0.LoewModel.A,outequi0.LoewModel.B,outequi0.LoewModel.C,0,outequi0.LoewModel.E);
sysEqui_Std1   = dss(outequi1.LoewModel.A,outequi1.LoewModel.B,outequi1.LoewModel.C,0,outequi1.LoewModel.E);
sysEqui_Std2   = dss(outequi2.LoewModel.A,outequi2.LoewModel.B,outequi2.LoewModel.C,0,outequi2.LoewModel.E);
sysEqui_Std3   = dss(outequi3.LoewModel.A,outequi3.LoewModel.B,outequi3.LoewModel.C,0,outequi3.LoewModel.E);

norm2_orig          = norm(sys_orig);

H2NormAdpt.Std0     = h2Norm(out0,origsys,Maxfreq,Minfreq);
H2NormAdpt.Std1     = h2Norm(out1,origsys,Maxfreq,Minfreq);
H2NormAdpt.Std2     = h2Norm(out2,origsys,Maxfreq,Minfreq);
H2NormAdpt.Std3     = h2Norm(out3,origsys,Maxfreq,Minfreq);

H2NormEqui.Std0     = h2Norm(outequi0,origsys,Maxfreq,Minfreq);
H2NormEqui.Std1     = h2Norm(outequi1,origsys,Maxfreq,Minfreq);
H2NormEqui.Std2     = h2Norm(outequi2,origsys,Maxfreq,Minfreq);
H2NormEqui.Std3     = h2Norm(outequi3,origsys,Maxfreq,Minfreq);
%}

%% Compare the transfer functions of the original and identified systems
%{
%IdenpHModel=H_Idennew;
%wt  =  logspace(-1,3,Ntotal);
wt  =  sort([w,snew]);
sH_orig = zeros(1,length(wt));
sH_Loew0 = zeros(1,length(wt));
sH_Loew1 = zeros(1,length(wt));
sH_Loew2 = zeros(1,length(wt));
sH_Loew3 = zeros(1,length(wt));
sH_equi = zeros(1,length(wt));

for j = 1:length(wt)
    sH_orig(j) = H_orig(wt(j)*1i);
    sH_Loew0(j) = out0.LoewModel.TF(wt(j)*1i);
    sH_Loew1(j) = out1.LoewModel.TF(wt(j)*1i);
    sH_Loew2(j) = out2.LoewModel.TF(wt(j)*1i);
    sH_Loew3(j) = out3.LoewModel.TF(wt(j)*1i);
    sH_equi(j)= outequi.LoewModel.TF(wt(j)*1i);
    
end

figure(1)
subplot(1,2,1)
loglog(wt,abs(sH_orig),'b',wt,abs(sH_equi),'-.r',wt,abs(sH_Loew0),'--c',wt,abs(sH_Loew1),'--g',wt,abs(sH_Loew2),'--k',wt,abs(sH_Loew3),'--m' )
legend('Original','Equidistant Loewner','Adaptive Loewner Std=0','Adaptive Loewner Std=10^-4','Adaptive Loewner Std=10^-5','Adaptive Loewner Std=10^-6','Location','Best')
subplot(1,2,2)
loglog(wt,abs(sH_orig-sH_equi),'-.r',wt,abs(sH_orig-sH_Loew0),'--c',wt,abs(sH_orig-sH_Loew1),'--g',wt,abs(sH_orig-sH_Loew2),'--k',wt,abs(sH_orig-sH_Loew3),'--m')
legend('Equidistant Loewner','Adaptive Loewner Std=0','Adaptive Loewner Std=10^-4','Adaptive Loewner Std=10^-5','Adaptive Loewner Std=10^-6','Location','Best')
%}


