%% Adaptive non-intrusive time domain pH realization: RLC example
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
% -------------------------------------------------------------------------

% Load Oseen example and set parameters
addpath(genpath('../Functions/'));

clearvars
clc
close all
warning off;


Hdata=[];
Select_Idx=[];
Tend        = 100;        % end time
deltaT      = 1e-02;      % time step
dim         = 100;%20;    % nominal dimension of discrete time learned model
kminPercent = 75;         % skip the first 25% of time step

t = 0:deltaT:Tend-deltaT; %time steps

stotal = 2*pi*1j*(-length(t)/2:(length(t)/2- 1))/length(t);
Ntotal = length(stotal);

xaxis = 1:Ntotal;

% Choosing the filter 'Fbell' as descirbed in the maunscript
Filter_type = 'Fbell';

bw = 0.001; % Bandwidth of the filter


load('rlc_serkan200.mat'); %load fom example as (A,B,C,E) system
E = eye(size(A));
origsys.A = A; origsys.B = B; origsys.C = C; origsys.E = E;

%compute the discrete time system
[ Orig.Ad, Orig.Bd, Orig.Cd, Orig.Ed ] = bwdEuler_sys( t, A, B, C, E );
Orig.H   = @(z) Cd*((z*Ed-Ad)\Bd); 


%% Comparison of the first two initial sets of intpoleration points

%expS=exp(stotal);
stotalarray = stotal(1:Ntotal/2);

% select interpolation points in this range
minFreq = 1e-02;
maxFreq = 1e03; %3;

minFreq_time = (1+minFreq*1i*deltaT/2)/(1-minFreq*1i*deltaT/2);
minrep       = repmat(real(minFreq_time), 1,Ntotal/2);
[minValue,position_min] = min(abs(real(exp(stotalarray))-minrep));


maxFreq_time = (1+maxFreq*1i*deltaT/2)/(1-maxFreq*1i*deltaT/2);
maxrep       = repmat(real(maxFreq_time), 1,Ntotal/2);
[maxValue,position_max] = min(abs(real(exp(stotalarray))-maxrep));


%set the inital interpolation points and their complex conjugate
%SelectedIndices1=[2 , length(stotal)-2+2 , length(stotal)/2 , length(stotal)-(length(stotal)/2)+2];
SelectedIndices1 = [position_max+1 , length(stotal)-(position_max+1)+2 , position_min-1 , length(stotal)-(position_min-1)+2];

index1 = position_max+double((uint16((length(stotal)-(position_min-position_max))/2*0.33)))+2;
index2 = position_max+double((uint16((length(stotal)-(position_min-position_max))/2*0.66)))+3;

%SelectedIndices2=[167,length(stotal)-167+2 333 length(stotal)-333+2];
SelectedIndices2 = [index1 , length(stotal)-index1+2 , index2 , length(stotal)-index2+2];



%compute the input the first set of initial interpolation points
u1 = getInputSignal2(t, stotal, length(SelectedIndices1), minFreq, maxFreq,SelectedIndices1);

uInF  = fftshift(fft(u1));
uTInF = zeros(size(uInF));
uTInF(SelectedIndices1) = uInF(SelectedIndices1);
uT    = ifft(ifftshift(uTInF)); % 'truncated' u

% compute the time domain simulation
yT = timeStep(Orig.Ad, Orig.Bd, Orig.Cd, Orig.Ed, u1, zeros(size(Orig.Ad,1))); 

% Compute a discrete time realization from time domain data
params.kminPercent = kminPercent;
params.dim = dim;

[tROM1,Hdata,Select_Idx] = time_DiscreteModel( uInF, yT, stotal, Select_Idx,SelectedIndices1,Hdata, params);

H_Iden1    = @(s) tROM1.Cr*((s*tROM1.Er-tROM1.Ar)\tROM1.Br) ; 

%compute the input the second set of initial interpolation points
u2 = getInputSignal2(t, stotal, length(SelectedIndices2), minFreq, maxFreq,SelectedIndices2);

uInF  = fftshift(fft(u2));
uTInF = zeros(size(uInF));
uTInF(SelectedIndices2) = uInF(SelectedIndices2);
uT    = ifft(ifftshift(uTInF)); % 'truncated' u

% Compute a discrete time realization from time domain data
yT = timeStep(Orig.Ad, Orig.Bd, Orig.Cd, Orig.Ed, u2, zeros(size(Orig.Ad,1))); 

% Compute a realization from time domain data
params.kminPercent  = kminPercent;
params.dim  = dim;


[tROM2,Hdata,Select_Idx] = time_DiscreteModel( uInF, yT, stotal, Select_Idx,SelectedIndices2,Hdata, params);

H_Iden2  = @(s) tROM2.Cr*((s*tROM2.Er-tROM2.Ar)\tROM2.Br) ; 


% compute the error between
sH_Iden1  = zeros(1,length(stotal));
sH_Iden2  = zeros(1,length(stotal));
stest     = exp(stotal);
for j = 1:length(stest)
    sH_Iden2(j) = tROM2.H(stest(j));
    sH_Iden1(j)  = tROM1.H(stest(j));
    err(j)       = abs( sH_Iden2(j)-sH_Iden1(j));
end

for j=1:position_max-1 
err(j)=0;
end

for j=position_min+1:length(stotal)-position_min+2 -1 
err(j)=0;
end

for j=length(stotal)-position_max+2+1 : length(stotal)
err(j)=0;
end

%err(1)=0; err(501)=0;
for j=1:length(SelectedIndices2)
    [ FiltVal ] = Filter_time( SelectedIndices2(j),bw,xaxis,Filter_type);
err=FiltVal.*err;
end

% compute  the maximum error and its complex conjugate
[MaxErr,indexMax] = max(err);

[ FiltVal ] = Filter_time( indexMax,bw,xaxis,Filter_type);
err         = FiltVal.*err;

wo=Ntotal-indexMax+2;
[ FiltVal ] = Filter_time( wo,bw,xaxis,Filter_type);
err         = FiltVal.*err;

% compute the second maximum error
[MaxErr2,indexMax2] = max(err);


sH_Idenorig = sH_Iden2;
tROMnew     = tROM2;


%% add new points and compute the error until the latter is small enough
while (MaxErr>10^-5 && length(Select_Idx)<100)
    
Added_Idx = [indexMax length(stotal)-indexMax+2 indexMax2 length(stotal)-indexMax2+2];

% compute the input with the new set of intrpolation points
unew  = getInputSignal2(t, stotal, length(Added_Idx), minFreq, maxFreq,Added_Idx);

uInF  = fftshift(fft(unew));
uTInF = zeros(size(uInF));
uTInF(Added_Idx) = uInF(Added_Idx);
uT    = ifft(ifftshift(uTInF)); % 'truncated' u

% compute the time domain simulation
yT = timeStep(Orig.Ad, Orig.Bd, Orig.Cd, Orig.Ed, unew, zeros(size(Orig.Ad,1))); 

% Compute a realization from time domain data
params.kminPercent = kminPercent;
params.dim = dim;

[tROMnew,Hdata,Select_Idx] = time_DiscreteModel( uInF, yT, stotal, Select_Idx,Added_Idx,Hdata, params);
H_Idennew                  = @(s) tROMnew.Cr*((s*tROMnew.Er-tROMnew.Ar)\tROMnew.Br) ; 


% compare the error between the old and the new systems
sH_Idennew  = zeros(1,length(stotal));
for j = 1:length(stest)
    sH_Idennew(j) = tROMnew.H(stest(j));
    err(j)        = abs(sH_Idenorig(j)-sH_Idennew(j));
end

for j=1:position_max-1 
err(j) = 0;
end

for j=position_min+1:length(stotal)-position_min+2 -1 
err(j) = 0;
end

for j=length(stotal)-position_max+2+1 : length(stotal)
err(j) = 0;
end


%filter out the new added points
%err(1)=0; err(501)=0;
for j = 1:length(Added_Idx)
[ FiltVal ] = Filter_time( Added_Idx(j),bw,xaxis,Filter_type);
err         = FiltVal.*err;
end

% compute  the maximum error and its complex conjugate
[MaxErr,indexMax]=max(err);

[ FiltVal ] = Filter_time( indexMax,bw,xaxis,Filter_type);
err         = FiltVal.*err;

wo          = Ntotal-indexMax+2;
[ FiltVal ] = Filter_time( wo,bw,xaxis,Filter_type);
err         = FiltVal.*err;

% compute the second maximum error
[MaxErr2,indexMax2] = max(err);

sH_Idenorig         = sH_Idennew;


end

% compute a continous time pH realization from the obtained discrete system
[pHModel] = Discrete2pH(tROMnew , deltaT);

%% Constrcut a realization with logarithmically equdistant interpolation points

load('RLCSerkan_24EquiPoints.mat');

%% Compare the transfer functions of the original and pH systems
tikzflag = 0 ;

H_sys   = @(s) origsys.C*((s*origsys.E-origsys.A)\origsys.B)+D; 
H_Iden  = @(s) pHModel.C*((s*pHModel.E-pHModel.A)\pHModel.B) + pHModel.D+D; 
H_equi  = @(s) pH_Loew.C*((s*pH_Loew.E-pH_Loew.A)\pH_Loew.B) + pH_Loew.D+D; 

w       = logspace(log10(minFreq),log10(maxFreq),300);
sH_sys  = zeros(1,length(w));
sH_Iden = zeros(1,length(w));
sH_equi = zeros(1,length(w));
for j = 1:length(w)
    sH_sys(j)  = H_sys(w(j)*1i);
    sH_Iden(j) = H_Iden(w(j)*1i);
    sH_equi(j) = H_equi(w(j)*1i);
end


%plot the bode plot for the original and identified pH model and the error
%between them.
figure(1)
subplot(1,2,1)
%loglog(w,abs(sH_sys),'*r',w,abs(sH_Iden),'b',w,abs(sH_equi),'-*m')
loglog(w,abs(sH_sys),'b',w,abs(sH_Iden),'--g',w,abs(sH_equi),'-.r')
legend('original system','adaptive','Equidistant','Location','Best')
subplot(1,2,2)
loglog(w,abs(sH_sys-sH_Iden),'--g',w,abs(sH_sys-sH_equi),'-.r')
legend('adaptive','Equidistant','Location','Best')

if tikzflag == 1
addpath('M2tikz/')
matlab2tikz('Exm_RLCSerkan_time_TF.tikz','width','\fwidth','height','\fheight')
end

%plot step reponse of the original and identified pH model
ttest = 0:0.01:4;
u   =  sin(2*ttest)+sin(20*ttest); %+sin(0.001*ttest);

OrigSys = dss(origsys.A,origsys.B,origsys.C,D,origsys.E);
IdenSys = dss(pHModel.A,pHModel.B,pHModel.C,pHModel.D+D,pHModel.E);
EquiSys = dss(pH_Loew.A,pH_Loew.B,pH_Loew.C,pH_Loew.D+D,pH_Loew.E);

y1 = lsim(OrigSys,u,ttest);
y2 = lsim(IdenSys,u,ttest);
y3 = lsim(EquiSys,u,ttest);


figure(2)
plot(ttest,y1,'b',ttest,y2,'--g',ttest,y3,'-.r');
legend('original system','adaptive','Equidistant','Location','Best')

if tikzflag == 1
matlab2tikz('Exm_RLCSerkan_time_step.tikz','width','\fwidth','height','\fheight')
end

