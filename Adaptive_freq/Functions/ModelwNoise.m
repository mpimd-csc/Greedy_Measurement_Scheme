function [sH_Idennew,out, snew] = ModelwNoise( Std,H_orig,D,tol,bw)
% -------------------------------------------------------------------------
% This function identifies a realization using frequency
% response data corrupted with Gaussian white noise. 
%
% -------------------------------------------------------------------------
% The function inputs are:
% Std              --  Standard deviation of the Gausisian white noise,
% H_orig           --  Transfer function handle of the original system,
% D                --  Direct feed-through term,
% tol              --  tolerance for the Loewner pencil that determines the order of
%           the identified state-space model. Its default value is 1e-8,
% bw               --  Bandwidth of the filter.
% -------------------------------------------------------------------------
%
% The function outputs are:
% sH_Idennew   -- Identified Loewner model, the system matrices and transfer
%                   function,
% Out          -- It stores additional information,
% snew         -- Set of inteprolation points used. 
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

Filter_type='Fbell'; % Choose the notch filter type


%% Generate the frequency response data by taking Np in a defined frequency range
Select_Idx = [];
Ntotal     = 1000; % total number of frequencies considered
Maxfreq    = 3 ;  % Maximum frequency
Minfreq    = -1 ; % minimum frequency

w = logspace(Minfreq ,Maxfreq ,Ntotal);

% Select Np initial interpolation points %Np=6;

add_index1  =  double((uint16(length(w)*0.20)));
add_index2  =  double((uint16(length(w)*0.40)));
add_index3  =  double((uint16(length(w)*0.60)));
add_index4  =  double((uint16(length(w)*0.80)));

w1      =  [w(1),w(add_index1),w(add_index2),w(add_index3),w(add_index4),w(end)];

F1_NoN = zeros(1,length(w1));
for j = 1:length(w1)
   F_NoN(j) = H_orig(1i*w1(j));  
end

Wn_real = Std.*randn(1,length(F_NoN)); % White noise with Standard deviation Std
Wn_imag = Std.*randn(1,length(F_NoN)); % White noise with Standard deviation Std

%crrupt the data with noise
F1 = (real(F_NoN) + Wn_real.*(abs(real(F_NoN)))) + 1i *(imag(F_NoN)+ Wn_imag.*(abs(imag(F_NoN))));


% select the next two points
for j = 2:length(w1)
   abs_err(j-1) = abs(F1(j)-F1(j-1));
end
[absMaxErr,AbsIndexMax] = max(abs_err);
MaxPos = (AbsIndexMax*2-1)*0.1;
abs_err(AbsIndexMax) = 0;

[absMaxErr2,AbsIndexMax2] = max(abs_err);
MaxPos2 = (AbsIndexMax2*2-1)*0.1;

index1  =  double((uint16(length(w)*MaxPos)));
index2  =  double((uint16(length(w)*MaxPos2)));

w2      =  [w(index1),w(index2)];

Select_Idx = [1,add_index1,add_index2,add_index3,add_index4, Ntotal, index1, index2 ];

% Evaluate the transfer function at the chosen interpolation points

F_NoN2 = zeros(1,length(w2));
snew   = sort([w1,w2]);

for j = 1:length(snew)
   F_NoN2(j) = H_orig(1i*snew(j));  
end

Wn_real2 = Std.*randn(1,length(F_NoN2)); % White noise with Standard deviation Std
Wn_imag2 = Std.*randn(1,length(F_NoN2)); % White noise with Standard deviation Std

F2 = (real(F_NoN2) + Wn_real2.*(abs(real(F_NoN2)))) + 1i *(imag(F_NoN2)+ Wn_imag2.*(abs(imag(F_NoN2))));

%% Identification of the inital models 
% Construct the two initial models
[IdenModel1, out1] = Construct_aModel(w1,F1,D,tol);
[IdenModel2, out2] = Construct_aModel(snew,F2,D,tol);

% compute the error betwwen the two initial models
sH_Loew1 = zeros(1,length(w));
sH_Loew2 = zeros(1,length(w));
for j = 1:length(w)
    sH_Loew1(j) = IdenModel1.TF(w(j)*1i);
    sH_Loew2(j) = IdenModel2.TF(w(j)*1i);
    err(j) = abs(sH_Loew2(j)-sH_Loew1(j));
end

% Apply the filter to avoid choosing the same interpolation points
[ FiltVal ] = Filter( w(index1),bw,w,Filter_type);
err         = FiltVal.*err;
[ FiltVal ] = Filter( w(index2),bw,w,Filter_type);
err         = FiltVal.*err;

% Compute the maximum error and filter the error
[MaxErr,indexMax] = max(err);

[ FiltVal ] = Filter( w(indexMax), bw, w, Filter_type);
err         = FiltVal.*err;

% Compute the second maximum error
[MaxErr2,indexMax2] = max(err);

% Store the new system
sH_Idenorig = sH_Loew2;

if Std==0
    threshold = 10^-6;
else
    threshold = Std*10^2;
end
%% loop until tolerance reached
while (MaxErr> threshold && length(Select_Idx)<30)
    
% Add the new interpolation points
snew       = sort([snew w(indexMax) w(indexMax2) ]);
Select_Idx = [Select_Idx indexMax  indexMax2];

% Compute the realization of the new system
for j = 1:length(snew)
    Fnew_NoN(j) = H_orig(snew(j)*1i);
end

Wn_real3 = Std.*randn(1,length(Fnew_NoN)); % White noise with Standard deviation Std
Wn_imag3 = Std.*randn(1,length(Fnew_NoN)); % White noise with Standard deviation Std

Fnew = (real(Fnew_NoN) + Wn_real3.*(abs(real(Fnew_NoN)))) + 1i *(imag(Fnew_NoN)+ Wn_imag3.*(abs(imag(Fnew_NoN))));

% Construct a model with the noisy data
[H_Idennew, out] = Construct_aModel(snew,Fnew,D,tol);


% Compare the between the new and the old system
sH_Idennew = zeros(1,length(w));
for j = 1:length(w)
    sH_Idennew(j) = H_Idennew.TF(w(j)*1i);
    err(j)        = abs(sH_Idenorig(j)-sH_Idennew(j));
end

% Apply a filter on all new interpolation points
[ FiltVal ] = Filter( w(indexMax),bw,w,Filter_type);
err         = FiltVal.*err;

[ FiltVal ] = Filter( w(indexMax2),bw,w,Filter_type);
err         = FiltVal.*err;

% Compute the maximum error
[MaxErr,indexMax] = max(err);

[ FiltVal ] = Filter( w(indexMax), bw, w, Filter_type);
err         = FiltVal.*err;

% Compute the second maximum error
[MaxErr2,indexMax2] = max(err);

% Store the new system
sH_Idenorig = sH_Idennew;

end








end
