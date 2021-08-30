function [IdenaModel, out, wt, outequi]=adaptive_freq(H_orig, Maxite, Maxfreq, Minfreq, params)
% -------------------------------------------------------------------------
% This function computes a data driven realization based on the adaptive
% algorithm presented in the paper:
% Cherifi K., Goyal, P., and Benner, P., Adaptive selection of 
% interpolation points for data driven modelling, 2020.
%
% -------------------------------------------------------------------------
% Inputs:
% H_orig               ---  TF of the original system,
% Maxite               ---  Maximum number if interations of the algorithm,
% Maxfreq              ---  Maximum frequency to be considered,
% Minfreq              ---  Minimum frequency to be considered,
% params.tol           ---  Tolerance of the Loewner framework truncation,
% params.bw            ---  Bandwidth of the filter,
% params.NormFlag      ---  A flag for norm comparison.
%
% -------------------------------------------------------------------------
%
% Outputs:
% IdenaModel           --- Identified model using the adaptive algorithm,
% out                  --- More information about the identified model,
% wt                   --- Set of interpolation points used,
% outequi              --- Model based on equidistant interpolation points.
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


addpath(genpath('../Functions/'))

%choose the filter 'notch' or 'Fbell'
Filter_type = 'Fbell'; %'notch';

tol = params.tol;
bw  = params.bw;

%% Generate the frequency response data by taking Np points in a defined frequency range
Ntotal     = 1000;


w = logspace(Minfreq ,Maxfreq ,Ntotal);
%select Np initial interpolation points Np=6

add_index1  =  double((uint16(length(w)*0.20)));
add_index2  =  double((uint16(length(w)*0.40)));
add_index3  =  double((uint16(length(w)*0.60)));
add_index4  =  double((uint16(length(w)*0.80)));

w1      =  [w(1),w(add_index1),w(add_index2),w(add_index3),w(add_index4),w(end)];

F1 = zeros(1,length(w1));
for j = 1:length(w1)
   F1(j) = H_orig(1i*w1(j));  
end

for j = 2:length(w1)
   abs_err(j-1)=abs(F1(j)-F1(j-1));
end
[absMaxErr,AbsIndexMax] = max(abs_err);

MaxPos                  =  (AbsIndexMax*2-1)*0.1;
abs_err(AbsIndexMax)    =  0;

[absMaxErr2,AbsIndexMax2] = max(abs_err);
MaxPos2                   = (AbsIndexMax2*2-1)*0.1;

index1  =  double((uint16(length(w)*MaxPos)));
index2  =  double((uint16(length(w)*MaxPos2)));


w2      =  [w(index1),w(index2)];

Select_Idx = [1,add_index1,add_index2,add_index3,add_index4, Ntotal, index1, index2 ];

% evaluate the transfer function at the chosen interpolation points

F2   = zeros(1,length(w2));
snew = sort([w1,w2]);

for j = 1:length(snew)
   F2(j) = H_orig(1i*snew(j));  
end

%% Identification of the inital models 
% construct the two initial models
[IdenModel1, out1] = Construct_aModel(w1,F1,0,tol);
[IdenModel2, out2] = Construct_aModel(snew,F2,0,tol);

% compute the error betwwen the two initial models
sH_Loew1 = zeros(1,length(w));
sH_Loew2 = zeros(1,length(w));
for j = 1:length(w)
    sH_Loew1(j) = IdenModel1.TF(w(j)*1i);
    sH_Loew2(j) = IdenModel2.TF(w(j)*1i);
    err(j)      = abs(sH_Loew2(j)-sH_Loew1(j));
end

[ FiltVal ] = Filter( w(index1),bw,w,Filter_type);
err = FiltVal.*err;
[ FiltVal ] = Filter( w(index2),bw,w,Filter_type);
err = FiltVal.*err;

%compute the maximum error and filter the error
[MaxErr,indexMax] = max(err);

[ FiltVal ] = Filter( w(indexMax), bw, w, Filter_type);
err         = FiltVal.*err;

%compute the second maximum error
[MaxErr2,indexMax2] = max(err);

%store the new system
sH_Idenorig = sH_Loew2;


%% loop until tolerance reached
while (MaxErr>10^-8 && length(Select_Idx)<Maxite)
%add the new interpolation points
snew       = sort([snew w(indexMax) w(indexMax2) ]);
Select_Idx = [Select_Idx indexMax  indexMax2];

%compute the realization of the new system
for j = 1:length(snew)
    Fnew(j) = H_orig(snew(j)*1i);
end

[H_Idennew, out] = Construct_aModel(snew,Fnew,0,tol);


% compare the between the new and the old system
sH_Idennew = zeros(1,length(w));
for j = 1:length(w)
    sH_Idennew(j) = H_Idennew.TF(w(j)*1i);
    err(j)        = abs(sH_Idenorig(j)-sH_Idennew(j));
end

% apply a filter on all new interpolation points
[ FiltVal ] = Filter( w(indexMax),bw,w,Filter_type);
err         = FiltVal.*err;


[ FiltVal ] = Filter( w(indexMax2),bw,w,Filter_type);
err         = FiltVal.*err;

%compute the maximum error
[MaxErr,indexMax] = max(err);

[ FiltVal ] = Filter( w(indexMax), bw, w, Filter_type);
err         = FiltVal.*err;

%compute the second maximum error
[MaxErr2,indexMax2] = max(err);

%store the new system
sH_Idenorig = sH_Idennew;

end

IdenaModel = H_Idennew;
wt         =  sort([w,snew]);

out.snew   = snew;


%% construct model with equidistant interpolation points
if params.NormFlag == 1
    w_equi = logspace(Minfreq ,Maxfreq ,Maxite);
else
    w_equi = logspace(Minfreq ,Maxfreq ,length(snew));
end

for j = 1:length(w_equi)
   F_equi(j) = H_orig(1i*w_equi(j));  
end

[H_equi, outequi] = Construct_aModel(w_equi,F_equi,0,tol);

outequi.w_equi = w_equi;


end
