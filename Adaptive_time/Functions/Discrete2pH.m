function [pHModel] = Discrete2pH( tROM , deltaT)
% -------------------------------------------------------------------------
% This function computes a pH realization from a discrete time state space
% system
%
% -------------------------------------------------------------------------
% The function inputs are:
% tROM          --  Discrete time state sapce system,
% deltaT        --  Sampling time.
% -------------------------------------------------------------------------
%
% The function outputs are:
% pHModel       --  pH realization.
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
% Copyright (C) 2020 Karim Cherifi, Pawan Goyal, Peter Benner
% Contact: Karim Cherifi, cherifi@mpi-magdeburg.mpg.de
% ------------------------------------------------------------------------- 

%% discrete to continous transformation

%{
clear A B C D E
sysd=dss(tROM.Ar,tROM.Br,-tROM.Cr,0,tROM.Er,deltaT);
sysc = d2c(sysd,'tustin');
[A,B,C,D,E] = dssdata(sysc);
%}


% based on Backward Euler
clear A B C D E
A=(tROM.Ar-tROM.Er)/deltaT;
B=tROM.Br/deltaT;
C=-tROM.Cr;
D=0;
E=tROM.Ar;
%}

%% Define the transfer function for the resulting system
n=size(A,1);
% Define the interpolated transfer function 
H_orig  = @(s) C*((s*E-A)\B) + D;

%define interpolation points used in Contruct_pHModel
Np      = 2*n;
w = logspace(-1,3,Np);
F = zeros(1,Np);
for j = 1:Np
   F(j) = H_orig(1i*w(j));  
end
D=H_orig(10^10);

%% pH realization
% compute spectral zeros, spectral directions and deduce the pH representation 
[pHModel, out] = Construct_pHModel(w,F,D);

%% Compute PH framework matrices
S = [-pHModel.A -pHModel.B; pHModel.C pHModel.D];

Vph =(S-S')/2;
Wph =(S+S')/2;

nph=size(pHModel.A,1);

pHModel.J = Vph(1:nph,1:nph);
pHModel.G = Vph(1:nph,nph+1:end);
pHModel.N = Vph(nph+1:end,nph+1:end);
pHModel.R = Wph(1:nph,1:nph);
pHModel.P = Wph(1:nph,nph+1:end);
pHModel.S = Wph(nph+1:end,nph+1:end);
pHModel.Q = inv(pHModel.E);

end
