function [tROM,Hdata2,sInterpIndx] = time_DiscreteModel( uInF, yT, s, sInterpIndx,newindices,Hdata, params)
% -------------------------------------------------------------------------
% This function computes a discrete time realization using time domain data
%
% -------------------------------------------------------------------------
% Inputs:
% uInF                 ---  input vector,
% yT                   ---  output vector,
% s                    ---  sampling points,
% sInterpIndx          ---  interpolation points,
% newindicies          ---  new interpolation points,
% Hdata                ---  transfer function evaluations,
% params.kminPercent   ---  percent of points used in reg problem,
% params.dim           ---  proposed dimension of the learned model.
%
% Remark: The order of the resulting pH system maybe less than params.dim if the order 
%  of the minimal realization of the pH system is lower than params.dim.
% -------------------------------------------------------------------------
%
% Outputs:
% tROM                  --- Discrete time realization
% Hdata2                --- updated set of transfer function evaluations
% sInterpIndx           ---  updated set of interpolation points,
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

% default parameters
if  nargin < 6
    params.kminPercent=75;
    params.dim=20;
end

%% interpolate frequency domain data from time domain data
% frequency domain data Hdata corresponding to interpolation points
% sInterpIndx is computed in a least squares problem. RegCond are the
% condition numbers.
[HdataNew, ~, regCond] = regProb(params.kminPercent, uInF, yT, s, newindices);

%% create Loewner (take care of z Transform (discrete!))
%params.dim is the dimension proposed by the user, however this value may
%change if the system is not minimal
Hdata2=[Hdata ;HdataNew];
sInterpIndx=[sInterpIndx newindices];

%sort
[Hdata2,sortIdx] = sort(Hdata2);
sInterpIndx = sInterpIndx(sortIdx);

tROM = createLoewner(Hdata2, exp(s(sInterpIndx)), params.dim);

tROM.H= @(z) tROM.Cr*((z*tROM.Er-tROM.Ar)\tROM.Br); 


% store auxiliary data
tROM.sInterpIndx = sInterpIndx;
tROM.sInterp = s(sInterpIndx);
tROM.uInFatInterp = uInF(sInterpIndx);
tROM.regCond = regCond;

end
