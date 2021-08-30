function [IdenAModel, out] = Construct_aModel(w,F,D, tol)
% -------------------------------------------------------------------------
% This function identifies a realization using frequency
% response data. 
%
% -------------------------------------------------------------------------
% The function inputs are:
% w                --  Frequencies on the j-omega axis, omega > 0,
% F                --  Transfer function values at the frequencies w,
% D                --  Direct feed-through term
% tol (optional)   --  tolerance for the Loewner pencil that determines the order of
%           the identified state-space model. Its default
%           value is 1e-8.
% -------------------------------------------------------------------------
%
% The function outputs are:
% IdenaModel  -- Identified Loewner model, the system matrices and transfer
%                   function
% Out         -- It stores additional information.
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

out = [];
if  nargin < 4
    tol = 1e-8;
end
if nargin < 3
    error('Error! The number of inputs should be 3. Please provide frequencies transfer function values at these frequencies, and the direct feed-through term.')
end
%
Np = length(w);
m  = size(D,1);
[L,sL,~,~,V,W,TL,TR]  = loewner(1i*w,F);

% Matching at infinity and transformation into real realization.
L   = real(TL*L*TR);             
sL  = real(TL*sL*TR - TL*ones(Np,1)*D*ones(1,Np)*TR);     
W    = real(W*TR - ones(1,Np)*D*TR);
V    = real(TL*V - TL*D*ones(Np,1));

% Compression of the realization and construction of minimal realization
[Y,svL, ~]   = svd([L sL]);
[~,~,X]      = svd([L; sL]);
svL          = diag(svL);
r            = sum((svL./svL(1))>tol);
Yr           = Y(:,1:r);
Xr           = X(:,1:r);

% State space model;
El = -Yr'*L*Xr;
Al = -Yr'*sL*Xr;
Bl = Yr'*V;
Cl = W*Xr;

out.LoewModel.E     = El; 
out.LoewModel.A     = Al; 
out.LoewModel.B     = Bl; 
out.LoewModel.C     = Cl;
out.LoewModel.D     = D;
out.LoewModel.TF    = @(s) Cl*((s*El-Al)\Bl) + D;

IdenAModel=out.LoewModel;
