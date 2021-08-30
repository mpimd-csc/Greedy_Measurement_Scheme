function H2NormVal=h2Norm(out,origsys,Maxfreq,Minfreq)
% -------------------------------------------------------------------------
% This function computes the H2 norm error between the orignal system and the idetnfied system
% as part of the paper:
% Cherifi K., Goyal, P., and Benner, P., Adaptive selection of 
% interpolation points for data driven modelling, 2020.
%
% -------------------------------------------------------------------------
% Inputs:
% out                  ---  Identified model,
% origsys              ---  Original system,
% Maxfreq              ---  Maximum frequency to be considered,
% Minfreq              ---  Minimum frequency to be considered,
%
% Remark: The order of the resulting pH system maybe less than params.dim if the order 
%  of the minimal realization of the pH system is lower than params.dim.
% -------------------------------------------------------------------------
%
% Outputs:
% H2NormVal            --- H2 norm error
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


Aequi1=out.LoewModel.E\out.LoewModel.A;
Bequi1=out.LoewModel.E\out.LoewModel.B;
Al=[origsys.A zeros(size(origsys.A,1),size(Aequi1,2));zeros(size(Aequi1,1),size(origsys.A,2)) Aequi1];
El=eye(size(Al));
Bl=[origsys.B;Bequi1];
Cl=[origsys.C -out.LoewModel.C];

%% Define options for the algorithm.
opts = struct( ...
    'FctUpdate' , 1, ...
    'Freqs'     , [10^Minfreq, 10^Maxfreq], ...
    'Info'      , 2, ...
    'MaxIter'   , 100, ...
    'MinIter'   , 10, ...
    'ModGramian', 0, ...
    'Npts'      , 1601, ...
    'Shifts'    , 'imaginary', ...
    'Solver'    , 'logm', ...
    'SolverFreq', 1, ...
    'StoreFacE' , 0, ...
    'StoreV'    , 0, ...
    'TolComp'   , eps, ...
    'TolLyap'   , 1.0e-12, ...
    'TolRHS'    , 1.0e-14, ...
    'TrueRes'   , 1, ...
    'L'         , [], ...
    'pL'        , []);


%% Compute the limited controllability Gramian.
fprintf(1, 'Solve the limited controllability Gramian:\n');
[ZC, YC, infoC] = freq_lyap_rksm(Al, Bl, El, opts);

F      = real(1i / pi * logm(full((Al + 1i * opts.Freqs(1) * El) ...
    \ (Al + 1i * opts.Freqs(2) * El)))) / El;
Bomega = El * F * Bl;
P      = ZC*YC*ZC';

H2NormVal=trace(Cl*P*Cl');


%% Compute the limited observability Gramian.
fprintf(1, 'Solve the limited observability Gramian:\n');
[ZO, YO, infoO] = freq_lyap_rksm(Al', Cl', El', opts);

F      = real(1i / pi * logm(full((Al + 1i * opts.Freqs(1) * El) ...
    \ (Al + 1i * opts.Freqs(2) * El)))) / El;
Comega = Cl * F * El;
Q      = ZO*YO*ZO';

H2norm2=trace(Bl'*Q*Bl); % this value should be equal to H2NormVal

end
