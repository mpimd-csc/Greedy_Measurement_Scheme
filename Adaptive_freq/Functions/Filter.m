function [ FiltVal ] = Filter(x0,bw,filtw,Filter_type)
% -------------------------------------------------------------------------
% This function outputs the discretized notch filter centered at x0  
%
% -------------------------------------------------------------------------
% The function inputs are:
% x0            --  Frequency at which the notch filter is centered,
% bw            --  Bandwidth of the filter,
% filtw         --  Array of frequencies where the filter is evaluated,
% Filter_type   --  chooe between two implementation of the notch filter:
%                  - 'Fbell': a bell shaped filter.
%                  - 'notch': based on built-in MATLAB function 'irrnotch'.
% -------------------------------------------------------------------------
%
% The function outputs are:
% FiltVal       -- Filter values in the range filtw
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


%%
if Filter_type == 'Fbell'
    esp = 10^-15;
% Construct the bell shaped filter
Filt  =  @(x) 1-exp(-bw*(log(abs(x)+esp)-log(abs(x0)+esp)).^2);

% evaluate the filter at the required frequencies 
FiltVal = Filt(filtw);


elseif Filter_type == 'notch'
% set the notch frequency
wo  = (x0-1)/Ntotal ; 

% construct the filter
[num,den]   = iirnotch(wo,bw);

% evaluate the filter at the required frequencies 
[h,ang_vec] = freqz(num,den,Ntotal);
FiltVal = abs(h)';

end

end

