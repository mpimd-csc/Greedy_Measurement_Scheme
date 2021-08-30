function [ y, X ] = timeStep( Ad, Bd, Cd, Ed, u, xInit )
%TIMESTEP Simulation of the discrete LTI system
% In
%   Ad, Bd, Cd, Ed      ...     Discrete system
%   u                   ...     input vector
%   xInit               ...     initial condition
% Out
%   X                   ...     Snapshots
%   y                   ...     output vector

%%
%start simulation
Tsteps = length(u);
y = zeros(1, Tsteps);
[L1, L2] = lu(Ed);

% prepare storing snapshots
if(nargout > 1)
    X = zeros(size(Ad, 1), Tsteps);
    X(:, 1) = xInit;
    indx = 1:Tsteps;
else
    X = xInit;
    indx = ones(1, Tsteps);
end

% simulate
y(1) = Cd*X(:, 1);
for i=2:Tsteps
    X(:, indx(i)) = L2\(L1\(Ad*X(:, indx(i - 1)) + Bd*u(i - 1)));
    y(i) = Cd*X(:, indx(i));
end

end

