function [ I ] = mylogspaceInteger( N, m )
% mylogspaceInteger Generates m (almost) log distributed integer points in 1:N
% In
%   N       ...     number of integer points
%   m       ...     number of points to select
% Out
%   I       ...     selected points

%%
assert(N >= m);

%intialize I
I = unique(ceil(logspace(0, log(N)/log(10), m)));
I(I > N) = N;
I = unique(I);

%extend I until it reaches length m
c = 1;
while(length(I) < m)
    I = unique(ceil(logspace(0, log(N)/log(10), m + c)));
    I(I > N) = N;
    I = unique(I);
    c = c + 1;
end

%ensure I has length m
I = I(1:m);

end

