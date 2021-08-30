function [ Ad, Bd, Cd, Ed ] = bwdEuler_sys( t, A, B, C, E )
% bwdEuler is a Discretization procedure based on Backward Euler
% In
%   A,B,C,E          ... continous time system
%   u                ... input vector
%   t                ... time samples
% Out
%   Ad,Bd,Cd,Ed      ... Discrete time system
%   y                ... output vector

%%
% initialization
deltaT = t(2) - t(1);

%{
sys=dss(A,B,C,0,E);
sysd = c2d(sys,deltaT,'tustin');
[Ad,Bd,Cd,Dd,Ed] = dssdata(sysd);
%}


%create discrete system (BwdEuler)
Ed = (E - deltaT*A);
Ad = E;
Bd = deltaT*B;
Cd = C;
%}

%{
% create discrete system (implicit midpoint)
Ed = (E - deltaT*A/2);
Ad = (E+deltaT*A/2);
Bd = deltaT*B;
Cd = C;
%}

% matched
%sysd = c2d(dss(A,B,C,0,E),deltaT,'matched');
%sysd=dss(A,B,C,0,E,deltaT);
%[Ad,Bd,Cd,Dd,Ed]=dssdata(sysd);


end

