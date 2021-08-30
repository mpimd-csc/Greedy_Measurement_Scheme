function [ ROM ] = createLoewner( Hdata, s, dim )
%createLoewner Creates a Loewner ROM of dimension dim
% In
%   Hdata           ... transfer function data
%   s               ... interpolation points
%   dim             ... dimension of ROM
% Out
%   L, Ls           ... Loewner matrices
%   Ar, Br, Cr, Er  ... Loewner ROM

%% sanity checks
assert(length(Hdata) == length(s));
%assert(2*dim <= length(s));
assert(mod(length(s), 2) == 0);
M = length(s)/2;

%% create initial Loewner ROM
L = zeros(M, M);
Ls = zeros(M, M);
Br = zeros(M, 1);
Cr = zeros(1, M);
for i=1:M
    for k=1:M
        HwY = Hdata(i);
        HwdimY = Hdata(M + k);
        L(i, k) = -(HwY - HwdimY)/(s(i) - s(M + k));
        Ls(i, k) = -(s(i)*HwY - s(M + k)*HwdimY)/(s(i) - s(M + k));
    end
    Br(i) = Hdata(i);
    Cr(1, i) = Hdata(M + i);
end
%
%real
%L=ones(N,1);
%R=L.';

%compute the initial system
J=(1/sqrt(2))*[1 -1i;1 1i];
Jb=blkdiag(kron(eye(floor(length(L)/2)),J));
L=real(Jb'*L*Jb);
Ls=real(Jb'*Ls*Jb);
Br=real(Jb'*Br);
Cr=real(Cr*Jb);
Er = L;
Ar = Ls;

%% truncate to order dim
[U1,S1,V1]=svd([L Ls]);
[U2,S2,V2]=svd([L;Ls]);


%determine if the realization is minimal with the proposed dimension dim
tol= 1e-14;%1e-8;
svL          = diag(S1);
minr         = sum((svL./svL(1))>tol);
if minr<dim
    dim=minr;
    fprintf('Warning: the order of the realization was changed to %d \n',dim);
end


% Compression of the realization and construction of minimal realization
Y=U1(:,1:dim);
X=V2(:,1:dim);
% Projected model
Er=-Y'*Er*X;
Ar=-Y'*Ar*X;
Br=Y'*Br;
Cr=Cr*X;

%% assemble ROM
ROM.Er = Er;
ROM.Ar = Ar;
ROM.Br = Br;
ROM.Cr = Cr;
ROM.dim = dim;
ROM.Hdata = Hdata;
ROM.s = s;


end

