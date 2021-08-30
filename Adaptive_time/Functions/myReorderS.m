function [ I ] = myReorderS( Isorted )
% myReorderS reorders the interpolation points such that
% complex conjugates follow each other:
% from
% a - ib
% c - id
% a + ib
% c + id
% to
% a + ib
% a - ib
% c + id
% c - id

%%
%intialization
m = length(Isorted);
I = zeros(size(Isorted));

% reorder
for i=1:m/2
    I(2*i - 1) = Isorted(m/2 + i);
    I(2*i) = Isorted(m/2 - i + 1);
end

end

