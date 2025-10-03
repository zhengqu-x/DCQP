function [X] = sMat(x,n)

% ======================================================================= %
% Construct a symmetric matrix X=(X_ij)_{i,j \in [n]} such that
% x=(X_11, sqrt(2)X_21, ..., sqrt(2)X_n1, X_22, sqrt(2)X_32,...,X_nn)'
%
% INPUT
% x                  A column vector of dimension n(n+1)/2
% n                  Dimension of X
%
% OUTPUT             
% X                  x=sVec(X)
% ======================================================================= %

X=zeros(n,n);

idx=1;
k=n-1;
for i = 1:n-1
    X(i,i)=x(idx);
    X(i+1:end,i)=x(idx+1:idx+k)/sqrt(2);
    X(i,i+1:end)=x(idx+1:idx+k)'/sqrt(2);
    idx=idx+k+1;
    k=k-1;
end

X(end,end)=x(end);



end