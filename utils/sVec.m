function [X_vec] = sVec(X)

% ======================================================================= %
% Vectorize a symmetric matrix X=(X_ij)_{i,j \in [n]} as
% (X_11, sqrt(2)X_21, ..., sqrt(2)X_n1, X_22, sqrt(2)X_32,...,X_nn)'
%
% INPUT
% X                  Symmetric matrix X
%
% OUTPUT             
% X_vec              (X_11, sqrt(2)X_21, ..., sqrt(2)X_n1, X_22,
%                     sqrt(2)X_32,...,X_nn)'
% ======================================================================= %


n=size(X,1);
idx=1;
k=n-1;
X_vec=zeros(0.5*n*(n+1),1);
for i = 1:n-1
    X_vec(idx)=X(i,i);
    X_vec(idx+1:idx+k)=X(i+1:i+k,i)*sqrt(2);
    idx=idx+k+1;
    k=k-1;
end

X_vec(end)=X(end,end);


end