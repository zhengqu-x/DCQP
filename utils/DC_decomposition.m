function [M,N] = DC_decomposition(Q,spn)

% ======================================================================= %
% Find a DC decomposition of matrix Q
%
% INPUT
% Q                Symmetric matrix to be decomposed
% spn              Nonnegative scalar parameter
%
% OUTPUT
% M N              Positive semidefinite matrices such that M - N = Q
% ======================================================================= %

n = size(Q, 1);

[V, D] = eig(Q);
D_v = diag(D);

M = V * diag(D_v.*(D_v>0)+spn*ones(n,1)) * V';
N = M-Q;

end