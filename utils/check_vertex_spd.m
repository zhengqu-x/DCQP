function [vertex_spd_flag,vertex_flag,I_x_c,B,BH] = check_vertex_spd(Q,A,b,Aeq,x,epsilon)

% ======================================================================= %
% Check whether the point x is a vertex, or whether Q|_H_I_x > 0
%
% INPUT
% Q              Matrix parameter in the quadratic objective function
% A, b           Constraint parameters in Ax <= b
% x              Candidate point to be tested
% eps            Numerical tolerance parameter (epsilon > 0)
%
% OUTPUT
% vertex_spd_flag   Indicator:
%                   0: x is not a vertex and Q|_{H_{I_x}} is not positive definite
%                   1: x is a vertex, or Q|_H_I_x > 0
% vertex_flag       Indicator:
%                   0: x is not a vertex
%                   1: x is a vertex
% I_x_c             Index set of inactive constraints {i : a_i^T x < b_i}
% B                 Null space basis of A_{I_x} null(A_I_x)
% BH                B^TQB
% ======================================================================= %

vertex_flag = 0;
spd_flag = 0;
BH=[];
B=[];

% Intermediate variable

Ax = A * x;
m  = size(A,1);
n=size(A,2);

I_x = find(abs(b - Ax) <= 1e-8); % {i: a_i^Tx = b_i}

I_x_c = setdiff(1:m, I_x)';


% Intermediate variable

A_I_x = [A(I_x,:);Aeq];

if isempty(A_I_x)
    B=eye(n);
else
    if rank(A_I_x) == size(A, 2)
        vertex_flag = 1;
    else
        % Check whether Q|_H_I_x > 0
        B = null(A_I_x);
        BH=B'*Q*B;
        eigBH=eig(BH);
        if min(eigBH)>-epsilon     
            spd_flag = 1;
            %if min(eigBH)<epsilon
                [V,D]=eig(BH);
                dd=diag(D);
                BH=V*diag(dd.*(dd>epsilon))*V';
            %end
        end
    end
end


vertex_spd_flag = vertex_flag || spd_flag;

end