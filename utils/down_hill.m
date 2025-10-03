function [x_bar] = down_hill(Q,d,B,A,b,x_hat,I_x_hat_c)

% ======================================================================= %
% Find a point x_bar \in {x: Ax <= b} such that 
% \Phi(x_bar) <= \Phi(x_hat),
% and I_x_hat is a proper subset of I_x_bar.
%
% INPUT
% Q, d             Parameters of the quadratic objective function
% B                Basis matrix
% A, b             Parameters of the inequality constraint A x <= b
% x_hat            x_hat from Algorithm 1
% I_x_hat_c        {i: a_i^Tx_hat < b_i}
%
% OUTPUT
% x_bar            A feasible point x_bar satisfying 
%                  \Phi(x_bar) <= \Phi(x_hat),
%                  and I_x_hat is a proper subset of I_x_bar.
% ======================================================================= %





[V, D] = eig(B'*Q*B);


h_x_hat = B*V(:, find(diag(D) <= 0, 1));



if (Q*x_hat + d)'*h_x_hat > 0 
    h_x_hat = -h_x_hat;
end


alpha = Inf;

for i = I_x_hat_c'
    denominator = A(i,:)*h_x_hat;
    if denominator > 0 
        alpha = min(alpha, (b(i)-A(i,:)*x_hat)/denominator);
    end
end

if alpha==Inf
    error('The feasible region may be unbounded.\n');
end

x_bar = x_hat + alpha * h_x_hat;




end