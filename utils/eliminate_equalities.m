function [Q_new, d_new, A_new, b_new, M, x0,  obj_constant] = eliminate_equalities(Q, d, A_bar, b_bar, Aeq, beq,n)
    
    
    % 1. QR with column pivoting for maximum numerical stability
    [~, R_qr, P] = qr(full(Aeq), 'vector');
    
    tol = max(size(Aeq)) * eps(norm(full(Aeq)));
    rank_Aeq = sum(abs(diag(R_qr)) > tol);
    
    if rank([full(Aeq), beq]) > rank_Aeq
        error('The equality constraints are structurally infeasible.');
    end
    
    % 2. Select variables
    idx_D = P(1:rank_Aeq);
    idx_I = P(rank_Aeq+1:end);
    
    num_ind = length(idx_I);
    
    % 3. Extract submatrices
    A_D = full(Aeq(:, idx_D));
    A_I = full(Aeq(:, idx_I));
    
    % 4. Compute T and v using robust left-division
    v = A_D \ beq;
    T = -(A_D \ A_I);
    
    % 5. Create the mapping x = M * x_I + x0
    M = zeros(n, num_ind);
    M(idx_I, :) = eye(num_ind);
    M(idx_D, :) = T;
    
    x0 = zeros(n, 1);
    x0(idx_D) = v;
    
    % 6. Transform Objective and Calculate Constant
    Q_new = M' * Q * M;
    d_new = M' * (Q * x0 + d);
    Q_new = 0.5 * (Q_new + Q_new'); % Symmetrize
    
    % The exact constant calculation you requested
    obj_constant = (x0' * Q * x0) + 2 * (d' * x0);
    
    % 7. Transform original inequalities 
    % (This automatically translates the [0,1] bounds since they are already in A_bar)
    if ~isempty(A_bar)
        A_new = A_bar * M;
        b_new = b_bar - A_bar * x0;
    else
        A_new = [];
        b_new = [];
    end
end