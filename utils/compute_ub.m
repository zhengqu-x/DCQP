
function [ub,sol]=compute_ub(Q,d,A_bar,b_bar,Aeq,beq,n,p,sol,nb_rounds)

% ======================================================================= %
% Compute the value of ub
%
% INPUT
% Q, d               Coefficients of the quadratic objective function
% A_bar, b_bar       Inequality constraints: A_bar x <= b_bar
% Aeq, beq           Equality constraints: Aeq x = beq
% n                  Number of decision variables
% p                  Struct of algorithmic parameters:
%                       p.eps_dcdecomposition : Threshold for DC decomposition
%                       p.eps_checkpsd        : Tolerance for PSD check
%                       p.metgqp              : QP solver method (for gurobiqp)
%                       p.tol_gqp             : Tolerance for QP solver
% sol               Initial feasible solution (can be empty)
% nb_rounds         Number of random initializations for refinement
%
% OUTPUT
% ub                Upper bound on the objective value over the feasible set
% sol               Best feasible solution corresponding to ub
% ======================================================================= %


    ub=Inf;
    % Perform DC decomposition
    [M,N] = DC_decomposition(Q,p.eps_dcdecomposition);
    
if ~isempty(sol)
 
    [x_kkt] = search_of_kkt_point(Q,d,A_bar,b_bar,Aeq,beq,M,N,sol,p.eps_checkpsd,p.metgqp,p.tol_gqp,n);
    ub=x_kkt'*Q*x_kkt+2*d'*x_kkt;
end
    
    for i=1:min(n,nb_rounds)
   
    vt=randn(n,1);
    vt=vt/norm(vt);
    x0 = gurobiqp(eye(n),-vt,A_bar,b_bar,Aeq,beq,p.metgqp,p.tol_gqp,n);

    [x_kkt] = search_of_kkt_point(Q,d,A_bar,b_bar,Aeq,beq,M,N,x0,p.eps_checkpsd,p.metgqp,p.tol_gqp,n);
    v_kkt=x_kkt'*Q*x_kkt+2*d'*x_kkt;
    if v_kkt<ub
        sol=x_kkt;
    end
    ub=min(ub,x_kkt'*Q*x_kkt+2*d'*x_kkt);
    end
end