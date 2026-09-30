
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
%                       p.dc_regularization   : Threshold for DC decomposition
%                       p.psd_check_tolerance : Tolerance for PSD check
%                       p.gurobi_qp_method    : QP solver method
%                       p.gurobi_qp_tolerance : QP solver tolerance
% sol               Initial feasible solution (can be empty)
% nb_rounds         Number of random initializations for refinement
%
% OUTPUT
% ub                Upper bound on the objective value over the feasible set
% sol               Best feasible solution corresponding to ub
% ======================================================================= %


    ub=Inf;
    % Perform DC decomposition
    [M,N] = DC_decomposition(Q,p.dc_regularization);
    
if ~isempty(sol)
 
    [x_kkt] = search_of_kkt_point(Q,d,A_bar,b_bar,Aeq,beq,M,N,sol, ...
        p.psd_check_tolerance,p.gurobi_qp_method,p.gurobi_qp_tolerance,n);
    ub=x_kkt'*Q*x_kkt+2*d'*x_kkt;
end
    
    for i=1:min(n,nb_rounds)
   
    vt=randn(n,1);
    vt=vt/norm(vt);
    x0 = gurobiqp(eye(n),-vt,A_bar,b_bar,Aeq,beq, ...
        p.gurobi_qp_method,p.gurobi_qp_tolerance,n);

    [x_kkt] = search_of_kkt_point(Q,d,A_bar,b_bar,Aeq,beq,M,N,x0, ...
        p.psd_check_tolerance,p.gurobi_qp_method,p.gurobi_qp_tolerance,n);
    v_kkt=x_kkt'*Q*x_kkt+2*d'*x_kkt;
    if v_kkt<ub
        sol=x_kkt;
    end
    ub=min(ub,x_kkt'*Q*x_kkt+2*d'*x_kkt);
    end
end
