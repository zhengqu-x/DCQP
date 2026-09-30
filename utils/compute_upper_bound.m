function [barx,v] = compute_upper_bound( ...
    Q,d,A,b,Aeq,beq,M,N,barz,best_sol,bestub,iteration,known_start,p)
%COMPUTE_UPPER_BOUND Reproduce the stored candidate-selection policy.

n=size(A,2);
x0=barz;
if max(A*barz-b)>1e-9 || ...
        (~isempty(beq) && norm(Aeq*barz-beq)>1e-9)
    x0=gurobiqp(eye(n),-barz,A,b,Aeq,beq, ...
        p.gurobi_qp_method,p.gurobi_qp_tolerance,n);
end

if iteration==1 && known_start
    barx=best_sol;
else
    barx=search_of_kkt_point(Q,d,A,b,Aeq,beq,M,N,x0, ...
        p.psd_check_tolerance,p.gurobi_qp_method, ...
        p.gurobi_qp_tolerance,n);
end
v=barx'*Q*barx+2*d'*barx;

if iteration==1 && ~known_start && v>bestub
    x0_best=best_sol;
    if max(A*x0_best-b)>1e-9 || ...
            (~isempty(beq) && norm(Aeq*x0_best-beq)>1e-9)
        x0_best=gurobiqp(eye(n),-x0_best,A,b,Aeq,beq, ...
            p.gurobi_qp_method,p.gurobi_qp_tolerance,n);
    end
    alternative=search_of_kkt_point(Q,d,A,b,Aeq,beq,M,N,x0_best, ...
        p.psd_check_tolerance,p.gurobi_qp_method, ...
        p.gurobi_qp_tolerance,n);
    alternative_value=alternative'*Q*alternative+2*d'*alternative;
    if alternative_value<=v
        barx=alternative;
        v=alternative_value;
    end
end
end
