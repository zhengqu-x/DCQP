function [corrected_lb,barz] = compute_lower_bound( ...
    Q,d,A,b,Aeq,beq,tolerance,m,n,tstar,verbose,mosek_quiet)
%COMPUTE_LOWER_BOUND Solve and correct the retained-region DNN bound.

[raw_lb,status,S,result]=lower_bound_dnn( ...
    Q,d,A,b,Aeq,beq,tolerance,m,n,mosek_quiet);
if status~=1
    error('qpsolver:DNNLowerBound', ...
        ['DNN lower bound was not solved successfully. Consider ', ...
         'increasing the MOSEK tolerance and retrying.']);
end

delta=min(eig(S));
if verbose
    fprintf('  raw SDP lb = %.4e, min_eig(S) = %.2e, tstar = %.2e\n', ...
        raw_lb,delta,tstar);
end
corrected_lb=raw_lb+delta+tstar^2*min(delta,0);
u=result.sol.itr.doty;
U=sMat(u,n+1);
barz=U(end,1:end-1)';
end
