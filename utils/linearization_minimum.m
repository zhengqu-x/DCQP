function rho = linearization_minimum(Q,d,barx,A,b,Aeq,beq,method,tolerance)
%LINEARIZATION_MINIMUM Signed minimum first-order variation on the region.
% rho=min_{x in P}(Q*barx+d)'*(x-barx), which is nonpositive for a
% feasible barx and zero at an exact first-order KKT point.

g=Q*barx+d;
[~,minimum_value,status]=gurobilp(full(g),A,b,Aeq,beq,[],[],method,tolerance);
if status~=1 || isempty(minimum_value) || ~isfinite(minimum_value)
    error('linearization_minimum:LPSolve', ...
        'Could not minimize the objective linearization over the current region.');
end
rho=minimum_value-g'*barx;
end
