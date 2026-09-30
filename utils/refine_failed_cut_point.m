function [x,info] = refine_failed_cut_point(Q,d,A,b,Aeq,beq,M,N,barx,beta,p,source,force_steps)
%REFINE_FAILED_CUT_POINT Failure-only recovery; normal KKT search is unchanged.
% rho=max(0,(Q*barx+d)'*barx-min_{y in P}(Q*barx+d)'*y).
% By default, refine only if rho>=beta (SDP failure recovery). Konno recovery
% passes source='Konno', force_steps=true: record KKT diagnostics, but run
% exactly 100 DC updates regardless of rho or a failed diagnostic LP.
% The caller must independently validate cuts at the refined anchor.
if nargin<12, source='SDP'; end
if nargin<13, force_steps=false; end
x=barx;
info=struct('attempted',false,'refined',false,'status','not_checked', ...
    'steps',0,'rho_before',NaN,'rho_after',NaN,'beta_before',beta, ...
    'value_before',barx'*Q*barx+2*d'*barx,'value_after',NaN, ...
    'anchor_before',barx,'anchor_after',[],'message','', ...
    'source',source,'force_steps',force_steps, ...
    'primal_violation_before',constraint_violation(barx,A,b,Aeq,beq), ...
    'primal_violation_after',NaN,'kkt_check_before_error','','kkt_check_after_error','');
try
    try
        info.rho_before=stationarity_gap(Q*barx+d,barx,A,b,Aeq,beq,p);
    catch ME
        info.kkt_check_before_error=ME.message;
        if ~force_steps, rethrow(ME); end
        if p.verbose, fprintf('%s KKT diagnostic failed: %s; continuing DC recovery.\n',source,ME.message); end
    end
    if p.verbose
        fprintf('%s recovery check: rho=%.2e, beta=%.2e, rho/beta=%.2e, primal violation=%.2e\n', ...
            source,info.rho_before,beta,info.rho_before/beta,info.primal_violation_before);
    end
    if ~force_steps && info.rho_before<beta
        info.status='rho_below_beta';
        info.message='No extra DC steps: rho is below beta.';
        return
    end
    info.attempted=true;
    if p.verbose
        fprintf('%s recovery: running 100 additional DC steps (unconditional=%d).\n',source,force_steps);
    end
    for k=1:100
        candidate=gurobiqp(0.5*M,0.5*(d-N*x),A,b,Aeq,beq, ...
            p.gurobi_qp_method,p.gurobi_qp_tolerance,numel(x));
        if numel(candidate)~=numel(x) || any(~isfinite(candidate))
            error('dcqp:RecoveryPoint','A DC recovery subproblem returned an invalid point.');
        end
        violation=constraint_violation(candidate,A,b,Aeq,beq);
        if violation>max(1e-8,10*p.gurobi_qp_tolerance)
            error('dcqp:RecoveryPoint','A DC recovery subproblem returned an infeasible point.');
        end
        x=candidate; info.steps=k;
    end
    info.primal_violation_after=constraint_violation(x,A,b,Aeq,beq);
    try
        info.rho_after=stationarity_gap(Q*x+d,x,A,b,Aeq,beq,p);
    catch ME
        info.kkt_check_after_error=ME.message;
        if ~force_steps, rethrow(ME); end
        if p.verbose, fprintf('%s KKT diagnostic failed after refinement: %s\n',source,ME.message); end
    end
    info.value_after=x'*Q*x+2*d'*x;
    info.anchor_after=x;
    info.refined=true; info.status='refined';
    if p.verbose
        fprintf('%s recovery complete: steps=%d, rho=%.2e, v=%.2e, movement=%.2e, primal violation=%.2e\n', ...
            source,info.steps,info.rho_after,info.value_after,norm(x-barx),info.primal_violation_after);
    end
catch ME
    x=barx; info.status='recovery_failed'; info.message=ME.message;
    if p.verbose, fprintf('%s recovery failed: %s\n',source,ME.message); end
end
end

function violation=constraint_violation(x,A,b,Aeq,beq)
violation=max([0;A*x-b]);
if ~isempty(Aeq), violation=max([violation;abs(Aeq*x-beq)]); end
end

function rho=stationarity_gap(g,x,A,b,Aeq,beq,p)
model=struct('obj',full(g(:)),'A',[sparse(A);sparse(Aeq)], ...
    'rhs',[b;beq],'sense',[repmat('<',size(A,1),1);repmat('=',size(Aeq,1),1)], ...
    'lb',-Inf(numel(x),1),'modelsense','min');
settings=struct('Method',p.gurobi_lp_method,'OutputFlag',0, ...
    'FeasibilityTol',p.gurobi_lp_tolerance, ...
    'OptimalityTol',p.gurobi_lp_tolerance);
result=gurobi(model,settings);
if strcmp(result.status,'INF_OR_UNBD')
    settings.DualReductions=0; result=gurobi(model,settings);
end
if ~strcmp(result.status,'OPTIMAL') || ~isfield(result,'objval')
    error('dcqp:RecoveryLP','Stationarity LP returned %s.',result.status);
end
rho=max(0,g'*x-result.objval);
end
