function [x,exitflag] = gurobiqp(Q,d,A,b,Aeq,beq,met_gqp,tol_gqp,n)

% ======================================================================= %
% Call Gurobi to solve the following convex quadratic program:
% min  x'*Q*x + 2*d'*x
% s.t. A*x   <= b
%      Aeq*x  = beq
%
% INPUT
%
% Q, d                    Parameters of the quadratic objective function
% A, b                    Inequality constraints: Ax <= b
% Aeq, beq                Equality constraints: Aeq x = beq
% met_gqp                 Gurobi method (-1=auto, 0=primal simplex,
%                         1=dual simplex, 2=barrier)
% tol_gqp                 Accuracy tolerance used in gurobiqp
% n                       Size(Q,1)
%
% OUTPUT
% x            Gurobi optimal solution
% exitflag     Status of Gurobi
% ======================================================================= %



% Build Gurobi model

model.Q = sparse(Q);
model.obj = 2*d;
model.modelsense = 'min';

model.A = [sparse(A);sparse(Aeq)];
model.rhs = [b;beq];
model.sense = [repmat('<',size(A,1),1); repmat('=',size(Aeq,1),1)];

model.lb = -Inf(n,1);
model.ub = -model.lb;

% Variable types: 'C' for continuous variable
model.vtype = 'C';

% Build Gurobi parameter
params.OutputFlag = 0; 
params.Method = met_gqp;

% Tolerance setting

params.OptimalityTol = tol_gqp; 
params.BarConvTol = tol_gqp; 
params.LPWarmStart = 2; 
% These auxiliary QPs must not silently enter nonconvex branch-and-bound.
params.NonConvex = 0;

% Solve model with Gurobi
original_model = model;
psd_shift = 0;
[result,model,psd_shift] = solve_convex_model(model,params,psd_shift);
if any(strcmp(result.status,{'SUBOPTIMAL','NUMERIC'}))
    % Retry once with simplex, retaining the requested tolerances.
    % If dual simplex was already selected, try primal simplex instead.
    first_status = result.status;
    params.Method = 1;
    if met_gqp == 1
        params.Method = 0;
    end
    fprintf('gurobiqp: %s with Method=%g; retrying with simplex Method=%g.\n', ...
        first_status,met_gqp,params.Method);
    [result,~,psd_shift] = solve_convex_model(model,params,psd_shift);
    fprintf('gurobiqp: simplex retry returned %s.\n',result.status);
end
x=[];


if strcmp(result.status,'OPTIMAL')
    if psd_shift>0 && (~isfield(result,'x') || ~isreal(result.x) || ...
            numel(result.x)~=n || any(~isfinite(result.x(:))))
        error('dcqp:gurobiqpPSDValidation', ...
            'PSD retry lacks a finite primal solution for original-QP validation.');
    end
    exitflag = 1; 
    x = result.x;
    v=x'*Q*x+2*d'*x;
    if isempty(Aeq) && isempty(beq)
        if v>0 && min(b)>=0
            x=zeros(n,1);
        end
    end
    if psd_shift>0
        % Validate the actual returned point, including the zero-point fallback.
        validate_psd_retry(original_model,result,x,tol_gqp,psd_shift);
    end
elseif strcmp(result.status,'INFEASIBLE')
    exitflag = -1; % Infeasible
elseif strcmp(result.status,'UNBOUNDED')
    exitflag = -2; % Unbounded
else
    error('dcqp:gurobiqpStatus', ...
        'Gurobi QP returned %s with Method=%g; no optimal solution is available.', ...
        result.status,params.Method);
end
end

function [result,model,psd_shift] = solve_convex_model(model,params,psd_shift)
try
    result = gurobi(model,params);
catch ME
    % One retry only, and only for Gurobi's explicit Q_NOT_PSD error.
    if psd_shift>0 || isempty(regexp(ME.message,'Gurobi error 10020([^0-9]|$)','once'))
        rethrow(ME);
    end
    H = full((model.Q+model.Q')/2);
    if ~isreal(H) || any(~isfinite(H(:)))
        rethrow(ME);
    end
    scale = max(1,norm(H,2));
    if ~isfinite(scale)
        rethrow(ME);
    end
    min_eigenvalue = min(eig(H));
    if ~isfinite(min_eigenvalue) || min_eigenvalue < -100*eps*scale
        rethrow(ME);
    end
    psd_shift = 1e-12*scale;
    model.Q = sparse(H+psd_shift*eye(size(H)));
    fprintf('gurobiqp: PSD retry, min_eig=%.2e, diagonal_shift=%.2e.\n', ...
        min_eigenvalue,psd_shift);
    result = gurobi(model,params);
end
end

function validate_psd_retry(model,result,x,tol,psd_shift)
% All variables in this wrapper are free, so there are no bound multipliers.
if ~isfield(result,'pi') || ~isreal(x) || any(~isfinite(x(:))) || ...
        numel(x)~=numel(model.obj) || ~isreal(result.pi) || ...
        numel(result.pi)~=numel(model.rhs) || any(~isfinite(result.pi(:)))
    error('dcqp:gurobiqpPSDValidation', ...
        'PSD retry lacks a finite primal/dual solution for original-QP validation.');
end
x=x(:); pi=result.pi(:);
Ax=model.A*x; slack=model.rhs-Ax;
ineq=model.sense=='<'; eq=model.sense=='=';
gradient=(model.Q+model.Q')*x+model.obj;
primal=max([0;-slack(ineq);abs(slack(eq))]);
stationarity=norm(gradient-model.A'*pi,Inf);
dual_sign=max([0;pi(ineq)]);
complementarity=max([0;abs(pi(ineq).*slack(ineq))]);
original_objective=x'*model.Q*x+model.obj'*x;

% Keep the requested tolerances, allowing only arithmetic roundoff in the
% residual calculations. Complementarity uses objective units.
constraint_scale=max([1;abs(model.A)*abs(x)+abs(model.rhs)]);
gradient_scale=max([1;abs(model.Q+model.Q')*abs(x)+abs(model.obj)+ ...
    abs(model.A')*abs(pi)]);
complementarity_scale=max([1;abs(pi).*(abs(model.A)*abs(x)+abs(model.rhs))]);
primal_tol=tol+100*eps*constraint_scale;
stationarity_tol=tol+100*eps*gradient_scale;
dual_tol=tol+100*eps*max([1;abs(pi)]);
complementarity_tol=tol*max(1,abs(original_objective))+100*eps*complementarity_scale;
residuals=[primal,stationarity,dual_sign,complementarity];
limits=[primal_tol,stationarity_tol,dual_tol,complementarity_tol];
if any(~isfinite([residuals,limits,original_objective])) || any(residuals>limits)
    error('dcqp:gurobiqpPSDValidation', ...
        ['PSD retry failed original-QP validation: primal %.2e (tol %.2e), ' ...
         'stationarity %.2e (tol %.2e), dual sign %.2e (tol %.2e), ' ...
         'complementarity %.2e (tol %.2e).'], ...
        primal,primal_tol,stationarity,stationarity_tol,dual_sign,dual_tol, ...
        complementarity,complementarity_tol);
end
fprintf(['gurobiqp: PSD retry validated on original QP: primal=%.2e, ' ...
    'stationarity=%.2e, dual_sign=%.2e, complementarity=%.2e, ' ...
    'objective_perturbation=%.2e.\n'], ...
    primal,stationarity,dual_sign,complementarity,psd_shift*norm(x)^2);
end
