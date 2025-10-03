function [x,fval,exitflag] = gurobilp(f,A,b,Aeq,beq,lb,ub,met_glp,tol_glp)

% ======================================================================= %
% Call Gurobi to solve the following linear program:
%     min_x   f'*x
%     s.t.    A*x <= b
%            Aeq*x = beq
%            lb <= x <= ub.
% Most of the codes are copied from linprog.m: 
% https://www.gurobi.com/documentation/current/examples/linprog_m.html
%
% INPUT
% f                       Coefficient vector of the linear objective
% A, b                    Inequality constraints Ax <= b
% Aeq, beq                Equality constraints Aeqx = beq
% lb, ub                  Variable bounds: lb <= x <= ub
% met_glp                 Gurobi LP solution method:
%                             0: Primal simplex
%                             1: Dual simplex
%                             2: Barrier
% tol_glp                 Numerical tolerance for Gurobi LP solver
%
% OUTPUT
% x                       Gurobi optimal solution x
% fval                    Gurobi optimal value
% exitflag                1: solved to optimality (subject to tolerances), 
%                             and an optimal solution is available
%                        -1: infeasible
%                        -2: unbounded
% ======================================================================= %

% Build Gurobi model
model.obj = f;
model.A = [sparse(A); sparse(Aeq)]; % A must be sparse
model.sense = [repmat('<',size(A,1),1); repmat('=',size(Aeq,1),1)];
model.rhs = [b; beq];
model.modelsense = 'min';

% Variable types: 'C' for continuous variable
model.vtype = 'C';

% Intermediate variable
if ~isempty(lb)
    model.lb = lb;
else
    n = size(model.A,2);
    model.lb = -Inf(n,1); % Default lb for MATLAB is -Inf
end

if ~isempty(ub)
    model.ub = ub;
else
    n = size(model.A,2);
    model.ub = Inf(n,1); % Default ub for MATLAB is Inf
end


% Build Gurobi parameter
params.OutputFlag = 0; 
params.Method = met_glp;

% Tolerance setting
params.OptimalityTol = tol_glp; 
params.LPWarmStart = 2; 

% Solve model with Gurobi
result = gurobi(model, params);

% Resolve model if status is INF_OR_UNBD
if strcmp(result.status,'INF_OR_UNBD')
    params.DualReductions = 0;
    result = gurobi(model,params);
end

% Collect results
x = [];

if isfield(result,'x')
x = result.x;
end

fval = [];

if isfield(result,'objval')
fval = result.objval;
end

if strcmp(result.status,'OPTIMAL')
    exitflag = 1; 
elseif strcmp(result.status,'INFEASIBLE')
    exitflag = -1; 
elseif strcmp(result.status,'UNBOUNDED')
    exitflag = -2; 
end

end