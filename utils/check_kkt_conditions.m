function [opt_val,vbasis,lambda] = check_kkt_conditions(Q,d,A,x,I_x,tol_glp,m,n)

% ======================================================================= %
% Check whether a given point x satisfies the KKT conditions for: 
%  min_x    x^T Q x + 2 d^T x
%  s.t.    A x <= b
%
% INPUT
% Q, d        Coefficients of the quadratic objective function
% A           Constraint matrix in A x <= b
% x           Candidate point to be tested
% I_x         Active set { i : a_i^T x = b_i }
% tol_glp     Tolerance parameter used in gurobilp
% m           Number of constraints, size(A,1)
% n           Number of variables, size(A,2)
%
% OUTPUT
% opt_val     Optimal objective value
% vbasis      Basis indices returned by solver
% lambda      Optimal Lagrange multipliers
% ======================================================================= %

v = ones(m, 1);
v(I_x) = 0;

% Build Gurobi model
model.obj = v;
model.A = sparse(A'); % A must be sparse
model.sense = repmat('=',n,1);
model.rhs = -(Q*x + d); % rhs must be dense
model.modelsense = 'min';

% Variable types: 'C' for continuous variable
model.vtype = 'C';
model.lb = zeros(m,1);
model.ub = Inf(m,1);

% Build Gurobi parameter
params.OutputFlag = 0; 
params.Method = 0; 
% Tolerance setting
params.FeasibilityTol = tol_glp; 
params.OptimalityTol = tol_glp; 
params.LPWarmStart = 2; 

% Solve model with Gurobi
result = gurobi(model, params);

% Resolve model if status is INF_OR_UNBD
if strcmp(result.status,'INF_OR_UNBD')
    params.DualReductions = 0;
    result = gurobi(model,params);
end

lambda = result.x;
opt_val = result.objval;
vbasis = result.vbasis;

end

