function [result] = gurobiqp_with_time_limit(Q,d,A,b,Aeq,beq,tol_gqp,n,Tlimit)

% ======================================================================= %
% Call Gurobi to solve the following convex quadratic program:
% min  x'*Q*x + 2*d'*x
% s.t. Ax <= b
%      Aeq x = beq
%
% INPUT
% Q, d             Parameters of the quadratic objective function
% A, b             Parameters of the inequality constraint A x <= b
% Aeq, beq         Equality constraints: Aeq x = beq
% tol_gqp          Accuracy tolerance used in Gurobi
% n                Size(Q,1)
% Tlimit           Time limit
%
% OUTPUT
% result           Gurobi optimal solution
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
params.OutputFlag = 1; 


% Tolerance setting
params.OptimalityTol = tol_gqp; 
params.BarConvTol = tol_gqp; 
params.LPWarmStart = 2; 
params.TimeLimit = Tlimit;

% Solve model with Gurobi
result = gurobi(model, params);

% Resolve model if status is INF_OR_UNBD
if strcmp(result.status,'INF_OR_UNBD')
   params.DualReductions = 0;
   result = gurobi(model,params);
end

end

