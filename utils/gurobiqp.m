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

% Tolerance setting

params.FeasibilityTol = tol_gqp; 
params.OptimalityTol = tol_gqp; 
params.BarConvTol = tol_gqp; 
params.LPWarmStart = 2; 

% Solve model with Gurobi
result = gurobi(model, params);
x=[];


if strcmp(result.status,'OPTIMAL')
    exitflag = 1; 
    x = result.x;
    v=x'*Q*x+2*d'*x;
    if isempty(Aeq) && isempty(beq)
        if v>0 && min(b)>=0
            x=zeros(n,1);
        end
    end
elseif strcmp(result.status,'INFEASIBLE')
    exitflag = -1; % Infeasible
elseif strcmp(result.status,'UNBOUNDED')
    exitflag = -2; % Unbounded
end





end

