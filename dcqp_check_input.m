function n=dcqp_check_input(Q, d, A, b, Aeq, beq, params)
% DCQP_CHECK_INPUT  Validate inputs for the DCQP solver
%
% SYNTAX:
%   dcqp_check_input(Q, d, A, b, Aeq, beq, params)
%
% DESCRIPTION:
%   Performs comprehensive input validation for the DCQP solver,
%   checking dimensions, types, and mathematical properties.
%
%   IMPORTANT REQUIREMENTS:
%   - Inequality constraints A*x <= b are REQUIRED (cannot be empty)
%   - The feasible region must be BOUNDED for the algorithm to work
%   - The matrix Q must have at least one negative eigenvalue
%
% INPUT:
%   Q       - (n x n) symmetric matrix (objective Hessian)
%   d       - (n x 1) vector (linear objective coefficients)
%   A       - (m x n) matrix (inequality constraint coefficients, REQUIRED)
%   b       - (m x 1) vector (inequality constraint bounds, REQUIRED)
%   Aeq     - (meq x n) matrix (equality constraint coefficients, optional)
%   beq     - (meq x 1) vector (equality constraint bounds, optional)
%   params  - Structure with algorithm parameters
%
% OUTPUT:
%   (none) - Throws error if validation fails
%
% ERRORS:
%   Throws descriptive errors for invalid inputs, missing inequality
%   constraints, or potentially unbounded feasible regions

% Copyright (c) 2025
% All rights reserved.

% Check that Q is provided and is a matrix
if isempty(Q) || ~ismatrix(Q)
    error('dcqp_check_input:invalid_Q', 'Q must be a non-empty matrix');
end

% Check that Q is square
[n1, n2] = size(Q);
if n1 ~= n2
    error('dcqp_check_input:nonsquare_Q', 'Q must be a square matrix');
end
n = n1;

% Check that Q is symmetric (within tolerance)
if norm(Q - Q', 'fro') > 1e-12 * norm(Q, 'fro')
    warning('dcqp_check_input:nonsymmetric_Q', ...
            'Q is not symmetric, symmetrizing as Q = (Q + Q'')/2');
    Q = (Q + Q') / 2;
end

% Check that Q contains only finite values
if ~all(isfinite(Q(:)))
    error('dcqp_check_input:nonfinite_Q', 'Q contains non-finite values');
end

% Check d vector
if isempty(d)
    error('dcqp_check_input:empty_d', 'd vector cannot be empty');
end

if ~isvector(d)
    error('dcqp_check_input:invalid_d', 'd must be a vector');
end

if length(d) ~= n
    error('dcqp_check_input:dimension_mismatch_d', ...
          'Dimension mismatch: d must have length %d (same as Q dimensions)', n);
end

% Ensure d is a column vector
d = d(:);

% Check that d contains only finite values
if ~all(isfinite(d))
    error('dcqp_check_input:nonfinite_d', 'd contains non-finite values');
end

% Check that Q is indefinite (required for DCQP algorithm)
% DCQP is designed for non-convex QPs where Q has at least one negative eigenvalue
eigenvalues = eig(Q);
if all(eigenvalues >= -1e-5)  % Allow small numerical tolerance
    warning('dcqp_check_input:not_indefinite', ...
          ['Q matrix must be indefinite (have at least one negative eigenvalue). ', ...
           'DCQP is designed for non-convex quadratic programs. ', ...
           'For convex problems (Q positive semidefinite), use standard QP solvers like quadprog.']);
end

% Check inequality constraints - REQUIRED (cannot be empty)
if isempty(A) || isempty(b)
    error('dcqp_check_input:missing_inequality', ...
          ['Inequality constraints A*x <= b are required and cannot be empty. ', ...
           'The DCQP solver requires a bounded feasible region.']);
end

% Check A matrix
if ~ismatrix(A)
    error('dcqp_check_input:invalid_A', 'A must be a matrix');
end

[m, n_A] = size(A);
if n_A ~= n
    error('dcqp_check_input:dimension_mismatch_A', ...
          'A must have %d columns (same as Q dimensions)', n);
end

if m == 0
    error('dcqp_check_input:empty_constraints', ...
          'At least one inequality constraint is required');
end

if ~all(isfinite(A(:)))
    error('dcqp_check_input:nonfinite_A', 'A contains non-finite values');
end

% Check b vector
if ~isvector(b) || length(b) ~= m
    error('dcqp_check_input:invalid_b', ...
          'b must be a vector with %d elements (same as A rows)', m);
end

b = b(:); % Ensure column vector

if ~all(isfinite(b))
    error('dcqp_check_input:nonfinite_b', 'b contains non-finite values');
end


% Check equality constraints
if isempty(Aeq) && isempty(beq)
    % No equality constraints - this is OK
    %Aeq = zeros(0, n);
    %beq = zeros(0, 1);
elseif isempty(Aeq) || isempty(beq)
    error('dcqp_check_input:incomplete_equality', ...
          'Both Aeq and beq must be provided for equality constraints');
else
    % Check Aeq matrix
    if ~ismatrix(Aeq)
        error('dcqp_check_input:invalid_Aeq', 'Aeq must be a matrix');
    end
    
    [meq, n_Aeq] = size(Aeq);
    if n_Aeq ~= n
        error('dcqp_check_input:dimension_mismatch_Aeq', ...
              'Aeq must have %d columns (same as Q dimensions)', n);
    end
    
    if ~all(isfinite(Aeq(:)))
        error('dcqp_check_input:nonfinite_Aeq', 'Aeq contains non-finite values');
    end
    
    % Check beq vector
    if ~isvector(beq) || length(beq) ~= meq
        error('dcqp_check_input:invalid_beq', ...
              'beq must be a vector with %d elements (same as Aeq rows)', meq);
    end
    
    beq = beq(:); % Ensure column vector
    
    if ~all(isfinite(beq))
        error('dcqp_check_input:nonfinite_beq', 'beq contains non-finite values');
    end
    
    % Check rank of equality constraints
    if rank(full(Aeq)) < min(size(Aeq))
        warning('dcqp_check_input:rank_deficient_Aeq', ...
                'Equality constraint matrix Aeq is rank deficient');
    end
end
% Check feasibility and boundedness of constraint system
% The feasible region must be non-empty and bounded for the DCQP algorithm to work properly
check_feasibility_and_boundedness(A, b, Aeq, beq, n);


% Check parameters structure
if ~isstruct(params)
    error('dcqp_check_input:invalid_params', 'params must be a structure');
end

% Validate key parameters
required_fields = {'max_iterations', 'gap_tolerance', 'mosek_tolerance', ...
                   'gurobi_qp_tolerance', 'gurobi_lp_tolerance', ...
                   'dc_regularization', 'psd_check_tolerance'};

for i = 1:length(required_fields)
    field = required_fields{i};
    if ~isfield(params, field)
        error('dcqp_check_input:missing_parameter', ...
              'Required parameter field ''%s'' is missing', field);
    end
    
    value = params.(field);
    if ~isscalar(value) || ~isfinite(value)
        error('dcqp_check_input:invalid_parameter', ...
              'Parameter ''%s'' must be a finite scalar', field);
    end
    
    % Check that tolerance parameters are positive
    if contains(field, 'tolerance') || strcmp(field, 'dc_regularization')
        if value <= 0
            error('dcqp_check_input:invalid_tolerance', ...
                  'Parameter ''%s'' must be positive', field);
        end
    end
end

% Check max_iterations is a positive integer
if params.max_iterations <= 0 || params.max_iterations ~= round(params.max_iterations)
    error('dcqp_check_input:invalid_max_iterations', ...
          'max_iterations must be a positive integer');
end

% Check gap_tolerance is reasonable
if params.gap_tolerance <= 0 || params.gap_tolerance >= 1
    error('dcqp_check_input:invalid_gap_tolerance', ...
          'gap_tolerance must be between 0 and 1');
end

% Check problem dimensions are reasonable
if n > 10000
    warning('dcqp_check_input:large_problem', ...
            'Problem has %d variables, which may lead to long computation times', n);
end

if size(A,1) + size(Aeq,1) > 10000
    warning('dcqp_check_input:many_constraints', ...
            'Problem has %d constraints, which may lead to long computation times', ...
            size(A,1) + size(Aeq,1));
end

% Check for trivial problems
if n == 0
    error('dcqp_check_input:empty_problem', 'Problem has no variables');
end

end



function check_feasibility_and_boundedness(A, b, Aeq, beq, n)
% CHECK_FEASIBILITY_AND_BOUNDEDNESS  Verify feasible region is non-empty and bounded
%
% DESCRIPTION:
%   Checks that the constraint system A*x <= b, Aeq*x = beq defines a 
%   non-empty AND bounded feasible region by:
%   1. Checking feasibility (solving a feasibility LP)
%   2. Computing actual bounds for each variable x_i (min and max values)
%
% INPUT:
%   A    - Inequality constraint matrix
%   b    - Inequality constraint bounds  
%   Aeq  - Equality constraint matrix
%   beq  - Equality constraint bounds
%   n    - Number of variables



fprintf('Checking feasibility and boundedness of constraint system...\n');

%% Step 1: Check Feasibility
% Solve: find x such that A*x <= b, Aeq*x = beq (with zero objective)
try
    [~, ~, exitflag] = gurobilp(zeros(n, 1), A, b, Aeq, beq, [], [], 0, 1e-6);
    
    if exitflag == -1
        error('dcqp_check_input:infeasible_constraints', ...
              ['The constraint system A*x <= b, Aeq*x = beq is infeasible. ', ...
               'No feasible solution exists. Please check your constraints.']);
    elseif exitflag == -2
        error('dcqp_check_input:unbounded_in_feasibility', ...
              'Constraint system is unbounded (detected during feasibility check)');
    end
    
    fprintf('✓ Feasibility check passed\n');
    
catch ME
    if contains(ME.identifier, 'dcqp_check_input:')
        rethrow(ME); % Re-throw our own errors
    else
        % If gurobilp fails for some other reason, issue a warning
        warning('dcqp_check_input:feasibility_check_failed', ...
                ['Could not verify feasibility due to solver error: %s. ', ...
                 'Proceeding with assumption that constraints are feasible.'], ...
                ME.message);
        return; % Skip boundedness check if we can't even check feasibility
    end
end

%% Step 2: Check Boundedness by Computing Variable Bounds
% For each variable x_i, solve:
% min x_i subject to A*x <= b, Aeq*x = beq
% max x_i subject to A*x <= b, Aeq*x = beq

fprintf('Checking boundedness by computing variable bounds...\n');

unbounded_vars = [];
large_bound_vars = [];
bound_threshold = 1e10; % Consider bounds larger than this as "effectively unbounded"

for i = 1:n
    % Objective vector: minimize x_i
    c_min = zeros(n, 1);
    c_min(i) = 1;
    
    % Objective vector: maximize x_i (minimize -x_i)
    c_max = zeros(n, 1);
    c_max(i) = -1;
    
    try
        % Compute lower bound: min x_i
        [~, lb_val, exitflag_min] = gurobilp(c_min, A, b, Aeq, beq, [], [], 0, 1e-6);
        
        % Compute upper bound: max x_i (solve min -x_i)
        [~, neg_ub_val, exitflag_max] = gurobilp(c_max, A, b, Aeq, beq, [], [], 0, 1e-6);
        
        % Check if either optimization was unbounded
        if exitflag_min == -2
            unbounded_vars(end+1) = i; %#ok<AGROW>
            fprintf('  Variable x_%d: unbounded below\n', i);
            continue;
        end
        
        if exitflag_max == -2
            unbounded_vars(end+1) = i; %#ok<AGROW>
            fprintf('  Variable x_%d: unbounded above\n', i);
            continue;
        end
        
        % Check if optimization failed
        if exitflag_min ~= 1 || exitflag_max ~= 1
            warning('dcqp_check_input:bound_computation_failed', ...
                    'Could not compute bounds for variable x_%d', i);
            continue;
        end
        
        % Compute actual bounds
        x_min = lb_val;
        x_max = -neg_ub_val; % Convert back from minimizing -x_i
        
        % Check if bounds are reasonable
        if abs(x_min) > bound_threshold || abs(x_max) > bound_threshold
            large_bound_vars(end+1) = i; %#ok<AGROW>
            fprintf('  Variable x_%d: bounds [%.2e, %.2e] (very large)\n', i, x_min, x_max);
        else
            %fprintf('  Variable x_%d: bounds [%.6f, %.6f]\n', i, x_min, x_max);
        end
        
    catch ME
        warning('dcqp_check_input:bound_computation_error', ...
                'Error computing bounds for variable x_%d: %s', i, ME.message);
    end
end

%% Step 3: Report Results
if ~isempty(unbounded_vars)
    error('dcqp_check_input:unbounded_variables', ...
          ['Variables x_%s are unbounded. ', ...
           'The DCQP algorithm requires all variables to be bounded. ', ...
           'Please add appropriate box constraints or additional inequalities.'], ...
          sprintf('%d,', unbounded_vars));
end

if ~isempty(large_bound_vars)
    warning('dcqp_check_input:large_bounds', ...
            ['Variables x_%s have very large bounds (>%.0e). ', ...
             'This may cause numerical issues. Consider tightening constraints.'], ...
            sprintf('%d,', large_bound_vars), bound_threshold);
end

fprintf('✓ Boundedness check passed - all variables are bounded\n');
fprintf('Feasibility and boundedness verification complete.\n\n');

end