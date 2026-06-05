function save_failed_instance(Q, d, A, b, Aeq, beq, params, error_info, solve_time)
% SAVE_FAILED_INSTANCE  Save failed DCQP instances for debugging
%
% SYNTAX:
%   save_failed_instance(Q, d, A, b, Aeq, beq, params, error_info, solve_time)
%
% DESCRIPTION:
%   Saves problem instances that failed to solve, along with detailed
%   error information, for later analysis and debugging.
%
% INPUT:
%   Q           - Quadratic coefficient matrix
%   d           - Linear coefficient vector  
%   A           - Inequality constraint matrix
%   b           - Inequality constraint RHS
%   Aeq         - Equality constraint matrix
%   beq         - Equality constraint RHS
%   params      - Solver parameters
%   error_info  - MATLAB exception object
%   solve_time  - Time elapsed before failure
%
% OUTPUT:
%   Saves a .mat file containing all problem data and error information

% Copyright (c) 2025
% All rights reserved.

try
    % Get the directory of this script and find the project directory
    script_dir = fileparts(mfilename('fullpath'));
    project_dir = script_dir;  % save_failed_instance.m is in the main directory
    
    % Create failed_instances directory if it doesn't exist
    failed_dir = fullfile(project_dir, 'failed_instances');
    if ~exist(failed_dir, 'dir')
        mkdir(failed_dir);
    end
    
    % Generate unique filename with timestamp
    timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
    filename = fullfile(failed_dir, sprintf('failed_instance_%s.mat', timestamp));
    
    % Collect problem data
    problem_data = struct();
    problem_data.Q = Q;
    problem_data.d = d;
    problem_data.A = A;
    problem_data.b = b;
    problem_data.Aeq = Aeq;
    problem_data.beq = beq;
    
    % Collect problem characteristics
    problem_info = struct();
    problem_info.n_vars = size(Q, 1);
    problem_info.n_ineq_constraints = size(A, 1);
    problem_info.n_eq_constraints = size(Aeq, 1);
    problem_info.Q_eigenvalues = eig(Q);
    problem_info.Q_condition_number = cond(Q);
    problem_info.Q_rank = rank(Q);
    problem_info.is_Q_indefinite = any(eig(Q) < -1e-8) && any(eig(Q) > 1e-8);
    problem_info.constraint_matrix_rank = rank([A; Aeq]);
    
    % Collect error information
    error_data = struct();
    error_data.message = error_info.message;
    error_data.identifier = error_info.identifier;
    error_data.stack = error_info.stack;
    error_data.solve_time = solve_time;
    error_data.timestamp = timestamp;
    
    % Collect solver parameters
    solver_params = params;
    
    % Save all data
    save(filename, 'problem_data', 'problem_info', 'error_data', 'solver_params', ...
         '-v7.3'); % Use v7.3 format for large files
    
    fprintf('Failed instance saved to: %s\n', filename);
    
catch save_error
    fprintf('Warning: Could not save failed instance: %s\n', save_error.message);
end

end
