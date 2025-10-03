function params = dcqp_default_params()
% DCQP_DEFAULT_PARAMS  Default parameters for the DC-QP solver
%
% SYNTAX:
%   params = dcqp_default_params()
%
% DESCRIPTION:
%   Returns a structure containing default algorithm parameters for the
%   DC decomposition based quadratic programming solver.
%
% OUTPUT:
%   params - Structure with the following fields:
%
%   Algorithm Control:
%     .max_iterations          - Maximum number of iterations (300)
%     .gap_tolerance          - Optimality gap tolerance (1e-4)
%     .max_time              - Maximum computation time in seconds (3600)
%     .verbose               - Display progress information (true)
%     .robust_mode           - Continue on errors when possible (false)
%
%   Solver Tolerances:
%     .mosek_tolerance       - MOSEK SDP solver tolerance (1e-9)
%     .gurobi_qp_tolerance   - Gurobi QP solver tolerance (1e-9)
%     .gurobi_lp_tolerance   - Gurobi LP solver tolerance (1e-9)
%
%   Algorithm Parameters:
%     .dc_regularization     - DC decomposition regularization (1e-5)
%     .psd_check_tolerance   - PSD checking tolerance (1e-8)
%     .initial_sampling_rounds - Random starts for upper bound (1)
%
%   Solver Methods:
%     .gurobi_lp_method      - Gurobi LP method (1=dual simplex)
%     .gurobi_qp_method      - Gurobi QP method (2=barrier)
%
%   Problem Bounds:
%     .lower_bounds          - Variable lower bounds ([] = -inf)
%     .upper_bounds          - Variable upper bounds ([] = +inf)
%
% EXAMPLES:
%   % Use default parameters
%   params = dcqp_default_params();
%
%   % Modify tolerance
%   params = dcqp_default_params();
%   params.gap_tolerance = 1e-6;
%
%   % Set variable bounds
%   params = dcqp_default_params();
%   params.lower_bounds = zeros(n, 1);  % Non-negative variables
%   params.upper_bounds = ones(n, 1);   % Upper bound of 1
%
% SEE ALSO: dcqp_solve, qpsolver

% Copyright (c) 2025, [Your Institution]
% All rights reserved.

params = struct();

% =================================================================
% Algorithm Control Parameters
% =================================================================
params.max_iterations = 300;           % Maximum iterations
params.gap_tolerance = 1e-4;           % Relative optimality gap tolerance
params.max_time = 3600;                % Maximum time in seconds (1 hour)
params.display_summary = true;                 % Display progress information
params.robust_mode = false;            % Continue on errors when possible
params.save_failed_instances = true;   % Save failed instances for debugging

% =================================================================
% Solver Tolerance Parameters
% =================================================================
params.mosek_tolerance = 1e-8;         % MOSEK SDP solver tolerance 
params.tol_mosek_cut =1e-8;            % MOSEK CUT SDP tolerance
params.gurobi_qp_tolerance = 1e-9;     % Gurobi QP solver tolerance  
params.gurobi_lp_tolerance = 1e-9;     % Gurobi LP solver tolerance

% =================================================================
% Algorithm-Specific Parameters
% =================================================================
params.dc_regularization = 1e-5;       % DC decomposition regularization (spn)
params.psd_check_tolerance = 1e-8;     % Tolerance for PSD checking (eps_checkpsd)
params.nb_rounds = 100;    % Number of random initializations
params.eta=0.9;

% =================================================================
% Solver Method Selection
% =================================================================
params.gurobi_lp_method = 1;           % 0=primal, 1=dual simplex, 2=barrier
params.gurobi_qp_method = 2;           % 0=auto, 1=simplex, 2=barrier

% =================================================================
% Variable Bounds (optional)
% =================================================================
params.lower_bounds = [];              % Variable lower bounds ([] = -inf)
params.upper_bounds = [];              % Variable upper bounds ([] = +inf)

% =================================================================
% Internal Parameters (usually not modified by users)
% =================================================================
params.scaling = 1;                    % Problem scaling factor
params.do_scaling=false;
params.filename = 'dcqp_result';       % Output filename prefix
params.verbose=false;                  % Conditional screen display

% Legacy parameter names for compatibility with existing code
params.max_N = params.max_iterations;
params.gap_tol = params.gap_tolerance;
params.eps_dcdecomposition = params.dc_regularization;
params.eps_checkpsd = params.psd_check_tolerance;
params.tol_mosek = params.mosek_tolerance;
params.tol_gqp = params.gurobi_qp_tolerance;
params.tol_glp = params.gurobi_lp_tolerance;
params.metglp = params.gurobi_lp_method;
params.metgqp = params.gurobi_qp_method;





end
