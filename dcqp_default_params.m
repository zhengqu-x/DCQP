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
%     .display_summary       - Display final solution summary (true)
%     .verbose               - Display iteration progress information (true)
%     .mosek_quiet           - Suppress MOSEK optimizer logs (true)
%     .robust_mode           - Continue on errors when possible (false)
%
%   Solver Tolerances:
%     .mosek_tolerance       - MOSEK SDP solver tolerance (1e-8)
%     .gurobi_qp_tolerance   - Gurobi QP solver tolerance (1e-9)
%     .gurobi_lp_tolerance   - Gurobi LP solver tolerance (1e-9)
%
%   Algorithm Parameters:
%     .dc_regularization     - DC decomposition regularization (1e-5)
%     .psd_check_tolerance   - PSD checking tolerance (1e-8)
%     .nb_rounds             - Random starts for upper bound (10)
%     .konnofirst            - Try an SDP-validated Konno cut before SDP cuts (false)
%     .accept_cut_below_threshold
%                            - Preserve legacy handling when cut_lb < nu (false)
%     .known_solution        - Known optimal point in input coordinates ([])
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

% Copyright (c) 2025
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
params.mosek_quiet = true;             % Suppress MOSEK optimizer logs

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
params.nb_rounds = 10;    % Number of random initializations
params.eta=0.9;
params.konnofirst=false;              % Try Konno first near the incumbent objective
params.accept_cut_below_threshold=false; % Throw below nu so Konno fallback can run
params.known_solution=[];             % Use directly; skip random initialization

% =================================================================
% Solver Method Selection
% =================================================================
params.gurobi_lp_method = 1;           % 0=primal, 1=dual simplex, 2=barrier
params.gurobi_qp_method = 2;           % -1=auto, 0=primal simplex, 1=dual simplex, 2=barrier

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
params.verbose=true;                  % Conditional screen display

end
