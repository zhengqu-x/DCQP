function dcqp_startup()
% DCQP_STARTUP  Initialize the DC-QP solver environment
%
% SYNTAX:
%   dcqp_startup()
%
% DESCRIPTION:
%   This function sets up the DCQP solver environment by:
%   - Adding necessary paths to MATLAB
%   - Checking for required dependencies
%   - Running basic functionality tests
%   - Displaying system information
%
% USAGE:
%   Run this function once after installing the DCQP solver:
%     dcqp_startup();
%
%   For automatic setup, add the following to your MATLAB startup.m file:
%     addpath('/path/to/DC-QP');
%     dcqp_startup();

% Copyright (c) 2025
% All rights reserved.

fprintf('=== DC-QP Solver Initialization ===\n\n');

% Get current directory
dcqp_root = fileparts(mfilename('fullpath'));
fprintf('DC-QP root directory: %s\n', dcqp_root);

% Add paths
fprintf('Setting up MATLAB paths...\n');
addpath(dcqp_root);

% Add legacy directory
legacy_dir = fullfile(dcqp_root, 'legacy');
if exist(legacy_dir, 'dir')
    addpath(legacy_dir);
    fprintf('  Added legacy directory\n');
end

% Add utils directory
utils_dir = fullfile(dcqp_root, 'utils');
if exist(utils_dir, 'dir')
    addpath(utils_dir);
    fprintf('  Added utils directory\n');
end

% Add plotstables directory (optional, if needed)
plotstables_dir = fullfile(dcqp_root, 'plotstables');
if exist(plotstables_dir, 'dir')
    addpath(plotstables_dir);
    fprintf('  Added plotstables directory\n');
end

% Add paper example scripts
paper_examples_dir = fullfile(dcqp_root, 'paper-examples');
if exist(paper_examples_dir, 'dir')
    addpath(paper_examples_dir);
    fprintf('  Added paper-examples directory\n');
end

% Check MATLAB version
matlab_version = version('-release');
matlab_year = str2double(matlab_version(1:4));
fprintf('\nMATLAB version: %s', matlab_version);
if matlab_year >= 2020
    fprintf(' OK\n');
else
    fprintf(' WARNING (Recommended: R2020a or later)\n');
end

% Check for required toolboxes
fprintf('\nChecking MATLAB toolboxes...\n');
if license('test', 'Optimization_Toolbox')
    fprintf('  Optimization Toolbox: OK\n');
else
    fprintf('  Optimization Toolbox: MISSING (Required)\n');
end

% Check for third-party solvers
fprintf('\nChecking third-party solvers...\n');

% Check MOSEK
try
    mosekopt('version');
    fprintf('  MOSEK: OK\n');
    mosek_available = true;
catch
    fprintf('  MOSEK: MISSING (Required for SDP computations)\n');
    mosek_available = false;
end

% Check Gurobi
try
    model = struct();
    model.A = sparse(1, 1);
    model.rhs = 0;
    model.sense = '=';
    model.obj = full(1);  % Ensure obj is dense
    params = struct();
    params.OutputFlag = 0;
    gurobi(model, params);
    fprintf('  Gurobi: OK\n');
    gurobi_available = true;
catch
    fprintf('  Gurobi: MISSING (Required for LP/QP subproblems)\n');
    gurobi_available = false;
end

% Summary of requirements
fprintf('\nDependency status:\n');
if mosek_available && gurobi_available
    fprintf('  All required solvers available\n');
    all_deps_ok = true;
else
    fprintf('  Missing required solvers\n');
    all_deps_ok = false;
    
    fprintf('\nTo install missing dependencies:\n');
    if ~mosek_available
        fprintf('  - MOSEK: Visit https://www.mosek.com/ for installation\n');
    end
    if ~gurobi_available
        fprintf('  - Gurobi: Visit https://www.gurobi.com/ for installation\n');
    end
end

% Run basic functionality test if all dependencies are available
if all_deps_ok
    fprintf('\nRunning basic functionality test...\n');
    try
        % Simple test problem with bounded feasible region - INDEFINITE Q
        Q = [-1, 0.5; 0.5, 1];  % Indefinite matrix (required for DCQP)
        d = [-1; -1];
        
        % Create bounded feasible region: 0 <= x <= 2
        A = [-1, 0;    % -x1 <= 0  -> x1 >= 0
             0, -1;    % -x2 <= 0  -> x2 >= 0
             1, 0;     %  x1 <= 2
             0, 1];    %  x2 <= 2
        b = [0; 0; 2; 2];
        
        % Validate the test problem before using it
        params = dcqp_default_params();
        params.verbose = false;
        
        try
            dcqp_check_input(Q, d, A, b, [], [], params);
        catch ME
            fprintf('  Test problem validation failed: %s\n', ME.message);
            fprintf('  Skipping functionality test.\n');
            return;
        end
        
        [x, ~, info] = dcqp_solve(Q, d, A, b, [], [], params);
        
        % Check if we got a reasonable result
        if info.gap<=params.gap_tolerance && all(isfinite(x)) && all(A*x <= b + 1e-6)
            fprintf('  Basic test: OK (solved in %.2f seconds)\n', info.time);
        else
            fprintf('  Basic test: WARNING (unexpected result: x=[%.3f,%.3f], status=%s)\n', ...
                    x(1), x(2), info.status);
        end
        
    catch ME
        fprintf('  Basic test: FAILED (error: %s)\n', ME.message);
    end
else
    fprintf('\nSkipping functionality test (missing dependencies)\n');
end

% Display usage information
fprintf('\n=== Quick Start ===\n');
fprintf('To solve a quadratic program:\n\n');
fprintf('  %% Define your problem\n');
fprintf('  Q = [1, -1; -1, -1];      %% Indefinite objective Hessian\n');
fprintf('  d = [1; 1];               %% Linear coefficients\n');
fprintf('  A = [1, 1; -1, 0; 0, -1]; %% Bounded inequality constraints\n');
fprintf('  b = [1; 0; 0];            %% Inequality bounds\n\n');
fprintf('  %% Solve\n');
fprintf('  [x, fval, info] = dcqp_solve(Q, d, A, b);\n\n');
fprintf('For examples: run dcqp_demo() or the scripts in paper-examples/\n');
fprintf('For help: type ''help dcqp_solve'' or ''doc dcqp_solve''\n');

% Display final status
fprintf('\n=== Initialization Complete ===\n');
if all_deps_ok
    fprintf('Status: Ready to use\n');
else
    fprintf('Status: Needs setup\n');
end

fprintf('DCQP solver initialized successfully!\n\n');

% Suggest next steps
if all_deps_ok
    fprintf('Suggested next steps:\n');
    fprintf('  1. Run demo: dcqp_demo()\n');
    fprintf('  2. Try existing test sets: solve_existing_testsets_with_dcqp(''qp20_10'')\n');
    fprintf('  3. Compare with Gurobi: solve_existing_testsets_with_gurobi(''qp20_10'')\n');
else
    fprintf('Next steps:\n');
    fprintf('  1. Install missing dependencies (see above)\n');
    fprintf('  2. Run dcqp_startup again to verify installation\n');
end

end
