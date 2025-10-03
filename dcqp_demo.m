%% EXAMPLE_BASIC - Basic usage of the DC-QP solver
%
% This example demonstrates the basic usage of the DC-QP solver on simple
% quadratic programming problems, including both convex and non-convex cases.

%% Clear workspace
clear; close all; clc;

diary('diaryfile_demo.txt');
diary on;

fprintf('=== DC-QP Basic Examples ===\n\n');

%% Example 1: Simple nonconvex QP 
fprintf('Example 1: Simple nonconvex QP\n');
fprintf('--------------------------------\n');

% Problem data - changed to nonconvex Q matrix  
n = 3;
Q1 = [1, -1, 0; -1, -2, -1; 0, -1, 1];  % nonconvex matrix
d1 = [-1; -2; -3];
A1 = [1, 1, 1; -1, 0, 0; 0, -1, 0; 0, 0, -1; 1, 0, 0; 0, 1, 0; 0, 0, 1];
b1 = [3; 0; 0; 0; 5; 5; 5];  % Box constraints: 0 <= x <= 5, sum(x) <= 3

% Solve with DC-QP
fprintf('Solving nonconvex QP with DC-QP...\n');
[x1, fval1, info1] = dcqp_solve(Q1, d1, A1, b1);

fprintf('Solution: x = [%.4f, %.4f, %.4f]\n', x1);
fprintf('Objective value: %.6f\n', fval1);
fprintf('Status: %s\n', info1.status);
fprintf('Optimality gap: %.2e\n', info1.gap);
fprintf('Computation time: %.2f seconds\n\n', info1.time);

%% Example 2: Nonconvex QP with better bounds
fprintf('Example 2: Nonconvex QP with box constraints\n');
fprintf('----------------------------------------------\n');

% Problem data - nonconvex Q matrix with proper bounds
Q2 = [1, -1 1 1 -1;-1 1 -1 1 1;1 -1 1 -1 1;1 1 -1 1 -1; -1 1 1 -1 1];  % nonconvex (has negative eigenvalue)
d2 = [0; 0;0;0; 0];
A2 = [-1 0 0 0 0; 0 -1 0 0 0; 0 0 -1 0 0; 0 0 0 -1 0; 0 0 0 0 -1; 1 1 1 1 1];  % x1+x2 <= 2, 0 <= x <= 3
b2 = [0; 0; 0; 0; 0; 1];

% Solve with DC-QP
fprintf('Solving nonconvex QP with DC-QP...\n');

params = dcqp_default_params();
params.filename='example2-Horn-matrix';
params.nb_rounds=1;
[x2, fval2, info2] = dcqp_solve(Q2, d2, A2, b2,[],[],params);

fprintf('Solution: x2 = [%.4f, %.4f,%.4f, %.4f, %.4f]\n', x2);
fprintf('Objective value: %.6f\n', fval2);
fprintf('Status: %s\n', info2.status);
fprintf('Optimality gap: %.2e\n', info2.gap);
fprintf('Computation time: %.2f seconds\n\n', info2.time);

%% Example 3: Nonconvex QP with equality constraints
fprintf('Example 3: QP with equality constraints\n');
fprintf('---------------------------------------\n');

% Problem data
Q3 = [4, -8; -8, 4];
d3 = [-1; -1];
A3 = [-1, 0; 0, -1];  % x >= 0
b3 = [0; 0];
Aeq3 = [1, 1];        % x1 + x2 = 1
beq3 = 1;

% Solve with DC-QP
fprintf('Solving nonconvex QP with equality constraints...\n');
[x3, fval3, info3] = dcqp_solve(Q3, d3, A3, b3, Aeq3, beq3);

fprintf('Solution: x = [%.4f, %.4f]\n', x3);
fprintf('Objective value: %.6f\n', fval3);
fprintf('Status: %s\n', info3.status);
fprintf('Constraint check: x1 + x2 = %.6f (should be 1.0)\n', sum(x3));
fprintf('Computation time: %.2f seconds\n\n', info3.time);

%% Example 4: Customizing solver parameters
fprintf('Example 4: Custom solver parameters\n');
fprintf('------------------------------------\n');

% Use the same problem as Example 2 but with tighter tolerance
params = dcqp_default_params();
params.gap_tolerance = 1e-6;  % Tighter tolerance
params.verbose = false;       % Suppress detailed output

fprintf('Solving with tighter tolerance (1e-6)...\n');
[x4, fval4, info4] = dcqp_solve(Q3, d3, A3, b3, Aeq3, beq3, params);

fprintf('Solution: x = [%.6f, %.6f]\n', x4);
fprintf('Objective value: %.8f\n', fval4);
fprintf('Status: %s\n', info4.status);
fprintf('Optimality gap: %.2e\n', info4.gap);
fprintf('Computation time: %.2f seconds\n\n', info4.time);



%% Summary
fprintf('\n=== Summary ===\n');
fprintf('All examples completed successfully!\n');
fprintf('The DC-QP solver can handle:\n');
fprintf('  - Non-convex (nonconvex) quadratic programs\n');
fprintf('  - Problems with equality and inequality constraints\n');
fprintf('  - Custom solver parameters\n');
diary off;