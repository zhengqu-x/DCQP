# DCQP: Doubly Nonnegative based Cutting Plane method for Quadratic Programming

[![MATLAB](https://img.shields.io/badge/MATLAB-R2020a+-blue.svg)](https://www.mathworks.com/products/matlab.html)
[![License](https://img.shields.io/badge/License-Academic-green.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-1.0.0-orange.svg)](https://github.com/zhengqu-x/DCQP)

## Overview

DCQP is a MATLAB package for solving **nonconvex quadratic programming** problems using a doubly nonnegative relaxation and cutting plane approach. The solver is specifically designed for problems where traditional convex optimization methods fail due to indefinite Hessian matrices.

### Problem Formulation

DCQP solves quadratic programming problems of the form:

```
minimize    x'*Q*x + 2*d'*x
subject to  A*x <= b
            Aeq*x = beq
```

where:
- **Q** is a symmetric matrix (possibly indefinite/nonconvex)
- **d** is the linear objective coefficient vector
- **A, b** define inequality constraints (**mandatory**)
- **Aeq, beq** define equality constraints (optional)

### Key Features

- ✅ **Nonconvex QP Solver**: Handles indefinite Hessian matrices with negative eigenvalues
- ✅ **Global Optimization**: Uses cutting plane methods for finding global optima
- ✅ **Bounded Feasible Region**: Requires bounded constraint sets
- ✅ **Multiple Solvers**: Integrates with Gurobi and MOSEK for subproblems
- ✅ **Comprehensive Examples**: Includes 64 benchmark and 140 synthetic test problems
- ✅ **Robust Implementation**: Error handling and debugging features

## Requirements

### Software Dependencies
- **MATLAB** R2020a or later
- **Gurobi Optimizer** (recommended version 12.0.1 or later)
- **MOSEK** (recommended version 11.0.12 or later)

### System Requirements
- Memory: At least 4GB RAM (8GB+ recommended for large problems)
- Operating System: Windows, macOS, or Linux

## Installation

1. **Clone or download** the DCQP package to your local machine
2. **Add to MATLAB path**:
   ```matlab
   addpath('/path/to/DCQP');
   dcqp_startup();
   ```
3. **Verify installation**:
   ```matlab
   dcqp_version();
   dcqp_demo();
   ```

### Automatic Setup
Add the following lines to your MATLAB `startup.m` file for automatic initialization:
```matlab
addpath('/path/to/DCQP');
dcqp_startup();
```

## Quick Start

### Basic Usage

```matlab
% Define a nonconvex QP problem
Q = [1 -1; -1 -1];    % Indefinite matrix (nonconvex)
d = [1; 1];           % Linear coefficients  
A = [1 1; -1 0; 0 -1]; % Inequality constraints
b = [1; 0; 0];        % Right-hand side

% Solve the problem
[x_opt, fval, info] = dcqp_solve(Q, d, A, b);

fprintf('Optimal solution: x = [%.4f, %.4f]\n', x_opt);
fprintf('Optimal value: %.6f\n', fval);
fprintf('Status: %s\n', info.status);
```

### With Equality Constraints

```matlab
% Add equality constraints
Aeq = [1 0];          % x1 = 0.5
beq = 0.5;

[x_opt, fval, info] = dcqp_solve(Q, d, A, b, Aeq, beq);
```

### Custom Parameters

```matlab
% Configure solver parameters
params = dcqp_default_params();
params.gap_tolerance = 1e-6;      % Tighter optimality tolerance
params.max_iterations = 500;      % More iterations
params.verbose = true;            % Display progress

[x_opt, fval, info] = dcqp_solve(Q, d, A, b, [], [], params);
```

## Algorithm Overview

DCQP employs a **Doubly nonnegative based Cutting plane** approach:

1. **Cutting Plane Method**: Iteratively adds linear cuts to the original QP problem
2. **Upper Bound Computation**: Uses local search methods for feasible solutions
3. **Lower Bound Computation**: Solves SDP relaxations for global lower bounds
4. **Convergence**: Terminates when the relative gap between lower and upper bounds is sufficiently small

### Key Algorithmic Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `gap_tolerance` | 1e-4 | Relative optimality gap tolerance |
| `max_iterations` | 300 | Maximum number of cutting plane iterations |
| `max_time` | 3600 | Maximum computation time (seconds) |
| `dc_regularization` | 1e-5 | DC decomposition regularization parameter |
| `nb_rounds` | 100 | Number of random initializations for upper bound |
| `do_scaling` | false | Scaling of the objective value |

## Datasets

### Benchmark Problems

The package includes 64 benchmark problems from the literature, organized in four groups:
- **qp20_10**: 20 variables, 10 constraints (16 instances)
- **qp30_15**: 30 variables, 15 constraints (16 instances)  
- **qp40_20**: 40 variables, 20 constraints (16 instances)
- **qp50_25**: 50 variables, 25 constraints (16 instances)



### Synthetic Problems

The synthetic problems were generated using `legacy/generateinstances_uniform.m` and `legacy/generateinstances_normal.m`. These create 140 nonconvex QPs organized in seven groups (20 instances each):

- **qp_n_0_1**: Normal distribution, density 0.1, no equality constraints
- **qp_n_0_3**: Normal distribution, density 0.3, no equality constraints  
- **qp_n_0_9**: Normal distribution, density 0.9, no equality constraints
- **qp_u_0_1**: Uniform distribution, density 0.1, no equality constraints
- **qp_u_0_3**: Uniform distribution, density 0.3, no equality constraints
- **qp_u_0_9**: Uniform distribution, density 0.9, no equality constraints
- **qp_u_25_1**: Uniform distribution, density 0.1, 25 equality constraints

All synthetic problems have:
- **Problem size**: 100 variables, 51 inequality constraints
- **Hessian construction**: Q = L₁ᵀL₁ - L₂ᵀL₂ where L₁ and L₂ are random sparse matrices (25×100 each) with entries drawn from uniform or normal distributions
- **Constraint generation**: 
  - First constraint: A₁ = [1, 1, ..., 1] (sum constraint)
  - Remaining constraints: A₂₋₅₁ are sparse random matrices with specified density
  - Right-hand side: b = A·x₀ + 0.1·rand() where x₀ is a feasible point
  - Equality constraints: Aeq (when present) generated similarly with beq = Aeq·x₀
- **Variable bounds**: 0 ≤ x ≤ 1
- **Feasible interior point**: x₀ constructed as normalized random vector to ensure all constraints are satisfiable




## Examples and Demos

### Run the Demo
```matlab
dcqp_demo();  % Runs basic examples with different problem types
```

### Benchmark Problems

```matlab
% Navigate to paper-examples directory first
cd('paper-examples/');

% Solve benchmark problems from the literature (16 instances in each group)
solve_benchmark_with_dcqp('qp20_10');  % 20 variables, 10 constraints
solve_benchmark_with_dcqp('qp30_15');  % 30 variables, 15 constraints
solve_benchmark_with_dcqp('qp40_20');  % 40 variables, 20 constraints
solve_benchmark_with_dcqp('qp50_25');  % 50 variables, 25 constraints
```

### Synthetic Problems

```matlab
% Navigate to paper-examples directory first
cd('paper-examples/');

% Test on randomly generated problems (20 instances in each group)
solve_synthetic_with_dcqp('qp_n_0_1');   % density 0.1, normal distribution, no equality constraint
solve_synthetic_with_dcqp('qp_n_0_3');   % density 0.3, normal distribution, no equality constraint
solve_synthetic_with_dcqp('qp_n_0_9');   % density 0.9, normal distribution, no equality constraint
solve_synthetic_with_dcqp('qp_u_0_1');   % density 0.1, uniform distribution, no equality constraint
solve_synthetic_with_dcqp('qp_u_0_3');   % density 0.3, uniform distribution, no equality constraint
solve_synthetic_with_dcqp('qp_u_0_9');   % density 0.9, uniform distribution, no equality constraint
solve_synthetic_with_dcqp('qp_u_25_1');  % density 0.1, uniform distribution, 25 equality constraints

% Test specific instance in a group
solve_synthetic_with_dcqp('qp_n_0_1', 5);  % Run only instance 5 from group
```

### Comparison With Gurobi

For performance comparison, the package includes Gurobi-based solvers that attempt to solve the same nonconvex QP problems:

```matlab
% Navigate to paper-examples directory first
cd('paper-examples/');

% Benchmark problems with Gurobi (with optional time limit)
solve_benchmark_with_gurobi('qp20_10');           % Default time limit (1 hour)
solve_benchmark_with_gurobi('qp30_15', 7200);     % 2 hours time limit

% Synthetic problems with Gurobi
solve_synthetic_with_gurobi('qp_n_0_1');          % All instances in group, defalut time limit (1 hour)
solve_synthetic_with_gurobi('qp_u_0_1', 5);       % Specific instance only, default time limit (1 hour)
solve_synthetic_with_gurobi('qp_u_0_1', 5, 1800); % With 30-minute time limit
```

### Output Results

Each benchmark and synthetic run generates detailed results:

**Console Output**: 
- **DCQP Solution Summary** for each instance displaying:
  - Instance name and problem dimensions (n variables, m inequality constraints, meq equality constraints)
  - Solution status message (e.g., "successfully reduced relative gap below 0.0001")
  - Best objective value (scientific notation, e.g., -3.000000e+01)
  - Relative optimality gap (scientific notation, e.g., 4.33e-11)
  - Computation time in seconds and total number of iterations
- **Real-time progress** (when `params.verbose = true`): iteration solver details including bounds

**Saved Files**:
- **Diary files**: `diaryfile-benchmark-dcqp.txt`, `diaryfile-synthetic-dcqp.txt`, etc. containing complete console logs
- **Individual result files**: Saved in `paper-examples/testresults/` directory with timestamps (e.g., `gurobi_qp20_10_1_1-2025-10-03_13-18-07.mat`) containing `bestsol` (optimal solution vector) and `info` (performance metrics for that instance)
- **Summary statistics**: Saved in `paper-examples/summary_results/` directory as `myrecord` matrices containing performance data for all solved instances. Each row represents one problem instance with 6 columns: [optimality_gap, max_constraint_violation, equality_constraint_violation, objective_value, lower_bound, computation_time]. For synthetic problems, this is only saved when all 20 instances in a group are solved together (not for individual instance runs).


## Function Reference

### Main Functions

- **`dcqp_solve(Q, d, A, b, Aeq, beq, params)`**: Main solver function
- **`dcqp_default_params()`**: Get default algorithm parameters
- **`dcqp_startup()`**: Initialize the DCQP environment
- **`dcqp_version()`**: Display version information
- **`dcqp_demo()`**: Run demonstration examples

### Utility Functions

- **`DC_decomposition(Q, spn)`**: Compute DC decomposition of matrix Q
- **`compute_ub(Q, d, A, b, Aeq, beq, n, params, sol, nb_rounds)`**: Compute upper bound
- **`generate_cut_dnn(Q, d, A, b, Aeq, beq, m, n, nuR, barx, x0, tol_mosek, beta)`**: Generate doubly nonnegative cutting plane
- **`lower_bound_dnn(Q, d, A, b, Aeq, beq, tol_mosek, m, n)`**: Compute doubly nonnegative relaxation lower bound
- **`check_kkt_conditions(x, Q, d, A, b, Aeq, beq)`**: Verify KKT conditions
- **`qpsolver(H, f, A, b, Aeq, beq, method, tolerance, n)`**: Generic QP solver interface

## Output Structure

The solver returns detailed information about the solution:

```matlab
info = struct(
    'status',      'successfully reduced relative gap below 0.0001',  % Solution status
    'gap',         4.33e-11,        % Relative optimality gap
    'iterations',  1,               % Number of iterations
    'time',        0.19,            % Total computation time (seconds)
    'upper_bound', -30.0,           % Best upper bound found
    'lower_bound', -30.0,           % Best lower bound achieved  
    'scaling',     1                % Problem scaling factor
);
```

### Status Codes
- `'successfully reduced relative gap below X'`: Solution found within gap tolerance X
- `'time_limit reached'`: Maximum computation time exceeded
- `'iteration_limit reached'`: Maximum number of iterations reached
- `'not solved'`: Algorithm terminated without meeting convergence criteria
- `'error'`: Solver encountered an error during execution

## Important Notes and Limitations

### Problem Requirements
1. **Negative eigenvalues**: Q must have at least one negative eigenvalue (nonconvex)
2. **Bounded feasible region**: The constraint set A*x ≤ b must be bounded
3. **Mandatory inequalities**: Inequality constraints A*x ≤ b cannot be empty
4. **No integer variables**: Continuous variables only

### Performance Considerations
- **Problem size**: Most efficient for problems with n ≤ 100 variables
- **Constraint density**: Performance degrades with very dense constraint matrices
- **Conditioning**: Ill-conditioned problems may require parameter tuning

### Troubleshooting
- **Unbounded problems**: Verify that the feasible region is bounded
- **Convergence issues**: Try increasing `gap_tolerance` or `max_iterations`
- **Memory errors**: Reduce `nb_rounds` for small problems
- **Solver failures**: Enable `save_failed_instances` for debugging

## File Structure

```
DCQP/
├── dcqp_solve.m              # Main solver function
├── dcqp_default_params.m     # Default parameters
├── dcqp_demo.m               # Demonstration examples
├── dcqp_startup.m            # Environment setup
├── dcqp_version.m            # Version information
├── dcqp_check_input.m        # Input validation
├── data/                     # Test problems
│   ├── benchmark/            # 64 literature benchmark problems
│   └── synthetic/            # Synthetic test instances (140 problems generated using legacy/generateinstances_uniform.m and legacy/generateinstances_normal.m)
├── utils/                    # Utility functions
│   ├── DC_decomposition.m    # DC decomposition
│   ├── compute_ub.m          # Upper bound computation
│   ├── qpsolver.m            # QP solver interface
│   └── ...                   # Other utilities  
├── paper-examples/           # Reproducible experiments
├── legacy/                   # Legacy functions
└── testresults/             # Saved results
```



## Authors

- **Zheng Qu, Defeng Sun, Jintao Xu** 

## Citation

If you use DCQP in your research, please cite the following paper:

```bibtex
@misc{qu2025progressiveboundstrengtheningdoubly,
      title={Progressive Bound Strengthening via Doubly Nonnegative Cutting Planes for Nonconvex Quadratic Programs}, 
      author={Zheng Qu and Defeng Sun and Jintao Xu},
      year={2025},
      eprint={2510.02948},
      archivePrefix={arXiv},
      primaryClass={math.OC},
      url={https://arxiv.org/abs/2510.02948}, 
}
```

**Plain text citation:**
> Zheng Qu, Defeng Sun, and Jintao Xu. "Progressive Bound Strengthening via Doubly Nonnegative Cutting Planes for Nonconvex Quadratic Programs." arXiv preprint arXiv:2510.02948, 2025. https://arxiv.org/abs/2510.02948

## License

This software is distributed under an Academic License. See `LICENSE` file for details.

## Support and Issues

- **Documentation**: See function help: `help dcqp_solve`
- **Examples**: Run `dcqp_demo()` for working examples
- **Issues**: Report bugs and feature requests on GitHub
- **Contact**: zhengqu@szu.edu.cn

## Version History

- **v1.0.0** (2025-10-03): Initial release
  - Core DCQP algorithm implementation
  - Benchmark and synthetic test examples
  - Comprehensive documentation

---

**Disclaimer**: This is research software. While extensively tested, use in production environments should be done with appropriate validation.
