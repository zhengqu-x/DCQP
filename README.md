# DCQP: Doubly Nonnegative based Cutting Plane method for Quadratic Programming

[![MATLAB](https://img.shields.io/badge/MATLAB-R2020a+-blue.svg)](https://www.mathworks.com/products/matlab.html)
[![License](https://img.shields.io/badge/License-Academic-green.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-1.0.1-orange.svg)](https://github.com/zhengqu-x/DCQP)

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
- ✅ **Comprehensive Examples**: Includes scripts to reproduce the experiments in the paper
- ✅ **Robust Implementation**: Error handling and debugging features

## Requirements

### Software Dependencies
- **MATLAB** R2020a or later
- **Gurobi Optimizer** (recommended version 12.0.1 or later)
- **MOSEK** 11.0.30 or earlier with the MATLAB interface configured. Newer MOSEK releases such as 11.2 are not currently supported because DCQP currently uses MOSEK's legacy `mosekopt` MATLAB toolbox interface and semidefinite-programming data structures.

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

## Direct use of `dcqp_solve`

### Model and function signature

The public solver call is

```matlab
[x_opt,fval,info] = dcqp_solve(Q,d,A,b,Aeq,beq,params);
```

and it solves

```text
minimize    x'*Q*x + 2*d'*x
subject to  A*x <= b
            Aeq*x = beq.
```

`A` and `b` are required and must describe a nonempty, bounded feasible
region. `Aeq`, `beq`, and `params` are optional. DCQP is intended for an
indefinite `Q`; use a convex QP solver when `Q` is positive semidefinite.


Without equality constraints, either omit them or pass empty matrices:

```matlab
[x_opt,fval,info] = dcqp_solve(Q,d,A,b);
% Equivalent:
[x_opt,fval,info] = dcqp_solve(Q,d,A,b,[],[]);
```

### Complete minimal example

```matlab
Q = [1 -1; -1 -1];
d = [1; 1];
A = [1 1; -1 0; 0 -1];
b = [1; 0; 0];

[x_opt,fval,info] = dcqp_solve(Q,d,A,b);

fprintf('x = [%.6g, %.6g]\n',x_opt);
fprintf('objective = %.8g\n',fval);
fprintf('relative gap = %.2e\n',info.gap);
fprintf('status = %s\n',info.status);
```


### Variable bounds

Variable bounds must currently be included in `A*x <= b`. For example,
to impose `lb <= x <= ub`:

```matlab
n = size(Q,1);
A = [A; eye(n); -eye(n)];
b = [b; ub(:); -lb(:)];
```

### Setting options

Always begin with the complete default structure, then override individual
fields. A partial structure such as `struct('max_time',600)` is not a valid
direct input to `dcqp_solve` because the solver requires the other fields too.

```matlab
params = dcqp_default_params();
params.gap_tolerance = 1e-5;
params.max_iterations = 500;
params.max_time = 1800;
params.verbose = true;
params.mosek_quiet = true;

[x_opt,fval,info] = dcqp_solve(Q,d,A,b,Aeq,beq,params);
```

## Algorithm Overview

DCQP employs a **Doubly nonnegative relaxation based Cutting plane** approach:

1. **Cutting Plane Method**: Iteratively adds linear cuts to the original QP problem by solving doubly nonnegative SDPs
2. **Upper Bound Computation**: Uses local search methods for feasible solutions
3. **Lower Bound Computation**: Solves doubly nonnegative relaxations for global lower bounds
4. **Convergence**: Terminates when the relative gap between lower and upper bounds is sufficiently small


### Complete option reference

The following are all fields returned by `dcqp_default_params()`.

#### Stopping and execution

| Parameter | Default | Description |
|---|---:|---|
| `max_iterations` | `300` | Maximum number of outer cutting-plane iterations. |
| `gap_tolerance` | `1e-4` | Target relative gap in the shifted solver objective. If the upper bound is close to zero, an absolute gap is used. |
| `max_time` | `3600` | Wall-clock time limit in seconds. |
| `robust_mode` | `false` | If `false`, rethrow solver errors. If `true`, return with `info.status='error'`. |
| `save_failed_instances` | `true` | Save problem data and diagnostics when an error is caught. |

#### Display and reporting

| Parameter | Default | Description |
|---|---:|---|
| `verbose` | `true` | Print initialization and per-iteration progress. |
| `display_summary` | `true` | Print the final solution or error summary. |
| `mosek_quiet` | `true` | Suppress MOSEK optimizer output while retaining DCQP's own messages. |
| `filename` | `'dcqp_result'` | Instance label shown in the final summary and used by debugging output. |

#### Algorithm controls

| Parameter | Default | Description |
|---|---:|---|
| `dc_regularization` | `1e-5` | Regularization used in the DC decomposition. |
| `psd_check_tolerance` | `1e-8` | Numerical tolerance used by second-order/KKT positive-semidefiniteness checks. |
| `nb_rounds` | `10` | Number of random starts used to obtain an initial upper bound, capped by the number of variables. |
| `eta` | `0.9` | Interpolation factor used to form the relaxed cut target `nuR`. |
| `konnofirst` | `false` | Try a validated generalized Konno cut before the DNN cut in each applicable iteration. |
| `accept_cut_below_threshold` | `false` | When `false`, reject a DNN cut whose validated discarded-region bound remains below `nu`, allowing the generalized Konno fallback to run. When `true`, preserve the legacy below-threshold handling; a DNN cut may be accepted unless a previously added generalized Konno cut makes it redundant. |
| `known_solution` | `[]` | Optional second-order KKT starting point in the original input coordinates. When supplied, DCQP checks feasibility, uses it as the initial incumbent, and skips the random upper-bound search. |
| `do_scaling` | `false` | Automatically rescale the objective from the initial upper bound before the cutting-plane iterations. Results are converted back before return. |

#### Subproblem tolerances and methods

| Parameter | Default | Description |
|---|---:|---|
| `mosek_tolerance` | `1e-8` | MOSEK tolerance used for DNN lower-bound and validation subproblems. |
| `tol_mosek_cut` | `1e-8` | Initial MOSEK tolerance used when generating a DNN cut; the cut routine may retry with adjusted tolerances. |
| `gurobi_qp_tolerance` | `1e-9` | Optimality tolerance passed to convex Gurobi QP subproblems. |
| `gurobi_lp_tolerance` | `1e-9` | Optimality tolerance passed to Gurobi LP subproblems. |
| `gurobi_lp_method` | `1` | Gurobi LP method: `0` primal simplex, `1` dual simplex, `2` barrier. |
| `gurobi_qp_method` | `2` | Gurobi QP method: `-1` automatic, `0` primal simplex, `1` dual simplex, `2` barrier. |

#### Internal or reserved fields

| Parameter | Default | Description |
|---|---:|---|
| `scaling` | `1` | Internal objective scale. Leave this at `1`; `dcqp_solve` updates it when `do_scaling=true`. |
| `lower_bounds` | `[]` | Reserved and currently unused. Encode lower bounds in `A,b`. |
| `upper_bounds` | `[]` | Reserved and currently unused. Encode upper bounds in `A,b`. |

## Datasets

### Existing Test Sets

Existing test-set `.mat` files are included under `data/existing_testsets/`. These problems are drawn from prior computational studies and are organized in four groups:
- **qp20_10**: 20 variables, 10 constraints (16 instances)
- **qp30_15**: 30 variables, 15 constraints (16 instances)  
- **qp40_20**: 40 variables, 20 constraints (16 instances)
- **qp50_25**: 50 variables, 25 constraints (16 instances)



### Newly Generated Synthetic Problems

Synthetic `.mat` files are included under `data/synthetic/`. All distributed
data generators are collected under `data/generators/`; the legacy synthetic
scripts are `generateinstances_uniform.m` and `generateinstances_normal.m`.
The standard synthetic collection contains 140 nonconvex QPs organized in
seven groups (20 instances each):

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

### Structured Paper Examples

The structured examples are stored under `data/structured/`:

- `inexact_stqp/`: 120 inexact StQP instances
- `boxqp/`: 120 structured BoxQP instances
- `many_local_minima/`: 15 StQPs with many local minima

See `data/structured/README.md` for construction details. All MATLAB
generation and preparation functions are collected under `data/generators/`.




## Examples and Demos

### Run the Demo
```matlab
dcqp_demo();  % Runs basic examples with different problem types
```

### Existing Test Sets

```matlab
% Navigate to paper-examples directory first
cd('paper-examples/');

% Solve existing QP test sets from prior computational studies (16 instances in each group)
solve_existing_testsets_with_dcqp('qp20_10');  % 20 variables, 10 constraints
solve_existing_testsets_with_dcqp('qp30_15');  % 30 variables, 15 constraints
solve_existing_testsets_with_dcqp('qp40_20');  % 40 variables, 20 constraints
[summary,result_folder] = solve_existing_testsets_with_dcqp('qp50_25');
```

### Newly Generated Synthetic Problems

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
[summary,result_folder] = solve_synthetic_with_dcqp('qp_u_25_1');

% Test specific instance in a group
[summary,result_folder] = solve_synthetic_with_dcqp('qp_n_0_1',5);
```

### Structured Paper Examples

Each function accepts `"all"`, numeric indices, or one or more instance IDs.
An optional structure controls the time limit and solver settings.

```matlab
solve_inexactStQP_with_dcqp("all");
solve_boxQP_with_dcqp(1:10);
solve_manyLocalMinima_with_dcqp("all");

solve_inexactStQP_with_gurobi(1, struct('max_time', 600));
solve_boxQP_with_gurobi("all");
solve_manyLocalMinima_with_gurobi("all");
```

The structured DCQP wrappers default to the settings used for the frozen
120-instance run: `konnofirst=true`, a known-solution start when available,
objective scaling, `gap_tolerance=1e-4`, and
`accept_cut_below_threshold=false`.

### Comparison With Gurobi

For performance comparison, the package includes Gurobi-based solvers that attempt to solve the same nonconvex QP problems:

```matlab
% Navigate to paper-examples directory first
cd('paper-examples/');

% Existing test sets with Gurobi (with optional time limit)
solve_existing_testsets_with_gurobi('qp20_10');
solve_existing_testsets_with_gurobi('qp30_15',"all",struct('max_time',7200));

% Newly generated synthetic problems with Gurobi
solve_synthetic_with_gurobi('qp_n_0_1');
solve_synthetic_with_gurobi('qp_u_0_1',5);
solve_synthetic_with_gurobi('qp_u_0_1',5,struct('max_time',1800));
```

### Output Results

Every paper-example runner uses the same result schema and writes below:

```text
paper-examples/results/
├── existing-tests/<group>/<solver>_<start-time>/
├── synthetic/<group>/<solver>_<start-time>/
└── structured/<family>/<solver>_<start-time>/
```

**Console Output**: 
- **DCQP Solution Summary** for each instance displaying:
  - Instance name and problem dimensions (n variables, m inequality constraints, meq equality constraints)
  - Solution status message (e.g., "successfully reduced relative gap below 0.0001")
  - Best objective value (scientific notation, e.g., -3.000000e+01)
  - Relative optimality gap (scientific notation, e.g., 4.33e-11)
  - Original-coordinate relative gap after adding back the shift constant
  - Variable shift norm and objective constant introduced by the internal shift
  - Computation time in seconds and total number of iterations
- **Real-time progress** (when `params.verbose = true`): iteration solver details including bounds

Each timestamped run contains `summary.csv`, `summary.mat`, `run.log`, and
one MAT file per instance. The common named columns include initial and final
relative gaps, time, bounds, iterations, number of added DCQP cuts, and
feasibility violation. Detailed DCQP iteration records are stored in
`raw.diagnostics` inside the per-instance MAT file; no separate `testresults/`
or `summary_results/` directory is created. See `paper-examples/README.md` for
the complete schema and calling conventions.


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
- **`qpsolver(Q, d, A, b, Aeq, beq, lb, ub, sol, parameters)`**: Internal cutting-plane solver used by `dcqp_solve`
- **`rescale_constraint_by_slack(a_row, b_value, A, b, Aeq, beq, met_glp, tol_glp)`**: Rescale an inequality row using its maximum feasible slack

## Output Structure

`dcqp_solve` returns the best point in `x_opt`, its objective value in
`fval`, and a diagnostic structure `info` with these fields:

| Field | Meaning |
|---|---|
| `status` | Termination description. |
| `gap` | Relative gap used by the solver in shifted coordinates (absolute gap when the upper bound is close to zero). |
| `original_gap` | Relative gap after restoring the objective constant; reported for reference. |
| `iterations` | Number of completed outer iterations. |
| `time` | Total wall-clock solution time in seconds. |
| `upper_bound`, `lower_bound` | Final bounds in the original objective coordinates. |
| `initial_relative_gap` | Relative gap after the first lower-bound computation. |
| `initial_upper_bound`, `initial_lower_bound` | Initial bounds in the original objective coordinates. |
| `added_cuts` | Total number of inequalities added by DCQP. |
| `dnn_cuts`, `konno_cuts` | Added-cut counts split by DNN and generalized Konno cuts. |
| `cut_records_match` | Whether the detailed cut records agree with the reported cut count. |
| `scaling` | Objective scaling factor used internally. |
| `variable_shift` | Shift `y = x - variable_shift` used internally. |
| `objective_constant` | Constant added when converting shifted bounds back to the original objective. |
| `initialization` | Either `'random'` or `'known_solution'`. |
| `initial_solution` | Initial incumbent in the original input coordinates. |
| `diagnostics` | Detailed per-iteration bounds, timings, cut values, and cut records. |

When `robust_mode=true` and an error occurs, `x_opt=[]`, `fval=Inf`,
`info.status='error'`, and `info.error_message` and `info.error_stack`
describe the failure.

### Status Codes
- `'successfully reduced relative gap below X'`: Solution found within relative gap tolerance X
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

Before solving, DCQP computes variable lower bounds and internally shifts the problem to `y = x - variable_shift`, so that the correction step is applied on a nonnegative feasible set. Returned solutions and objective values are converted back to the original variables.

### Performance Considerations
- **Problem size**: Most efficient for problems with n ≤ 100 variables
- **Constraint density**: Performance degrades with very dense constraint matrices
- **Conditioning**: Ill-conditioned problems may require parameter tuning

### Troubleshooting
- **Unbounded problems**: Verify that the feasible region is bounded
- **Convergence issues**: Try increasing `gap_tolerance` or `max_iterations`
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
├── data/                     # Packaged datasets and their generators
│   ├── existing_testsets/    # Existing test-set .mat files
│   ├── synthetic/            # Random synthetic .mat files
│   ├── structured/           # Structured paper-example .mat files
│   └── generators/           # MATLAB generation and preparation functions
├── utils/                    # Utility functions
│   ├── DC_decomposition.m    # DC decomposition
│   ├── compute_ub.m          # Upper bound computation
│   ├── qpsolver.m            # Internal cutting-plane solver
│   ├── rescale_constraint_by_slack.m # Constraint row rescaling helper
│   └── ...                   # Other utilities  
├── paper-examples/           # Reproducible experiments
├── legacy/                   # Legacy functions
└── tests/                    # Frozen-runtime integrity check
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

This software is distributed under an Academic License for academic research use only. Commercial use is prohibited without explicit written permission from the copyright holders. See `LICENSE` for details.

## Support and Issues

- **Documentation**: See function help: `help dcqp_solve`
- **Examples**: Run `dcqp_demo()` or use the scripts in `paper-examples/`
- **Issues**: Report bugs and feature requests on GitHub
- **Contact**: zhengqu@szu.edu.cn

## Version History

- **v1.1.0** (2026-09-30): Reproducibility and structured-example release
  - Added the structured BoxQP, inexact-StQP, and many-local-minima suites
  - Unified all paper-example runners and result summaries
  - Added generalized Konno-cut support and validated DNN-cut recovery
  - Added MOSEK quiet mode and correct added-cut reporting
  - Added frozen-runtime provenance and integrity verification

- **v1.0.1** (2026-06-03): Maintenance update
  - Added internal variable shifting for nonnegative correction coordinates
  - Added constraint-row rescaling helper for cutting-plane subproblems
  - Fixed solution rescaling logic in the internal solver
  - Updated documentation and repository metadata

- **v1.0.0** (2025-10-03): Initial release
  - Core DCQP algorithm implementation
  - Existing test-set and synthetic test examples
  - Comprehensive documentation

---

**Disclaimer**: This is research software. While extensively tested, use in production environments should be done with appropriate validation.
