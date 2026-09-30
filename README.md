# DCQP: Doubly Nonnegative based Cutting Plane method for Quadratic Programming

[![MATLAB](https://img.shields.io/badge/MATLAB-R2020a+-blue.svg)](https://www.mathworks.com/products/matlab.html)
[![License](https://img.shields.io/badge/License-Academic-green.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-1.0.1-orange.svg)](https://github.com/zhengqu-x/DCQP)

## Overview

DCQP is a MATLAB package for solving **nonconvex quadratic programming** problems using a doubly nonnegative relaxation and cutting plane approach.

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
- ✅ **Global Optimality Certificate**: Provides a certificate of global optimality within a specified tolerance using cutting-plane methods
- ✅ **Bounded Feasible Region**: Requires bounded constraint sets
- ✅ **Multiple Solvers**: Integrates with Gurobi and MOSEK for subproblems
- ✅ **Comprehensive Examples**: Includes scripts to reproduce the experiments in the associated paper


## Requirements

### Software Dependencies
- **MATLAB** R2020a or later
- **Gurobi Optimizer** (recommended version 12.0.1 or later)
- **MOSEK** 11.0.30 or earlier with the MATLAB interface configured. Newer MOSEK releases such as 11.2 are not currently supported because DCQP currently uses MOSEK's legacy `mosekopt` MATLAB toolbox interface and semidefinite-programming data structures.


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

To reproduce experiments in the associated paper, including
DCQP/Gurobi comparisons and result-file documentation, see
[`paper-examples/README.md`](paper-examples/README.md).

### Variable bounds

Variable bounds must be included in `A*x <= b`. For example,
to impose `lb <= x <= ub`:

```matlab
n = size(Q,1);
A = [A; eye(n); -eye(n)];
b = [b; ub(:); -lb(:)];
```

### Setting options

Always begin with the complete default structure, then override individual
fields.

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
