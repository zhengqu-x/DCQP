# Paper examples

The experiment entry points in this directory reproduce the three numerical
families in the paper:

- existing test sets (Section 6.2);
- random synthetic instances (Sections 6.3-6.4);
- structured instances: BoxQP, inexact StQP, and many-local-minima StQP
  (Section 7).

## Datasets

### Existing test sets

Existing test-set `.mat` files are included under `../data/existing_testsets/`.
These problems are drawn from prior computational studies and are organized
in four groups:

- **qp20_10**: 20 variables, 10 constraints (16 instances)
- **qp30_15**: 30 variables, 15 constraints (16 instances)
- **qp40_20**: 40 variables, 20 constraints (16 instances)
- **qp50_25**: 50 variables, 25 constraints (16 instances)

### Newly generated synthetic problems

Synthetic `.mat` files are included under `../data/synthetic/`. All
distributed data generators are collected under `../data/generators/`; the
legacy synthetic scripts are `generateinstances_uniform.m` and
`generateinstances_normal.m`. The standard synthetic collection contains 140
nonconvex QPs organized in seven groups (20 instances each):

- **qp_n_0_1**: Normal distribution, density 0.1, no equality constraints
- **qp_n_0_3**: Normal distribution, density 0.3, no equality constraints
- **qp_n_0_9**: Normal distribution, density 0.9, no equality constraints
- **qp_u_0_1**: Uniform distribution, density 0.1, no equality constraints
- **qp_u_0_3**: Uniform distribution, density 0.3, no equality constraints
- **qp_u_0_9**: Uniform distribution, density 0.9, no equality constraints
- **qp_u_25_1**: Uniform distribution, density 0.1, 25 equality constraints

All synthetic problems have:

- **Problem size**: 100 variables, 51 inequality constraints
- **Hessian construction**: Q = L₁ᵀL₁ - L₂ᵀL₂, where L₁ and L₂ are random
  sparse matrices (25×100 each) with entries drawn from uniform or normal
  distributions
- **Constraint generation**:
  - First constraint: A₁ = [1, 1, ..., 1] (sum constraint)
  - Remaining constraints: A₂₋₅₁ are sparse random matrices with the specified
    density
  - Right-hand side: b = A·x₀ + 0.1·rand(), where x₀ is a feasible point
  - Equality constraints: Aeq (when present) is generated similarly, with
    beq = Aeq·x₀
- **Variable bounds**: 0 ≤ x ≤ 1
- **Feasible interior point**: x₀ is a normalized random vector, ensuring that
  all constraints are satisfiable

### Structured paper examples

The structured examples are stored under `../data/structured/`:

- `inexact_stqp/`: 120 inexact StQP instances
- `boxqp/`: 120 structured BoxQP instances
- `many_local_minima/`: 15 StQPs with many local minima

See [`../data/structured/README.md`](../data/structured/README.md) for
construction details. All MATLAB generation and preparation functions are
collected under `../data/generators/`.

## Unified result layout

Every runner writes the same named summary columns and the same per-instance
MAT-file structure below a single `results/` tree:

```text
results/
├── existing-tests/
│   └── <group>/<solver>_<start-time>/
├── synthetic/
│   └── <group>/<solver>_<start-time>/
└── structured/
    └── <family>/<solver>_<start-time>/
```

Each timestamped run directory contains:

- `summary.csv`: named, one-row-per-instance results;
- `summary.mat`: the same table, rows, options, and run identity;
- `<instance-id>.mat`: the solution, returned solver information, detailed
  DCQP diagnostics, options, and instance metadata;
- `run.log`: console output for the complete run.

The common summary fields are `id`, `category`, `family`, `solver`,
`original_n`, `seed`, `status`, `solved`, `time`, `initial_relative_gap`,
`relative_gap`, `iterations`, `added_cuts`, `upper_bound`, `lower_bound`,
`feasibility_violation`, and `error_message`. Fields unavailable for a method
(for example, DCQP cut counts for Gurobi) are `NaN`.

The shared runner automatically applies the appropriate settings for each
experiment category. Structured instances use `konnofirst=true`; synthetic
and existing-test instances use `konnofirst=false`.



## Structured examples

```matlab
[summary,result_folder] = solve_inexactStQP_with_dcqp("all");
[summary,result_folder] = solve_manyLocalMinima_with_dcqp("all");
```

Structured data are stored by family under `../data/structured/`. The first
argument of each structured solver selects which instances to run:

```matlab
% Run every inexact-StQP instance.
[summary,result_folder] = solve_inexactStQP_with_dcqp("all");

% Run the first, fifth, and tenth BoxQP files in natural filename order.
[summary,result_folder] = solve_boxQP_with_dcqp([1 5 10]);

% Run one instance by its filename without the .mat extension.
[summary,result_folder] = solve_boxQP_with_dcqp( ...
    "boxqp_inexact_n25_seed1");

% Run several explicitly named instances.
ids = ["boxqp_inexact_n25_seed1", "boxqp_inexact_n25_seed10"];
[summary,result_folder] = solve_boxQP_with_dcqp(ids);
```

Numeric positions refer only to the selected family and use natural ordering,
so `seed2` comes before `seed10`. Instance IDs are the corresponding MAT-file
names without `.mat`. 

## Synthetic examples

```matlab
% Run all 20 instances in the qp_n_0_3 group.
[summary,result_folder] = solve_synthetic_with_dcqp( ...
    'qp_n_0_3',"all");

% Run instances 1, 5, and 10 in that group.
[summary,result_folder] = solve_synthetic_with_dcqp( ...
    'qp_n_0_3',[1 5 10]);

% Run one instance by its filename without the .mat extension.
[summary,result_folder] = solve_synthetic_with_gurobi( ...
    'qp_n_0_3',"qp_n_0_3_7");

% Run several explicitly named instances.
ids = ["qp_n_0_3_2", "qp_n_0_3_11"];
[summary,result_folder] = solve_synthetic_with_dcqp( ...
    'qp_n_0_3',ids);
```

Synthetic groups are `qp_n_0_1`, `qp_n_0_3`, `qp_n_0_9`, `qp_u_0_1`,
`qp_u_0_3`, `qp_u_0_9`, and `qp_u_25_1`. The equality-elimination diagnostic
is available through `solve_synthetic_with_dcqp_eliminate`.

## Existing test sets

```matlab
% Run all 16 instances in the qp20_10 group.
[summary,result_folder] = solve_existing_testsets_with_dcqp( ...
    'qp20_10',"all");

% Run instances 1, 5, and 10 in that group.
[summary,result_folder] = solve_existing_testsets_with_dcqp( ...
    'qp20_10',[1 5 10]);

% Run one instance by its filename without the .mat extension.
[summary,result_folder] = solve_existing_testsets_with_gurobi( ...
    'qp20_10',"qp20_10_2_3");

% Run several explicitly named instances.
ids = ["qp20_10_1_1", "qp20_10_4_4"];
[summary,result_folder] = solve_existing_testsets_with_dcqp( ...
    'qp20_10',ids);
```

Existing-test groups are `qp20_10`, `qp30_15`, `qp40_20`, and `qp50_25`.
For both synthetic and existing-test runners, the second argument follows the
same selection rules as the structured examples: `"all"`, numeric positions,
one instance ID, or an array of instance IDs.

All MATLAB data generators are under `../data/generators/`. Generated data
belong below `../data/`; experiment outputs belong below `results/`.
