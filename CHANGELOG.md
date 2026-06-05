# Changelog

All notable changes to the DCQP project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.1] - 2026-06-03

### Added
- Added `utils/rescale_constraint_by_slack.m` to rescale inequality rows using maximum feasible slack.
- Added automatic lower-bound variable shifting before solving so correction steps operate on nonnegative transformed variables while returned solutions remain in the original coordinates.

### Fixed
- Corrected internal solution rescaling in `qpsolver.m` when the initial scaling factor is not applied.
- Improved first-iteration cut selection in `qpsolver.m` by reseeding KKT search from the incumbent solution when the initial KKT point is worse than the current best upper bound.
- Computed reported optimality gap in shifted solver coordinates to reduce sensitivity to additive objective constants while preserving the paper's relative-gap formula.

### Documentation
- Updated README metadata, utility function descriptions, dataset notes, and license wording.

## [1.0.0] - 2025-10-03

### Added
- Initial release of DCQP solver
- Core algorithm implementation for nonconvex quadratic programming
- Main solver function `dcqp_solve()` with comprehensive parameter control
- Utility functions for DC decomposition, bound computation, and cutting plane generation
- 64 existing test-set problems from prior computational studies (4 groups: qp20_10, qp30_15, qp40_20, qp50_25)
- 140 newly generated synthetic test problems (7 groups with different sparsity and distribution patterns)
- Comparison framework with Gurobi optimizer
- Comprehensive demonstration examples in `dcqp_demo()`
- Automatic environment setup with `dcqp_startup()`
- Result logging and analysis capabilities
- Error handling and debugging features with failed instance saving

### Features
- Doubly nonnegative relaxation-based cutting plane method
- Global optimization for indefinite quadratic forms
- Integration with Gurobi and MOSEK solvers for subproblems
- Automatic problem scaling and conditioning
- Verbose output options for algorithm monitoring
- Batch processing capabilities for computational test sets

### Documentation
- Complete README with installation, usage, and API documentation
- Function help documentation for all main functions
- Example code for basic and advanced usage scenarios
- Detailed explanation of algorithm parameters and output structure

### Dependencies
- MATLAB R2020a or later
- Gurobi Optimizer (recommended version 12.0.1+)
- MOSEK 11.0.30 or earlier with the MATLAB interface configured

### Files Structure
- Main solver functions in root directory
- Utility functions in `utils/` directory
- Test datasets in `data/existing_testsets/` and `data/synthetic/`
- Reproducible experiments in `paper-examples/`
- Legacy generation functions in `legacy/`
