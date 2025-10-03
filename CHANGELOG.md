# Changelog

All notable changes to the DCQP project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2025-10-03

### Added
- Initial release of DCQP solver
- Core algorithm implementation for nonconvex quadratic programming
- Main solver function `dcqp_solve()` with comprehensive parameter control
- Utility functions for DC decomposition, bound computation, and cutting plane generation
- 64 benchmark problems from the literature (4 groups: qp20_10, qp30_15, qp40_20, qp50_25)
- 140 synthetic test problems (7 groups with different sparsity and distribution patterns)
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
- Batch processing capabilities for benchmark evaluation

### Documentation
- Complete README with installation, usage, and API documentation
- Function help documentation for all main functions
- Example code for basic and advanced usage scenarios
- Detailed explanation of algorithm parameters and output structure

### Dependencies
- MATLAB R2020a or later
- Gurobi Optimizer (recommended version 12.0.1+)
- MOSEK (recommended version 11.0.12+)

### Files Structure
- Main solver functions in root directory
- Utility functions in `utils/` directory
- Test datasets in `data/benchmark/` and `data/synthetic/`
- Reproducible experiments in `paper-examples/`
- Legacy generation functions in `legacy/`
