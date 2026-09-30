# Data generators

This folder contains all MATLAB data-generation functions distributed with
the package.

- `generateinstances_normal.m` and `generateinstances_uniform.m` generate the
  original random synthetic instances in `data/synthetic/`.
- `generate_stqp_inexact.m`, `generate_boxqp_inexact.m`, and
  `generate_stqp_many_local.m` implement the three structured constructions.
- `generate_structured_suite.m` and `prepare_structured_pilot.m` build custom
  structured suites below `data/structured/generated/`.
- The remaining preparation functions generate, reduce, or validate complete
  structured suites using MATLAB.

Generated data should be written below `data/`; experiment results belong
under `paper-examples/results/` and are ignored by Git.
