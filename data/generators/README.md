# Data generators

This folder contains all MATLAB data-generation functions distributed with
the package.

- `generateinstances_normal.m` and `generateinstances_uniform.m` generate the
  original random synthetic instances in `data/synthetic/`.
- `generate_stqp_inexact.m` is an independent MATLAB implementation of the
  ordinary StQP construction for custom suites. It does not reproduce the
  authors' Julia random stream. The committed 120-instance `inexact_stqp/`
  paper suite was generated with the authors' Julia COP code from
  [Zenodo](https://doi.org/10.5281/zenodo.15113944), without adding the
  cardinality constraint.
- `generate_boxqp_inexact.m` and `generate_stqp_many_local.m` implement the
  other two structured constructions in MATLAB.
- `generate_structured_suite.m` and `prepare_structured_pilot.m` build custom
  structured suites below `data/structured/generated/`.
- The remaining preparation functions generate, reduce, or validate complete
  structured suites using MATLAB.

Generated data should be written below `data/`; experiment results belong
under `paper-examples/results/` and are ignored by Git.
