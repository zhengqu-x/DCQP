# Structured paper examples

This folder contains the structured instances used to illustrate the strengths
and limitations of DCQP beyond the original random synthetic examples.

- `inexact_stqp/`: 120 StQPs with guaranteed inexact DNN relaxations—60 at
  original dimension 25 and 60 at original dimension 50. These instances were
  generated with the authors' Julia COP code archived on
  [Zenodo](https://doi.org/10.5281/zenodo.15113944), without its subsequent
  cardinality constraint, and then reduced by eliminating one simplex
  variable.
- `boxqp/`: 120 BoxQPs with guaranteed inexact RLT relaxations—60 each at
  dimensions 25 and 50.
- `many_local_minima/`: 1 saved StQP with a certified many-local-minima
  construction at original dimension 25. The saved model has dimension 24
  after eliminating one simplex variable.

All generation functions are stored together in `../generators/`. Public
solver entry points are in `paper-examples/`.
