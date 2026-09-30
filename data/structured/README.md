# Structured paper examples

This folder contains the structured instances used to illustrate the strengths
and limitations of DCQP beyond the original random synthetic examples.

- `inexact_stqp/`: 120 StQPs with guaranteed inexact DNN relaxations—60 at
  original dimension 25 and 60 at original dimension 50.
- `boxqp/`: 120 BoxQPs with guaranteed inexact RLT relaxations—60 each at
  dimensions 25 and 50.
- `many_local_minima/`: 15 saved StQPs with certified many-local-minima
  constructions at original dimensions 12, 24, and 48.

All generation functions are stored together in `../generators/`. Public
solver entry points are in `paper-examples/`.
