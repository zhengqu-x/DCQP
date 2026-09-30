# Algorithm provenance

The numerical runtime in this package is copied byte-for-byte from the
frozen source snapshot executed by the 120-instance inexact-StQP run that
started on September 27, 2026 at 20:10 and finished at 21:51.

The frozen runtime consists of:

- the six root-level `dcqp_*.m` files;
- every MATLAB file under `utils/`;
- `legacy/save_failed_instance.m`.

No solver or numerical-helper source was taken from the later refactored
`git-prepare` implementation. The approved post-snapshot changes do not alter
the optimization model or solver tolerances:

- the display-only `mosek_quiet` option selects MOSEK `echo(0)` instead of
  `echo(5)`;
- objective bounds and cut targets use four digits after the decimal point
  in scientific notation, while lower-level numerical diagnostics use two;
- each outer iteration is displayed as one compact block, without repeating
  fixed solver tolerances;
- reporting records the number of DNN and generalized Konno inequalities
  actually appended to the retained region.
- `qpsolver` returns its iteration diagnostics to `dcqp_solve` instead of
  writing timestamped MAT files relative to the current working directory;
  this changes result storage only, not any numerical branch or tolerance.
- `dcqp_version` contains the package release number and date only; updating
  this metadata does not affect the numerical runtime.

The cut-count reporting also corrects `alpha_record` to zero when independent
validation terminates on an empty retained region, because no inequality is
appended in that branch. This is bookkeeping only and does not change the
algorithmic state or numerical result.

All paper-example entry points use one result writer under
`paper-examples/results/`. The legacy synthetic and existing-test Gurobi
settings are retained (`OptimalityTol`, `BarConvTol`, `LPWarmStart`, and the
same time limit), so the directory and schema unification does not change the
solver model used by those scripts.
