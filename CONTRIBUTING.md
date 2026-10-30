# Contributing

The public repository is
[TangoRomeo80/OFA-Cycle](https://github.com/TangoRomeo80/OFA-Cycle).

Keep changes focused and explain their effect on the mathematical model,
algorithm, or reproducibility. Include the command and input that expose a bug
when reporting one. For simulation findings, include the Python version,
method, actual seed, and number of runs.

## Local workflow

Use Python 3.12 or newer from the repository root:

```sh
python -m unittest discover -s tests -v
python -m cycle_ofa validate
python -m cycle_ofa validate --exhaustive
python -m cycle_ofa reproduce --profile quick --output results/local-check
python -m cycle_ofa camera-ready --profile quick --output results/local-camera-check
```

The package should retain standard-library-only runtime dependencies. Add tests
for new mathematical behavior and meaningful regressions. Exact results should
use rational equality. Statistical tests should use documented tolerances and
fixed seeds, and should not require every confidence interval to cover the
truth or every sample variance ratio to exceed one.

Set `CYCLE_OFA_EXHAUSTIVE=1` when running the unit suite to include exhaustive
ternary canonicalization through length eight. The installed
`validate --exhaustive` command also exercises this mathematical check.

## Scientific invariants

- Gap lengths sum to the cycle circumference whenever a facility is active.
- Selection weights sum to `2*n`, including the one-facility self-gap case.
- A transition decrements exactly one residual capacity.
- Dihedral transformations preserve the exact remaining-cost value.
- A fresh uniform-capacity reduced DP visits the bracelet-count number of
  states, including the terminal state.
- Removing the penultimate active facility creates a single gap of length
  `n`; removing the final facility terminates the trajectory.
- The conditional estimator preserves the mean, rather than the direct
  trajectory's complete distribution.
- The offline optimum cannot exceed greedy's cost on the same arrival
  sequence, allowing for floating-point tolerance. Verify changes against
  the independent small-instance assignment oracle.
- The total-variance identity is `Var(C)=Var(Y)+E[S]`, where `S` includes the
  singleton's positional variance `n*n/48` for every remaining arrival.
- Diagnostic trajectories preserve the original conditional cost and RNG
  consumption for the same input seed.

Offline comparisons use paired arrivals. Report their ratio of sample means
and use the paired covariance when computing its delta-method interval; do not
substitute a mean of per-run ratios. Keep capacity fixed when studying how
variance contributions vary with facility count. A finite-grid trend does not
prove an asymptotic limit.

Changes to seed derivation, draw order, or sampling algorithms can change
simulation rows. Document such changes and version the seed scheme when its
definition changes. Never edit generated values manually to make them agree
with a previous table. Generate results with the intended configuration and
retain the corresponding manifest.
The current SHA-256 stream labels are `direct`, `conditional`,
`offline-comparison`, and `variance-diagnostic`; v0.2.0 preserves the first
two streams and their sampling order.

Keep regenerated local experiments outside the checked-in release results
unless updating those results intentionally. Before proposing new release
tables, run the full paper profile and explain any changes in exact or
simulation values. Also run
`python -m cycle_ofa camera-ready --runs 10000 --seed 20260822 --output results/camera-ready`
when updating the offline or variance tables. Commit the corresponding CSVs,
generated table file, and manifest together.

Contributions are made under the project's MIT License. The separately
distributed paper and publisher's files are outside that software license.
