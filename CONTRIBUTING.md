# Contributing

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
```

The package should retain standard-library-only runtime dependencies. Add tests
for new mathematical behavior and meaningful regressions. Exact results should
use rational equality. Statistical tests should use documented tolerances and
fixed seeds, and should not require every confidence interval to cover the
truth or every sample variance ratio to exceed one.

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

Changes to seed derivation, draw order, or sampling algorithms can change
simulation rows. Document such changes and version the seed scheme when its
definition changes. Never edit generated values manually to make them agree
with a previous table. Generate results with the intended configuration and
retain the corresponding manifest.

Keep regenerated local experiments outside the checked-in release results
unless updating those results intentionally. Before proposing new release
tables, run the full paper profile and explain any changes in exact or
simulation values.

Contributions are made under the project's MIT License. The separately
distributed paper and publisher's files are outside that software license.
