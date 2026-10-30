# cycle-ofa

Research software for **Multi-Capacity Greedy Online Facility Assignment on
Cycles: Exact Gap Recurrences and Symmetry-Reduced Algorithms**, by Md. Rawha
Siddiqi Riad, Md. Tanzeem Rahat, and Md. Manzurul Hasan, accepted at ICAA 2027.

The package implements the paper's exact rational dynamic program, dihedral
symmetry reduction, bracelet counts, direct continuous-arrival simulation, and
conditional gap simulation. It also exports experiment data and LaTeX tables.
The runtime uses only the Python standard library.

## Requirements and quick start

Use Python 3.12 or newer. Run these commands from this directory:

```sh
python -m cycle_ofa exact 3 2
python -m cycle_ofa simulate 8 2 --runs 10000 --method both
python -m cycle_ofa validate
python -m cycle_ofa reproduce --profile quick
```

The exact computation returns the total expected assignment cost. For `(3,2)`
the exact value is `67/27`; the per-arrival value is `67/162`. Commands print
JSON to standard output. Experiment commands additionally write the files
described below.

On Windows, `py -3.12` can replace `python` if that is how Python is installed.

Installation is optional when running from the source directory. To install
the package and its `cycle-ofa` command:

```sh
python -m pip install .
cycle-ofa exact 3 2
```

Installation uses setuptools as a build dependency; the installed software has
no third-party runtime dependencies.

## Model and implementation

There are `n >= 3` facilities at the vertices of a unit-edge cycle of
circumference `n`. Each initially has integer capacity `m >= 1`. Exactly `n*m`
independent arrivals are uniform on the continuous cycle. Greedy assigns each
arrival to a nearest facility with positive remaining capacity. A fixed tie
rule is immaterial for these continuous expectations.

For the clockwise gaps `g` between active facilities, the immediate expected
cost is `sum(g_i**2)/(4*n)`. Facility `i` is selected with probability
`(g_(i-1)+g_i)/(2*n)`. Its conditional expected cost is
`(g_(i-1)**2+g_i**2)/(4*(g_(i-1)+g_i))`. A single active facility has one gap
of length `n` and expected cost `n/4` per remaining assignment.

- **Exact DP** uses `fractions.Fraction`, so its reported fractions are exact.
  An explicit stack avoids dependence on Python's recursion limit. The default
  keys are the least rotations of a residual-capacity vector or its reversal,
  computed using Duval's linear-time word algorithm. `--no-symmetry` provides
  an independent state-space comparison.
- **Direct simulation** samples continuous arrival coordinates and finds a
  nearest active facility using an ordered list.
- **Conditional simulation** samples the selected facility from integer gap
  weights, then sums conditional expected costs. A Fenwick tree and circular
  neighbor links implement the paper's `O(n*m*log(n))` trajectory algorithm.

The conditional estimator has the same expectation as the direct total cost
and no greater population variance per trajectory. Individual samples can
show different variance ordering. The estimators have different distributions,
and variance reduction alone does not establish a wall-clock speed advantage.

Exact computation remains exponential. The CLI defaults to a limit of
1,000,000 states; `--max-states` changes that limit. Arbitrary-precision
numerators and denominators also grow with the problem size.

## Commands

```sh
python -m cycle_ofa exact 8 2
python -m cycle_ofa exact 5 2 --no-symmetry --max-states 10000
python -m cycle_ofa simulate 20 5 --method conditional --runs 10000 --seed 20260822
python -m cycle_ofa simulate 20 5 --method both --runs 10000 --seed 20260822
python -m cycle_ofa reproduce --profile quick --output results/quick
python -m cycle_ofa reproduce --profile paper --output results/paper
python -m cycle_ofa validate --exhaustive
```

`simulate` defaults to both methods, 10,000 runs, and seed `20260822`.
With `--method both`, this is a base seed used to derive separate random
streams. With `--method direct` or `--method conditional`, the supplied seed
is passed directly to that method's random generator.
The reported variances are unbiased sample variances, and confidence interval
half-widths use the paper's normal approximation `1.96*s/sqrt(runs)`.

Use `python -m cycle_ofa --help` or append `--help` to a subcommand for its
arguments. `reproduce --runs` overrides the profile's simulation run count;
the resulting output is then a custom-run experiment.

## Reproducing the experiments

The **quick** profile is a smaller demonstration. It computes capacities
`m=1,2` for `n=3..6`, plus `(8,2)`, and simulates `(6,1)` and `(8,2)` with
1,000 runs per method. Its tables cover this smaller grid and do not reproduce
the full manuscript.

The **paper** profile computes the manuscript's exact grid:

| Capacity | Facility counts |
| --- | --- |
| `m=1,2` | `n=3..12` |
| `m=3` | `n=3..9` |
| `m=4` | `n=3..8` |

It simulates `(8,2)`, `(8,4)`, `(20,5)`, `(50,5)`, and `(200,10)`, with
10,000 runs per method. This profile performs substantially more work than
the quick profile. The expected-cost table displays `n=3..8` and `m=1..4`,
while the raw exact data retain the complete grid.

Each reproduction writes:

| File | Contents |
| --- | --- |
| `exact.csv` | Exact total costs, per-arrival costs, and state counts |
| `simulation.csv` | Means, sample variances, interval widths, and actual seeds |
| `paper_tables.tex` | LaTeX tables generated from this run |
| `manifest.json` | Configuration and environment information for the run |

The LaTeX tables use the manuscript's table labels and require `booktabs`.
The software writes to its selected results directory; it does not replace
the separate manuscript's existing tables. To use a regenerated table file,
review it and copy it into the manuscript's `generated` directory.

### Random-seed provenance

The supplied manuscript states base seed `20260822` but does not specify its
original seed derivation or provide the original Python program. For
`reproduce` and `simulate --method both`, this package therefore defines an
explicit, versioned seed scheme. For each `(n,m,method)`, the actual seed is
the big-endian integer represented by the first eight bytes of SHA-256 of the
ASCII string:

```text
cycle-ofa-v1:{base_seed}:{n}:{m}:{method}
```

Here `method` is `direct` or `conditional`; braces are replaced by the decimal
values and method name. The actual seeds are saved with the results. Each
method uses a local `random.Random` instance, and Python's process-randomized
`hash()` is not used.

The exact tables should agree with the manuscript. Newly generated simulation
rows are reproducible results of this implementation; they are not claimed to
match the manuscript's old rows digit for digit. Record the Python version
alongside seeds because pseudorandom sampling details can differ across Python
versions. Report regenerated simulation rows when citing this implementation.

## Python API

```python
from fractions import Fraction

from cycle_ofa.exact import ExactSolver
from cycle_ofa.simulation import simulate
from cycle_ofa.symmetry import bracelet_count, canonical

solver = ExactSolver(n=3, max_states=100_000)
assert solver.solve(m=2) == Fraction(67, 27)
assert solver.state_count == bracelet_count(3, 3)

# Expected remaining cost from a particular residual-capacity vector.
remaining_cost = solver.value((2, 0, 1))
representative = canonical((2, 0, 1))
statistics = simulate(8, 2, runs=10000, seed=12345, method="conditional")
```

The API's `simulate(..., seed=...)` and single-method CLI commands accept the
actual random-generator seed. `simulate --method both` and `reproduce` derive
per-method seeds from a base seed as documented above.
`ExactSolver.state_count` includes the terminal state and
counts the solver's accumulated memo entries. A fresh solver starting from a
uniform capacity visits exactly `bracelet_count(n, m+1)` states with symmetry,
or `(m+1)**n` without symmetry.

## Validation

```sh
python -m cycle_ofa validate
python -m cycle_ofa validate --exhaustive
python -m unittest discover -s tests -v
```

The installed `validate` command checks rational benchmark values, raw versus
reduced DP, bracelet counts, canonicalization, and the constant conditional
trajectory for `(3,1)`. The exhaustive option additionally compares
canonicalization with its direct definition for all ternary words through
length eight and randomized larger words. The separate source test suite also
checks gap boundary cases, Fenwick sampling, and simulation behavior.
Statistical confidence intervals are evidence with sampling uncertainty,
rather than deterministic correctness certificates.

The central regression values from the paper are:

```text
G(3,1) = 17/12
G(4,1) = 71/32
G(5,1) = 1569/500
G(3,2) = 67/27
```

## Citation and license

The software citation is in [CITATION.cff](CITATION.cff). Add the public
repository URL and an archived-release DOI there once they exist; neither is
invented in this source package. Cite the associated paper as well when using
its mathematical results.

The software is released under the [MIT License](LICENSE). This license does
not grant rights to the separate manuscript or publisher's materials.
