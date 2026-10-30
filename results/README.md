# Included experiment results

These results accompany software v0.2.0 at
[TangoRomeo80/OFA-Cycle](https://github.com/TangoRomeo80/OFA-Cycle).

`paper/` contains a completed full-profile run on Python 3.12.14:

- 33 exact rational instances covering the manuscript's entire exact grid.
- All 24 displayed per-arrival costs agree with the supplied manuscript to
  its four-decimal precision.
- All computed memo-table sizes agree with the bracelet counts, including
  the six instances displayed in the manuscript's state-count table.
- Five simulation instances with 10,000 trajectories per method, using
  base seed `20260822` and the documented `cycle-ofa-v1` seed derivation.

`quick/` contains the smaller demonstration profile, with 1,000 trajectories
per method. Its tables cover only that profile's subset of instances.

`camera-ready/` contains a completed 10,000-run paper-profile experiment:

- `offline.csv` compares greedy and the offline optimum on the same arrival
  sequences for `(8,2)`, `(8,4)`, `(20,5)`, `(50,5)`, and `(200,10)`. It records
  paired means, covariance, excess cost, and confidence intervals. The reported
  ratio divides the two sample means; it is neither the mean of per-run ratios
  nor a worst-case competitive-ratio guarantee.
- `variance.csv` uses fixed `m=5` and `n=20,50,100,200`. It records conditional
  sample variance, mean removed positional variance and its confidence interval,
  implied direct variance, and the diagnostic point estimate
  `1+mean_removed_variance/conditional_variance`. No confidence interval or
  asymptotic limit is asserted for this ratio.
- `reviewer_tables.tex` contains the generated offline and variance tables.

Each results directory contains CSV data, generated LaTeX tables, and a
manifest with the run configuration, environment, actual seeds (also in the
Monte Carlo CSVs), and SHA-256 checksums of its three data/table output files.
Repeating a run with the same
code and Python version reproduces the data; the manifest timestamp changes.
Exact results do not depend on random seeds. Monte Carlo samples depend on
the implementation and Python version as well as the recorded seeds.

Seeds are the first eight SHA-256 bytes, interpreted as a big-endian integer,
of `cycle-ofa-v1:{base_seed}:{n}:{m}:{method}`. Stream labels are `direct`,
`conditional`, `offline-comparison`, and `variance-diagnostic`.

From the repository root, regenerate with:

```sh
python -m cycle_ofa reproduce --profile paper --output results/paper
python -m cycle_ofa reproduce --profile quick --output results/quick
python -m cycle_ofa camera-ready --runs 10000 --seed 20260822 --output results/camera-ready
```

For a short check of the new experiments, use
`python -m cycle_ofa camera-ready --profile quick --output results/local-camera-check`.
It uses `(8,2)` and `(12,2)` with 100 runs per case.

Each run is staged before replacing existing results. A failed computation
preserves a previously completed run. A manifest is published last, so a
partially published output does not carry an obsolete success manifest.
