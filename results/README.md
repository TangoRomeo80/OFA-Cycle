# Included experiment results

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

The simulation rows are newly generated results. The original manuscript's
per-row seed derivation and program were unavailable, so its old Monte Carlo
rows are not reconstructed digit for digit. Exact results are unaffected by
the choice of random seed. The separate manuscript files have not been edited.

Each directory contains CSV data, generated LaTeX tables, and a manifest with
the run configuration, environment, actual seeds (in the simulation CSV),
and SHA-256 checksums of its three output files. Repeating a run with the same
code and Python version reproduces the data; the manifest timestamp changes.

From the repository root, regenerate with:

```sh
python -m cycle_ofa reproduce --profile paper --output results/paper
python -m cycle_ofa reproduce --profile quick --output results/quick
```

Each run is staged before replacing existing results. A failed computation
preserves a previously completed run. A manifest is published last, so a
partially published output does not carry an obsolete success manifest.
