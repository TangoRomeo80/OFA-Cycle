"""Exercise reviewer exports, seed provenance, and failure-safe publication."""

from contextlib import redirect_stderr, redirect_stdout
import csv
import hashlib
from io import StringIO
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from cycle_ofa.camera_ready import camera_ready
from cycle_ofa.cli import main
from cycle_ofa.experiments import derive_seed


class CameraReadyTests(unittest.TestCase):
    def test_computed_reproducible_exports_and_seeds(self):
        with TemporaryDirectory() as first, TemporaryDirectory() as second:
            a = camera_ready(first, profile="quick", runs=12, seed=51)
            b = camera_ready(second, profile="quick", runs=12, seed=51)
            self.assertEqual(a["output_sha256"], b["output_sha256"])
            for name, expected in a["output_sha256"].items():
                self.assertEqual(hashlib.sha256((Path(first) / name).read_bytes()).hexdigest(), expected)
            with (Path(first) / "offline.csv").open(encoding="utf-8", newline="") as stream:
                offline_rows = list(csv.DictReader(stream))
            with (Path(first) / "variance.csv").open(encoding="utf-8", newline="") as stream:
                variance_rows = list(csv.DictReader(stream))
            self.assertEqual(len(offline_rows), 2)
            self.assertEqual(len(variance_rows), 2)
            for row in offline_rows:
                n, m = int(row["n"]), int(row["m"])
                self.assertEqual(int(row["seed"]), derive_seed(51, n, m, "offline-comparison"))
                self.assertLessEqual(float(row["offline_mean"]), float(row["greedy_mean"]))
                self.assertGreaterEqual(float(row["ratio_of_means"]), 1)
            for row in variance_rows:
                n, m = int(row["n"]), int(row["m"])
                self.assertEqual(int(row["seed"]), derive_seed(51, n, m, "variance-diagnostic"))
                self.assertAlmostEqual(float(row["implied_direct_variance"]), float(row["conditional_variance"]) + float(row["mean_removed_variance"]))
            tables = (Path(first) / "reviewer_tables.tex").read_text(encoding="utf-8")
            self.assertEqual(tables.count(r"\begin{table}"), 2)
            self.assertIn(r"\label{tab:offline}", tables)
            self.assertIn(r"\label{tab:variance}", tables)
            self.assertIn("12 runs", tables)
            self.assertIn("$m=2$", tables)

    def test_cli_json_and_default_quick_runs(self):
        with TemporaryDirectory() as output:
            stdout, stderr = StringIO(), StringIO()
            with redirect_stdout(stdout), redirect_stderr(stderr):
                self.assertEqual(main(["camera-ready", "--profile", "quick", "--output", output]), 0)
            result = json.loads(stdout.getvalue())
            self.assertEqual(result["runs_per_case"], 100)
            self.assertEqual(result["profile"], "quick")
            self.assertIn("offline comparison", stderr.getvalue())

    def test_invalid_inputs_leave_no_directory(self):
        with TemporaryDirectory() as parent:
            output = Path(parent) / "unused"
            for kwargs in ({"profile": "typo"}, {"runs": 1}, {"seed": True}):
                with self.assertRaises(ValueError):
                    camera_ready(output, **kwargs)
                self.assertFalse(output.exists())

    def test_failed_rerun_preserves_published_results(self):
        with TemporaryDirectory() as output:
            camera_ready(output, profile="quick", runs=5)
            old_files = {p.name: p.read_bytes() for p in Path(output).iterdir()}
            with patch("cycle_ofa.camera_ready.diagnose_variance", side_effect=RuntimeError("interrupted")):
                with self.assertRaisesRegex(RuntimeError, "interrupted"):
                    camera_ready(output, profile="quick", runs=6)
            self.assertEqual(old_files, {p.name: p.read_bytes() for p in Path(output).iterdir()})


if __name__ == "__main__":
    unittest.main()
