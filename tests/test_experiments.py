"""Check public commands, exports, seed provenance, and reproducible files."""

from contextlib import redirect_stderr, redirect_stdout
import csv
import hashlib
from io import StringIO
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from cycle_ofa.cli import main
from cycle_ofa.experiments import derive_seed, reproduce


class ExperimentTests(unittest.TestCase):
    def test_seed_derivation(self):
        message = b"cycle-ofa-v1:20260822:8:2:direct"
        expected = int.from_bytes(hashlib.sha256(message).digest()[:8], "big")
        self.assertEqual(derive_seed(20260822, 8, 2, "direct"), expected)
        self.assertNotEqual(expected, derive_seed(20260822, 8, 2, "conditional"))

    def test_exports_are_computed_and_reproducible(self):
        with TemporaryDirectory() as first, TemporaryDirectory() as second:
            a = reproduce(first, runs=10)
            b = reproduce(second, runs=10)
            self.assertEqual(a["output_sha256"], b["output_sha256"])
            for name, expected_hash in a["output_sha256"].items():
                data = (Path(first) / name).read_bytes()
                self.assertEqual(hashlib.sha256(data).hexdigest(), expected_hash)
            with (Path(first) / "exact.csv").open(encoding="utf-8", newline="") as stream:
                rows = list(csv.DictReader(stream))
            triangle = next(row for row in rows if row["n"] == "3" and row["m"] == "2")
            self.assertEqual((triangle["numerator"], triangle["denominator"]), ("67", "27"))
            self.assertTrue(all(row["stored_states"] == row["bracelet_count"] for row in rows))
            tables = (Path(first) / "paper_tables.tex").read_text(encoding="utf-8")
            self.assertEqual(tables.count(r"\begin{table}"), 3)
            self.assertIn("10 runs", tables)
            self.assertIn("0.4136", tables)

    def test_cli_json_and_error(self):
        out = StringIO()
        with redirect_stdout(out):
            self.assertEqual(main(["exact", "3", "2"]), 0)
        self.assertEqual(json.loads(out.getvalue())["fraction"], "67/27")
        out = StringIO()
        with redirect_stdout(out):
            main(["simulate", "3", "1", "--runs", "5", "--method", "both"])
        stats = json.loads(out.getvalue())
        self.assertNotEqual(stats[0]["seed"], stats[1]["seed"])
        with redirect_stderr(StringIO()), self.assertRaises(SystemExit) as caught:
            main(["simulate", "3", "1", "--runs", "1"])
        self.assertEqual(caught.exception.code, 2)

    def test_bad_profile_does_not_create_output(self):
        with TemporaryDirectory() as parent:
            output = Path(parent) / "unused"
            with self.assertRaises(ValueError):
                reproduce(output, profile="typo")
            self.assertFalse(output.exists())

    def test_failed_rerun_preserves_published_results(self):
        with TemporaryDirectory() as output:
            reproduce(output, runs=4)
            old_files = {p.name: p.read_bytes() for p in Path(output).iterdir()}
            with patch("cycle_ofa.experiments.simulate", side_effect=RuntimeError("interrupted")):
                with self.assertRaisesRegex(RuntimeError, "interrupted"):
                    reproduce(output, runs=5)
            self.assertEqual(old_files, {p.name: p.read_bytes() for p in Path(output).iterdir()})


if __name__ == "__main__":
    unittest.main()
