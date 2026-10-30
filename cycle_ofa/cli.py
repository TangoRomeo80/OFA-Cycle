"""Command-line entry points; JSON to stdout, progress/errors to stderr."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import sys

from . import __version__
from .camera_ready import camera_ready
from .exact import ExactSolver
from .experiments import BASE_SEED, derive_seed, reproduce
from .simulation import simulate
from .symmetry import bracelet_count
from .validation import run_validation


def _positive(value: str) -> int:
    try:
        number = int(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("expected a positive integer") from error
    if number < 1:
        raise argparse.ArgumentTypeError("expected a positive integer")
    return number


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version=f"cycle-ofa {__version__}")
    commands = parser.add_subparsers(dest="command", required=True)
    exact = commands.add_parser("exact", help="compute a Fraction-valued DP result")
    exact.add_argument("n", type=_positive)
    exact.add_argument("m", type=_positive)
    exact.add_argument("--no-symmetry", action="store_true", help="use the raw residual-state DP")
    exact.add_argument("--max-states", type=_positive, default=1_000_000)
    simulation = commands.add_parser("simulate", help="estimate total cost by Monte Carlo")
    simulation.add_argument("n", type=_positive)
    simulation.add_argument("m", type=_positive)
    simulation.add_argument("--runs", type=_positive, default=10_000)
    simulation.add_argument("--seed", type=int, default=BASE_SEED)
    simulation.add_argument("--method", choices=("direct", "conditional", "both"), default="both")
    reproduction = commands.add_parser("reproduce", help="generate experiment CSVs and LaTeX tables")
    reproduction.add_argument("--profile", choices=("quick", "paper"), default="quick")
    reproduction.add_argument("--runs", type=_positive, default=None)
    reproduction.add_argument("--seed", type=int, default=BASE_SEED)
    reproduction.add_argument("--output", type=Path, default=None)
    reproduction.add_argument("--max-states", type=_positive, default=1_000_000)
    reviewer = commands.add_parser("camera-ready", help="generate offline baselines and variance diagnostics")
    reviewer.add_argument("--profile", choices=("quick", "paper"), default="paper")
    reviewer.add_argument("--runs", type=_positive, default=None)
    reviewer.add_argument("--seed", type=int, default=BASE_SEED)
    reviewer.add_argument("--output", type=Path, default=Path("results") / "camera-ready")
    validation = commands.add_parser("validate", help="run packaged mathematical self-checks")
    validation.add_argument("--exhaustive", action="store_true", help="check ternary words through length 8")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "exact":
            solver = ExactSolver(args.n, symmetry=not args.no_symmetry, max_states=args.max_states)
            value = solver.solve(args.m)
            result = dict(
                n=args.n, m=args.m, fraction=str(value), numerator=value.numerator,
                denominator=value.denominator, total_cost=float(value),
                per_arrival=float(value / (args.n * args.m)), symmetry=solver.symmetry,
                stored_states=solver.state_count, raw_states=(args.m + 1) ** args.n,
                bracelet_count=bracelet_count(args.n, args.m + 1),
            )
        elif args.command == "simulate":
            if args.method == "both":
                result = [asdict(simulate(
                    args.n, args.m, runs=args.runs,
                    seed=derive_seed(args.seed, args.n, args.m, method), method=method,
                )) for method in ("direct", "conditional")]
            else:
                result = asdict(simulate(args.n, args.m, runs=args.runs,
                                        seed=args.seed, method=args.method))
        elif args.command == "reproduce":
            destination = args.output or Path("results") / args.profile
            manifest = reproduce(
                destination, profile=args.profile, runs=args.runs, seed=args.seed,
                max_states=args.max_states, progress=lambda message: print(message, file=sys.stderr),
            )
            result = dict(output=str(destination.resolve()), **manifest)
        elif args.command == "camera-ready":
            manifest = camera_ready(
                args.output, profile=args.profile, runs=args.runs, seed=args.seed,
                progress=lambda message: print(message, file=sys.stderr),
            )
            result = dict(output=str(args.output.resolve()), **manifest)
        else:
            result = run_validation(exhaustive=args.exhaustive)
    except (ValueError, RuntimeError, OSError) as error:
        parser.exit(2, f"cycle-ofa: {error}\n")
    print(json.dumps(result, indent=2, allow_nan=False))
    return 0
