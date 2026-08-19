#!/usr/bin/env python3
"""Locate the annotation richness at which caching starts returning time.

A fast annotator on a small input can lose time under VCFcache because the work
caching avoids is smaller than the caching layer's own cost. Adding supplementary
annotation databases increases per-variant work and should eventually cross that
boundary. This script measures where.

Each database set is timed on two inputs from the same genome. The panel input is
small enough that its runtime is essentially the annotator's start-up with that
database set loaded, so subtracting it isolates the per-variant work. That
distinction matters: loading more databases raises start-up too, and caching
cannot return start-up.
"""

from __future__ import annotations

import argparse
import csv
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Any

from benchmarks.run_strategy_comparison import write_tsv


def load(path: Path) -> dict[tuple[str, str], list[float]]:
    """Group wall times by database set and assay."""
    grouped: dict[tuple[str, str], list[float]] = defaultdict(list)
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            grouped[(row["arm"], row["assay"])].append(float(row["wall_seconds"]))
    return grouped


def analyse(
    grouped: dict[tuple[str, str], list[float]],
    vcfcache_cost: float,
    hit_rate: float,
) -> list[dict[str, Any]]:
    """Derive per-variant work and the predicted cached outcome for each arm."""
    arms = sorted({key[0] for key in grouped}, key=lambda a: (a != "core", a))
    required = vcfcache_cost / hit_rate
    core_variable = None
    rows = []
    for arm in arms:
        panel = grouped.get((arm, "panel"))
        wes = grouped.get((arm, "wes"))
        if not panel or not wes:
            continue
        startup = statistics.median(panel)
        direct = statistics.median(wes)
        variable = direct - startup
        if core_variable is None:
            core_variable = variable
        # Caching leaves the start-up and the misses, and adds its own cost.
        cached = startup + (1 - hit_rate) * variable + vcfcache_cost
        rows.append(
            {
                "arm": arm,
                "repeats": len(wes),
                "startup_seconds": startup,
                "startup_cv_percent": (
                    statistics.stdev(panel) / statistics.mean(panel) * 100
                    if len(panel) > 1
                    else 0.0
                ),
                "direct_wes_seconds": direct,
                "direct_cv_percent": (
                    statistics.stdev(wes) / statistics.mean(wes) * 100
                    if len(wes) > 1
                    else 0.0
                ),
                "variable_work_seconds": variable,
                "multiple_of_core": variable / core_variable if core_variable else float("nan"),
                "required_variable_seconds": required,
                "predicted_cached_seconds": cached,
                "predicted_speedup": direct / cached if cached > 0 else float("nan"),
                "predicted_win": direct > cached,
            }
        )
    return rows


def main() -> None:
    """Print the richness ladder and the first arm that breaks even."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metrics", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--vcfcache-cost",
        type=float,
        default=26.5,
        help="Measured VCFcache cost for this input size, in seconds",
    )
    parser.add_argument("--hit-rate", type=float, default=0.932)
    args = parser.parse_args()

    rows = analyse(load(args.metrics), args.vcfcache_cost, args.hit_rate)
    write_tsv(args.output, rows)

    print(
        f"VCFcache cost {args.vcfcache_cost:.1f}s at hit rate {args.hit_rate:.3f} "
        f"-> needs {args.vcfcache_cost / args.hit_rate:.1f}s of per-variant work\n"
    )
    print(
        f"{'arm':26} {'start-up':>9} {'direct':>8} {'variable':>9} "
        f"{'xcore':>6} {'pred cached':>11} {'pred':>7}"
    )
    print("-" * 82)
    for row in rows:
        print(
            f"{row['arm']:26} {row['startup_seconds']:8.2f}s {row['direct_wes_seconds']:7.2f}s "
            f"{row['variable_work_seconds']:8.2f}s {row['multiple_of_core']:6.2f} "
            f"{row['predicted_cached_seconds']:10.2f}s {row['predicted_speedup']:6.2f}x"
        )
    winners = [r for r in rows if r["predicted_win"]]
    if winners:
        first = winners[0]
        print(
            f"\nFirst arm predicted to return time: {first['arm']} "
            f"at {first['multiple_of_core']:.2f}x core per-variant work "
            f"({first['predicted_speedup']:.2f}x)"
        )
    else:
        print("\nNo arm reached break-even; more annotation work is required.")


if __name__ == "__main__":
    main()
