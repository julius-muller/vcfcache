#!/usr/bin/env python3
"""Fit and report the VCFcache break-even boundary across assays and annotators.

The matched assay campaign shows caching losing time for small inputs and for a
fast annotator on an exome, while winning decisively for whole genomes. This
module derives the boundary that separates those regimes so it can be stated as
a rule rather than a per-assay verdict.

Decomposing an annotator into a fixed start-up cost and per-variant work,

    T_direct = s_ann + N * t_ann
    T_cached = s_ann + (1 - f) * N * t_ann + C_vcfcache(N)

so the annotator start-up cancels and caching returns time exactly when

    f * N * t_ann > C_vcfcache(N)

Equivalently, caching wins when the observed hit rate exceeds the break-even hit
rate C_vcfcache divided by the recoverable per-variant work. The annotator's
start-up cost is never recoverable, so it does not enter that threshold.

One caveat is visible in the residuals. The model treats a cache miss as costing
the average per-variant time of the input, but misses are enriched for rare and
coding variants, which are more expensive to annotate for a transcript-heavy
annotator. That bias inflates the apparent VCFcache overhead at exome scale for
VEP, while fastVEP, whose cost is dominated by payload handling rather than
transcript evaluation, fits cleanly across three orders of magnitude of input
size.
"""

from __future__ import annotations

import argparse
import csv
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Any, Sequence

from benchmarks.run_strategy_comparison import write_tsv

ASSAY_ORDER = ("panel", "wes", "wgs")
ASSAY_LABELS = {"panel": "Panel", "wes": "WES", "wgs": "WGS"}
TOOL_LABELS = {"vep": "VEP", "fastvep": "fastVEP"}
# The smallest assay carries too few variants for per-variant work to matter, so
# its direct runtime estimates the annotator start-up cost.
STARTUP_ASSAY = "panel"


def read_rows(path: Path) -> list[dict[str, str]]:
    """Read semantically valid paired assay rows."""
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    selected = [row for row in rows if row["assay"] in ASSAY_ORDER]
    if not selected:
        raise RuntimeError("No panel, WES, or WGS rows were supplied")
    if any(row["semantic_pass"].lower() != "true" for row in selected):
        raise RuntimeError("Break-even input contains a failed semantic comparison")
    return selected


def median_cells(rows: Sequence[dict[str, str]]) -> dict[tuple[str, str], dict[str, float]]:
    """Collapse per-sample rows to one median observation per tool and assay."""
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[(row["tool"], row["assay"])].append(row)

    cells = {}
    for key, group in grouped.items():
        cells[key] = {
            "samples": float(len(group)),
            "input_records": statistics.median(float(r["input_records"]) for r in group),
            "cache_hit_rate": statistics.median(float(r["cache_hit_rate"]) for r in group),
            "uncached_wall_seconds": statistics.median(
                float(r["uncached_wall_seconds"]) for r in group
            ),
            "cached_wall_seconds": statistics.median(
                float(r["cached_wall_seconds"]) for r in group
            ),
            "speedup": statistics.median(float(r["speedup"]) for r in group),
        }
    return cells


def annotator_startup(cells: dict[tuple[str, str], dict[str, float]], tool: str) -> float:
    """Estimate an annotator's fixed start-up cost from its smallest assay."""
    cell = cells.get((tool, STARTUP_ASSAY))
    if cell is None:
        raise RuntimeError(f"No {STARTUP_ASSAY} cell for {tool}; cannot estimate start-up")
    return cell["uncached_wall_seconds"]


def build_table(cells: dict[tuple[str, str], dict[str, float]]) -> list[dict[str, Any]]:
    """Derive per-cell break-even quantities."""
    table = []
    for tool in sorted({key[0] for key in cells}):
        startup = annotator_startup(cells, tool)
        for assay in ASSAY_ORDER:
            cell = cells.get((tool, assay))
            if cell is None:
                continue
            direct = cell["uncached_wall_seconds"]
            cached = cell["cached_wall_seconds"]
            hit_rate = cell["cache_hit_rate"]
            records = cell["input_records"]

            # Work the annotator does per variant, net of its fixed start-up.
            variable_work = max(direct - startup, 0.0)
            # Caching can only ever return the fraction of that work it skips.
            recoverable = hit_rate * variable_work
            # Whatever cached runtime is left over is VCFcache's own cost.
            overhead = cached - startup - (1 - hit_rate) * variable_work
            # Caching returns time when f * variable_work exceeds that cost, so
            # the threshold is measured against the recoverable work alone. An
            # assay with no meaningful per-variant work can never break even.
            break_even = overhead / variable_work if variable_work > 0 else float("inf")

            table.append(
                {
                    "tool": tool,
                    "tool_label": TOOL_LABELS.get(tool, tool),
                    "assay": assay,
                    "assay_label": ASSAY_LABELS[assay],
                    "samples": int(cell["samples"]),
                    "input_records": round(records),
                    "cache_hit_rate": hit_rate,
                    "annotator_startup_seconds": startup,
                    "direct_seconds": direct,
                    "cached_seconds": cached,
                    "variable_work_seconds": variable_work,
                    "recoverable_seconds": recoverable,
                    "vcfcache_overhead_seconds": overhead,
                    "vcfcache_overhead_per_variant_us": (
                        overhead / records * 1e6 if records else float("nan")
                    ),
                    "break_even_hit_rate": break_even,
                    "predicted_win": recoverable > overhead,
                    "measured_speedup": cell["speedup"],
                    "measured_win": cell["speedup"] > 1.0,
                }
            )
    return table


def fit_overhead_model(rows: Sequence[dict[str, Any]]) -> tuple[float, float]:
    """Least-squares fit of C_vcfcache(N) = c0 + c1 * N over the given cells."""
    xs = [float(row["input_records"]) for row in rows]
    ys = [row["vcfcache_overhead_seconds"] for row in rows]
    n = len(xs)
    if n < 2:
        raise RuntimeError("Need at least two cells to fit the overhead model")
    mean_x = sum(xs) / n
    mean_y = sum(ys) / n
    variance = sum((x - mean_x) ** 2 for x in xs)
    if variance == 0:
        return mean_y, 0.0
    slope = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys, strict=True)) / variance
    return mean_y - slope * mean_x, slope


def leave_one_out_residuals(table: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    """Predict each cell from an overhead model fitted without that cell.

    The annotator start-up and per-variant work are properties of the annotator
    and input, measured independently of VCFcache. What the model must predict
    is VCFcache's own cost, so the held-out check refits C_vcfcache(N) on the
    five retained cells and asks whether the sixth cell's outcome follows.
    """
    residuals = []
    for held_out in table:
        others = [row for row in table if row is not held_out]
        intercept, slope = fit_overhead_model(others)
        predicted_overhead = intercept + slope * held_out["input_records"]
        predicted_cached = (
            held_out["annotator_startup_seconds"]
            + (1 - held_out["cache_hit_rate"]) * held_out["variable_work_seconds"]
            + predicted_overhead
        )
        predicted_speedup = (
            held_out["direct_seconds"] / predicted_cached
            if predicted_cached > 0
            else float("nan")
        )
        residuals.append(
            {
                "tool": held_out["tool"],
                "assay": held_out["assay"],
                "measured_overhead_seconds": held_out["vcfcache_overhead_seconds"],
                "predicted_overhead_seconds": predicted_overhead,
                "measured_speedup": held_out["measured_speedup"],
                "predicted_speedup": predicted_speedup,
                "measured_win": held_out["measured_win"],
                "predicted_win": predicted_speedup > 1.0,
                "outcome_agrees": (predicted_speedup > 1.0) == held_out["measured_win"],
            }
        )
    return residuals


def project_required_richness(
    table: Sequence[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Ask how much extra annotation work each losing cell would need to win.

    For a losing cell the recipe is too cheap: the work caching could skip is
    smaller than VCFcache's own cost. Holding the hit rate and VCFcache cost
    fixed, the recipe must do `C_vcfcache / f` seconds of per-variant work to
    break even, which is a multiple of what the measured recipe does. The
    multiple is a lower bound, because a richer recipe also enlarges the
    annotation payload that VCFcache copies, recompresses and indexes.
    """
    projections = []
    for row in table:
        if row["measured_win"]:
            continue
        required = row["vcfcache_overhead_seconds"] / row["cache_hit_rate"]
        current = row["variable_work_seconds"]
        projections.append(
            {
                "tool": row["tool"],
                "assay": row["assay"],
                "current_variable_work_seconds": current,
                "required_variable_work_seconds": required,
                "required_multiple_of_current": (
                    required / current if current > 0 else float("inf")
                ),
            }
        )
    return projections


def main() -> None:
    """Write the break-even source table and print a readable summary."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--assay-metrics",
        type=Path,
        required=True,
        help="Frozen matched assay metrics TSV",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Destination TSV for the derived break-even table",
    )
    args = parser.parse_args()

    cells = median_cells(read_rows(args.assay_metrics))
    table = build_table(cells)
    write_tsv(args.output, table)

    header = (
        f"{'cell':16} {'N':>9} {'f':>6} {'s_ann':>7} {'var work':>9} "
        f"{'recover':>9} {'overhead':>9} {'break-even f':>12} {'pred':>5} {'meas':>7}"
    )
    print(header)
    print("-" * len(header))
    for row in table:
        threshold = row["break_even_hit_rate"]
        threshold_text = (
            f"{threshold * 100:11.1f}% " if threshold != float("inf") else f"{'never':>11}  "
        )
        print(
            f"{row['tool'] + '/' + row['assay']:16} "
            f"{row['input_records']:9d} "
            f"{row['cache_hit_rate'] * 100:5.1f}% "
            f"{row['annotator_startup_seconds']:7.1f} "
            f"{row['variable_work_seconds']:9.1f} "
            f"{row['recoverable_seconds']:9.1f} "
            f"{row['vcfcache_overhead_seconds']:9.1f} "
            f"{threshold_text}"
            f"{'WIN' if row['predicted_win'] else 'LOSS':>5} "
            f"{row['measured_speedup']:6.2f}x"
        )

    print("\nFitted VCFcache cost per annotator, C(N) = c0 + c1 * N:")
    for tool in sorted({row["tool"] for row in table}):
        tool_rows = [row for row in table if row["tool"] == tool]
        intercept, slope = fit_overhead_model(tool_rows)
        print(f"  {TOOL_LABELS.get(tool, tool):8} {intercept:6.1f} s + {slope * 1e6:5.1f} us/variant")
        for row in tool_rows:
            fitted = intercept + slope * row["input_records"]
            print(
                f"    {row['assay']:6} measured {row['vcfcache_overhead_seconds']:7.1f} s  "
                f"fitted {fitted:7.1f} s  residual {fitted - row['vcfcache_overhead_seconds']:+7.1f} s"
            )

    print(
        "\nExtra annotation work each losing cell would need to break even "
        "(lower bound, holding VCFcache cost fixed):"
    )
    for row in project_required_richness(table):
        multiple = row["required_multiple_of_current"]
        multiple_text = f"{multiple:.1f}x" if multiple != float("inf") else "unreachable"
        print(
            f"  {row['tool'] + '/' + row['assay']:16} "
            f"current {row['current_variable_work_seconds']:7.1f} s  "
            f"required {row['required_variable_work_seconds']:7.1f} s  "
            f"= {multiple_text} the measured recipe"
        )

    residuals = leave_one_out_residuals(table)
    agree = sum(1 for r in residuals if r["outcome_agrees"])
    print(f"Leave-one-out outcome agreement: {agree}/{len(residuals)} cells")
    for row in residuals:
        status = "ok" if row["outcome_agrees"] else "MISMATCH"
        print(
            f"  {row['tool'] + '/' + row['assay']:16} "
            f"overhead measured {row['measured_overhead_seconds']:7.1f} s  "
            f"predicted {row['predicted_overhead_seconds']:7.1f} s  "
            f"speedup {row['measured_speedup']:6.2f}x vs {row['predicted_speedup']:6.2f}x  {status}"
        )


if __name__ == "__main__":
    main()
