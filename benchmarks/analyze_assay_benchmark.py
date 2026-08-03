#!/usr/bin/env python3
"""Summarize paired panel, WES, and WGS runs and render manuscript evidence."""

from __future__ import annotations

import argparse
import csv
import random
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Any, Callable, Sequence

from benchmarks.run_strategy_comparison import write_tsv

ASSAY_ORDER = ("panel", "wes", "wgs")
ASSAY_LABELS = {"panel": "Panel", "wes": "WES", "wgs": "WGS"}
COLORS = {"panel": "#59A14F", "wes": "#F28E2B", "wgs": "#4C78A8"}
BOOTSTRAP_SEED = 20260803
BOOTSTRAP_REPLICATES = 10_000


def read_rows(paths: Sequence[Path]) -> list[dict[str, str]]:
    """Read semantically valid paired rows from one or more collectors."""
    rows = []
    for path in paths:
        with path.open(newline="") as handle:
            rows.extend(csv.DictReader(handle, delimiter="\t"))
    selected = [row for row in rows if row["assay"] in ASSAY_ORDER]
    if not selected:
        raise RuntimeError("No panel, WES, or WGS rows were supplied")
    if any(row["semantic_pass"].lower() != "true" for row in selected):
        raise RuntimeError("Assay input contains a failed semantic comparison")
    return selected


def percentile(values: Sequence[float], fraction: float) -> float:
    """Return a linearly interpolated percentile."""
    ordered = sorted(values)
    if not ordered:
        raise ValueError("Cannot calculate a percentile of no values")
    position = (len(ordered) - 1) * fraction
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    weight = position - lower
    return ordered[lower] * (1 - weight) + ordered[upper] * weight


def bootstrap_interval(
    values: Sequence[float],
    statistic: Callable[[Sequence[float]], float] = statistics.median,
    *,
    seed: int = BOOTSTRAP_SEED,
    replicates: int = BOOTSTRAP_REPLICATES,
) -> tuple[float, float]:
    """Return a deterministic sample-level percentile bootstrap interval."""
    generator = random.Random(seed)
    estimates = [
        statistic(generator.choices(values, k=len(values))) for _ in range(replicates)
    ]
    return percentile(estimates, 0.025), percentile(estimates, 0.975)


def summarize(rows: Sequence[dict[str, str]]) -> list[dict[str, Any]]:
    """Report descriptive sample-level statistics by assay."""
    groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        groups[row["assay"]].append(row)
    output = []
    for index, assay in enumerate(ASSAY_ORDER):
        group = groups.get(assay, [])
        if not group:
            continue
        relative = [float(row["relative_runtime"]) for row in group]
        speedup = [float(row["speedup"]) for row in group]
        hit_rate = [float(row["cache_hit_rate"]) for row in group]
        low, high = bootstrap_interval(relative, seed=BOOTSTRAP_SEED + index)
        output.append(
            {
                "assay": assay,
                "assay_label": ASSAY_LABELS[assay],
                "samples": len({row["sample"] for row in group}),
                "observations": len(group),
                "median_input_records": statistics.median(
                    int(row["input_records"]) for row in group
                ),
                "median_hit_rate": statistics.median(hit_rate),
                "hit_rate_q1": percentile(hit_rate, 0.25),
                "hit_rate_q3": percentile(hit_rate, 0.75),
                "median_relative_runtime": statistics.median(relative),
                "relative_runtime_q1": percentile(relative, 0.25),
                "relative_runtime_q3": percentile(relative, 0.75),
                "relative_runtime_bootstrap_low": low,
                "relative_runtime_bootstrap_high": high,
                "median_speedup": statistics.median(speedup),
                "total_wall_hours_saved": sum(
                    float(row["wall_seconds_saved"]) for row in group
                )
                / 3600,
                "semantic_passes": sum(
                    row["semantic_pass"].lower() == "true" for row in group
                ),
                "unexpected_annotation_mismatches": sum(
                    int(row["annotation_mismatches"]) for row in group
                ),
                "ignored_known_vep_mismatches": sum(
                    int(row["ignored_annotation_mismatches"]) for row in group
                ),
                "annotation_order_only": sum(
                    int(row["annotation_order_only"]) for row in group
                ),
            }
        )
    return output


def _scale(value: float, low: float, high: float, start: float, end: float) -> float:
    if high == low:
        return (start + end) / 2
    return start + (value - low) * (end - start) / (high - low)


def render_svg(
    path: Path, rows: Sequence[dict[str, str]], summaries: Sequence[dict[str, Any]]
) -> None:
    """Render a three-panel assay-size, mechanism, and correctness figure."""
    width, height = 1500, 590
    panels = ((70, 470), (550, 1010), (1090, 1430))
    top, bottom = 80, 440
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        "<style>text{font-family:Arial,sans-serif;fill:#222}.t{font-size:19px;font-weight:700}.a{font-size:13px}.big{font-size:34px;font-weight:700}.small{font-size:11px}</style>",
    ]
    for left, right in panels[:2]:
        svg.extend(
            [
                f'<line x1="{left}" y1="{bottom}" x2="{right}" y2="{bottom}" stroke="#333"/>',
                f'<line x1="{left}" y1="{top}" x2="{left}" y2="{bottom}" stroke="#333"/>',
            ]
        )

    left, right = panels[0]
    svg.append(f'<text class="t" x="{left}" y="35">A  Time remaining by assay</text>')
    for assay_index, assay in enumerate(ASSAY_ORDER):
        group = [row for row in rows if row["assay"] == assay]
        x = _scale(assay_index, 0, 2, left + 65, right - 65)
        for index, row in enumerate(group):
            jitter = ((index * 37) % 31 - 15) * 0.8
            y = _scale(float(row["relative_runtime"]), 0, 1, bottom, top)
            svg.append(
                f'<circle cx="{x + jitter:.1f}" cy="{y:.1f}" r="3" fill="{COLORS[assay]}" opacity="0.45"/>'
            )
        summary = next(value for value in summaries if value["assay"] == assay)
        median_y = _scale(float(summary["median_relative_runtime"]), 0, 1, bottom, top)
        svg.append(
            f'<line x1="{x - 35:.1f}" y1="{median_y:.1f}" x2="{x + 35:.1f}" y2="{median_y:.1f}" stroke="#111" stroke-width="4"/>'
        )
        svg.append(
            f'<text class="a" x="{x:.1f}" y="466" text-anchor="middle">{ASSAY_LABELS[assay]}</text>'
        )
    svg.append(
        f'<text class="a" transform="translate({left - 43} 300) rotate(-90)">Cached / uncached wall time</text>'
    )

    left, right = panels[1]
    svg.append(
        f'<text class="t" x="{left}" y="35">B  Individual impact follows hit rate</text>'
    )
    max_speed = max(2.0, max(float(row["speedup"]) for row in rows) * 1.08)
    for row in rows:
        x = _scale(float(row["cache_hit_rate"]), 0, 1, left, right)
        y = _scale(float(row["speedup"]), 1, max_speed, bottom, top)
        svg.append(
            f'<circle cx="{x:.1f}" cy="{y:.1f}" r="3.5" fill="{COLORS[row["assay"]]}" opacity="0.55"/>'
        )
    svg.extend(
        [
            f'<text class="a" x="{(left + right) / 2}" y="470" text-anchor="middle">Observed cache hit rate</text>',
            f'<text class="a" transform="translate({left - 43} 310) rotate(-90)">Measured speedup</text>',
        ]
    )

    left, right = panels[2]
    total = sum(int(value["observations"]) for value in summaries)
    unexpected = sum(
        int(value["unexpected_annotation_mismatches"]) for value in summaries
    )
    ignored = sum(int(value["ignored_known_vep_mismatches"]) for value in summaries)
    order_only = sum(int(value["annotation_order_only"]) for value in summaries)
    svg.extend(
        [
            f'<text class="t" x="{left}" y="35">C  Output correctness</text>',
            f'<rect x="{left}" y="90" width="{right-left}" height="310" rx="22" fill="#EAF6EE"/>',
            f'<text class="big" x="{(left + right) / 2}" y="170" text-anchor="middle" fill="#2E8B57">{total}/{total} PASS</text>',
            f'<text class="a" x="{left + 28}" y="225">Unexpected annotation mismatches: {unexpected}</text>',
            f'<text class="a" x="{left + 28}" y="265">Known VEP HGNC_ID differences: {ignored}</text>',
            f'<text class="a" x="{left + 28}" y="305">CSQ member-order differences: {order_only}</text>',
            f'<text class="small" x="{left + 28}" y="350">Known upstream differences are counted and reported;</text>',
            f'<text class="small" x="{left + 28}" y="370">all other key, header, or CSQ changes fail a task.</text>',
        ]
    )
    svg.append(
        '<text class="small" x="70" y="535">Bars show medians; points are biological samples. Timing includes lookup, VEP, output assembly, compression, indexing, and accounting.</text>'
    )
    svg.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(svg) + "\n")


def run(args: argparse.Namespace) -> None:
    """Write assay summaries and the manuscript SVG."""
    rows = read_rows(args.metrics)
    summaries = summarize(rows)
    missing = set(ASSAY_ORDER) - {row["assay"] for row in summaries}
    if missing:
        raise RuntimeError(f"Missing publication assays: {sorted(missing)}")
    args.output.mkdir(parents=True, exist_ok=True)
    write_tsv(args.output / "assay_summary.tsv", summaries)
    render_svg(args.output / "assay_benchmark.svg", rows, summaries)
    print(f"Wrote assay publication outputs under {args.output}")


def parser() -> argparse.ArgumentParser:
    """Build the assay analysis CLI."""
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--metrics", type=Path, nargs="+", required=True)
    result.add_argument("--output", type=Path, required=True)
    return result


def main() -> None:
    """Run assay analysis."""
    run(parser().parse_args())


if __name__ == "__main__":
    main()
