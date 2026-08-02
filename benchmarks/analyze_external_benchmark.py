#!/usr/bin/env python3
"""Aggregate external-WGS measurements and render the scaling-model figure."""

from __future__ import annotations

import argparse
import csv
import json
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Any, Sequence

from benchmarks.run_strategy_comparison import write_tsv, xml_escape

STRATEGY_ORDER = ("gnomad_af_0.1", "gnomad_af_0.01", "cohort_3_genomes")
STRATEGY_LABELS = {
    "gnomad_af_0.1": "Bundled gnomAD AF >= 10%",
    "gnomad_af_0.01": "Bundled gnomAD AF >= 1%",
    "cohort_3_genomes": "Three-genome cohort cache",
}
COLORS = {"kpgp": "#4C78A8", "sgdp": "#59A14F", "pgp": "#E15759"}


def read_rows(path: Path) -> list[dict[str, str]]:
    """Read and validate measured publication-source rows."""
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise RuntimeError(f"No external benchmark rows in {path}")
    if any(row["semantic_pass"].lower() != "true" for row in rows):
        raise RuntimeError("Figure input contains a failed semantic comparison")
    return rows


def strategy_summary(rows: Sequence[dict[str, str]]) -> list[dict[str, Any]]:
    """Summarize paired sample replicates without pooling cohorts."""
    groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        groups[(row["cohort"], row["strategy"])].append(row)
    values = []
    for (cohort, strategy), group in sorted(groups.items()):
        hit_rates = [float(row["cache_hit_rate"]) for row in group]
        cached = [float(row["cached_wall_seconds"]) for row in group]
        uncached = [float(row["uncached_wall_seconds"]) for row in group]
        speedups = [float(row["speedup"]) for row in group]
        values.append(
            {
                "cohort": cohort,
                "strategy": strategy,
                "samples": len({row["sample"] for row in group}),
                "replicates": len(group),
                "median_hit_rate": statistics.median(hit_rates),
                "median_cached_seconds": statistics.median(cached),
                "median_uncached_seconds": statistics.median(uncached),
                "median_speedup": statistics.median(speedups),
                "median_relative_runtime": statistics.median(
                    float(row["relative_runtime"]) for row in group
                ),
            }
        )
    return values


def build_costs(strategy_manifest: Path) -> dict[str, float]:
    """Read per-cohort custom-cache build costs from the frozen strategy manifest."""
    value = json.loads(strategy_manifest.read_text())
    return {
        cohort: float(strategy["build_wall_seconds"])
        for cohort, strategy in value["cohort_strategies"].items()
    }


def scaling_rows(
    summaries: Sequence[dict[str, Any]], costs: dict[str, float]
) -> list[dict[str, Any]]:
    """Evaluate T_eff = T_cached + T_build/S for cohort sizes 1 through 1000."""
    output = []
    for row in summaries:
        build = costs[row["cohort"]] if row["strategy"] == "cohort_3_genomes" else 0.0
        for cohort_size in range(1, 1001):
            effective = float(row["median_cached_seconds"]) + build / cohort_size
            output.append(
                {
                    "cohort": row["cohort"],
                    "strategy": row["strategy"],
                    "cohort_size": cohort_size,
                    "cache_build_seconds": build,
                    "effective_seconds_per_sample": effective,
                    "effective_speedup": float(row["median_uncached_seconds"])
                    / effective,
                }
            )
    return output


def _scale(value: float, low: float, high: float, start: float, end: float) -> float:
    if high == low:
        return (start + end) / 2
    return start + (value - low) * (end - start) / (high - low)


def make_svg(
    path: Path,
    rows: Sequence[dict[str, str]],
    summaries: Sequence[dict[str, Any]],
    modeled: Sequence[dict[str, Any]],
) -> None:
    """Render a dependency-free three-panel publication SVG."""
    width, height = 1500, 560
    panels = ((70, 450), (550, 930), (1030, 1450))
    top, bottom = 65, 410
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        "<style>text{font-family:Arial,sans-serif;fill:#222}.t{font-size:18px;font-weight:700}.a{font-size:13px}.l{font-size:12px}</style>",
    ]
    for left, right in panels:
        svg.extend(
            [
                f'<line x1="{left}" y1="{bottom}" x2="{right}" y2="{bottom}" stroke="#222"/>',
                f'<line x1="{left}" y1="{top}" x2="{left}" y2="{bottom}" stroke="#222"/>',
            ]
        )

    left, right = panels[0]
    svg.append(f'<text class="t" x="{left}" y="30">A  Relative runtime</text>')
    for index, strategy in enumerate(STRATEGY_ORDER):
        x = _scale(index, 0, 2, left + 70, right - 70)
        group = [row for row in summaries if row["strategy"] == strategy]
        for offset, value in enumerate(group):
            y = _scale(float(value["median_relative_runtime"]), 0, 1, bottom, top)
            svg.append(
                f'<circle cx="{x + (offset - 1) * 12:.1f}" cy="{y:.1f}" r="6" fill="{COLORS[value["cohort"]]}"/>'
            )
        label = STRATEGY_LABELS[strategy]
        svg.append(
            f'<text class="l" text-anchor="middle" x="{x:.1f}" y="440">{xml_escape(label)}</text>'
        )

    left, right = panels[1]
    svg.append(f'<text class="t" x="{left}" y="30">B  Speedup versus hit rate</text>')
    max_speed = max(float(row["speedup"]) for row in rows) * 1.05
    for row in rows:
        x = _scale(float(row["cache_hit_rate"]), 0, 1, left, right)
        y = _scale(float(row["speedup"]), 1, max_speed, bottom, top)
        svg.append(
            f'<circle cx="{x:.1f}" cy="{y:.1f}" r="2.4" fill="{COLORS[row["cohort"]]}" opacity="0.45"/>'
        )
    svg.append(
        f'<text class="a" x="{(left + right) / 2}" y="455" text-anchor="middle">Cache hit rate</text>'
    )

    left, right = panels[2]
    svg.append(f'<text class="t" x="{left}" y="30">C  Build-cost amortization</text>')
    custom = [row for row in modeled if row["strategy"] == "cohort_3_genomes"]
    max_effective = max(float(row["effective_speedup"]) for row in custom) * 1.05
    for cohort in sorted({row["cohort"] for row in custom}):
        points = [row for row in custom if row["cohort"] == cohort]
        coordinates = " ".join(
            f'{_scale(float(row["cohort_size"]), 1, 1000, left, right):.1f},'
            f'{_scale(float(row["effective_speedup"]), 0, max_effective, bottom, top):.1f}'
            for row in points
        )
        svg.append(
            f'<polyline points="{coordinates}" fill="none" stroke="{COLORS[cohort]}" stroke-width="2"/>'
        )
    for marker in (5, 10, 100):
        x = _scale(marker, 1, 1000, left, right)
        svg.append(
            f'<line x1="{x:.1f}" y1="{top}" x2="{x:.1f}" y2="{bottom}" stroke="#aaa" stroke-dasharray="3,3"/>'
        )
        svg.append(
            f'<text class="l" x="{x:.1f}" y="430" text-anchor="middle">S={marker}</text>'
        )
    svg.append(
        f'<text class="a" x="{(left + right) / 2}" y="455" text-anchor="middle">Cohort size S</text>'
    )

    legend_x = 585
    for index, cohort in enumerate(sorted(COLORS)):
        svg.append(
            f'<circle cx="{legend_x + index * 95}" cy="510" r="6" fill="{COLORS[cohort]}"/>'
            f'<text class="a" x="{legend_x + 10 + index * 95}" y="515">{cohort.upper()}</text>'
        )
    svg.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(svg) + "\n")


def run(args: argparse.Namespace) -> None:
    """Write aggregate tables, model grid, and SVG figure."""
    rows = read_rows(args.metrics)
    summaries = strategy_summary(rows)
    modeled = scaling_rows(summaries, build_costs(args.strategies))
    args.output.mkdir(parents=True, exist_ok=True)
    write_tsv(args.output / "external_wgs_summary.tsv", summaries)
    write_tsv(args.output / "external_wgs_scaling_model.tsv", modeled)
    make_svg(args.output / "external_wgs_scaling.svg", rows, summaries, modeled)
    print(f"Wrote external-WGS publication outputs under {args.output}")


def parser() -> argparse.ArgumentParser:
    """Build the analysis CLI."""
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--metrics", type=Path, required=True)
    result.add_argument("--strategies", type=Path, required=True)
    result.add_argument("--output", type=Path, required=True)
    return result


def main() -> None:
    """Run publication analysis."""
    args = parser().parse_args()
    run(args)


if __name__ == "__main__":
    main()
