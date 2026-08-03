#!/usr/bin/env python3
"""Render a simple user-facing VCFcache impact calculator as SVG and TSV."""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


@dataclass(frozen=True)
class Scenario:
    """One hit-rate situation shown in the time-left bars."""

    label: str
    hit_rate: float


DEFAULT_SCENARIOS = (
    Scenario("Distant WGS", 0.50),
    Scenario("Related WGS", 0.80),
    Scenario("WES", 0.90),
    Scenario("Very high reuse", 0.95),
)


def cached_seconds(
    baseline_seconds: float, hit_rate: float, overhead_seconds: float
) -> float:
    """Return modeled cached wall time, capped at the uncached baseline."""
    if baseline_seconds <= 0:
        raise ValueError("Baseline runtime must be positive")
    if not 0 <= hit_rate <= 1:
        raise ValueError("Hit rate must be between zero and one")
    if overhead_seconds < 0:
        raise ValueError("Lookup overhead cannot be negative")
    return min(baseline_seconds, overhead_seconds + (1 - hit_rate) * baseline_seconds)


def impact_rows(
    scenarios: Sequence[Scenario],
    pipeline_minutes: Sequence[float],
    sample_counts: Sequence[int],
    overhead_seconds: float,
    focus_hit_rate: float,
    build_minutes: float,
) -> list[dict[str, str | int | float]]:
    """Create the machine-readable source values used by the SVG."""
    rows: list[dict[str, str | int | float]] = []

    def seconds(value: float) -> float:
        return round(value, 6)

    normalized = 100.0
    for scenario in scenarios:
        after = cached_seconds(normalized, scenario.hit_rate, overhead_seconds=0)
        rows.append(
            {
                "panel": "time_left",
                "label": scenario.label,
                "hit_rate": scenario.hit_rate,
                "baseline_seconds": normalized,
                "cached_seconds": seconds(after),
                "saved_seconds": seconds(normalized - after),
                "sample_count": 1,
                "build_seconds": 0,
            }
        )
    for minutes in pipeline_minutes:
        baseline = minutes * 60
        after = cached_seconds(baseline, focus_hit_rate, overhead_seconds)
        rows.append(
            {
                "panel": "pipeline_cost",
                "label": f"{minutes:g} minute pipeline",
                "hit_rate": focus_hit_rate,
                "baseline_seconds": baseline,
                "cached_seconds": seconds(after),
                "saved_seconds": seconds(baseline - after),
                "sample_count": 1,
                "build_seconds": 0,
            }
        )
    scale_baseline = pipeline_minutes[len(pipeline_minutes) // 2] * 60
    scale_cached = cached_seconds(scale_baseline, focus_hit_rate, overhead_seconds)
    build_seconds = build_minutes * 60
    for count in sample_counts:
        before = count * scale_baseline
        after = build_seconds + count * scale_cached
        rows.append(
            {
                "panel": "sample_scale",
                "label": f"{count:,} samples",
                "hit_rate": focus_hit_rate,
                "baseline_seconds": before,
                "cached_seconds": seconds(after),
                "saved_seconds": seconds(before - after),
                "sample_count": count,
                "build_seconds": build_seconds,
            }
        )
    return rows


def human_time(seconds: float) -> str:
    """Format seconds compactly for a nontechnical reader."""
    if seconds < 60:
        return f"{seconds:.0f} sec"
    if seconds < 3600:
        minutes = seconds / 60
        if math.isclose(minutes, round(minutes), abs_tol=0.05):
            return f"{minutes:.0f} min"
        return f"{minutes:.1f} min"
    if seconds < 48 * 3600:
        return f"{seconds / 3600:.1f} h"
    return f"{seconds / 86400:.1f} days"


def _write_tsv(path: Path, rows: Sequence[dict[str, str | int | float]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _escape(value: object) -> str:
    return str(value).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def render_svg(
    path: Path,
    rows: Sequence[dict[str, str | int | float]],
    focus_hit_rate: float,
    build_minutes: float,
) -> None:
    """Render a dependency-free, presentation-ready impact graphic."""
    width, height = 1500, 940
    dark, saved, ink, muted, trust = (
        "#28536B",
        "#A7D8B7",
        "#18323F",
        "#60757F",
        "#2E8B57",
    )
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#FCFDFB"/>',
        "<style>text{font-family:Arial,sans-serif}.title{font-size:34px;font-weight:700}.h{font-size:21px;font-weight:700}.b{font-size:17px}.small{font-size:14px}</style>",
        f'<text class="title" x="70" y="58" fill="{ink}">How much annotation time do you get back?</text>',
        f'<text class="b" x="70" y="92" fill="{muted}">New time ≈ lookup + (1 − hit rate) × your current time</text>',
    ]

    time_rows = [row for row in rows if row["panel"] == "time_left"]
    bar_x, bar_width = 300, 690
    svg.append(
        f'<text class="h" x="70" y="145" fill="{ink}">Find your expected cache reuse</text>'
    )
    for index, row in enumerate(time_rows):
        y = 190 + index * 82
        fraction = float(row["cached_seconds"]) / float(row["baseline_seconds"])
        remaining = bar_width * fraction
        svg.extend(
            [
                f'<text class="b" x="70" y="{y + 25}" fill="{ink}">{_escape(row["label"])}</text>',
                f'<rect x="{bar_x}" y="{y}" width="{bar_width}" height="38" rx="12" fill="{saved}"/>',
                f'<rect x="{bar_x}" y="{y}" width="{remaining:.1f}" height="38" rx="12" fill="{dark}"/>',
                f'<text class="small" x="{bar_x + remaining + 18:.1f}" y="{y + 25}" fill="{trust}">← time returned</text>',
                f'<text class="b" x="{bar_x + bar_width + 24}" y="{y + 25}" fill="{ink}">{float(row["hit_rate"]):.0%} hit  ·  {fraction:.0%} time left  ·  {1 / fraction:.1f}× faster</text>',
            ]
        )
    svg.extend(
        [
            '<rect x="70" y="535" width="860" height="180" rx="22" fill="#EDF5F8"/>',
            f'<text class="h" x="100" y="575" fill="{ink}">Fast or slow pipeline? The same rule applies.</text>',
        ]
    )
    cost_rows = [row for row in rows if row["panel"] == "pipeline_cost"]
    for index, row in enumerate(cost_rows):
        x = 100 + index * 270
        svg.extend(
            [
                f'<text class="small" x="{x}" y="615" fill="{muted}">Today: {human_time(float(row["baseline_seconds"]))}</text>',
                f'<text class="h" x="{x}" y="650" fill="{dark}">With cache: {human_time(float(row["cached_seconds"]))}</text>',
                f'<text class="b" x="{x}" y="682" fill="{trust}">Get back {human_time(float(row["saved_seconds"]))}</text>',
            ]
        )

    svg.extend(
        [
            '<rect x="970" y="535" width="460" height="180" rx="22" fill="#EAF6EE"/>',
            f'<text class="h" x="1000" y="578" fill="{trust}">✓ Same annotations</text>',
            f'<text class="b" x="1000" y="618" fill="{ink}">Cached and uncached outputs are compared.</text>',
            f'<text class="b" x="1000" y="650" fill="{ink}">Any unexpected annotation change fails.</text>',
            f'<text class="small" x="1000" y="689" fill="{muted}">Known VEP HGNC_ID nondeterminism is reported separately.</text>',
            f'<text class="h" x="70" y="770" fill="{ink}">From 2 samples to 1,000: savings add up</text>',
        ]
    )
    scale_rows = [row for row in rows if row["panel"] == "sample_scale"]
    for index, row in enumerate(scale_rows):
        x = 70 + index * 350
        saved_seconds = float(row["saved_seconds"])
        saving_label = (
            f"Get back {human_time(saved_seconds)}"
            if saved_seconds >= 0
            else f"Build cost still exceeds savings by {human_time(-saved_seconds)}"
        )
        saving_color = trust if saved_seconds >= 0 else muted
        svg.extend(
            [
                f'<rect x="{x}" y="800" width="310" height="105" rx="18" fill="white" stroke="#D8E3E7" stroke-width="2"/>',
                f'<text class="h" x="{x + 22}" y="835" fill="{ink}">{_escape(row["label"])}</text>',
                f'<text class="b" x="{x + 22}" y="868" fill="{saving_color}">{saving_label}</text>',
                f'<text class="small" x="{x + 22}" y="894" fill="{muted}">{focus_hit_rate:.0%} hit · {"bundled cache" if build_minutes == 0 else "build cost included"}</text>',
            ]
        )
    svg.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(svg) + "\n")


def parser() -> argparse.ArgumentParser:
    """Build the command-line parser."""
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--output", type=Path, required=True)
    result.add_argument("--lookup-overhead-seconds", type=float, default=0)
    result.add_argument(
        "--model-json",
        type=Path,
        help="Use measured overhead from analyze_controlled_runtime.py",
    )
    result.add_argument("--model-pipeline", default="everything")
    result.add_argument("--focus-hit-rate", type=float, default=0.80)
    result.add_argument(
        "--pipeline-minutes", type=float, nargs=3, default=(10, 60, 600)
    )
    result.add_argument(
        "--sample-counts", type=int, nargs=4, default=(2, 10, 100, 1000)
    )
    result.add_argument("--build-minutes", type=float, default=0)
    return result


def main() -> None:
    """Render the user-facing figure and its machine-readable source table."""
    args = parser().parse_args()
    if any(value <= 0 for value in args.pipeline_minutes):
        raise ValueError("Pipeline runtimes must be positive")
    if any(value <= 0 for value in args.sample_counts):
        raise ValueError("Sample counts must be positive")
    overhead = args.lookup_overhead_seconds
    if args.model_json:
        model = json.loads(args.model_json.read_text())
        try:
            measured = model["pipelines"][args.model_pipeline]
        except KeyError as error:
            raise ValueError(
                f"Pipeline {args.model_pipeline!r} is absent from {args.model_json}"
            ) from error
        overhead = float(measured["lookup_preprocessing_overhead_seconds"])
    rows = impact_rows(
        DEFAULT_SCENARIOS,
        args.pipeline_minutes,
        args.sample_counts,
        overhead,
        args.focus_hit_rate,
        args.build_minutes,
    )
    args.output.mkdir(parents=True, exist_ok=True)
    _write_tsv(args.output / "user_impact_source.tsv", rows)
    render_svg(
        args.output / "user_impact.svg", rows, args.focus_hit_rate, args.build_minutes
    )
    saved_per_sample = float(
        next(row for row in rows if row["panel"] == "pipeline_cost")["saved_seconds"]
    )
    if saved_per_sample > 0 and args.build_minutes > 0:
        break_even = math.ceil(args.build_minutes * 60 / saved_per_sample)
        print(f"Custom-cache break-even: {break_even} samples")
    print(f"Wrote user-facing impact figure under {args.output}")


if __name__ == "__main__":
    main()
