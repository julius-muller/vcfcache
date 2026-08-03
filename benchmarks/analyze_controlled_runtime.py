#!/usr/bin/env python3
"""Fit the runtime scaling rule and render the controlled manuscript panel."""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence

from benchmarks.run_pilot import write_json_atomic
from benchmarks.run_strategy_comparison import write_tsv, xml_escape

PIPELINE_ORDER = ("vanilla", "everything", "delay_5ms", "delay_20ms")
PIPELINE_LABELS = {
    "vanilla": "Vanilla VEP",
    "everything": "VEP --everything",
    "delay_5ms": "Vanilla + 5 ms/consequence",
    "delay_20ms": "Vanilla + 20 ms/consequence",
}
PIPELINE_COLORS = {
    "vanilla": "#4C78A8",
    "everything": "#F28E2B",
    "delay_5ms": "#59A14F",
    "delay_20ms": "#B279A2",
}


def read_rows(path: Path) -> list[dict[str, str]]:
    """Read a complete, semantically validated controlled matrix."""
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    expected = len(PIPELINE_ORDER) * 5
    if len(rows) != expected:
        raise RuntimeError(f"Expected {expected} controlled rows, found {len(rows)}")
    if any(row["semantic_pass"].lower() != "true" for row in rows):
        raise RuntimeError("Controlled input contains a failed semantic comparison")
    if {row["pipeline"] for row in rows} != set(PIPELINE_ORDER):
        raise RuntimeError("Controlled input does not contain the frozen pipelines")
    return rows


def fit_models(
    rows: Sequence[dict[str, str]],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Fit the prespecified unit-slope model using a robust overhead intercept."""
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["pipeline"]].append(row)
    modeled: list[dict[str, Any]] = []
    models: dict[str, Any] = {}
    for pipeline in PIPELINE_ORDER:
        group = sorted(
            grouped[pipeline], key=lambda row: float(row["observed_hit_rate"])
        )
        uncached_values = {
            round(float(row["uncached_wall_seconds"]), 6) for row in group
        }
        if len(uncached_values) != 1:
            raise RuntimeError(
                f"Pipeline {pipeline} has more than one baseline runtime"
            )
        uncached = float(group[0]["uncached_wall_seconds"])
        residuals = [
            float(row["cached_wall_seconds"])
            - (1 - float(row["observed_hit_rate"])) * uncached
            for row in group
        ]
        overhead = max(0.0, statistics.median(residuals))
        errors = []
        for row in group:
            hit_rate = float(row["observed_hit_rate"])
            cached = float(row["cached_wall_seconds"])
            predicted = overhead + (1 - hit_rate) * uncached
            errors.append(cached - predicted)
            modeled.append(
                {
                    "pipeline": pipeline,
                    "pipeline_label": PIPELINE_LABELS[pipeline],
                    "delay_us": int(row["delay_us"]),
                    "target_hit_rate": float(row["target_hit_rate"]),
                    "observed_hit_rate": hit_rate,
                    "uncached_wall_seconds": uncached,
                    "cached_wall_seconds": cached,
                    "measured_relative_runtime": cached / uncached,
                    "predicted_cached_seconds": predicted,
                    "predicted_relative_runtime": predicted / uncached,
                    "residual_seconds": cached - predicted,
                    "speedup": uncached / cached,
                    "semantic_pass": True,
                }
            )
        rmse = math.sqrt(statistics.mean(error * error for error in errors))
        models[pipeline] = {
            "label": PIPELINE_LABELS[pipeline],
            "uncached_wall_seconds": uncached,
            "lookup_preprocessing_overhead_seconds": overhead,
            "lookup_preprocessing_overhead_fraction": overhead / uncached,
            "rmse_seconds": rmse,
            "observations": len(group),
            "formula": "T_cached = overhead + (1 - hit_rate) * T_uncached",
        }
    document = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "model": "T_cached = T_overhead + (1 - f) * T_uncached",
        "fit": "unit slope fixed by the prespecified model; overhead is the nonnegative median residual",
        "pipelines": models,
        "recommended_user_figure_pipeline": "everything",
    }
    return modeled, document


def _scale(value: float, low: float, high: float, start: float, end: float) -> float:
    return start + (value - low) * (end - start) / (high - low)


def render_svg(
    path: Path, rows: Sequence[dict[str, Any]], model: dict[str, Any]
) -> None:
    """Render four small multiples with observed points and model curves."""
    width, height = 1300, 780
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        "<style>text{font-family:Arial,sans-serif;fill:#222}.title{font-size:25px;font-weight:700}.h{font-size:17px;font-weight:700}.a{font-size:13px}.small{font-size:11px}</style>",
        '<text class="title" x="65" y="42">Measured runtime follows cache misses across pipeline costs</text>',
        '<text class="a" x="65" y="70">Points: observed · line: T cached = measured overhead + (1 − hit rate) × T uncached · one real WES input, no technical repeats</text>',
    ]
    panel_width, panel_height = 535, 255
    positions = ((65, 105), (690, 105), (65, 430), (690, 430))
    for pipeline, (left, top) in zip(PIPELINE_ORDER, positions, strict=True):
        right, bottom = left + panel_width, top + panel_height
        plot_left, plot_right = left + 55, right - 20
        plot_top, plot_bottom = top + 38, bottom - 42
        group = sorted(
            (row for row in rows if row["pipeline"] == pipeline),
            key=lambda row: float(row["observed_hit_rate"]),
        )
        maximum = max(
            0.6, max(float(row["measured_relative_runtime"]) for row in group) * 1.12
        )
        svg.extend(
            [
                f'<rect x="{left}" y="{top}" width="{panel_width}" height="{panel_height}" rx="12" fill="#FAFBFC" stroke="#D7DEE2"/>',
                f'<text class="h" x="{left + 18}" y="{top + 27}">{xml_escape(PIPELINE_LABELS[pipeline])}</text>',
                f'<line x1="{plot_left}" y1="{plot_bottom}" x2="{plot_right}" y2="{plot_bottom}" stroke="#333"/>',
                f'<line x1="{plot_left}" y1="{plot_top}" x2="{plot_left}" y2="{plot_bottom}" stroke="#333"/>',
            ]
        )
        for tick in (0.5, 0.8, 0.9, 1.0):
            x = _scale(tick, 0.5, 1.0, plot_left, plot_right)
            svg.append(
                f'<text class="small" x="{x:.1f}" y="{plot_bottom + 18}" text-anchor="middle">{tick:.0%}</text>'
            )
        coordinates = " ".join(
            f'{_scale(float(row["observed_hit_rate"]), 0.5, 1.0, plot_left, plot_right):.1f},'
            f'{_scale(float(row["predicted_relative_runtime"]), 0, maximum, plot_bottom, plot_top):.1f}'
            for row in group
        )
        svg.append(
            f'<polyline points="{coordinates}" fill="none" stroke="{PIPELINE_COLORS[pipeline]}" stroke-width="3"/>'
        )
        for row in group:
            x = _scale(float(row["observed_hit_rate"]), 0.5, 1.0, plot_left, plot_right)
            y = _scale(
                float(row["measured_relative_runtime"]),
                0,
                maximum,
                plot_bottom,
                plot_top,
            )
            svg.append(
                f'<circle cx="{x:.1f}" cy="{y:.1f}" r="5" fill="white" stroke="{PIPELINE_COLORS[pipeline]}" stroke-width="3"/>'
            )
        details = model["pipelines"][pipeline]
        svg.append(
            f'<text class="small" x="{plot_right}" y="{top + 27}" text-anchor="end">uncached {float(details["uncached_wall_seconds"]):.1f}s · overhead {float(details["lookup_preprocessing_overhead_seconds"]):.1f}s</text>'
        )
        svg.append(
            f'<text class="small" transform="translate({left + 16} {plot_bottom}) rotate(-90)">Relative runtime</text>'
        )
        svg.append(
            f'<text class="a" x="{(plot_left + plot_right) / 2}" y="{bottom - 8}" text-anchor="middle">Observed cache hit rate</text>'
        )
    svg.append(
        '<text class="small" x="65" y="755">Synthetic delays emit no annotation fields and run only for variants sent to VEP; semantic equivalence is required for every point.</text>'
    )
    svg.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(svg) + "\n")


def run(args: argparse.Namespace) -> None:
    """Write model source data, a machine-readable fit, and the figure."""
    rows = read_rows(args.metrics)
    modeled, model = fit_models(rows)
    args.output.mkdir(parents=True, exist_ok=True)
    write_tsv(args.output / "controlled_runtime_model.tsv", modeled)
    write_json_atomic(args.output / "runtime_model.json", model)
    render_svg(args.output / "controlled_runtime.svg", modeled, model)
    print(f"Wrote controlled-runtime publication outputs under {args.output}")


def parser() -> argparse.ArgumentParser:
    """Build the analysis CLI."""
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--metrics", type=Path, required=True)
    result.add_argument("--output", type=Path, required=True)
    return result


def main() -> None:
    """Run controlled-runtime analysis."""
    args = parser().parse_args()
    run(args)


if __name__ == "__main__":
    main()
