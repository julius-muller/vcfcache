#!/usr/bin/env python3
"""Run one paired independent Panel/WES annotator benchmark task."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from benchmarks.run_pilot import (
    PilotConfig,
    preflight,
    run_one,
    semantic_compare,
    strict_semantic_compare,
    write_json_atomic,
)


def task_row(path: Path, task_id: int) -> dict[str, str]:
    """Read one position-stable task row."""
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not 0 <= task_id < len(rows) or int(rows[task_id]["task_id"]) != task_id:
        raise RuntimeError(f"Task {task_id} is absent or unstable in {path}")
    return rows[task_id]


def config_for(row: dict[str, str], run_root: Path) -> PilotConfig:
    """Resolve a task row to the immutable pilot configuration."""
    return PilotConfig(
        data_root=run_root,
        input_vcf=Path(row["input_vcf"]),
        cache_dir=Path(row["cache_dir"]),
        params_file=Path(row["params_file"]),
        replicate=int(row["replicate"]),
    )


def warmup_config_for(row: dict[str, str], run_root: Path) -> PilotConfig:
    """Resolve the separate, excluded calibration input used before timing."""
    return PilotConfig(
        data_root=run_root / "untimed_warmup",
        input_vcf=Path(row["warmup_input_vcf"]),
        cache_dir=Path(row["cache_dir"]),
        params_file=Path(row["params_file"]),
        replicate=1,
    )


def finalize(row: dict[str, str], config: PilotConfig, run_root: Path) -> None:
    """Compare paired outputs strictly and freeze the task summary."""
    uncached = json.loads((config.run_dir("uncached") / "metrics.json").read_text())
    cached = json.loads((config.run_dir("cached") / "metrics.json").read_text())
    comparator = strict_semantic_compare if row["tool"] == "fastvep" else semantic_compare
    comparison = comparator(
        config.run_dir("cached") / "output.bcf",
        config.run_dir("uncached") / "output.bcf",
    )
    write_json_atomic(config.comparison_path, comparison)
    if comparison.get("semantic_pass") is not True:
        raise RuntimeError(f"Semantic comparison failed: {config.comparison_path}")
    uncached_seconds = float(uncached["wall_seconds"])
    cached_seconds = float(cached["wall_seconds"])
    summary = {
        **{key: row[key] for key in (
            "task_id", "tool", "assay", "cohort", "sample", "assembly",
            "cache_strategy", "cache_kind", "cache_alias", "condition_order",
        )},
        "replicate": int(row["replicate"]),
        "input_vcf": row["input_vcf"],
        "untimed_warmup_input_vcf": row["warmup_input_vcf"],
        "untimed_warmup_excluded": True,
        "input_records": int(uncached["output_records"]),
        "cache_hit_rate": cached.get("cache_hit_rate"),
        "uncached_wall_seconds": uncached_seconds,
        "cached_wall_seconds": cached_seconds,
        "relative_runtime": cached_seconds / uncached_seconds,
        "speedup": uncached_seconds / cached_seconds,
        "wall_seconds_saved": uncached_seconds - cached_seconds,
        "semantic_pass": True,
        "semantic_comparator": comparison.get("comparator"),
        "records_compared": comparison.get("records_compared"),
        "raw_annotation_mismatches": comparison.get("raw_annotation_mismatches", 0),
        "ignored_annotation_mismatches": comparison.get(
            "ignored_annotation_mismatches", 0
        ),
    }
    write_json_atomic(run_root / "external_assay_summary.json", summary)
    print(json.dumps(summary, indent=2, sort_keys=True))


def main() -> None:
    """Run, compare, or report one task condition."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task-manifest", type=Path, required=True)
    parser.add_argument("--task-id", type=int, required=True)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument(
        "--condition", choices=("warmup", "uncached", "cached", "finalize")
    )
    parser.add_argument("--print-order", action="store_true")
    args = parser.parse_args()
    row = task_row(args.task_manifest, args.task_id)
    if args.print_order:
        print(row["condition_order"].replace(",", "\n"))
        return
    if args.condition is None:
        parser.error("--condition is required unless --print-order is used")
    config = config_for(row, args.run_root)
    if args.condition == "finalize":
        finalize(row, config, args.run_root)
        return
    if args.condition == "warmup":
        warmup = warmup_config_for(row, args.run_root)
        preflight(warmup)
        run_one(warmup, "cached")
        return
    preflight(config)
    run_one(config, args.condition)


if __name__ == "__main__":
    main()
