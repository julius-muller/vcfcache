#!/usr/bin/env python3
"""Run one baseline or cached controlled-runtime condition."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

from benchmarks.run_pilot import (
    PilotConfig,
    preflight,
    run_one,
    semantic_compare,
    sha256sum,
    write_json_atomic,
)


def read_task(path: Path, task_id: int) -> dict[str, str]:
    """Read one stable task row."""
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not 0 <= task_id < len(rows) or int(rows[task_id]["task_id"]) != task_id:
        raise RuntimeError(f"Task {task_id} is absent or unstable in {path}")
    return rows[task_id]


def _symlink(source: Path, destination: Path, root: Path) -> None:
    destination.symlink_to(source.relative_to(root))


def execute(args: argparse.Namespace) -> dict[str, Any]:
    """Execute and validate one controlled-runtime task."""
    task = read_task(args.task_manifest, args.task_id)
    input_vcf = Path(task["input_vcf"])
    cache = Path(task["cache_dir"])
    params = Path(task["params_file"])
    if sha256sum(input_vcf) != task["input_sha256"]:
        raise RuntimeError(f"Controlled input checksum mismatch: {input_vcf}")
    metadata = json.loads((cache / "controlled_cache.json").read_text())
    if metadata["pipeline"] != task["pipeline"]:
        raise RuntimeError(f"Controlled cache pipeline mismatch: {cache}")
    if float(metadata["target_hit_rate"]) != float(task["target_hit_rate"]):
        raise RuntimeError(f"Controlled cache hit-rate mismatch: {cache}")
    config = PilotConfig(
        data_root=args.run_root / "run",
        input_vcf=input_vcf,
        cache_dir=cache,
        params_file=params,
        replicate=1,
    )
    preflight(config)
    mode = task["mode"]
    metrics = run_one(config, mode)
    output = config.run_dir(mode) / "output.bcf"
    summary: dict[str, Any] = {
        "task": task,
        "metrics": metrics,
        "pipeline": task["pipeline"],
        "mode": mode,
        "target_hit_rate": float(task["target_hit_rate"]),
        "delay_us": int(task["delay_us"]),
    }
    if mode == "uncached":
        _symlink(output, args.run_root / "baseline_output.bcf", args.run_root)
        _symlink(
            Path(f"{output}.csi"),
            args.run_root / "baseline_output.bcf.csi",
            args.run_root,
        )
        summary.update(
            {
                "uncached_wall_seconds": float(metrics["wall_seconds"]),
                "semantic_pass": None,
            }
        )
    else:
        baseline_root = Path(task["baseline_result"])
        baseline_output = baseline_root / "baseline_output.bcf"
        baseline_summary = json.loads(
            (baseline_root / "controlled_summary.json").read_text()
        )
        comparison = semantic_compare(output, baseline_output)
        write_json_atomic(args.run_root / "semantic_comparison.json", comparison)
        if comparison.get("semantic_pass") is not True:
            raise RuntimeError(f"Controlled semantic comparison failed: {comparison}")
        cached_seconds = float(metrics["wall_seconds"])
        uncached_seconds = float(baseline_summary["uncached_wall_seconds"])
        summary.update(
            {
                "observed_hit_rate": metrics.get("cache_hit_rate"),
                "cached_wall_seconds": cached_seconds,
                "uncached_wall_seconds": uncached_seconds,
                "relative_runtime": cached_seconds / uncached_seconds,
                "speedup": uncached_seconds / cached_seconds,
                "semantic_pass": True,
                "semantic_comparison": comparison,
            }
        )
    write_json_atomic(args.run_root / "controlled_summary.json", summary)
    return summary


def parser() -> argparse.ArgumentParser:
    """Build the controlled-task CLI."""
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--task-manifest", type=Path, required=True)
    result.add_argument("--task-id", type=int, required=True)
    result.add_argument("--run-root", type=Path, required=True)
    return result


def main() -> None:
    """Run one controlled task."""
    args = parser().parse_args()
    args.run_root.mkdir(parents=True, exist_ok=True)
    execute(args)


if __name__ == "__main__":
    main()
