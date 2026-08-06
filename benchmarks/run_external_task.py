#!/usr/bin/env python3
"""Execute one four-condition external-WGS benchmark task."""

from __future__ import annotations

import argparse
import csv
import json
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from benchmarks.run_pilot import (
    PilotConfig,
    git_output,
    preflight,
    run_one,
    semantic_compare,
    sha256sum,
    write_json_atomic,
)


def read_task(path: Path, task_id: int) -> dict[str, str]:
    """Read and validate one stable task-manifest row."""
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not 0 <= task_id < len(rows):
        raise IndexError(f"Task {task_id} is absent from {path}")
    row = rows[task_id]
    if int(row["task_id"]) != task_id:
        raise RuntimeError(f"Task manifest is not position stable at {task_id}")
    return row


def load_strategies(
    path: Path, cohort: str, assembly: str
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Load common bundled strategies plus the cohort-specific custom cache."""
    document = json.loads(path.read_text())
    if document["cohort_assemblies"].get(cohort) != assembly:
        raise RuntimeError(f"Strategy assembly mismatch for {cohort}: {assembly}")
    strategies = list(document["bundled_strategies_by_assembly"][assembly])
    strategies.append(document["cohort_strategies"][cohort])
    names = [item["name"] for item in strategies]
    if names != ["gnomad_af_0.1", "gnomad_af_0.01", "cohort_3_genomes"]:
        raise RuntimeError(f"Unexpected strategy order: {names}")
    for strategy in strategies:
        if strategy.get("assembly") != assembly:
            raise RuntimeError(
                f"Strategy {strategy['name']} is not built for {assembly}"
            )
        cache = Path(strategy["cache_dir"])
        required = (
            cache / "annotation.yaml",
            cache / "params.snapshot.yaml",
            cache / "vcfcache_annotated.bcf",
            cache / "vcfcache_annotated.bcf.csi",
        )
        missing = [str(item) for item in required if not item.exists()]
        if missing:
            raise FileNotFoundError(
                f"Incomplete strategy {strategy['name']}: {missing}"
            )
        observed = sha256sum(cache / "annotation.yaml")
        if observed != strategy["annotation_yaml_sha256"]:
            raise RuntimeError(f"Annotation recipe changed for {strategy['name']}")
        provenance_path = Path(strategy["provenance_path"])
        provenance = json.loads(provenance_path.read_text())
        for key, expected in strategy["provenance_expected"].items():
            if provenance.get(key) != expected:
                raise RuntimeError(
                    f"Provenance mismatch for {strategy['name']} at {key}: "
                    f"{provenance.get(key)!r} != {expected!r}"
                )
    hashes = {item["annotation_yaml_sha256"] for item in strategies}
    if len(hashes) != 1:
        raise RuntimeError(
            "All cache strategies must use the identical annotation.yaml"
        )
    return document, strategies


def _config(
    run_root: Path,
    input_vcf: Path,
    strategy_name: str,
    cache: Path,
    params_file: Path,
    replicate: int,
) -> PilotConfig:
    return PilotConfig(
        data_root=run_root / "runs" / strategy_name,
        input_vcf=input_vcf,
        cache_dir=cache,
        params_file=params_file,
        replicate=replicate,
    )


def summarize_completed_runs(
    *,
    task: dict[str, str],
    document: dict[str, Any],
    strategies: list[dict[str, Any]],
    execution_order: list[str],
    run_dirs: dict[str, Path],
    run_root: Path,
    comparison_workers: int = 1,
) -> dict[str, Any]:
    """Validate four completed conditions and write their external summary."""
    baseline_dir = run_dirs["uncached"]
    baseline = baseline_dir / "output.bcf"
    baseline_metrics = json.loads((baseline_dir / "metrics.json").read_text())
    if comparison_workers < 1:
        raise ValueError("comparison_workers must be at least one")
    cached_outputs = {
        strategy["name"]: run_dirs[strategy["name"]] / "output.bcf"
        for strategy in strategies
    }
    if comparison_workers == 1:
        comparisons = {
            name: semantic_compare(output, baseline)
            for name, output in cached_outputs.items()
        }
    else:
        with ProcessPoolExecutor(max_workers=comparison_workers) as executor:
            futures = {
                name: executor.submit(semantic_compare, output, baseline)
                for name, output in cached_outputs.items()
            }
            comparisons = {name: future.result() for name, future in futures.items()}

    rows = []
    for strategy in strategies:
        name = strategy["name"]
        run_dir = run_dirs[name]
        metrics = json.loads((run_dir / "metrics.json").read_text())
        comparison = comparisons[name]
        comparison_path = run_root / f"semantic_{name}.json"
        original_path = run_root / f"semantic_{name}.original_failed.json"
        if comparison_path.exists() and not original_path.exists():
            original_path.write_bytes(comparison_path.read_bytes())
        write_json_atomic(comparison_path, comparison)
        if comparison.get("semantic_pass") is not True:
            raise RuntimeError(f"Semantic comparison failed: {comparison_path}")
        cached_seconds = float(metrics["wall_seconds"])
        uncached_seconds = float(baseline_metrics["wall_seconds"])
        rows.append(
            {
                "cohort": task["cohort"],
                "tool": task.get("tool", document.get("tool", "vep")),
                "assembly": task["assembly"],
                "sample": task["sample"],
                "phase": task["phase"],
                "replicate": int(task["replicate"]),
                "strategy": name,
                "strategy_kind": strategy["kind"],
                "bundled_alias": strategy.get("alias", ""),
                "zenodo_doi": strategy.get("doi", ""),
                "input_records": int(task["input_records"]),
                "cache_hit_rate": metrics.get("cache_hit_rate"),
                "cached_wall_seconds": cached_seconds,
                "uncached_wall_seconds": uncached_seconds,
                "speedup": uncached_seconds / cached_seconds,
                "relative_runtime": cached_seconds / uncached_seconds,
                "semantic_pass": True,
                "raw_annotation_mismatches": comparison.get(
                    "raw_annotation_mismatches", 0
                ),
                "ignored_annotation_mismatches": comparison.get(
                    "ignored_annotation_mismatches", 0
                ),
            }
        )
    summary = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "campaign_commit": document["commit"],
        "comparison_commit": git_output("rev-parse", "HEAD"),
        "task": task,
        "execution_order": execution_order,
        "rows": rows,
    }
    write_json_atomic(run_root / "external_summary.json", summary)
    return summary


def execute(args: argparse.Namespace) -> dict[str, Any]:
    """Run one uncached baseline and three cached conditions exactly once."""
    task = read_task(args.task_manifest, args.task_id)
    document, strategies = load_strategies(
        args.strategies, task["cohort"], task["assembly"]
    )
    if task.get("tool", "vep") != document.get("tool", "vep"):
        raise RuntimeError("Task annotator does not match strategy bundle")
    input_vcf = Path(task["input_vcf"])
    if sha256sum(input_vcf) != task["input_sha256"]:
        raise RuntimeError(f"Input checksum mismatch: {input_vcf}")
    replicate = int(task["replicate"])
    runtime_params = document["runtime_params_by_assembly"][task["assembly"]]
    params_file = Path(runtime_params["path"])
    if not params_file.exists():
        raise FileNotFoundError(f"Runtime params are missing: {params_file}")
    if sha256sum(params_file) != runtime_params["sha256"]:
        raise RuntimeError(f"Runtime params changed: {params_file}")
    order = task["strategy_order"].split(",")
    expected_order = {"uncached", *(item["name"] for item in strategies)}
    if len(order) != 4 or set(order) != expected_order:
        raise RuntimeError(f"Invalid four-condition order: {order}")
    by_name = {item["name"]: item for item in strategies}
    baseline_strategy = by_name["gnomad_af_0.01"]
    completed: dict[str, Path] = {}
    args.run_root.mkdir(parents=True, exist_ok=True)
    for name in order:
        strategy = baseline_strategy if name == "uncached" else by_name[name]
        cache = Path(strategy["cache_dir"])
        config = _config(args.run_root, input_vcf, name, cache, params_file, replicate)
        preflight(config)
        run_one(config, "uncached" if name == "uncached" else "cached")
        completed[name] = config.run_dir("uncached" if name == "uncached" else "cached")

    return summarize_completed_runs(
        task=task,
        document=document,
        strategies=strategies,
        execution_order=order,
        run_dirs=completed,
        run_root=args.run_root,
    )


def parser() -> argparse.ArgumentParser:
    """Build the task-runner CLI."""
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--task-manifest", type=Path, required=True)
    result.add_argument("--strategies", type=Path, required=True)
    result.add_argument("--task-id", type=int, required=True)
    result.add_argument("--run-root", type=Path, required=True)
    return result


def main() -> None:
    """Run one Slurm task."""
    args = parser().parse_args()
    execute(args)


if __name__ == "__main__":
    main()
