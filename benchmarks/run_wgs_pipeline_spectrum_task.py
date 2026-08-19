#!/usr/bin/env python3
"""Run or compare one condition of a paired real-WGS virtual pipeline load."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

from benchmarks.run_pilot import (
    KNOWN_VEP_IGNORED_CSQ_FIELDS,
    KNOWN_VEP_UNORDERED_CSQ_FIELDS,
    PUBLICATION_STATISTICS_MODE,
    PilotConfig,
    _canonical_record,
    _csq_fields,
    _csq_header,
    _locus_groups,
    _query_process,
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


def config_for(task: dict[str, str], run_root: Path) -> PilotConfig:
    """Build the immutable pilot configuration for this task."""
    return PilotConfig(
        data_root=run_root / "run",
        input_vcf=Path(task["input_vcf"]),
        cache_dir=Path(task["cache_dir"]),
        params_file=Path(task["params_file"]),
        replicate=1,
    )


def execute_mode(args: argparse.Namespace, task: dict[str, str]) -> None:
    """Run exactly one separately-accounted cached or uncached cell."""
    if args.mode not in {"cached", "uncached"}:
        raise ValueError(args.mode)
    input_vcf = Path(task["input_vcf"])
    if sha256sum(input_vcf) != task["input_sha256"]:
        raise RuntimeError(f"WGS spectrum input checksum mismatch: {input_vcf}")
    config = config_for(task, args.run_root)
    preflight(config)
    metrics = run_one(config, args.mode)
    write_json_atomic(args.run_root / f"{args.mode}_metrics.json", metrics)


def semantic_fingerprint(path: Path) -> dict[str, Any]:
    """Return the same canonical VEP digest used by the paired comparator."""
    header = _csq_header(path)
    fields = _csq_fields(header)
    ignored_indices = tuple(
        fields.index(field)
        for field in KNOWN_VEP_IGNORED_CSQ_FIELDS
        if field in fields
    )
    unordered_indices = tuple(
        fields.index(field)
        for field in KNOWN_VEP_UNORDERED_CSQ_FIELDS
        if field in fields
    )
    process = _query_process(path)
    assert process.stdout is not None
    digest = hashlib.sha256()
    records = 0
    for _, group in _locus_groups(process.stdout):
        group.sort(key=lambda value: tuple(value[:8]))
        for value in group:
            key, canonical = _canonical_record(
                value, ignored_indices, unordered_indices
            )
            digest.update(("\t".join(key) + "\t" + canonical + "\n").encode())
            records += 1
    stderr = process.communicate()[1]
    if process.returncode:
        raise RuntimeError(f"bcftools fingerprint query failed: {stderr}")
    return {
        "algorithm": "vep_semantic_fingerprint_v1",
        "semantic_sha256": digest.hexdigest(),
        "csq_header_sha256": hashlib.sha256(header.encode()).hexdigest(),
        "records": records,
        "ignored_csq_fields": [
            field for field in KNOWN_VEP_IGNORED_CSQ_FIELDS if field in fields
        ],
        "unordered_csq_fields": [
            field for field in KNOWN_VEP_UNORDERED_CSQ_FIELDS if field in fields
        ],
    }


def summarize_cached_extension(
    args: argparse.Namespace, task: dict[str, str]
) -> dict[str, Any]:
    """Summarize a new cached cell against a frozen direct baseline."""
    config = config_for(task, args.run_root)
    metrics = json.loads((args.run_root / "cached_metrics.json").read_text())
    fingerprint = semantic_fingerprint(config.run_dir("cached") / "output.bcf")
    cached = float(metrics["wall_seconds"])
    uncached = float(task["baseline_uncached_wall_seconds"])
    hit_rate = float(metrics["cache_hit_rate"])
    expected_hit_rate = float(task["expected_hit_rate"])
    if abs(hit_rate - expected_hit_rate) > 1e-12:
        raise RuntimeError(
            f"Unexpected hit rate for task {task['task_id']}: "
            f"{hit_rate} != {expected_hit_rate}"
        )
    if int(metrics["output_records"]) != int(task["input_records"]):
        raise RuntimeError("Cached extension output record count changed")
    summary: dict[str, Any] = {
        "sample": task["sample"],
        "cohort": task["cohort"],
        "assembly": task["assembly"],
        "pipeline": task["pipeline"],
        "delay_us": int(task["delay_us"]),
        "input_records": int(metrics["output_records"]),
        "cache_hit_rate": hit_rate,
        "uncached_wall_seconds": uncached,
        "cached_wall_seconds": cached,
        "relative_runtime": cached / uncached,
        "speedup": uncached / cached,
        "time_saved_seconds": uncached - cached,
        "statistics_mode": PUBLICATION_STATISTICS_MODE,
        "baseline_reused": True,
        "baseline_campaign_id": task["baseline_campaign_id"],
        "baseline_summary_sha256": task["baseline_summary_sha256"],
        "semantic_validation": "pending_composed_validation",
        "source_semantic_pass": task["source_semantic_pass"].lower() == "true",
        "source_semantic_campaign_id": task["source_semantic_campaign_id"],
        "source_semantic_comparator": task["source_semantic_comparator"],
        "semantic_fingerprint": fingerprint,
        "task": task,
    }
    write_json_atomic(
        args.run_root / "wgs_pipeline_cached_summary.json", summary
    )
    return summary


def compare(args: argparse.Namespace, task: dict[str, str]) -> dict[str, Any]:
    """Validate the paired outputs and write the final task summary."""
    config = config_for(task, args.run_root)
    cached_metrics = json.loads((args.run_root / "cached_metrics.json").read_text())
    uncached_metrics = json.loads(
        (args.run_root / "uncached_metrics.json").read_text()
    )
    comparison = semantic_compare(
        config.run_dir("cached") / "output.bcf",
        config.run_dir("uncached") / "output.bcf",
    )
    write_json_atomic(args.run_root / "semantic_comparison.json", comparison)
    if comparison.get("semantic_pass") is not True:
        raise RuntimeError(f"WGS pipeline semantic comparison failed: {comparison}")
    cached = float(cached_metrics["wall_seconds"])
    uncached = float(uncached_metrics["wall_seconds"])
    summary: dict[str, Any] = {
        "sample": "KPGP-00319",
        "cohort": "KPGP",
        "assembly": "GRCh38",
        "pipeline": task["pipeline"],
        "delay_us": int(task["delay_us"]),
        "input_records": int(cached_metrics["output_records"]),
        "cache_hit_rate": float(cached_metrics["cache_hit_rate"]),
        "uncached_wall_seconds": uncached,
        "cached_wall_seconds": cached,
        "relative_runtime": cached / uncached,
        "speedup": uncached / cached,
        "time_saved_seconds": uncached - cached,
        "semantic_pass": True,
        "semantic_comparator": comparison["comparator"],
        "statistics_mode": PUBLICATION_STATISTICS_MODE,
        "first_mode": task["first_mode"],
        "second_mode": task["second_mode"],
        "task": task,
    }
    write_json_atomic(args.run_root / "wgs_pipeline_summary.json", summary)
    return summary


def main() -> None:
    """Run one mode or compare the completed pair."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task-manifest", type=Path, required=True)
    parser.add_argument("--task-id", type=int, required=True)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument(
        "--mode",
        choices=("cached", "uncached", "compare", "cached-extension-summary"),
        required=True,
    )
    args = parser.parse_args()
    args.run_root.mkdir(parents=True, exist_ok=True)
    task = read_task(args.task_manifest, args.task_id)
    if args.mode == "compare":
        compare(args, task)
    elif args.mode == "cached-extension-summary":
        summarize_cached_extension(args, task)
    else:
        execute_mode(args, task)


if __name__ == "__main__":
    main()
