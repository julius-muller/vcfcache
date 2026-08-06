#!/usr/bin/env python3
"""Collect compact, provenance-rich benchmark tables for figure development."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


def read_json(path: Path) -> dict[str, Any]:
    """Read one JSON object."""
    return json.loads(path.read_text())


def sha256(path: Path) -> str:
    """Return a file SHA-256 digest."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    """Write a nonempty TSV atomically."""
    if not rows:
        raise RuntimeError(f"Refusing to write an empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_suffix(path.suffix + ".partial")
    with partial.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    partial.replace(path)


def assay_label(input_vcf: str) -> str:
    """Infer the assay from a frozen prepared-data path."""
    if "/panel_acmg_sf_v3.3/" in input_vcf:
        return "panel"
    if "/wes_twist_core_targets/" in input_vcf:
        return "wes_target_control"
    if "/wes_twist_core/" in input_vcf:
        return "wes"
    if "/hprc_r2/" in input_vcf:
        return "hprc_wgs"
    return "wgs"


def paired_rows(campaign_root: Path, phase: str) -> list[dict[str, Any]]:
    """Collect every finalized, semantically valid paired task in a phase."""
    campaign = read_json(campaign_root / "campaign.json")
    task_roots = sorted(
        (campaign_root / phase / "tasks").glob("task-*"),
        key=lambda path: int(path.name.removeprefix("task-")),
    )
    rows: list[dict[str, Any]] = []
    for root in task_roots:
        metadata = read_json(root / "slurm-task.json")
        summaries = list(root.glob("**/summary_r*.json"))
        if len(summaries) != 1:
            raise RuntimeError(f"Expected one summary under {root}, found {summaries}")
        summary = read_json(summaries[0])
        comparison = summary.get("semantic_comparison", {})
        if comparison.get("semantic_pass") is not True:
            raise RuntimeError(f"Paired task is not semantically valid: {root}")
        uncached = float(summary["uncached_wall_seconds"])
        cached = float(summary["cached_wall_seconds"])
        rows.append(
            {
                "campaign_id": campaign["campaign_id"],
                "campaign_commit": campaign.get("commit", ""),
                "run_commit": summary["commit"],
                "phase": phase,
                "task_id": metadata["task_id"],
                "sample": summary["sample"],
                "assay": assay_label(metadata["input_vcf"]),
                "population": metadata.get("population", ""),
                "superpopulation": metadata.get("superpopulation", ""),
                "input_records": summary["input_records"],
                "replicate": summary["replicate"],
                "first_mode": metadata["first_mode"],
                "cache_hit_rate": summary["cache_hit_rate"],
                "cached_wall_seconds": cached,
                "uncached_wall_seconds": uncached,
                "relative_runtime": cached / uncached,
                "speedup": uncached / cached,
                "wall_seconds_saved": uncached - cached,
                "semantic_pass": True,
                "validation_status": "semantically_validated",
                "records_compared": comparison.get("records_compared", 0),
                "key_mismatches": comparison.get("key_mismatches", 0),
                "annotation_mismatches": comparison.get("annotation_mismatches", 0),
                "raw_annotation_mismatches": comparison.get(
                    "raw_annotation_mismatches", 0
                ),
                "ignored_annotation_mismatches": comparison.get(
                    "ignored_annotation_mismatches", 0
                ),
                "annotation_order_only": comparison.get("annotation_order_only", 0),
                "record_order_only_loci": comparison.get("record_order_only_loci", 0),
                "slurm_job_id": metadata["slurm_job_id"],
                "slurm_node": metadata["slurm_node"],
            }
        )
    return rows


def manifest_tasks(path: Path) -> dict[int, dict[str, str]]:
    """Index a frozen external task manifest by task ID."""
    with path.open(newline="") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        return {int(row["task_id"]): row for row in rows}


def external_validated_rows(task_root: Path) -> list[dict[str, Any]]:
    """Read rows produced by a completed external task."""
    summary = read_json(task_root / "external_summary.json")
    return [
        {
            **row,
            "task_id": int(summary["task"]["task_id"]),
            "validation_status": "semantically_validated",
            "slurm_job_id": read_json(task_root / "slurm-task.json")["slurm_job_id"],
            "slurm_node": read_json(task_root / "slurm-task.json")["slurm_node"],
        }
        for row in summary["rows"]
    ]


def latest_attempt(task_attempts: Path) -> Path | None:
    """Return the most recent archived job attempt for one task."""
    jobs = list(task_attempts.glob("job-*"))
    if not jobs:
        return None
    return max(jobs, key=lambda path: int(path.name.removeprefix("job-")))


def one_metrics(root: Path, mode: str) -> dict[str, Any]:
    """Read exactly one external strategy metrics file."""
    pattern = (
        "**/uncached_r01/metrics.json"
        if mode == "uncached"
        else "**/cached_r01/metrics.json"
    )
    matches = list(root.glob(pattern))
    if len(matches) != 1:
        raise RuntimeError(
            f"Expected one {mode} metrics file under {root}, found {matches}"
        )
    return read_json(matches[0])


def external_provisional_rows(
    attempt: Path, task: dict[str, str], task_id: int
) -> list[dict[str, Any]]:
    """Extract runtime metrics whose semantic comparison awaits fixed post-processing."""
    run_root = attempt / "run" / "runs"
    uncached = one_metrics(run_root / "uncached", "uncached")
    uncached_wall = float(uncached["wall_seconds"])
    marker = "success" if (attempt / "success.json").exists() else "failure"
    rows: list[dict[str, Any]] = []
    for strategy in ("gnomad_af_0.1", "gnomad_af_0.01", "cohort_3_genomes"):
        metrics = one_metrics(run_root / strategy, "cached")
        cached_wall = float(metrics["wall_seconds"])
        rows.append(
            {
                "assembly": task["assembly"],
                "bundled_alias": "",
                "cache_hit_rate": metrics["cache_hit_rate"],
                "cached_wall_seconds": cached_wall,
                "cohort": task["cohort"],
                "ignored_annotation_mismatches": "",
                "input_records": task["input_records"],
                "phase": task["phase"],
                "raw_annotation_mismatches": "",
                "relative_runtime": cached_wall / uncached_wall,
                "replicate": task["replicate"],
                "sample": task["sample"],
                "semantic_pass": "",
                "speedup": uncached_wall / cached_wall,
                "strategy": strategy,
                "strategy_kind": (
                    "custom_cohort"
                    if strategy == "cohort_3_genomes"
                    else "bundled_zenodo"
                ),
                "uncached_wall_seconds": uncached_wall,
                "zenodo_doi": "",
                "task_id": task_id,
                "validation_status": f"comparator_pending_{marker}",
                "slurm_job_id": attempt.name.removeprefix("job-"),
                "slurm_node": "",
            }
        )
    return rows


def external_rows(campaign_root: Path, phase: str) -> list[dict[str, Any]]:
    """Collect validated and comparator-pending completed external tasks."""
    tasks = manifest_tasks(campaign_root / "manifests" / f"{phase}.tsv")
    phase_root = campaign_root / phase
    finalized = {
        int(path.name.removeprefix("task-")): path
        for path in (phase_root / "tasks").glob("task-*")
    }
    attempts = {
        int(path.name.removeprefix("task-")): latest_attempt(path)
        for path in (phase_root / "attempts").glob("task-*")
    }
    rows: list[dict[str, Any]] = []
    for task_id in sorted(tasks):
        if task_id in finalized:
            rows.extend(external_validated_rows(finalized[task_id]))
            continue
        attempt = attempts.get(task_id)
        if attempt is None or not any(
            (attempt / name).exists() for name in ("success.json", "failure.json")
        ):
            continue
        rows.extend(external_provisional_rows(attempt, tasks[task_id], task_id))
    return rows


def count_values(rows: Iterable[dict[str, Any]], key: str) -> dict[str, int]:
    """Count values for snapshot provenance."""
    counts: dict[str, int] = {}
    for row in rows:
        value = str(row[key])
        counts[value] = counts.get(value, 0) + 1
    return dict(sorted(counts.items()))


def run(args: argparse.Namespace) -> None:
    """Write compact tables and a machine-readable snapshot manifest."""
    args.output.mkdir(parents=True, exist_ok=True)
    primary = paired_rows(args.primary, "warmup")
    assay = paired_rows(args.assay, "measured")
    external = external_rows(args.external, "warmup")
    expected_external = len(manifest_tasks(args.external / "manifests" / "warmup.tsv"))
    external_tasks = {int(row["task_id"]) for row in external}
    validated_external_tasks = {
        int(row["task_id"])
        for row in external
        if row["validation_status"] == "semantically_validated"
    }
    external_complete = (
        len(external_tasks) == expected_external
        and validated_external_tasks == external_tasks
        and len(external) == expected_external * 3
    )
    status = "FINAL" if external_complete else "PRELIMINARY"
    paths = {
        "primary_wgs_metrics.tsv": primary,
        "assay_metrics.tsv": assay,
        "external_wgs_metrics.tsv": external,
    }
    for name, rows in paths.items():
        write_tsv(args.output / name, rows)

    strategies_source = args.external / "manifests" / "strategies.json"
    strategies_target = args.output / "external_strategies.json"
    shutil.copyfile(strategies_source, strategies_target)
    files = [args.output / name for name in paths] + [strategies_target]
    snapshot = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "created_on_host": platform.node(),
        "status": status,
        "campaigns": {
            "primary": read_json(args.primary / "campaign.json")["campaign_id"],
            "assay": read_json(args.assay / "campaign.json")["campaign_id"],
            "external": read_json(args.external / "campaign.json")["campaign_id"],
        },
        "row_counts": {name: len(rows) for name, rows in paths.items()},
        "sample_counts": {
            "primary_wgs": len({row["sample"] for row in primary}),
            "assay_by_type": count_values(assay, "assay"),
            "external_completed": len(external_tasks),
            "external_semantically_validated": len(validated_external_tasks),
            "external_expected": expected_external,
        },
        "external_validation_status": count_values(external, "validation_status"),
        "notes": [
            "Primary WGS contains the 49 successful samples; HG02888 is the documented cgroup-OOM exclusion.",
            (
                "All external WGS tasks and all three cache strategies passed corrected semantic validation."
                if external_complete
                else "Comparator-pending external rows are exploratory only until corrected semantic post-processing passes."
            ),
            "No BCF or other bulky run output is included in this snapshot.",
        ],
        "files": {
            path.name: {"bytes": path.stat().st_size, "sha256": sha256(path)}
            for path in files
        },
    }
    snapshot_path = args.output / "SNAPSHOT.json"
    snapshot_path.write_text(json.dumps(snapshot, indent=2, sort_keys=True) + "\n")
    print(json.dumps(snapshot, indent=2, sort_keys=True))


def parser() -> argparse.ArgumentParser:
    """Build the command-line parser."""
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--primary", type=Path, required=True)
    result.add_argument("--assay", type=Path, required=True)
    result.add_argument("--external", type=Path, required=True)
    result.add_argument("--output", type=Path, required=True)
    return result


if __name__ == "__main__":
    run(parser().parse_args())
