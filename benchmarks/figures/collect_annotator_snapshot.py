#!/usr/bin/env python3
"""Freeze fastVEP, VEP, and statistics-calibration figure source data."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

STRATEGIES = ("gnomad_af_0.1", "gnomad_af_0.01", "cohort_3_genomes")


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


def read_tsv(path: Path) -> list[dict[str, str]]:
    """Read a TSV into dictionaries."""
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    """Write a non-empty TSV atomically."""
    if not rows:
        raise RuntimeError(f"Refusing to write an empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_suffix(path.suffix + ".partial")
    with partial.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    partial.replace(path)


def fastvep_rows(campaign_root: Path) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Collect and validate every finalized fastVEP task."""
    campaign = read_json(campaign_root / "campaign.json")
    if campaign.get("tool") != "fastvep":
        raise RuntimeError("The supplied campaign is not marked as fastVEP")
    if campaign.get("statistics_mode") != "light":
        raise RuntimeError("The fastVEP campaign did not use --statistics light")

    phase = next(iter(campaign["phases"]))
    expected = int(campaign["phases"][phase]["tasks"])
    task_roots = sorted(
        (campaign_root / phase / "tasks").glob("task-*"),
        key=lambda path: int(path.name.removeprefix("task-")),
    )
    if len(task_roots) != expected:
        raise RuntimeError(
            f"Expected {expected} finalized fastVEP tasks, found {len(task_roots)}"
        )

    rows: list[dict[str, Any]] = []
    comparator_versions: set[str] = set()
    comparison_commits: set[str] = set()
    for root in task_roots:
        summary = read_json(root / "external_summary.json")
        slurm = read_json(root / "slurm-task.json")
        task = summary["task"]
        task_id = int(task["task_id"])
        if task_id != int(root.name.removeprefix("task-")):
            raise RuntimeError(f"Task ID mismatch under {root}")
        if task.get("tool") != "fastvep" or summary.get("statistics_mode") != "light":
            raise RuntimeError(f"Unexpected tool/statistics mode under {root}")
        if len(summary["rows"]) != len(STRATEGIES):
            raise RuntimeError(f"Expected three cache strategies under {root}")

        observed_strategies = {row["strategy"] for row in summary["rows"]}
        if observed_strategies != set(STRATEGIES):
            raise RuntimeError(f"Unexpected strategies under {root}: {observed_strategies}")
        comparator = str(summary.get("semantic_comparator", ""))
        if not comparator.startswith("fastvep_complete_record_and_header_"):
            raise RuntimeError(f"Non-strict fastVEP comparator under {root}: {comparator}")
        comparator_versions.add(comparator)
        comparison_commits.add(str(summary.get("comparison_commit", "")))

        for row in summary["rows"]:
            if row.get("semantic_pass") is not True:
                raise RuntimeError(f"Semantic validation failed under {root}")
            rows.append(
                {
                    **row,
                    "task_id": task_id,
                    "provider": task.get("provider", ""),
                    "validation_status": "strict_complete_record_and_header_validated",
                    "slurm_job_id": slurm["slurm_job_id"],
                    "slurm_node": slurm["slurm_node"],
                }
            )

    metadata = {
        "campaign_id": campaign["campaign_id"],
        "campaign_commit": campaign["commit"],
        "phase": phase,
        "expected_tasks": expected,
        "completed_tasks": len(task_roots),
        "semantic_comparators": sorted(comparator_versions),
        "comparison_commits": sorted(comparison_commits),
        "statistics_mode": campaign["statistics_mode"],
    }
    return rows, metadata


def vep_light_rows(campaign_root: Path) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Collect and validate a completed VEP light-statistics campaign."""
    campaign = read_json(campaign_root / "campaign.json")
    if campaign.get("tool") != "vep" or campaign.get("statistics_mode") != "light":
        raise RuntimeError("The calibration campaign is not VEP --statistics light")

    phase = next(iter(campaign["phases"]))
    expected = int(campaign["phases"][phase]["tasks"])
    task_roots = sorted(
        (campaign_root / phase / "tasks").glob("task-*"),
        key=lambda path: int(path.name.removeprefix("task-")),
    )
    if len(task_roots) != expected:
        raise RuntimeError(
            f"Expected {expected} finalized VEP light tasks, found {len(task_roots)}"
        )

    rows: list[dict[str, Any]] = []
    comparator_versions: set[str] = set()
    comparison_commits: set[str] = set()
    for root in task_roots:
        summary = read_json(root / "external_summary.json")
        slurm = read_json(root / "slurm-task.json")
        task = summary["task"]
        task_id = int(task["task_id"])
        if task_id != int(root.name.removeprefix("task-")):
            raise RuntimeError(f"Task ID mismatch under {root}")
        if task.get("tool") != "vep" or summary.get("statistics_mode") != "light":
            raise RuntimeError(f"Unexpected tool/statistics mode under {root}")
        if len(summary["rows"]) != len(STRATEGIES):
            raise RuntimeError(f"Expected three cache strategies under {root}")
        if {row["strategy"] for row in summary["rows"]} != set(STRATEGIES):
            raise RuntimeError(f"Unexpected strategies under {root}")

        comparator = str(summary.get("semantic_comparator", ""))
        if comparator != "vep_semantic_with_documented_exceptions_v1":
            raise RuntimeError(f"Unexpected VEP semantic comparator under {root}")
        comparator_versions.add(comparator)
        comparison_commits.add(str(summary.get("comparison_commit", "")))

        for row in summary["rows"]:
            if row.get("semantic_pass") is not True:
                raise RuntimeError(f"Semantic validation failed under {root}")
            rows.append(
                {
                    **row,
                    "task_id": task_id,
                    "provider": task.get("provider", ""),
                    "validation_status": "semantically_validated",
                    "slurm_job_id": slurm["slurm_job_id"],
                    "slurm_node": slurm["slurm_node"],
                }
            )

    metadata = {
        "campaign_id": campaign["campaign_id"],
        "campaign_commit": campaign["commit"],
        "phase": phase,
        "expected_tasks": expected,
        "completed_tasks": len(task_roots),
        "semantic_comparators": sorted(comparator_versions),
        "comparison_commits": sorted(comparison_commits),
        "statistics_mode": campaign["statistics_mode"],
    }
    return rows, metadata


def normalize_vep_rows(path: Path) -> list[dict[str, Any]]:
    """Add explicit annotator/method labels to the frozen VEP rows."""
    rows = read_tsv(path)
    if len(rows) != 52 * len(STRATEGIES):
        raise RuntimeError(f"Expected 156 VEP rows, found {len(rows)}")
    if {row["strategy"] for row in rows} != set(STRATEGIES):
        raise RuntimeError("The VEP snapshot does not contain the expected strategies")
    if any(row["validation_status"] != "semantically_validated" for row in rows):
        raise RuntimeError("The VEP snapshot contains unvalidated rows")
    for row in rows:
        row.update(
            {
                "tool": "vep",
                "statistics_mode": "full",
                "semantic_comparator": "vep_semantic_comparator",
            }
        )
    return rows


def comparison_rows(
    fastvep: list[dict[str, Any]], vep: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    """Create a common, paired table for R plotting."""
    columns = (
        "tool",
        "statistics_mode",
        "sample",
        "cohort",
        "assembly",
        "strategy",
        "strategy_kind",
        "input_records",
        "cache_hit_rate",
        "uncached_wall_seconds",
        "cached_wall_seconds",
        "relative_runtime",
        "speedup",
        "semantic_pass",
        "semantic_comparator",
        "validation_status",
    )
    combined: list[dict[str, Any]] = []
    for source in (vep, fastvep):
        for row in source:
            combined.append({column: row.get(column, "") for column in columns})

    keys_by_tool: dict[str, set[tuple[str, str, str]]] = {}
    for row in combined:
        keys_by_tool.setdefault(str(row["tool"]), set()).add(
            (str(row["sample"]), str(row["cohort"]), str(row["strategy"]))
        )
    if keys_by_tool.get("vep") != keys_by_tool.get("fastvep"):
        raise RuntimeError("VEP and fastVEP rows are not exactly paired")
    return combined


def statistics_calibration_rows(
    light: list[dict[str, Any]], full: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    """Pair light-statistics measurements with their original full runs."""
    full_by_key = {
        (str(row["sample"]), str(row["cohort"]), str(row["strategy"])): row
        for row in full
    }
    paired: list[dict[str, Any]] = []
    for light_row in light:
        key = (
            str(light_row["sample"]),
            str(light_row["cohort"]),
            str(light_row["strategy"]),
        )
        if key not in full_by_key:
            raise RuntimeError(f"VEP calibration row is absent from full snapshot: {key}")
        full_row = full_by_key[key]
        full_hit_rate = float(full_row["cache_hit_rate"])
        light_hit_rate = float(light_row["cache_hit_rate"])
        if abs(full_hit_rate - light_hit_rate) > 1e-12:
            raise RuntimeError(f"Cache hit rate changed between statistics modes: {key}")

        full_uncached = float(full_row["uncached_wall_seconds"])
        light_uncached = float(light_row["uncached_wall_seconds"])
        full_cached = float(full_row["cached_wall_seconds"])
        light_cached = float(light_row["cached_wall_seconds"])
        full_speedup = float(full_row["speedup"])
        light_speedup = float(light_row["speedup"])
        paired.append(
            {
                "sample": key[0],
                "cohort": key[1],
                "assembly": light_row["assembly"],
                "strategy": key[2],
                "cache_hit_rate": light_hit_rate,
                "full_uncached_wall_seconds": full_uncached,
                "light_uncached_wall_seconds": light_uncached,
                "full_cached_wall_seconds": full_cached,
                "light_cached_wall_seconds": light_cached,
                "full_relative_runtime": float(full_row["relative_runtime"]),
                "light_relative_runtime": float(light_row["relative_runtime"]),
                "full_speedup": full_speedup,
                "light_speedup": light_speedup,
                "speedup_ratio_light_over_full": light_speedup / full_speedup,
                "uncached_runtime_ratio_light_over_full": light_uncached / full_uncached,
                "cached_runtime_ratio_light_over_full": light_cached / full_cached,
                "semantic_pass": light_row["semantic_pass"],
            }
        )
    expected = len({(row["sample"], row["cohort"]) for row in light}) * len(
        STRATEGIES
    )
    if len(paired) != expected:
        raise RuntimeError(
            f"Expected {expected} paired statistics rows, found {len(paired)}"
        )
    return paired


def run(args: argparse.Namespace) -> None:
    """Build the frozen source-data bundle."""
    args.output.mkdir(parents=True, exist_ok=False)
    fastvep, fastvep_metadata = fastvep_rows(args.fastvep_campaign)
    vep_full_statistics = normalize_vep_rows(args.vep_metrics)

    vep_light_full: list[dict[str, Any]] = []
    vep_light_full_pairs: list[dict[str, Any]] = []
    vep_light_full_metadata: dict[str, Any] | None = None
    if args.vep_light_full_campaign is not None:
        vep_light_full, vep_light_full_metadata = vep_light_rows(
            args.vep_light_full_campaign
        )
        if len({row["sample"] for row in vep_light_full}) != 52:
            raise RuntimeError("The full VEP light campaign does not contain 52 samples")
        vep_light_full_pairs = statistics_calibration_rows(
            vep_light_full, vep_full_statistics
        )

    comparison_vep = vep_light_full or vep_full_statistics
    combined = comparison_rows(fastvep, comparison_vep)

    calibration: list[dict[str, Any]] = []
    calibration_pairs: list[dict[str, Any]] = []
    calibration_metadata: dict[str, Any] | None = None
    if args.vep_light_campaign is not None:
        calibration, calibration_metadata = vep_light_rows(args.vep_light_campaign)
        calibration_pairs = statistics_calibration_rows(
            calibration, vep_full_statistics
        )

    paths = {
        "fastvep_external_wgs_metrics.tsv": fastvep,
        "vep_external_wgs_full_statistics_metrics.tsv": vep_full_statistics,
        "annotator_external_wgs_metrics.tsv": combined,
    }
    if vep_light_full:
        paths.update(
            {
                "vep_external_wgs_light_statistics_metrics.tsv": vep_light_full,
                "vep_statistics_full_cohort_comparison.tsv": vep_light_full_pairs,
            }
        )
    if calibration:
        paths.update(
            {
                "vep_light_calibration_metrics.tsv": calibration,
                "vep_statistics_calibration.tsv": calibration_pairs,
            }
        )
    for name, rows in paths.items():
        write_tsv(args.output / name, rows)

    source_vep_snapshot = read_json(args.vep_snapshot)
    snapshot = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "created_on_host": platform.node(),
        "status": (
            "ANNOTATOR_LIGHT_MATCHED_FINAL"
            if vep_light_full_metadata is not None
            else (
                "ANNOTATOR_COMPARISON_FINAL"
                if calibration_metadata is not None
                else "FASTVEP_FINAL_VEP_COMPARISON_PROVISIONAL"
            )
        ),
        "fastvep": fastvep_metadata,
        "vep": {
            "source_snapshot": str(args.vep_snapshot),
            "source_snapshot_status": source_vep_snapshot["status"],
            "campaign_id": (
                vep_light_full_metadata["campaign_id"]
                if vep_light_full_metadata is not None
                else source_vep_snapshot["campaigns"]["external"]
            ),
            "statistics_mode": "light" if vep_light_full else "full",
            "completed_tasks": 52,
        },
        "vep_light_full": vep_light_full_metadata,
        "vep_light_calibration": calibration_metadata,
        "row_counts": {name: len(rows) for name, rows in paths.items()},
        "sample_counts": {
            "fastvep": len({row["sample"] for row in fastvep}),
            "vep": len({row["sample"] for row in comparison_vep}),
            "vep_light_full": len({row["sample"] for row in vep_light_full}),
            "vep_light_calibration": len({row["sample"] for row in calibration}),
        },
        "notes": [
            "The standalone fastVEP dataset is final: all 52 genomes and 156 cached outputs passed strict complete-record and relevant-header comparison.",
            (
                "The paired annotator view uses matched --statistics light measurements for all 52 VEP and 52 fastVEP genomes."
                if vep_light_full_metadata is not None
                else (
                    "The paired annotator view retains the measured VEP --statistics full values and accompanies them with a completed six-genome --statistics light sensitivity calibration."
                    if calibration_metadata is not None
                    else "The paired annotator view is provisional because VEP used --statistics full whereas fastVEP used --statistics light."
                )
            ),
            "The same genomes, variants, cache blueprints, and three cache strategies are paired across annotators.",
            "The VEP statistics sensitivity comparison is paired within sample and strategy; it is not used to rescale any observation.",
            "The full VEP light campaign is analyzed as measured and paired with the historical full-statistics campaign only for sensitivity reporting.",
            "No technical replicate is included; genomes are the inferential units.",
        ],
    }
    snapshot["files"] = {
        name: {
            "bytes": (args.output / name).stat().st_size,
            "sha256": sha256(args.output / name),
        }
        for name in paths
    }
    (args.output / "SNAPSHOT.json").write_text(
        json.dumps(snapshot, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(snapshot, indent=2, sort_keys=True))


def parser() -> argparse.ArgumentParser:
    """Build the command-line parser."""
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--fastvep-campaign", type=Path, required=True)
    result.add_argument("--vep-light-campaign", type=Path)
    result.add_argument("--vep-light-full-campaign", type=Path)
    result.add_argument("--vep-metrics", type=Path, required=True)
    result.add_argument("--vep-snapshot", type=Path, required=True)
    result.add_argument("--output", type=Path, required=True)
    return result


if __name__ == "__main__":
    run(parser().parse_args())
