#!/usr/bin/env python3
"""Collect a complete paired Slurm phase into publication-source tables."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any


def assay_label(input_vcf: str) -> str:
    """Infer the frozen assay label from its prepared-data path."""
    if "/panel_acmg_sf_v3.3/" in input_vcf:
        return "panel"
    if "/wes_twist_core_targets/" in input_vcf:
        return "wes_target_control"
    if "/wes_twist_core/" in input_vcf:
        return "wes"
    if "/hprc_r2/" in input_vcf:
        return "hprc_wgs"
    return "wgs"


def collect_rows(campaign_root: Path, phase: str) -> list[dict[str, Any]]:
    """Read every completed task and require the frozen phase to be complete."""
    campaign = json.loads((campaign_root / "campaign.json").read_text())
    expected = int(campaign["phases"][phase]["tasks"])
    task_roots = sorted(
        (campaign_root / phase / "tasks").glob("task-*"),
        key=lambda path: int(path.name.removeprefix("task-")),
    )
    if len(task_roots) != expected:
        raise RuntimeError(
            f"Expected {expected} {phase} tasks, found {len(task_roots)}"
        )
    rows = []
    for root in task_roots:
        metadata = json.loads((root / "slurm-task.json").read_text())
        summaries = list(root.glob("**/summary_r*.json"))
        if len(summaries) != 1:
            raise RuntimeError(
                f"Expected one paired summary under {root}, found {summaries}"
            )
        summary = json.loads(summaries[0].read_text())
        comparison = summary.get("semantic_comparison", {})
        if comparison.get("semantic_pass") is not True:
            raise RuntimeError(f"Task did not pass semantic comparison: {root}")
        uncached = float(summary["uncached_wall_seconds"])
        cached = float(summary["cached_wall_seconds"])
        rows.append(
            {
                "campaign_id": campaign["campaign_id"],
                "commit": summary["commit"],
                "phase": phase,
                "task_id": metadata["task_id"],
                "sample": summary["sample"],
                "assay": assay_label(metadata["input_vcf"]),
                "population": metadata["population"],
                "superpopulation": metadata["superpopulation"],
                "input_vcf": metadata["input_vcf"],
                "input_sha256": metadata["input_sha256"],
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


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    """Atomically write a nonempty publication TSV."""
    if not rows:
        raise ValueError("Publication table must not be empty")
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_suffix(path.suffix + ".partial")
    with partial.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    partial.replace(path)


def run(args: argparse.Namespace) -> Path:
    """Collect and write one campaign phase."""
    rows = collect_rows(args.campaign_root, args.phase)
    output = args.output or args.campaign_root / "publication/paired_metrics.tsv"
    write_tsv(output, rows)
    print(f"Collected {len(rows)} paired tasks -> {output}")
    return output


def parser() -> argparse.ArgumentParser:
    """Build the collector CLI."""
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--campaign-root", type=Path, required=True)
    result.add_argument("--phase", default="measured")
    result.add_argument("--output", type=Path)
    return result


def main() -> None:
    """Collect paired benchmark output."""
    run(parser().parse_args())


if __name__ == "__main__":
    main()
