#!/usr/bin/env python3
"""Collect the independent Panel/WES extension with matched frozen WGS rows."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_WGS = (
    REPO_ROOT
    / "benchmarks/figures/source_data/final"
    / "2026-08-09T223718Z-light-matched-final"
    / "annotator_external_wgs_metrics.tsv"
)
DEFAULT_OUTPUT = (
    REPO_ROOT
    / "benchmarks/external_assay_v3/source_data/external_assay_v3_metrics.tsv"
)
FIELDS = (
    "tool",
    "assay",
    "cohort",
    "sample",
    "assembly",
    "input_records",
    "cache_strategy",
    "cache_kind",
    "cache_hit_rate",
    "uncached_wall_seconds",
    "cached_wall_seconds",
    "relative_runtime",
    "speedup",
    "wall_seconds_saved",
    "semantic_pass",
    "semantic_comparator",
    "records_compared",
    "raw_annotation_mismatches",
    "ignored_annotation_mismatches",
    "source_campaign",
)


def sha256sum(path: Path) -> str:
    """Return a streaming SHA-256 digest."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def assay_rows(campaign_root: Path) -> list[dict[str, object]]:
    """Read and validate the 48 paired Panel/WES task summaries."""
    paths = sorted(campaign_root.glob("tasks/task-*/external_assay_summary.json"))
    rows = [json.loads(path.read_text()) for path in paths]
    if len(rows) != 48:
        raise RuntimeError(f"Expected 48 completed Panel/WES tasks, found {len(rows)}")
    task_ids = sorted(int(row["task_id"]) for row in rows)
    if task_ids != list(range(48)):
        raise RuntimeError("Panel/WES task IDs are incomplete or duplicated")
    if not all(row.get("semantic_pass") is True for row in rows):
        raise RuntimeError("At least one Panel/WES semantic comparison failed")
    return rows


def matched_wgs_rows(
    wgs_path: Path, samples: set[tuple[str, str]]
) -> list[dict[str, object]]:
    """Extract AF >=1% WGS observations for exactly the selected GRCh38 genomes."""
    with wgs_path.open(encoding="utf-8", newline="") as handle:
        candidates = list(csv.DictReader(handle, delimiter="\t"))
    rows: list[dict[str, object]] = []
    for row in candidates:
        key = (row["cohort"], row["sample"])
        if (
            key in samples
            and row["strategy"] == "gnomad_af_0.01"
            and row["assembly"] == "GRCh38"
            and row["tool"] in {"vep", "fastvep"}
        ):
            rows.append(
                {
                    "tool": row["tool"],
                    "assay": "wgs",
                    "cohort": row["cohort"],
                    "sample": row["sample"],
                    "assembly": row["assembly"],
                    "input_records": row["input_records"],
                    "cache_strategy": row["strategy"],
                    "cache_kind": row["strategy_kind"],
                    "cache_hit_rate": row["cache_hit_rate"],
                    "uncached_wall_seconds": row["uncached_wall_seconds"],
                    "cached_wall_seconds": row["cached_wall_seconds"],
                    "relative_runtime": row["relative_runtime"],
                    "speedup": row["speedup"],
                    "wall_seconds_saved": (
                        float(row["uncached_wall_seconds"])
                        - float(row["cached_wall_seconds"])
                    ),
                    "semantic_pass": row["semantic_pass"],
                    "semantic_comparator": row["semantic_comparator"],
                    "records_compared": row["input_records"],
                    "raw_annotation_mismatches": row.get(
                        "raw_annotation_mismatches", ""
                    ),
                    "ignored_annotation_mismatches": row.get(
                        "ignored_annotation_mismatches", ""
                    ),
                    "source_campaign": "external-wgs-light-matched-final",
                }
            )
    if len(rows) != 24:
        raise RuntimeError(f"Expected 24 matched WGS rows, found {len(rows)}")
    return rows


def normalize_assay(row: dict[str, object]) -> dict[str, object]:
    """Select stable publication columns from one Panel/WES task summary."""
    result = {key: row.get(key, "") for key in FIELDS}
    result["source_campaign"] = "external-assay-v3-independent-af001-warm"
    return result


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    """Write the combined source table atomically with Unix line endings."""
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_suffix(path.suffix + ".partial")
    with partial.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    partial.replace(path)


def collect(args: argparse.Namespace) -> None:
    """Combine and validate a complete matched assay-by-annotator table."""
    assay = assay_rows(args.campaign_root)
    samples = {(str(row["cohort"]), str(row["sample"])) for row in assay}
    if len(samples) != 12:
        raise RuntimeError(f"Expected 12 matched genomes, found {len(samples)}")
    combined = [normalize_assay(row) for row in assay]
    combined.extend(matched_wgs_rows(args.wgs_metrics, samples))
    combined.sort(
        key=lambda row: (
            str(row["tool"]), str(row["assay"]),
            str(row["cohort"]), str(row["sample"])
        )
    )
    expected = {
        (tool, assay_name, cohort, sample)
        for cohort, sample in samples
        for tool in ("vep", "fastvep")
        for assay_name in ("panel", "wes", "wgs")
    }
    observed = {
        (str(row["tool"]), str(row["assay"]), str(row["cohort"]), str(row["sample"]))
        for row in combined
    }
    if len(combined) != 72 or observed != expected:
        raise RuntimeError("The 2 tools x 3 assays x 12 genomes matrix is incomplete")
    if not all(str(row["semantic_pass"]).lower() == "true" for row in combined):
        raise RuntimeError("At least one combined semantic gate is not true")
    write_tsv(args.output, combined)
    manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "campaign_root": str(args.campaign_root),
        "wgs_metrics": str(args.wgs_metrics),
        "wgs_metrics_sha256": sha256sum(args.wgs_metrics),
        "output": str(args.output),
        "output_sha256": sha256sum(args.output),
        "rows": len(combined),
        "samples": len(samples),
        "tools": 2,
        "assays": 3,
        "technical_repeats": 0,
        "all_semantic_pass": True,
    }
    manifest_path = args.output.with_suffix(".manifest.json")
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


def parser() -> argparse.ArgumentParser:
    """Build the command-line parser."""
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--campaign-root", type=Path, required=True)
    result.add_argument("--wgs-metrics", type=Path, default=DEFAULT_WGS)
    result.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return result


if __name__ == "__main__":
    collect(parser().parse_args())
