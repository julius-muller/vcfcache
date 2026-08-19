#!/usr/bin/env python3
"""Assemble the complete VEP figure snapshot with light-statistics WGS data."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def read_json(path: Path) -> dict[str, Any]:
    """Read one JSON object."""
    return json.loads(path.read_text())


def read_tsv(path: Path) -> list[dict[str, str]]:
    """Read a TSV file."""
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    """Return a file SHA-256 digest."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def copy_file(source: Path, target: Path) -> None:
    """Copy one frozen source file."""
    if not source.is_file():
        raise RuntimeError(f"Missing source file: {source}")
    shutil.copyfile(source, target)


def run(args: argparse.Namespace) -> None:
    """Create a complete immutable source bundle for the R renderer."""
    base_manifest = read_json(args.base / "SNAPSHOT.json")
    annotator_manifest = read_json(args.annotator / "SNAPSHOT.json")
    if base_manifest.get("status") != "FINAL":
        raise RuntimeError("The base VEP snapshot is not final")
    if annotator_manifest.get("status") != "ANNOTATOR_LIGHT_MATCHED_FINAL":
        raise RuntimeError("The annotator snapshot is not light-mode matched and final")

    light_source = args.annotator / "vep_external_wgs_light_statistics_metrics.tsv"
    light_rows = read_tsv(light_source)
    if len(light_rows) != 156:
        raise RuntimeError(f"Expected 156 VEP light rows, found {len(light_rows)}")
    if len({row["sample"] for row in light_rows}) != 52:
        raise RuntimeError("Expected 52 unique VEP light samples")
    if any(row["statistics_mode"] != "light" for row in light_rows):
        raise RuntimeError("The VEP replacement contains non-light rows")
    if any(row["validation_status"] != "semantically_validated" for row in light_rows):
        raise RuntimeError("The VEP replacement contains unvalidated rows")

    args.output.mkdir(parents=True, exist_ok=False)
    files = {
        "primary_wgs_metrics.tsv": args.base / "primary_wgs_metrics.tsv",
        "assay_metrics.tsv": args.base / "assay_metrics.tsv",
        "external_wgs_metrics.tsv": light_source,
        "external_strategies.json": args.base / "external_strategies.json",
    }
    for name, source in files.items():
        copy_file(source, args.output / name)

    manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "created_on_host": platform.node(),
        "status": "FINAL",
        "campaigns": {
            "primary": base_manifest["campaigns"]["primary"],
            "assay": base_manifest["campaigns"]["assay"],
            "external": annotator_manifest["vep_light_full"]["campaign_id"],
        },
        "row_counts": {
            "primary_wgs_metrics.tsv": len(
                read_tsv(args.output / "primary_wgs_metrics.tsv")
            ),
            "assay_metrics.tsv": len(read_tsv(args.output / "assay_metrics.tsv")),
            "external_wgs_metrics.tsv": len(light_rows),
        },
        "sample_counts": {
            **base_manifest["sample_counts"],
            "external_completed": 52,
            "external_semantically_validated": 52,
            "external_expected": 52,
        },
        "external_validation_status": {"semantically_validated": 156},
        "statistics_mode": {"external_wgs": "light"},
        "source_snapshots": {
            "base_vep_assay": str(args.base),
            "annotator_light_matched": str(args.annotator),
        },
        "notes": [
            "Primary WGS and assay tables are byte-identical to the original final VEP snapshot.",
            "External WGS contains all 52 VEP --statistics light reruns and 156 semantically validated cache-strategy rows.",
            "The historical full-rescan external WGS measurements remain frozen in the base snapshot and are not copied into this bundle.",
            "No BCF or other bulky run output is included in this snapshot.",
        ],
    }
    manifest["files"] = {
        name: {
            "bytes": (args.output / name).stat().st_size,
            "sha256": sha256(args.output / name),
        }
        for name in files
    }
    (args.output / "SNAPSHOT.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


def parser() -> argparse.ArgumentParser:
    """Build the CLI parser."""
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--base", type=Path, required=True)
    result.add_argument("--annotator", type=Path, required=True)
    result.add_argument("--output", type=Path, required=True)
    return result


if __name__ == "__main__":
    run(parser().parse_args())
