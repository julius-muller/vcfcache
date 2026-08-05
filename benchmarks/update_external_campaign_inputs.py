#!/usr/bin/env python3
"""Apply a verified prepared-input repair to a frozen external campaign."""

from __future__ import annotations

import argparse
import csv
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from benchmarks.run_pilot import git_output, sha256sum, write_json_atomic


def _write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    partial = path.with_name(path.name + ".partial")
    with partial.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    partial.replace(path)


def update_campaign(
    campaign_root: Path, repair_report: Path, qc_path: Path
) -> dict[str, Any]:
    """Update affected task rows and campaign provenance atomically."""
    report = json.loads(repair_report.read_text())
    repairs = {item["sample"]: item for item in report["repairs"]}
    campaign_path = campaign_root / "campaign.json"
    campaign = json.loads(campaign_path.read_text())
    phase_updates: dict[str, list[int]] = {}
    for phase, metadata in campaign["phases"].items():
        manifest = campaign_root / "manifests" / f"{phase}.tsv"
        with manifest.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = list(reader.fieldnames or ())
            rows = list(reader)
        changed = []
        for row in rows:
            repair = repairs.get(row["sample"])
            if repair is None:
                continue
            row.update(
                {
                    "input_vcf": repair["destination"],
                    "input_records": str(repair["repaired"]["records"]),
                    "input_sha256": repair["destination_sha256"],
                }
            )
            changed.append(int(row["task_id"]))
        if changed:
            backup = manifest.with_name(manifest.name + ".before_duplicate_repair")
            if not backup.exists():
                shutil.copy2(manifest, backup)
            _write_tsv(manifest, rows, fields)
        metadata["manifest"] = str(manifest)
        metadata["manifest_sha256"] = sha256sum(manifest)
        phase_updates[phase] = changed

    qc_destination = campaign_root / "manifests/external_wgs_qc.tsv"
    report_destination = campaign_root / "manifests/external_wgs_duplicate_repair.json"
    shutil.copy2(qc_path, qc_destination)
    shutil.copy2(repair_report, report_destination)
    campaign["qc_source"] = str(qc_path)
    campaign["qc"] = str(qc_destination)
    campaign["qc_sha256"] = sha256sum(qc_destination)
    campaign["input_repair"] = {
        "applied_at": datetime.now(timezone.utc).isoformat(),
        "applied_with_commit": git_output("rev-parse", "HEAD"),
        "method": report["method"],
        "report": str(report_destination),
        "report_sha256": sha256sum(report_destination),
        "phase_task_ids": phase_updates,
        "cache_rebuild_required": report["cache_rebuild_required"],
        "cache_rebuild_rationale": report["cache_rebuild_rationale"],
    }
    write_json_atomic(campaign_path, campaign)
    return campaign["input_repair"]


def parser() -> argparse.ArgumentParser:
    """Build the command-line parser."""
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--campaign-root", required=True, type=Path)
    result.add_argument("--repair-report", required=True, type=Path)
    result.add_argument("--qc", required=True, type=Path)
    return result


def main() -> None:
    """Apply the repair and print the campaign provenance block."""
    args = parser().parse_args()
    print(
        json.dumps(
            update_campaign(
                args.campaign_root,
                args.repair_report,
                args.qc,
            ),
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
