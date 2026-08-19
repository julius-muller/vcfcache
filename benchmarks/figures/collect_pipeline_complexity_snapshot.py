#!/usr/bin/env python3
"""Freeze the real-WGS virtual pipeline campaign with annotator anchors."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

EXPECTED_DELAYS_US = {500, 1_000, 2_000, 4_000, 7_000, 10_000}


def sha256sum(path: Path) -> str:
    """Return the SHA-256 digest of a file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    """Read one tab-separated table into string-valued rows."""
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    """Validate and freeze the controlled and annotator source tables."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wgs-spectrum-campaign-root", type=Path, required=True)
    parser.add_argument("--annotator-snapshot", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    campaign_root = args.wgs_spectrum_campaign_root.resolve()
    campaign_path = campaign_root / "campaign.json"
    metrics_path = campaign_root / "publication/wgs_pipeline_spectrum.tsv"
    annotator_path = (
        args.annotator_snapshot.resolve() / "annotator_external_wgs_metrics.tsv"
    )
    inputs = {
        "wgs_pipeline_campaign.json": campaign_path,
        "wgs_pipeline_spectrum.tsv": metrics_path,
        "annotator_external_wgs_metrics.tsv": annotator_path,
    }
    for path in inputs.values():
        if not path.is_file():
            raise FileNotFoundError(path)

    campaign = json.loads(campaign_path.read_text())
    spectrum = read_tsv(metrics_path)
    annotators = read_tsv(annotator_path)
    if campaign.get("statistics_mode") != "light":
        raise RuntimeError("WGS spectrum campaign did not use light statistics")
    if campaign.get("task_count") != 6 or campaign.get("technical_repeats") != 0:
        raise RuntimeError("WGS spectrum campaign is not the frozen six-load design")
    if len(spectrum) != 6:
        raise RuntimeError(f"Expected 6 WGS spectrum rows, found {len(spectrum)}")
    if {int(row["delay_us"]) for row in spectrum} != EXPECTED_DELAYS_US:
        raise RuntimeError("WGS virtual pipeline loads are incomplete")
    if any(row["sample"] != "KPGP-00319" for row in spectrum):
        raise RuntimeError("WGS spectrum does not use the frozen KPGP genome")
    if any(row["semantic_pass"].lower() != "true" for row in spectrum):
        raise RuntimeError("A WGS spectrum output failed semantic comparison")
    if any(row.get("statistics_mode") != "light" for row in spectrum):
        raise RuntimeError("A WGS spectrum timed cell did not use light statistics")
    if len(annotators) != 312 or {row["tool"] for row in annotators} != {
        "vep",
        "fastvep",
    }:
        raise RuntimeError("Matched VEP/fastVEP anchor table is incomplete")
    matched = [
        row
        for row in annotators
        if row["sample"] == "KPGP-00319"
        and row["strategy"] == "gnomad_af_0.1"
    ]
    if len(matched) != 2 or {row["tool"] for row in matched} != {"vep", "fastvep"}:
        raise RuntimeError("Matched fastVEP/VEP KPGP-00319 anchors are absent")

    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(output)
    output.mkdir(parents=True)
    hashes: dict[str, str] = {}
    provenance: dict[str, str] = {}
    for name, source in inputs.items():
        destination = output / name
        shutil.copy2(source, destination)
        hashes[name] = sha256sum(destination)
        provenance[name] = str(source)
    snapshot = {
        "status": "PIPELINE_COMPLEXITY_WGS_FINAL",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "campaign_id": campaign["campaign_id"],
        "campaign_commit": campaign["commit"],
        "statistics_mode": "light",
        "sample": "KPGP-00319",
        "input_sha256": campaign["input_sha256"],
        "wgs_spectrum_design": {
            "samples": 1,
            "virtual_pipeline_loads": 6,
            "delay_us_per_transcript_consequence": sorted(EXPECTED_DELAYS_US),
            "prior_observed_hit_rate": campaign["prior_observed_hit_rate"],
            "source_alias": campaign["source_alias"],
            "source_doi": campaign["source_doi"],
            "technical_repeats": 0,
            "timed_cells": 12,
            "validated_cached_outputs": 6,
        },
        "files": hashes,
        "provenance": provenance,
    }
    (output / "SNAPSHOT.json").write_text(
        json.dumps(snapshot, indent=2, sort_keys=True) + "\n"
    )
    print(output)


if __name__ == "__main__":
    main()
