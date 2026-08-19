#!/usr/bin/env python3
"""Freeze matched AF >=10% and AF >=1% real-WGS pipeline spectra."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

EXPECTED_DELAYS_US = {500, 1_000, 2_000, 4_000, 7_000, 10_000}
AF1_ALIAS = "cache-gnomad-v4.1-GRCh38-joint-af001-vep115.2-e"
AF1_DOI = "10.5281/zenodo.18190046"


def sha256sum(path: Path) -> str:
    """Return the SHA-256 digest of a file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    """Read a tab-separated table."""
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def validate_spectrum(rows: list[dict[str, str]], *, reused: bool) -> None:
    """Validate one complete six-load spectrum."""
    if len(rows) != 6:
        raise RuntimeError(f"Expected six spectrum rows, found {len(rows)}")
    if {int(row["delay_us"]) for row in rows} != EXPECTED_DELAYS_US:
        raise RuntimeError("Pipeline loads are incomplete")
    if any(row["sample"] != "KPGP-00319" for row in rows):
        raise RuntimeError("Pipeline spectrum uses an unexpected sample")
    if any(row["statistics_mode"] != "light" for row in rows):
        raise RuntimeError("A timed cell did not use light statistics")
    if any(row["semantic_pass"].lower() != "true" for row in rows):
        raise RuntimeError("A cached output failed semantic validation")
    if reused and any(row["baseline_reused"].lower() != "true" for row in rows):
        raise RuntimeError("An AF >=1% row did not reuse the frozen direct baseline")


def main() -> None:
    """Validate and freeze the matched two-cache source tables."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--af10-snapshot", type=Path, required=True)
    parser.add_argument("--af1-campaign-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    af10_root = args.af10_snapshot.resolve()
    af1_root = args.af1_campaign_root.resolve()
    af10_snapshot_path = af10_root / "SNAPSHOT.json"
    af10_metrics_path = af10_root / "wgs_pipeline_spectrum.tsv"
    annotator_path = af10_root / "annotator_external_wgs_metrics.tsv"
    af1_campaign_path = af1_root / "campaign.json"
    af1_metrics_path = af1_root / "publication/wgs_pipeline_spectrum_af001.tsv"
    af1_validation_path = af1_root / "publication/semantic_validation.json"
    inputs = {
        "af10_snapshot.json": af10_snapshot_path,
        "wgs_pipeline_spectrum_af10.tsv": af10_metrics_path,
        "annotator_external_wgs_metrics.tsv": annotator_path,
        "af1_campaign.json": af1_campaign_path,
        "wgs_pipeline_spectrum_af1.tsv": af1_metrics_path,
        "af1_semantic_validation.json": af1_validation_path,
    }
    for path in inputs.values():
        if not path.is_file():
            raise FileNotFoundError(path)

    af10_snapshot = json.loads(af10_snapshot_path.read_text())
    af1_campaign = json.loads(af1_campaign_path.read_text())
    validation = json.loads(af1_validation_path.read_text())
    af10 = read_tsv(af10_metrics_path)
    af1 = read_tsv(af1_metrics_path)
    annotators = read_tsv(annotator_path)
    if af10_snapshot.get("status") != "PIPELINE_COMPLEXITY_WGS_FINAL":
        raise RuntimeError("AF >=10% source snapshot is not final")
    validate_spectrum(af10, reused=False)
    validate_spectrum(af1, reused=True)
    if validation.get("semantic_pass") is not True:
        raise RuntimeError("AF >=1% composed semantic validation did not pass")
    if validation.get("new_cached_outputs") != 6:
        raise RuntimeError("AF >=1% semantic validation is incomplete")
    if (
        af1_campaign.get("source_alias") != AF1_ALIAS
        or af1_campaign.get("source_doi") != AF1_DOI
        or af1_campaign.get("new_timed_cells") != 6
        or af1_campaign.get("reused_direct_cells") != 6
        or af1_campaign.get("technical_repeats") != 0
    ):
        raise RuntimeError("AF >=1% campaign metadata is not the frozen design")
    by_delay_af10 = {int(row["delay_us"]): row for row in af10}
    for row in af1:
        baseline = by_delay_af10[int(row["delay_us"])]
        if abs(
            float(row["uncached_wall_seconds"])
            - float(baseline["uncached_wall_seconds"])
        ) > 1e-9:
            raise RuntimeError("AF >=1% row does not use the matched direct baseline")
        if row["baseline_campaign_id"] != af10_snapshot["campaign_id"]:
            raise RuntimeError("AF >=1% baseline campaign ID is inconsistent")
    if len(annotators) != 312 or {row["tool"] for row in annotators} != {
        "vep",
        "fastvep",
    }:
        raise RuntimeError("Matched annotator anchor table is incomplete")

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
        "status": "PIPELINE_COMPLEXITY_DUAL_CACHE_FINAL",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "sample": "KPGP-00319",
        "input_sha256": af10_snapshot["input_sha256"],
        "statistics_mode": "light",
        "technical_repeats": 0,
        "direct_baselines": 6,
        "cached_measurements": 12,
        "semantic_passes": 12,
        "cache_strategies": [
            {
                "label": "Bundled gnomAD any-stratum AF >=10%",
                "alias": af10_snapshot["wgs_spectrum_design"]["source_alias"],
                "doi": af10_snapshot["wgs_spectrum_design"]["source_doi"],
                "hit_rate": float(af10[0]["cache_hit_rate"]),
                "measurement_kind": "paired_cached_and_direct",
            },
            {
                "label": "Bundled gnomAD any-stratum AF >=1%",
                "alias": af1_campaign["source_alias"],
                "doi": af1_campaign["source_doi"],
                "hit_rate": float(af1[0]["cache_hit_rate"]),
                "measurement_kind": "new_cached_with_frozen_matched_direct",
            },
        ],
        "af10_campaign_id": af10_snapshot["campaign_id"],
        "af1_campaign_id": af1_campaign["campaign_id"],
        "af1_semantic_comparator": validation["comparator"],
        "files": hashes,
        "provenance": provenance,
    }
    (output / "SNAPSHOT.json").write_text(
        json.dumps(snapshot, indent=2, sort_keys=True) + "\n"
    )
    print(output)


if __name__ == "__main__":
    main()
