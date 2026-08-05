#!/usr/bin/env python3
"""Repair legacy prepared external-WGS inputs containing duplicate allele keys."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

from benchmarks.run_pilot import sha256sum, write_json_atomic

DEFAULT_SAMPLES = (
    "hu365511",
    "hu6664F3",
    "hu9AF7CC",
    "huC8266A",
    "huDEFDD1",
    "huCDD5EE",
    "huCF305F",
)


def variant_key_summary(path: Path) -> dict[str, Any]:
    """Stream normalized keys once and return content/QC statistics."""
    command = [
        "bcftools",
        "query",
        "--format",
        r"%CHROM\t%POS\t%REF\t%ALT\t%TYPE\n",
        str(path),
    ]
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if process.stdout is None:
        raise RuntimeError("bcftools query did not expose stdout")
    digest = hashlib.sha256()
    previous = ""
    records = duplicates = snps = indels = multiallelic = 0
    contigs: set[str] = set()
    for line in process.stdout:
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 5:
            process.kill()
            raise RuntimeError(f"Malformed bcftools query output for {path}")
        chrom, pos, ref, alt, variant_type = fields
        key = "\t".join((chrom, pos, ref, alt))
        records += 1
        duplicates += key == previous
        if key != previous:
            digest.update(key.encode())
            digest.update(b"\n")
        previous = key
        contigs.add(chrom)
        multiallelic += "," in alt
        snps += variant_type == "SNP"
        indels += variant_type == "INDEL"
    stderr = process.stderr.read() if process.stderr is not None else ""
    status = process.wait()
    if status:
        raise RuntimeError(
            f"bcftools query failed for {path} with {status}: {stderr.strip()}"
        )
    return {
        "records": records,
        "duplicate_keys": duplicates,
        "unique_key_sha256": digest.hexdigest(),
        "snps": snps,
        "indels": indels,
        "multiallelic": multiallelic,
        "contigs": sorted(contigs),
    }


def repair_vcf(source: Path, destination: Path) -> dict[str, Any]:
    """Deduplicate one VCF atomically and prove its unique-key set is unchanged."""
    original = variant_key_summary(source)
    if original["duplicate_keys"] < 1:
        raise RuntimeError(f"Expected duplicate keys in repair target: {source}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_name(destination.name + ".partial")
    partial_index = Path(f"{partial}.tbi")
    partial.unlink(missing_ok=True)
    partial_index.unlink(missing_ok=True)
    try:
        subprocess.run(
            [
                "bcftools",
                "norm",
                "--rm-dup",
                "exact",
                "--output-type",
                "z",
                "--output",
                str(partial),
                str(source),
            ],
            check=True,
        )
        subprocess.run(
            ["bcftools", "index", "--tbi", "--force", str(partial)], check=True
        )
        repaired = variant_key_summary(partial)
        if repaired["duplicate_keys"]:
            raise RuntimeError(f"Duplicate keys remain after repair: {source}")
        if repaired["unique_key_sha256"] != original["unique_key_sha256"]:
            raise RuntimeError(f"Unique variant-key set changed during repair: {source}")
        if original["records"] - repaired["records"] != original["duplicate_keys"]:
            raise RuntimeError(f"Unexpected record loss during repair: {source}")
        partial.replace(destination)
        partial_index.replace(Path(f"{destination}.tbi"))
    except BaseException:
        partial.unlink(missing_ok=True)
        partial_index.unlink(missing_ok=True)
        raise
    return {
        "source": str(source),
        "source_sha256": sha256sum(source),
        "source_bytes": source.stat().st_size,
        "destination": str(destination),
        "destination_sha256": sha256sum(destination),
        "destination_bytes": destination.stat().st_size,
        "original": original,
        "repaired": repaired,
    }


def _write_qc(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    partial = path.with_name(path.name + ".partial")
    with partial.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    partial.replace(path)


def repair(root: Path, samples: Iterable[str]) -> dict[str, Any]:
    """Repair selected PGP inputs and atomically update the preparation QC."""
    qc_path = root / "qc/external_wgs_qc.tsv"
    with qc_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or ())
        rows = list(reader)
    if "duplicate_keys" not in fields:
        fields.append("duplicate_keys")
    requested = set(samples)
    selected = {row["sample"]: row for row in rows if row["sample"] in requested}
    missing = requested - selected.keys()
    if missing:
        raise RuntimeError(f"Repair samples absent from QC: {sorted(missing)}")
    repairs = []
    for sample in sorted(requested):
        row = selected[sample]
        if row["cohort"] != "pgp" or row["assembly"] != "GRCh37":
            raise RuntimeError(f"Unexpected repair target: {row}")
        source = Path(row["path"])
        destination = (
            root
            / "samples/GRCh37/pgp_deduplicated"
            / f"{sample}.GRCh37.small_variants.vcf.gz"
        )
        if destination.exists():
            destination.unlink()
            Path(f"{destination}.tbi").unlink(missing_ok=True)
        result = repair_vcf(source, destination)
        repaired = result["repaired"]
        row.update(
            {
                "path": str(destination),
                "records": str(repaired["records"]),
                "snps": str(repaired["snps"]),
                "indels": str(repaired["indels"]),
                "multiallelic": str(repaired["multiallelic"]),
                "contigs": ",".join(repaired["contigs"]),
                "bytes": str(result["destination_bytes"]),
                "sha256": result["destination_sha256"],
                "status": "PASS",
                "errors": "",
                "duplicate_keys": "0",
            }
        )
        result.update({"sample": sample, "role": row["role"]})
        repairs.append(result)
    for row in rows:
        row.setdefault("duplicate_keys", "0")
    backup = qc_path.with_name(qc_path.name + ".before_duplicate_repair")
    if not backup.exists():
        shutil.copy2(qc_path, backup)
    _write_qc(qc_path, rows, fields)
    report = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "method": "bcftools norm --rm-dup exact",
        "policy": "original prepared inputs retained; repaired paths are new files",
        "cache_rebuild_required": False,
        "cache_rebuild_rationale": (
            "The affected training input has the identical unique variant-key digest "
            "before and after repair; cache membership is unchanged."
        ),
        "qc": str(qc_path),
        "qc_sha256": sha256sum(qc_path),
        "repairs": repairs,
    }
    report_path = root / "manifests/external_wgs_duplicate_repair.json"
    write_json_atomic(report_path, report)
    return report


def parser() -> argparse.ArgumentParser:
    """Build the command-line parser."""
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--root", type=Path, required=True)
    result.add_argument("--samples", nargs="+", default=list(DEFAULT_SAMPLES))
    return result


def main() -> None:
    """Repair requested inputs and print the provenance report."""
    print(json.dumps(repair(**vars(parser().parse_args())), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
