"""Expose frozen Rust-adaptive timings to the same renderer, without relabelling."""

import csv
import json
import math
import tarfile
from pathlib import Path

from benchmarks.python_adaptive_results.collect import table
from benchmarks.run_pilot import sha256sum, write_json_atomic


def main():
    """Audit the original archive and map units/labels, never measurement values."""
    root = Path(__file__).parent / "final-20261103"
    original = root / "source/original"
    snapshot = json.loads((original / "source/SNAPSHOT.json").read_text())
    archive = original / "source/terminal-evidence.tar.gz"
    assert sha256sum(archive) == snapshot["archive_sha256"]
    python_root = Path(__file__).parents[1] / "python_adaptive_results/final-20261103"
    with (python_root / "tables/cells.tsv").open() as stream:
        python = {r["condition_id"]: r for r in csv.DictReader(stream, delimiter="\t")}

    def convert(row, evidence):
        # The old exported rows identify exactly which archived attempt was used.
        prefix = row["source_archive_member"]
        metrics = json.load(evidence.extractfile(prefix + "/resources.json"))
        assert math.isclose(
            metrics["wall_seconds"], float(row["cached_seconds"]), abs_tol=1e-8
        )
        assert (
            metrics["exit_code"] == 0
            and len(metrics["affinity"]) == int(row["cpus"]) == 8
        )
        assert row["engine_requested"] == row["engine_effective"] == "adaptive"
        assert row["helper_version"].startswith("vcfcache-engine ")
        assert row["input_records"] == row["output_records"]
        if row["condition_id"] == row["measurement_id"]:
            p = python[row["condition_id"]]
            assert math.isclose(
                float(row["cached_seconds"]),
                float(p["previous_rust_cached_seconds"]),
                abs_tol=1e-8,
            )
            assert math.isclose(
                float(row["direct_seconds"]), float(p["direct_seconds"]), abs_tol=1e-8
            )
            assert row["cache_payload_sha256"] == p["cache_sha256"]
            assert row["input_sha256"] == p["input_sha256"]
        row = dict(row)
        row["cached_peak_gib"] = float(row["cached_peak_bytes"]) / 1024**3
        row["direct_peak_gib"] = (
            float(row["direct_peak_bytes"]) / 1024**3
            if row["direct_peak_bytes"]
            else None
        )
        row["evidence_class"] = {
            "validated": "validated_with_recorded_scope",
            "veprs_csq_timing_only": "veprs_clinical_timing_only",
            "snpeff_warning_failed_strict": "snpeff_warning_excluded",
        }[row["evidence_class"]]
        # No previous-Rust self-comparison is manufactured for this arm.
        row["previous_rust_cached_seconds"] = None
        return row

    with tarfile.open(archive) as evidence:
        for filename in ["cells.tsv", "plot_views.tsv"]:
            with (original / "tables" / filename).open() as stream:
                rows = [
                    convert(r, evidence) for r in csv.DictReader(stream, delimiter="\t")
                ]
            assert len(rows) == (600 if filename == "cells.tsv" else 612)
            table(root / "tables" / filename, rows)
    write_json_atomic(
        root / "source/SNAPSHOT.json",
        dict(
            implementation="Rust-adaptive (compiled fused router with adaptive access); NOT resident/native",
            original_snapshot=snapshot,
            archive_sha256=sha256sum(archive),
            notes="Historical measured campaign; no jobs rerun. All 600 canonical timings match the Rust controls used by the Python replay, with identical direct denominators, input and cache hashes. Unit/label mapping only.",
        ),
    )
    print(
        "PASS: 600 Rust-adaptive measurements and 612 views bound to original archive and matching Python replay controls"
    )


if __name__ == "__main__":
    main()
