"""Inventory historical manuscript cells without silently approving baseline reuse."""

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

from benchmarks.run_pilot import REPO_ROOT, sha256sum, write_json_atomic

NATIVE = Path(
    "benchmarks/figures/output/final/2026-08-25T224623Z-native-annotation-final-v3/source"
)
EXTERNAL = Path(
    "benchmarks/figures/source_data/final/2026-08-09T223718Z-light-matched-final/annotator_external_wgs_metrics.tsv"
)


def read_tsv(path):
    """Read an explicit source table with preserved textual provenance."""
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def inventory(repo):
    """Retain all experiment families and known exclusions; never join unlike recipes."""
    records = []

    def add(
        family, path, row, index, tool, level, assay, strategy, direct, cached, accepted
    ):
        record = dict(
            family=family,
            source=str(path),
            source_sha256=sha256sum(repo / path),
            source_row=index,
            sample=row["sample"],
            cohort=row["cohort"],
            assembly=row["assembly"],
            tool=tool,
            level=level,
            assay=assay,
            strategy=strategy,
            input_records=row.get("input_records"),
            historical_direct_seconds=float(direct),
            historical_cached_seconds=float(cached),
            historical_acceptance=accepted,
            reference_class="historical_timing",
            reuse_decision="pending_provenance_audit",
            direct_rerun_authorized=False,
        )
        identity = {
            k: record[k]
            for k in [
                "family",
                "source",
                "source_row",
                "sample",
                "tool",
                "level",
                "assay",
                "strategy",
            ]
        }
        record["id"] = hashlib.sha256(
            json.dumps(identity, sort_keys=True).encode()
        ).hexdigest()[:16]
        records.append(record)

    path = NATIVE / "campaign_metrics.tsv"
    native = read_tsv(repo / path)
    direct = {
        (r["sample"], r["tool"], r["level"], r["assembly"], r["assay"]): r
        for r in native
        if r["condition"] == "direct"
    }
    for index, row in enumerate(native, 2):
        if row["condition"] != "cached":
            continue
        d = direct[
            (row["sample"], row["tool"], row["level"], row["assembly"], row["assay"])
        ]
        add(
            "native_clinical_assays",
            path,
            row,
            index,
            row["tool"],
            row["level"],
            row["assay"],
            "gnomad_af_0.01",
            d["wall_seconds"],
            row["wall_seconds"],
            row["concordant"],
        )
        records[-1]["historical_direct_digest"] = d["digest"]
        records[-1]["output_evidence"] = (
            "annotator-specific historical digest; full direct output retention must be checked"
        )
    path = NATIVE / "external_assay_v3_metrics.tsv"
    for index, row in enumerate(read_tsv(repo / path), 2):
        if row["tool"].lower() != "vep":
            continue
        add(
            "vep_matched_assays",
            path,
            row,
            index,
            "vep",
            "everything",
            row["assay"],
            row["cache_strategy"],
            row["uncached_wall_seconds"],
            row["cached_wall_seconds"],
            row["semantic_pass"],
        )
    for index, row in enumerate(read_tsv(repo / EXTERNAL), 2):
        add(
            "full52_wgs_strategies",
            EXTERNAL,
            row,
            index,
            row["tool"].lower(),
            "original_external_recipe",
            "wgs",
            row["strategy"],
            row["uncached_wall_seconds"],
            row["cached_wall_seconds"],
            row["semantic_pass"],
        )
    path = NATIVE / "dual_cache_pipeline_spectrum.tsv"
    for index, row in enumerate(read_tsv(repo / path), 2):
        add(
            "synthetic_pipeline_spectrum",
            path,
            row,
            index,
            "vep",
            row["pipeline"],
            "wgs",
            row["cache_strategy"],
            row["uncached_wall_seconds"],
            row["cached_wall_seconds"],
            row["semantic_pass"],
        )
    assert len({r["id"] for r in records}) == len(records)
    return records


def reference_hashes(repo):
    """Bind shared original genomes, recipe snapshots and established GRCh38 regions."""
    expected = {}
    tasks = (
        repo
        / "benchmarks/results/raw/campaigns/external-vep-light-full52-bbef41d7db8e/manifests/measured.tsv"
    )
    for row in read_tsv(tasks):
        expected[row["input_vcf"]] = row["input_sha256"]
    snapshots = repo / "benchmarks/reproduction/final_clinical_campaign/cache_recipes"
    for snapshot in sorted(snapshots.glob("*.annotation.yaml")):
        tool, level, assembly = snapshot.name.split(".", 1)[0].split("_")
        name = f"{tool}_{level}_{assembly.lower()}_af001"
        prefix = (
            "/results/richness/caches"
            if tool == "fastvep"
            else "/results/snpeff/caches"
        )
        root = f"{prefix}/{name}/cache/{name}"
        expected[root + "/annotation.yaml"] = sha256sum(snapshot)
        expected[root + "/params.snapshot.yaml"] = sha256sum(
            snapshots / f"{tool}_{level}_{assembly}.params.yaml"
        )
    for row in read_tsv(repo / "benchmarks/manifests/assay_regions.tsv"):
        path = row["bed"]
        digest = row.get("sha256")
        if path and digest:
            expected["/results/richness/regions/" + Path(path).name] = digest
    return expected


def main():
    """Write a source-fingerprinted planning inventory, not a submitted campaign."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    args.output.mkdir(parents=True)
    records = inventory(REPO_ROOT)
    expected = reference_hashes(REPO_ROOT)
    (args.output / "inputs_and_recipes.sha256").write_text(
        "".join(f"{digest}  {path}\n" for path, digest in sorted(expected.items()))
    )
    fields = list(dict.fromkeys(k for r in records for k in r))
    with (args.output / "baseline_reuse_ledger.tsv").open("w") as stream:
        writer = csv.DictWriter(stream, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)
    summary = dict(
        status="planning_inventory_not_submission",
        families=dict(Counter(r["family"] for r in records)),
        rows=len(records),
        cached_node_hours=sum(r["historical_cached_seconds"] for r in records) / 3600,
        note="Possible cross-family duplicates are not deduplicated until full command/input fingerprints agree; direct columns repeat across strategies and must not be summed as unique work.",
        approved_direct_reuses=0,
        rerun_direct_by_default=False,
        veprs_extension="Pinned core and clinical, same twelve GRCh38 assay genomes; prepared after candidate freeze, no upstream repair.",
        deadline="2026-11-02T18:00:00+01:00",
    )
    write_json_atomic(args.output / "INVENTORY.json", summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
