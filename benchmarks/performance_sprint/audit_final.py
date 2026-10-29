"""Derive local provenance and validation audits from the frozen evidence only.

No cluster access, annotation, cache mutation or mismatch exemptions are added.
Build costs are accepted only when provenance binds the executed payload hash.
"""

import argparse
import csv
import hashlib
import json
import math
import statistics
import tarfile
from collections import Counter
from pathlib import Path

import yaml

from benchmarks.performance_sprint.snapshot_final import table


def read_rows(path):
    """Read a tabular export without coercing identifiers or missing values."""
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def build_threshold(build_seconds, direct_seconds, cached_seconds):
    """Minimum integer cohort size giving strictly positive total time savings."""
    saving = direct_seconds - cached_seconds
    return math.floor(build_seconds / saving) + 1 if saving > 0 else None


def seed_audit(strategies, rows):
    """Held-out means disjoint from this cache's seeds, not from gnomAD donors."""
    result = []
    for cohort, strategy in strategies["cohort_strategies"].items():
        targets = {
            r["sample"]
            for r in rows
            if r["family"] == "full52_wgs_strategies" and r["cohort"] == cohort
        }
        seeds = set(strategy["training_samples"])
        if not targets or targets & seeds:
            raise ValueError("Empty targets or cohort-cache training leakage")
        result.append(
            dict(
                cohort=cohort,
                assembly=strategy["assembly"],
                training_samples=";".join(sorted(seeds)),
                n_seeds=len(seeds),
                n_test_genomes=len(targets),
                n_overlapping_ids=len(targets & seeds),
                definition="sample-ID disjoint from three cache seeds; not a gnomAD donor audit",
            )
        )
    return result


def audit(root):
    """Write payload-bound provenance and comparator audits for a final snapshot."""
    rows = read_rows(root / "tables/cells.tsv")
    views = read_rows(root / "tables/plot_views.tsv")
    snapshot = json.loads((root / "source/SNAPSHOT.json").read_text())
    archive_path = root / "source/terminal-evidence.tar.gz"
    archive_sha256 = hashlib.sha256(archive_path.read_bytes()).hexdigest()
    if archive_sha256 != snapshot["archive_sha256"]:
        raise ValueError("Frozen evidence archive SHA256 differs from snapshot")
    with tarfile.open(archive_path) as archive:

        def data(name):
            return json.load(archive.extractfile(name))

        main = data("main/manifest.json")
        repair = data("repair/manifest.json")
        metadata = data("cache_metadata/INDEX.json")
        provenance = {}
        recipes_by_sha256 = {}
        for entry in metadata:
            content = archive.extractfile(entry["archive_path"]).read()
            if hashlib.sha256(content).hexdigest() != entry["sha256"]:
                raise ValueError("Archived metadata hash differs")
            if entry["path"].endswith("/annotation.yaml"):
                recipes_by_sha256[entry["sha256"]] = yaml.safe_load(content)
            if entry["path"].endswith("provenance.json"):
                value = json.loads(content)
                if value.get("cache_sha256"):
                    provenance[value["cache_sha256"]] = (entry, value)

        strategies = [data(f"cache_metadata/strategies-{i}.json") for i in range(2)]
        vep = next(s for s in strategies if s["tool"] == "vep")
        seeds = seed_audit(vep, rows)
        table(root / "tables/cohort_seed_audit.tsv", seeds)
        attestations = main["freeze"]["asset_attestations"] | repair["freeze"].get(
            "asset_attestations", {}
        )
        inventories = {}
        comparisons = []
        for row in rows:
            prefix = row["source_archive_member"]
            command = data(prefix + "/command.json")
            task = command["task"]
            comparison = data(prefix + "/comparison.json")
            cache = task["cache"]
            payload = cache["files"]["vcfcache_annotated.bcf"]
            recipe_hash = cache["files"]["annotation.yaml"]
            key = payload, recipe_hash
            if key not in inventories:
                original = cache.get("original_directory", cache["directory"])
                canonical = cache.get("canonical_directory", cache["directory"])
                sizes = [
                    a["bytes"]
                    for path, a in attestations.items()
                    if path
                    in {d + "/vcfcache_annotated.bcf" for d in (original, canonical)}
                    and a.get("sha256") == payload
                ]
                recipe = recipes_by_sha256[recipe_hash]
                entry, prov = provenance.get(payload, ({}, {}))
                inventories[key] = dict(
                    tool=row["tool"],
                    level=row["level"],
                    assembly=row["assembly"],
                    strategy=row["strategy"],
                    payload_sha256=payload,
                    recipe_sha256=recipe_hash,
                    cache_bytes=sizes[0] if sizes else None,
                    original_directory=original,
                    annotation_command=recipe.get("annotation_cmd", "")
                    .strip()
                    .replace("\n", " "),
                    required_tool_version=recipe.get("required_tool_version"),
                    provenance_archive_member=entry.get("archive_path"),
                    payload_bound_build_seconds=prov.get("wall_seconds"),
                    build_scope=(
                        "historical local build; not a new matched build or download measurement"
                        if prov.get("wall_seconds") is not None
                        else "not measured locally for this payload"
                    ),
                    provenance_kind=prov.get("kind", prov.get("source", "")),
                    n_canonical_conditions=0,
                )
            inventories[key]["n_canonical_conditions"] += 1
            comparisons.append(
                dict(
                    condition_id=row["condition_id"],
                    tool=row["tool"],
                    level=row["level"],
                    assay=row["assay"],
                    sample=row["sample"],
                    evidence_class=row["evidence_class"],
                    scope=row["validation_scope"],
                    record_mismatches=row["mismatched_records"],
                    raw_vep_direct_differences=row["vep_raw_annotation_differences"],
                    headers_equal=row["headers_equal"],
                    non_csq_pass=comparison.get("non_csq_comparison", {}).get(
                        "semantic_pass"
                    ),
                    ignored_fields=json.dumps(comparison.get("ignored_fields", [])),
                    direct_comparator=comparison.get("against_direct", {}).get(
                        "comparator", comparison.get("comparator", "")
                    ),
                    direct_ignored_csq_fields=json.dumps(
                        comparison.get("against_direct", {}).get(
                            "ignored_csq_fields", []
                        )
                    ),
                    direct_unordered_csq_fields=json.dumps(
                        comparison.get("against_direct", {}).get(
                            "unordered_csq_fields", []
                        )
                    ),
                    cached_reference_strict_pass=comparison.get(
                        "against_cached", {}
                    ).get("semantic_pass"),
                    documented_cached_semantic_exception_pass=comparison.get(
                        "documented_vep_cached", {}
                    ).get("semantic_pass"),
                    example_loci=json.dumps(
                        [e.get("locus") for e in comparison.get("examples", [])]
                    ),
                    source_archive_member=prefix + "/comparison.json",
                )
            )

        table(root / "tables/cache_inventory.tsv", list(inventories.values()))
        table(root / "tables/comparator_audit.tsv", comparisons)
        amortization = []
        for strategy_set in strategies:
            for cohort, local_strategy in strategy_set["cohort_strategies"].items():
                selected = [
                    local_strategy,
                    *strategy_set["bundled_strategies_by_assembly"][
                        local_strategy["assembly"]
                    ],
                ]
                for strategy in selected:
                    conditions = [
                        r
                        for r in views
                        if r["family"] == "full52_wgs_strategies"
                        and r["tool"] == strategy_set["tool"]
                        and r["cohort"] == cohort
                        and r["strategy"] == strategy["name"]
                    ]
                    hashes = {r["cache_payload_sha256"] for r in conditions}
                    if len(hashes) != 1:
                        raise ValueError(
                            "Strategy view is not one frozen cache payload"
                        )
                    payload = next(iter(hashes))
                    entry, prov = provenance.get(payload, ({}, {}))
                    build = prov.get("wall_seconds")
                    direct = statistics.median(
                        float(r["direct_seconds"]) for r in conditions
                    )
                    cached = statistics.median(
                        float(r["cached_seconds"]) for r in conditions
                    )
                    saving = statistics.median(
                        float(r["saved_seconds"]) for r in conditions
                    )
                    amortization.append(
                        dict(
                            tool=strategy_set["tool"],
                            cohort=cohort,
                            assembly=strategy["assembly"],
                            strategy=strategy["name"],
                            kind=strategy["kind"],
                            source_doi=strategy.get(
                                "source_doi", strategy.get("doi", "")
                            ),
                            source_alias=strategy.get(
                                "source_alias", strategy.get("alias", "")
                            ),
                            payload_sha256=payload,
                            payload_bound_historical_build_seconds=build,
                            strategy_manifest_build_seconds=strategy.get(
                                "build_wall_seconds"
                            ),
                            recorded_cost_difference_seconds=(
                                build - strategy["build_wall_seconds"]
                                if build is not None
                                and "build_wall_seconds" in strategy
                                else None
                            ),
                            n_test_genomes=len(conditions),
                            median_direct_seconds=direct,
                            median_cached_seconds=cached,
                            median_paired_saved_seconds=saving,
                            minimum_samples_to_repay_recorded_build=(
                                build_threshold(build, saving, 0)
                                if build is not None
                                else None
                            ),
                            provenance_archive_member=entry.get("archive_path"),
                            cost_scope=(
                                "historical annotation build from available blueprint; download and blueprint preparation excluded"
                                if strategy_set["tool"] == "fastvep"
                                else (
                                    "historical cohort-cache build; not matched to current runtime hardware/settings"
                                    if build is not None
                                    else "ready-made Zenodo cache; local build not required, download time not measured"
                                )
                            ),
                        )
                    )
        table(root / "tables/recorded_build_amortization.tsv", amortization)

    summary = dict(
        canonical_conditions=len(rows),
        plotting_views=len(views),
        distinct_genomes=len({r["sample"] for r in rows}),
        accepted_under_declared_comparators=sum(r["validated"] == "True" for r in rows),
        evidence_classes=dict(Counter(r["evidence_class"] for r in rows)),
        comparison_scopes=dict(Counter(r["validation_scope"] for r in rows)),
        baseline_origins=dict(Counter(r["direct_baseline_origin"] for r in rows)),
        executed_builds=dict(Counter(r["build_commit"] for r in rows)),
        cache_payload_recipe_combinations=len(inventories),
        payload_bound_build_costs=len(provenance),
        terminal_archive_sha256=archive_sha256,
        snapshot_created_utc=snapshot["created_utc"],
        cohort_seed_id_audit_passed=True,
        gnomad_donor_disjointness="not established",
    )
    (root / "tables/FINAL_AUDIT.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("snapshot_directory", type=Path)
    args = parser.parse_args()
    print(json.dumps(audit(args.snapshot_directory), indent=2))
