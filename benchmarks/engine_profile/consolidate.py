"""Reproduce backend-aware consolidation without rerunning annotation or audits."""

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

from benchmarks.engine_profile.review import union_seconds
from benchmarks.engine_profile.run import verify
from benchmarks.run_pilot import sha256sum, write_json_atomic


def write_table(path, rows):
    """Emit stable, rectangular TSVs with explicit NA values."""
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w") as stream:
        writer = csv.DictWriter(stream, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {key: "NA" if row.get(key) is None else row.get(key) for key in fields}
            )


def consolidate(root):
    """Keep original/recovery configurations distinct and use latest explicit attempts."""
    rows, attempts, phase_rows, provenance = [], [], [], []
    for campaign in ("original", "recovery"):
        source = root / "source" / campaign
        manifest = verify(source)
        if campaign == "recovery" and manifest["recovery"][
            "parent_manifest_sha256"
        ] != sha256sum(root / "source/original/manifest.json"):
            raise ValueError("recovery does not belong to the frozen original campaign")
        for file in sorted(source.rglob("*")):
            if file.is_file():
                provenance.append(
                    dict(path=str(file.relative_to(root)), sha256=sha256sum(file))
                )
        for group in manifest["groups"]:
            for cell in group["cells"]:
                files = sorted((source / "attempts").glob(f"*/{cell['id']}/state.json"))
                for file in files:
                    state = json.loads(file.read_text())
                    attempts.append(
                        dict(
                            campaign=campaign,
                            cell=cell["id"],
                            attempt=file.parent.parent.name,
                            state=state["state"],
                            error=state.get("error"),
                            latest=file == files[-1],
                        )
                    )
                file = files[-1] if files else None
                state = json.loads(file.read_text()) if file else {}
                evidence = file.parent if file else None
                stats_file = evidence / "statistics.json" if evidence else None
                stats = (
                    json.loads(stats_file.read_text())
                    if stats_file and stats_file.exists()
                    else {}
                )
                metrics = state.get("metrics", {})
                if not metrics and evidence and (evidence / "resources.json").exists():
                    metrics = json.loads((evidence / "resources.json").read_text())
                counts = stats.get("variant_counts", {})
                n, missing = counts.get("input_variants"), counts.get(
                    "missing_variants"
                )
                events = (
                    [
                        json.loads(line)
                        for trace in sorted(evidence.glob("trace/*.jsonl"))
                        for line in trace.read_text().splitlines()
                    ]
                    if evidence
                    else []
                )
                spans = defaultdict(list)
                for event in events:
                    start = event["start_unix_ns"]
                    spans[event["name"]].append((start, start + event["duration_ns"]))
                phases = {
                    name: union_seconds(intervals) for name, intervals in spans.items()
                }
                for name, duration in phases.items():
                    phase_rows.append(
                        dict(
                            campaign=campaign,
                            cell=cell["id"],
                            phase=name,
                            union_seconds=duration,
                        )
                    )
                annotation = (
                    sum(
                        x["wall_seconds"]
                        for x in stats.get("stage_timings", [])
                        if x["stage"] == "annotate-misses"
                    )
                    or None
                )
                wall = metrics.get("wall_seconds")
                completed = state.get("state") in {
                    "complete",
                    "complete_with_direct_difference",
                }
                comparison_file = evidence / "against_direct.json" if evidence else None
                comparison = (
                    json.loads(comparison_file.read_text())
                    if comparison_file and comparison_file.exists()
                    else {}
                )
                row = dict(
                    campaign=campaign,
                    cell=cell["id"],
                    group=group["id"],
                    sample=cell["sample"],
                    recipe=cell["recipe"],
                    assay=cell["assay"],
                    role=cell["role"],
                    purpose=cell["purpose"],
                    state=state.get("state", "missing"),
                    attempt=file.parent.parent.name if file else None,
                    requested=cell["engine"],
                    effective=state.get("effective_engine"),
                    cpus=cell["cpus"],
                    configured_workers=cell["workers"],
                    executed_workers=stats.get("workers"),
                    memory_mb=cell["memory_mb"],
                    profile=cell["profile"],
                    form=cell["form"],
                    skip_split=cell["skip_split"],
                    wall_seconds=wall,
                    cpu_seconds=metrics.get("cpu_seconds"),
                    peak_gib=(
                        (metrics.get("cgroup_after") or {}).get("memory_peak_bytes", 0)
                        / 1024**3
                        if metrics
                        else None
                    ),
                    records=n,
                    misses=missing,
                    hit_fraction=1 - missing / n if n and missing is not None else None,
                    annotation_seconds=annotation,
                    cache_side_residual_seconds=(
                        wall - annotation if wall and annotation else None
                    ),
                    annotation_fraction=(
                        annotation / wall if wall and annotation else None
                    ),
                    ideal_zero_cache_overhead_gain=(
                        wall / annotation if completed and wall and annotation else None
                    ),
                    direct_pass=state.get("direct_pass"),
                    instrumentation_equivalent=state.get("instrumentation_equivalent"),
                    record_mismatches=comparison.get("record_mismatches"),
                    fallback=stats.get("fallback_reason"),
                    candidate_records=stats.get("candidate_records"),
                    lookup_queries=stats.get("lookup_queries"),
                    buffer_bytes=stats.get("buffer_bytes"),
                    native_attempt_seconds=phases.get("native.execute"),
                    vep_forks=(
                        cell.get("annotation_overrides", {}).get(
                            "vep_forks", cell["cpus"]
                        )
                        if cell["recipe"] == "vep"
                        else None
                    ),
                    vep_buffer=cell.get("annotation_overrides", {}).get(
                        "vep_buffer",
                        manifest["recipes"][cell["recipe"]]["params"].get("vep_buffer"),
                    ),
                )
                for name in (
                    "pipeline.route",
                    "pipeline.assemble",
                    "pipeline.preprocess",
                    "pipeline.index",
                    "native.lookup",
                    "native.construct_encode_output",
                ):
                    row[name] = phases.get(name)
                rows.append(row)
    references = {
        (r["campaign"], r["sample"], r["recipe"], r["assay"]): r
        for r in rows
        if r["role"] == "direct" and r["state"] == "complete"
    }
    matched = []
    for row in rows:
        reference = references.get(
            (row["campaign"], row["sample"], row["recipe"], row["assay"])
        )
        if (
            reference
            and row["role"] != "direct"
            and row["state"] in {"complete", "complete_with_direct_difference"}
        ):
            row["direct_seconds"] = reference["wall_seconds"]
            row["speedup"] = reference["wall_seconds"] / row["wall_seconds"]
            row["time_saved_seconds"] = reference["wall_seconds"] - row["wall_seconds"]
            if row["purpose"] in {"primary", "reference"}:
                matched.append(row)
    pairs = []
    calibration = defaultdict(dict)
    for row in rows:
        if row["purpose"].startswith("instrumentation-pair-"):
            calibration[row["purpose"]][row["profile"]] = row
    for pair, values in calibration.items():
        on, off = values[True], values[False]
        pairs.append(
            dict(
                pair=pair,
                on_seconds=on["wall_seconds"],
                off_seconds=off["wall_seconds"],
                relative_delta=on["wall_seconds"] / off["wall_seconds"] - 1,
            )
        )
    write_table(root / "cells.tsv", rows)
    recovery = {r["cell"]: r for r in rows if r["campaign"] == "recovery"}
    originals = {r["cell"]: r for r in rows if r["campaign"] == "original"}
    for identifier, replacement in recovery.items():
        parent = originals[identifier]
        if any(
            parent[key] != replacement[key]
            for key in ("sample", "recipe", "assay", "role", "requested")
        ):
            raise ValueError(f"incompatible recovery cell: {identifier}")
    canonical = [
        recovery.get(r["cell"], r) for r in rows if r["campaign"] == "original"
    ]
    write_table(root / "canonical_conditions.tsv", canonical)
    write_table(root / "attempts.tsv", attempts)
    write_table(root / "matched.tsv", matched)
    write_table(root / "phases.tsv", phase_rows)
    write_table(root / "instrumentation.tsv", pairs)
    targets = []
    for row in matched:
        if row["purpose"] != "primary" or row["annotation_seconds"] is None:
            continue
        for target in (1, 2, 3):
            allowance = row["direct_seconds"] / target - row["annotation_seconds"]
            residual = row["cache_side_residual_seconds"]
            targets.append(
                dict(
                    campaign=row["campaign"],
                    cell=row["cell"],
                    recipe=row["recipe"],
                    assay=row["assay"],
                    effective=row["effective"],
                    target_speedup=target,
                    allowed_cache_side_seconds=allowance,
                    required_cache_side_reduction_seconds=max(0, residual - allowance),
                    required_cache_side_reduction_fraction=max(0, residual - allowance)
                    / residual,
                    possible_with_unchanged_annotation=allowance >= 0,
                    assumption="observed miss stage unchanged; ignores measurement uncertainty",
                )
            )
    write_table(root / "targets.tsv", targets)
    write_json_atomic(
        root / "SOURCE_MANIFEST.json",
        dict(
            files=provenance,
            digest=hashlib.sha256(
                json.dumps(provenance, sort_keys=True).encode()
            ).hexdigest(),
            notes="Latest explicit attempts; no failed timings treated as successful; separate original/recovery baselines; no pooled annotator rankings.",
        ),
    )
    return rows


def main():
    """Create reproducible tables inside a new consolidation directory."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    consolidate(args.directory)


if __name__ == "__main__":
    main()
