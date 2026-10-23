#!/usr/bin/env python3
"""Enrich archived cells with observed backends, preserving superseded evidence."""

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

import yaml


def effective_backend(arm, stats):
    """Use observed stats, not an arm's requested engine, to label execution."""
    if arm == "direct":
        return "uncached"
    if arm == "original-fused":
        return "original_fused"
    actual = stats.get("effective_engine", stats.get("engine"))
    if actual == "fused":
        return "adaptive"
    if actual:
        return actual
    # Historical main and legacy stats predate backend instrumentation.
    return "legacy" if arm in {"legacy", "main-legacy"} else "unrecorded"


def collect(root):
    """Keep every timed row; mark supersession explicitly from the job manifest."""
    jobs = {
        str(job["id"]): job
        for job in json.loads((root / "RUN_MANIFEST.json").read_text())["jobs"]
    }
    rows = []
    for measurements in sorted(root.glob("*/benchmark/measurements.json")):
        condition = measurements.parent.parent
        fields = dict(
            line.split("=", 1)
            for line in (condition / "CONDITIONS.txt").read_text().splitlines()
            if "=" in line
        )
        job_id = fields["job"]
        status = jobs[job_id]["status"]
        primary = status == "complete"
        recipe = fields.get("annotator", fields.get("pipeline"))
        for cell in json.loads(measurements.read_text()):
            scale = Path(cell["input"]).stem
            directory = measurements.parent / f"{scale}-{cell['arm']}-{cell['rep']}"
            stats_files = list(directory.glob("*_vcstats/compare_stats.yaml"))
            if len(stats_files) > 1:
                raise RuntimeError(f"ambiguous stats: {directory}")
            stats = yaml.safe_load(stats_files[0].read_text()) if stats_files else {}
            comparison = json.loads((directory / "comparison.json").read_text())
            counts = stats.get("variant_counts", {})
            total = counts.get("input_variants")
            missing = counts.get("missing_variants")
            timings = stats.get("stage_timings", [])
            miss_seconds = (
                sum(
                    row["wall_seconds"]
                    for row in timings
                    if row["stage"] == "annotate-misses"
                )
                if timings
                else None
            )
            rows.append(
                dict(
                    recipe=recipe,
                    job=job_id,
                    job_status=status,
                    primary=primary,
                    family="heavy" if "pipeline" in fields else "matrix",
                    scale=scale,
                    requested_arm=cell["arm"],
                    repetition=cell["rep"],
                    requested_engine=stats.get("requested_engine", "legacy"),
                    effective_backend=effective_backend(cell["arm"], stats),
                    fallback_reason=stats.get("fallback_reason", ""),
                    wall_seconds=cell["wall_seconds"],
                    cpus=fields["cpus"],
                    records=total,
                    misses=missing,
                    hit_fraction=(
                        (1 - missing / total) if total and missing is not None else None
                    ),
                    miss_annotation_seconds=miss_seconds,
                    semantic_pass=cell["semantic_pass"],
                    audit_role=(
                        "baseline_self_check"
                        if cell["arm"] == "direct"
                        else "against_direct"
                    ),
                    mismatched_records=comparison.get("record_mismatches"),
                    comparison=str((directory / "comparison.json").relative_to(root)),
                    stats=str(stats_files[0].relative_to(root)) if stats_files else "",
                )
            )
    keys = [(r["job"], r["scale"], r["requested_arm"], r["repetition"]) for r in rows]
    if len(keys) != len(set(keys)):
        raise RuntimeError("duplicate cells: refusing to silently overwrite")
    return rows


def main():
    """Write a new review directory, never modify the original archived tables."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = collect(args.root)
    args.output.mkdir(parents=True, exist_ok=False)
    with (args.output / "cells_backend_audit.tsv").open("w") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    counts = Counter(
        (r["requested_arm"], r["effective_backend"]) for r in rows if r["primary"]
    )
    lines = [
        "# Backend-aware review",
        "",
        f"{len(rows)} archived timed rows retained; "
        f"{sum(r['primary'] and r['family'] == 'matrix' for r in rows)} primary matrix cells; "
        f"{sum(r['primary'] and r['family'] == 'heavy' for r in rows)} primary heavy cells; "
        f"{sum(not r['primary'] for r in rows)} superseded cells.",
        "",
        "Requested arm and observed backend are distinct. Historical arm `auto` "
        "means adaptive lookup, not the new availability-based CLI engine `auto`.",
        "",
        "| Requested arm | Effective backend | Primary cells |",
        "| --- | --- | ---: |",
    ]
    lines += [
        f"| {arm} | {backend} | {count} |"
        for (arm, backend), count in sorted(counts.items())
    ]
    lines += [
        "",
        "## Native requests that fell back",
        "",
        "| Recipe | Scale | Job | Observed backend | Reason |",
        "| --- | --- | --- | --- | --- |",
    ]
    lines += [
        f"| {r['recipe']} | {r['scale']} | {r['job']} | {r['effective_backend']} | {r['fallback_reason']} |"
        for r in rows
        if r["primary"]
        and r["requested_arm"] == "native"
        and r["effective_backend"] != "native"
    ]
    lines += [
        "",
        "Audits are unchanged. A baseline self-check is not independent "
        "validation. Clinical mismatches remain failures, not ignored fields. "
        "Small timing differences and slopes spanning fallback backends do "
        "not identify an engine advantage.",
        "",
    ]
    (args.output / "BACKEND_AUDIT.md").write_text("\n".join(lines))


if __name__ == "__main__":
    main()
