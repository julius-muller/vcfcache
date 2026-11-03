"""Read-only terminal evidence collection for the Python-adaptive manuscript replay."""

import argparse
import csv
import json
import subprocess
import tarfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from benchmarks.performance_sprint.repair_headers import direct_binding
from benchmarks.run_pilot import sha256sum, write_json_atomic

HOST = "ubuntu@10.133.241.128"
REMOTE = (
    "/mnt/data/slurm-results/python-adaptive-manuscript-20261101-ee2162473/campaign"
)
REQUEST = "/mnt/data/slurm-results/deployment/manuscript-final-20261027T191111Z-fe197f34d/benchmarks/performance_sprint/reference-request-v2-20261027.json"


def read(path):
    """Read a frozen JSON document."""
    return json.loads(path.read_text())


def table(path, rows):
    """Export missing values as missing, never as invented zero measurements."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w") as stream:
        writer = csv.DictWriter(stream, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def cpu(metrics):
    """Prefer the aggregate independent Slurm-step CPU delta."""
    before = metrics.get("cgroup_before", {}).get("cpu_stat", {}).get("usage_usec")
    after = metrics.get("cgroup_after", {}).get("cpu_stat", {}).get("usage_usec")
    return (after - before) / 1e6 if before is not None and after is not None else None


def flatten(task, directory, manifest):
    """Validate bindings, terminal state, counts and comparison before exporting."""
    state = read(directory / "state.json")
    command = read(directory / "command.json")
    metrics = read(directory / "resources.json")
    stats = read(directory / "statistics.json")
    compare = read(directory / "comparison.json")
    assert direct_binding(task) == direct_binding(command["task"])
    assert task["cache"]["files"] == command["task"]["cache"]["files"]
    assert state["source_manifest_sha256"] == sha256sum(manifest)
    assert task["engine"] == "adaptive" and state["build_commit"] == "ee2162473"
    assert metrics["exit_code"] == 0 and len(metrics["affinity"]) == task["cpus"] == 8
    assert metrics["cgroup_after"] and metrics["wall_seconds"] > 0
    counts = stats["variant_counts"]
    assert (
        counts["input_variants"] == counts["total_output"] == task["input"]["records"]
    )
    assert counts["dropped_variants"] == 0
    assert 0 <= counts["missing_variants"] <= counts["input_variants"]
    status = state["state"]
    if status == "complete_with_direct_difference":
        assert (task["tool"], task["level"]) == ("veprs", "clinical")
        assert (
            compare["non_csq_comparison"]["semantic_pass"]
            and compare["relevant_headers_equal"]
        )
    elif status == "failed":
        assert state["error"] == "new strict paired output differs"
        assert task["tool"] == "snpeff" and task["sample"] == "hu6664F3"
        assert compare["record_mismatches"] == 1 and compare["relevant_headers_equal"]
        assert compare["examples"][0]["locus"] == ["chr2", "11273459"]
    else:
        assert status == "complete"
        if task["baseline"]["kind"] == "new_paired_reference":
            assert compare["semantic_pass"]
        else:
            assert compare.get("against_direct", {}).get("semantic_pass") or task[
                "baseline"
            ].get("inherited_direct_difference")
            if task["baseline"]["kind"] == "composed_vep_reference":
                assert compare["composed_cached_reference"]["semantic_pass"]
            else:
                assert compare["against_cached"]["semantic_pass"] or (
                    compare.get("documented_vep_cached", {}).get("semantic_pass")
                    and compare.get("non_csq_cached", {}).get("semantic_pass")
                    and compare["against_cached"]["relevant_headers_equal"]
                )
    direct = (
        read(directory / "direct/resources.json")
        if (directory / "direct/resources.json").exists()
        else {}
    )
    if direct:
        assert direct["exit_code"] == 0 and len(direct["affinity"]) == task["cpus"]
        assert state.get("reused_direct"), "No new direct runs were authorized"
    seconds = metrics["wall_seconds"]
    baseline = state.get("direct_seconds", task["baseline"].get("direct_seconds"))
    assert baseline > 0
    return dict(
        condition_id=task["id"],
        measurement_id=task["id"],
        is_alias=False,
        **{
            key: task[key]
            for key in [
                "family",
                "tool",
                "level",
                "assay",
                "cohort",
                "assembly",
                "sample",
                "strategy",
            ]
        },
        state=status,
        validated=status == "complete",
        evidence_class={
            "complete": "validated_with_recorded_scope",
            "complete_with_direct_difference": "veprs_clinical_timing_only",
            "failed": "snpeff_warning_excluded",
        }[status],
        validation_scope=compare["validation_scope"],
        input_records=counts["input_variants"],
        output_records=counts["total_output"],
        missing_variants=counts["missing_variants"],
        cache_hit_rate=1 - counts["missing_variants"] / counts["input_variants"],
        direct_seconds=baseline,
        cached_seconds=seconds,
        speedup=baseline / seconds,
        saved_seconds=baseline - seconds,
        saved_fraction=1 - seconds / baseline,
        samples_per_day_direct=86400 / baseline,
        samples_per_day_cached=86400 / seconds,
        previous_rust_cached_seconds=task["historical_cached_seconds"],
        previous_rust_to_python_ratio=task["historical_cached_seconds"] / seconds,
        direct_kind=task["baseline"]["kind"],
        direct_reused=bool(direct),
        cpus=task["cpus"],
        cached_cpu_seconds=cpu(metrics),
        direct_cpu_seconds=cpu(direct),
        cached_peak_gib=metrics["cgroup_after"]["memory_peak_bytes"] / 1024**3,
        direct_peak_gib=(
            direct.get("cgroup_after", {}).get("memory_peak_bytes", 0) / 1024**3
            if direct
            else None
        ),
        job_id=state["job_id"],
        node=state["node"],
        build_commit=state["build_commit"],
        input_sha256=task["input"]["sha256"],
        cache_sha256=task["cache"]["files"]["vcfcache_annotated.bcf"],
        recipe_sha256=task["cache"]["files"]["annotation.yaml"],
        raw_record_mismatches=compare.get("record_mismatches"),
        raw_vep_annotation_differences=compare.get("against_direct", {}).get(
            "raw_annotation_mismatches"
        ),
        source_evidence=str(directory.relative_to(manifest.parent)),
    )


def collect(output):
    """Download bounded metadata only; retain immutable raw evidence in an archive."""
    output.mkdir(parents=True, exist_ok=False)
    source = output / "source"
    source.mkdir()
    raw = source / "campaign"
    filters = [
        "--exclude=cli/",
        "--exclude=logs/",
        "--exclude=resources.jsonl",
        "--include=*/",
    ]
    filters += [
        "--include=" + name
        for name in [
            "manifest.json",
            "READY",
            "submissions.json",
            "START.json",
            "AUTO_DISPATCH.json",
            "state.json",
            "command.json",
            "comparison.json",
            "statistics.json",
            "resources.json",
            "REUSED.json",
            "review-*.json",
        ]
    ]
    subprocess.run(
        ["rsync", "-a", *filters, "--exclude=*", f"{HOST}:{REMOTE}/", str(raw) + "/"],
        check=True,
    )
    subprocess.run(
        ["rsync", "-a", f"{HOST}:{REQUEST}", str(source / "reference-request.json")],
        check=True,
    )
    manifest = raw / "manifest.json"
    assert sha256sum(manifest) == (raw / "READY").read_text().strip()
    data = read(manifest)
    latest = {}
    attempts = []
    for path in sorted((raw / "attempts").glob("*/*/state.json")):
        latest[path.parent.name] = path.parent
        attempts.append(
            dict(
                condition_id=path.parent.name,
                attempt=path.parent.parent.name,
                **read(path),
            )
        )
    assert len(latest) == len(data["tasks"]) == 600
    rows = [flatten(t, latest[t["id"]], manifest) for t in data["tasks"]]
    assert Counter(r["state"] for r in rows) == {
        "complete": 572,
        "complete_with_direct_difference": 24,
        "failed": 4,
    }
    request = read(source / "reference-request.json")
    expected = {
        x["sha256"]
        for t in data["tasks"]
        for x in t["baseline"]["approval_evidence"]
        if x["path"].endswith("reference-request-v2-20261027.json")
    }
    assert expected == {sha256sum(source / "reference-request.json")}
    request_rows = {r["id"]: r for r in request["rows"]}
    lookup = {r["condition_id"]: r for r in rows}
    views = list(rows)
    for alias, canonical in data["freeze"]["aliases"].items():
        row = dict(lookup[canonical])
        row.update(
            {
                k: request_rows[alias][k]
                for k in [
                    "family",
                    "tool",
                    "level",
                    "assay",
                    "cohort",
                    "assembly",
                    "sample",
                    "strategy",
                ]
            }
        )
        row.update(condition_id=alias, is_alias=True)
        views.append(row)
    assert len(views) == 612 and sum(r["is_alias"] for r in views) == 12
    assert len([r for r in views if r["family"] == "full52_wgs_strategies"]) == 312
    table(output / "tables/cells.tsv", rows)
    table(output / "tables/plot_views.tsv", views)
    table(output / "tables/exclusions.tsv", [r for r in rows if not r["validated"]])
    # Flatten only scalar attempt metadata; complete documents remain archived.
    table(
        output / "tables/attempts.tsv",
        [
            {
                k: r.get(k)
                for k in [
                    "condition_id",
                    "attempt",
                    "state",
                    "job_id",
                    "node",
                    "started_unix",
                    "finished_unix",
                    "error",
                ]
            }
            for r in attempts
        ],
    )
    archive = source / "terminal-evidence.tar.gz"
    with tarfile.open(archive, "w:gz") as stream:
        stream.add(raw, arcname="campaign")
        stream.add(source / "reference-request.json", arcname="reference-request.json")
    write_json_atomic(
        source / "SNAPSHOT.json",
        dict(
            collected_utc=datetime.now(timezone.utc).isoformat(),
            remote=REMOTE,
            manifest_sha256=sha256sum(manifest),
            archive_sha256=sha256sum(archive),
            states=dict(Counter(r["state"] for r in rows)),
            measurements=600,
            plot_views=612,
            reused_paired_direct_controls=sum(r["direct_reused"] for r in rows),
            new_direct_runs=0,
            product_commit=data["freeze"]["build"]["commit"],
            implementation="Python/bcftools adaptive; no compiled helper",
            notes="Aliases are plotting views, not extra measurements. Missing phase timings are not inferred. Previous cached timings are Rust-adaptive historical controls, not main-branch comparisons.",
        ),
    )
    print(json.dumps(read(source / "SNAPSHOT.json"), indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    collect(parser.parse_args().output)
