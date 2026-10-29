"""Archive terminal manuscript evidence and export auditable, nonduplicated tables.

Read-only against the cluster. No annotation, retries, scheduling or cache writes.
Dense sampler streams and large BCFs are deliberately not downloaded.
"""

import argparse
import csv
import json
import shlex
import subprocess
import tarfile
import tempfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from benchmarks.engine_profile.deploy import HOST
from benchmarks.performance_sprint.repair_headers import direct_binding
from benchmarks.run_pilot import sha256sum, write_json_atomic

MAIN = "/results/manuscript-final/20261027T191111Z-fe197f34d/campaign"
REPAIR = "/results/manuscript-final/20261028T210624Z-787a3081a/campaign"


def load(path):
    """Read one frozen JSON evidence file."""
    return json.loads(path.read_text())


def table(path, rows):
    """Export a typed-by-column TSV; missing evidence stays empty, not zero."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with path.open("w") as stream:
        writer = csv.DictWriter(stream, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: "" if v is None else v for k, v in row.items()})


def copy_remote(source, destination, filters=()):
    """Copy only explicitly selected evidence, preserving all remote originals."""
    mapped = source.replace("/results/", "/mnt/data/slurm-results/", 1)
    subprocess.run(
        [
            "rsync",
            "-a",
            "--copy-links",
            "--rsh=ssh -o LogLevel=ERROR",
            "--rsync-path=sudo -u appuser rsync",
            *filters,
            f"{HOST}:{mapped}",
            str(destination),
        ],
        check=True,
    )


def controller_paths(paths):
    """Resolve worker-prefix symlinks in the controller's shared mount view."""
    script = """import json,os,sys
from pathlib import Path
def mapped(value): return Path(str(value).replace("/results/","/mnt/data/slurm-results/",1))
result={}
for original in json.load(sys.stdin):
 p=mapped(original)
 for _ in range(12):
  if not p.is_symlink(): break
  link=Path(os.readlink(p));p=mapped(link) if link.is_absolute() else p.parent/link
 if not p.is_file(): raise RuntimeError("unavailable metadata: "+original+" -> "+str(p))
 result[original]=str(p)
print(json.dumps(result))"""
    result = subprocess.run(
        ["ssh", "-o", "LogLevel=ERROR", HOST, "python3 -c " + shlex.quote(script)],
        input=json.dumps(paths),
        text=True,
        capture_output=True,
        check=True,
    )
    return json.loads(result.stdout)


def cpu(metrics):
    """Use aggregate step CPU delta rather than one child's RSS/CPU."""
    before = metrics.get("cgroup_before", {}).get("cpu_stat", {}).get("usage_usec")
    after = metrics.get("cgroup_after", {}).get("cpu_stat", {}).get("usage_usec")
    return (
        (after - before) / 1e6
        if before is not None and after is not None
        else metrics.get("cpu_seconds")
    )


def peak(metrics):
    """Peak is valid here only because each condition used a fresh Slurm step."""
    return metrics.get("cgroup_after", {}).get("memory_peak_bytes")


def reference_metrics(directory, state):
    """Prefer the original explicit direct measurement, including reused controls."""
    file = directory / "direct/resources.json"
    return load(file) if file.exists() else state.get("direct_metrics", {})


def audit_cell(task, directory, state, primary, repair, carried_builds):
    """Fail closed on wrong attempts, broken metrics, CPUs or validation evidence."""
    metrics = load(directory / "resources.json")
    stats = load(directory / "statistics.json")
    compare = load(directory / "comparison.json")
    command = load(directory / "command.json")
    actual_task = command["task"]
    if direct_binding(task) != direct_binding(actual_task):
        raise ValueError(f"condition binding differs: {task['id']}")
    if task["cache"]["files"] != actual_task["cache"]["files"]:
        raise ValueError("cache recipe/payload/index changed")
    if metrics["exit_code"] or metrics["wall_seconds"] <= 0:
        raise ValueError("no successful cached timing")
    if len(metrics["affinity"]) != task["cpus"] or not metrics.get("cgroup_after"):
        raise ValueError("CPU/isolation metrics missing")
    counts = stats["variant_counts"]
    if (
        counts["input_variants"] != task["input"]["records"]
        or counts["total_output"] != task["input"]["records"]
    ):
        raise ValueError("input/output counts differ")
    if state["state"] not in {"complete", "complete_with_direct_difference", "failed"}:
        raise ValueError("campaign is not terminal")
    if (
        state["state"] == "failed"
        and state.get("error") != "new strict paired output differs"
    ):
        raise ValueError(f"unexpected operational failure: {state.get('error')}")
    direct = reference_metrics(directory, state)
    if direct and (
        direct.get("exit_code") != 0 or len(direct.get("affinity", [])) != task["cpus"]
    ):
        raise ValueError("direct resource metrics invalid")
    baseline_kind = actual_task["baseline"]["kind"]
    if state["state"] == "complete_with_direct_difference":
        if (
            task["tool"] != "veprs"
            or task["level"] != "clinical"
            or not compare["non_csq_comparison"]["semantic_pass"]
            or not compare["relevant_headers_equal"]
        ):
            raise ValueError("unapproved timing-only exception")
    elif state["state"] == "failed":
        if (
            task["tool"] != "snpeff"
            or task["baseline"].get("historical_acceptance") != "FAIL_CONTENT"
            or compare["record_mismatches"] != 1
            or not compare["relevant_headers_equal"]
        ):
            raise ValueError("unexpected strict comparison failure")
    elif baseline_kind == "new_paired_reference":
        if not compare["semantic_pass"]:
            raise ValueError("paired validation did not pass")
    else:
        if not compare.get("against_direct", {}).get(
            "semantic_pass"
        ) and not actual_task["baseline"].get("inherited_direct_difference"):
            raise ValueError("direct comparator did not pass")
        if baseline_kind == "composed_vep_reference":
            if not compare["composed_cached_reference"]["semantic_pass"]:
                raise ValueError("composed reference did not pass")
        elif not compare["against_cached"]["semantic_pass"]:
            if not (
                compare.get("documented_vep_cached", {}).get("semantic_pass")
                and compare.get("non_csq_cached", {}).get("semantic_pass")
                and compare["against_cached"]["relevant_headers_equal"]
            ):
                raise ValueError(
                    "cached reference differs beyond documented equivalence"
                )
    build = state.get("build_commit") or (
        primary["freeze"]["build"]["commit"]
        if int(directory.parent.name[1:]) < 3
        else repair["freeze"]["build"]["commit"]
    )
    if state.get("carried_forward"):
        build = carried_builds[state["carried_forward"]["campaign"]]
    return metrics, stats, compare, direct, actual_task, build


def flatten(task, directory, state, primary, repair, carried_builds=None):
    """One row per executed condition, with absolute savings and phase headroom."""
    metrics, stats, compare, direct, actual_task, build = audit_cell(
        task, directory, state, primary, repair, carried_builds or {}
    )
    measured = metrics["wall_seconds"]
    direct_seconds = state.get("direct_seconds") or task["baseline"].get(
        "direct_seconds"
    )
    if not direct_seconds or direct_seconds <= 0:
        raise ValueError("missing direct denominator")
    engine = stats["engine_wall_seconds"]
    phases = {p["phase"]: p["wall_seconds"] for p in stats.get("phase_timings", [])}
    annotate = sum(
        p["wall_seconds"]
        for p in stats.get("stage_timings", [])
        if p["stage"] == "annotate-misses"
    )
    residual = measured - annotate
    counts = stats["variant_counts"]
    n, misses = counts["input_variants"], counts["missing_variants"]
    f = 1 - misses / n
    quality = (
        "validated"
        if state["state"] == "complete"
        else (
            "veprs_csq_timing_only"
            if state["state"] == "complete_with_direct_difference"
            else "snpeff_warning_failed_strict"
        )
    )
    scope = compare["validation_scope"]
    raw_vep = compare.get("against_direct", {}).get("raw_annotation_mismatches", 0)
    reused = bool(state.get("reused_direct"))
    carried = bool(state.get("carried_forward"))
    row = dict(
        condition_id=task["id"],
        measurement_id=task["id"],
        is_alias=False,
        family=task["family"],
        tool=task["tool"],
        level=task["level"],
        assay=task["assay"],
        cohort=task["cohort"],
        assembly=task["assembly"],
        sample=task["sample"],
        strategy=task["strategy"],
        state=state["state"],
        evidence_class=quality,
        validated=quality == "validated",
        validation_scope=scope,
        input_records=n,
        missing_variants=misses,
        output_records=counts["total_output"],
        cache_hit_rate=f,
        direct_seconds=direct_seconds,
        cached_seconds=measured,
        speedup=direct_seconds / measured,
        saved_seconds=direct_seconds - measured,
        saved_fraction=1 - measured / direct_seconds,
        samples_per_day_direct=86400 / direct_seconds,
        samples_per_day_cached=86400 / measured,
        historical_cached_seconds=task["historical_cached_seconds"],
        historical_to_new_speedup=task["historical_cached_seconds"] / measured,
        historical_cached_reduction=1 - measured / task["historical_cached_seconds"],
        direct_kind=actual_task["baseline"]["kind"],
        direct_reused_from_failed_attempt=reused,
        direct_baseline_origin=(
            "reused_completed_new_control"
            if reused
            else ("fresh_paired" if direct else "historical_reuse")
        ),
        build_commit=build,
        engine_requested=task["engine"],
        engine_effective=stats.get("effective_engine"),
        helper_version=stats.get("native_version"),
        cpus=task["cpus"],
        cached_cpu_seconds=cpu(metrics),
        direct_cpu_seconds=cpu(direct) if direct else None,
        cached_peak_bytes=peak(metrics),
        direct_peak_bytes=peak(direct) if direct else None,
        cached_step=metrics["cgroup_after"].get("path"),
        direct_step=direct.get("cgroup_after", {}).get("path"),
        timed_start_unix=metrics["start_unix_ns"] / 1e9,
        direct_start_unix=direct.get("start_unix_ns", 0) / 1e9 if direct else None,
        job_seconds=state["finished_unix"] - state["started_unix"],
        untimed_job_seconds=state["finished_unix"]
        - state["started_unix"]
        - measured
        - (direct_seconds if direct and not reused else 0),
        job_id=state["job_id"],
        node=state["node"],
        attempt=directory.parent.name,
        carried_forward=carried,
        source_manifest=state.get("source_manifest"),
        source_evidence=MAIN + "/attempts/" + directory.parent.name + "/" + task["id"],
        source_archive_member="main/attempts/"
        + directory.parent.name
        + "/"
        + task["id"],
        input_path=task["input"]["path"],
        input_sha256=task["input"]["sha256"],
        cache_path=task["cache"]["directory"],
        cache_payload_sha256=task["cache"]["files"]["vcfcache_annotated.bcf"],
        cache_recipe_sha256=task["cache"]["files"]["annotation.yaml"],
        mismatched_records=compare.get("record_mismatches", 0),
        vep_raw_annotation_differences=raw_vep,
        headers_equal=compare.get(
            "relevant_headers_equal",
            compare.get("against_cached", {}).get("relevant_headers_equal"),
        ),
        annotate_misses_seconds=annotate,
        cache_side_residual_seconds=residual,
        engine_seconds=engine,
        cli_other_seconds=measured - engine,
        preprocess_seconds=phases.get("preprocess", 0),
        route_seconds=phases.get("route", 0),
        assemble_seconds=phases.get("assemble", 0),
        engine_other_seconds=engine - annotate - sum(phases.values()),
        annotation_fraction=annotate / measured,
        non_annotation_fraction=residual / measured,
        hypothetical_zero_non_annotation_gain=measured / annotate if annotate else None,
        equal_cost_zero_overhead_speedup=1 / (1 - f) if f < 1 else None,
        effective_direct_us_per_record=direct_seconds / n * 1e6,
        effective_miss_us_per_record=annotate / misses * 1e6 if misses else None,
        experimental_flags_active=any(
            stats.get(k, False)
            for k in [
                "experimental_indexed_misses",
                "experimental_native_route",
                "experimental_trust_routed_hits",
            ]
        ),
    )
    if min(row["engine_other_seconds"], row["cli_other_seconds"], residual) < -0.05:
        raise ValueError("phase wall times do not form a nonnegative decomposition")
    if row["experimental_flags_active"] or row["engine_effective"] != "adaptive":
        raise ValueError("unapproved engine configuration")
    return row


def collect(output, resume_preparation=False):
    """Download a bounded terminal snapshot, compress originals, export analysis data."""
    if output.exists():
        if (
            not resume_preparation
            or (output / "source/SNAPSHOT.json").exists()
            or (output / "figures").exists()
        ):
            raise FileExistsError(
                "completed snapshots/figures are immutable; choose a new output"
            )
    else:
        output.mkdir(parents=True)
    (output / "source").mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="vcfcache-final-snapshot-") as temporary:
        root = Path(temporary)
        main = root / "main"
        repair = root / "repair"
        main.mkdir()
        repair.mkdir()
        for remote, local in [(MAIN, main), (REPAIR, repair)]:
            for name in ["manifest.json", "READY", "submissions.json"]:
                copy_remote(remote + "/" + name, local / (name))
        for name in [
            "AUTO_DISPATCH.json",
            "FINAL_PARENT_COLLECTOR.json",
            "DIRECT_REUSE_PREFLIGHT.json",
            "CONTROLLER_DISPATCH.json",
        ]:
            copy_remote(REPAIR + "/" + name, repair / name)
        filters = [
            "--exclude=cli/",
            "--exclude=resources.jsonl",
            "--include=*/",
            *[
                "--include=" + name
                for name in [
                    "state.json",
                    "command.json",
                    "comparison.json",
                    "statistics.json",
                    "resources.json",
                    "REUSED.json",
                ]
            ],
            "--exclude=*",
        ]
        (main / "attempts").mkdir()
        copy_remote(MAIN + "/attempts/", str(main / "attempts") + "/", filters)
        primary = load(main / "manifest.json")
        fixed = load(repair / "manifest.json")
        for local in [main, repair]:
            if (
                sha256sum(local / "manifest.json")
                != (local / "READY").read_text().strip()
            ):
                raise ValueError("downloaded manifest does not match its READY digest")
        tasks = {t["id"]: t for t in primary["tasks"]}
        if len(tasks) != 600:
            raise ValueError("manuscript inventory changed")
        latest = {}
        attempts = []
        for file in sorted((main / "attempts").glob("*/*/state.json")):
            state = load(file)
            directory = file.parent
            identifier = directory.name
            if identifier not in tasks:
                raise ValueError("unknown condition in attempt tree")
            latest[identifier] = (state, directory)
            attempts.append(
                dict(
                    condition_id=identifier,
                    attempt=directory.parent.name,
                    state=state["state"],
                    job_id=state.get("job_id"),
                    node=state.get("node"),
                    started_unix=state.get("started_unix"),
                    finished_unix=state.get("finished_unix"),
                    error=state.get("error"),
                    latest=False,
                )
            )
        if set(latest) != set(tasks):
            raise ValueError("missing final condition")
        carried_builds = {}
        carried_folder = root / "carried_manifests"
        carried_folder.mkdir()
        for state, _directory in latest.values():
            if state.get("carried_forward"):
                campaign = state["carried_forward"]["campaign"]
                if campaign not in carried_builds:
                    target = carried_folder / (Path(campaign).parent.name + ".json")
                    copy_remote(campaign + "/manifest.json", target)
                    carried_manifest = load(target)

                    def product(manifest):
                        return {
                            k: v
                            for k, v in manifest["freeze"]["build"]["files"].items()
                            if k.startswith(("vcfcache/", "native/", "resources/"))
                        }

                    if product(primary) != product(carried_manifest):
                        raise ValueError("carried product files differ")
                    carried_builds[campaign] = carried_manifest["freeze"]["build"][
                        "commit"
                    ]
        rows = [
            flatten(
                task,
                latest[identifier][1],
                latest[identifier][0],
                primary,
                fixed,
                carried_builds,
            )
            for identifier, task in tasks.items()
        ]
        # Resolve the twelve duplicate WGS plotting views without creating measurements.
        request = load(Path(__file__).with_name("reference-request-v2-20261027.json"))
        request_rows = {r["id"]: r for r in request["rows"]}
        aliases = primary["freeze"]["aliases"]
        views = [dict(r) for r in rows]
        lookup = {r["condition_id"]: r for r in rows}
        for identifier, canonical in aliases.items():
            view = dict(lookup[canonical])
            source = request_rows[identifier]
            for field in [
                "family",
                "tool",
                "level",
                "assay",
                "cohort",
                "assembly",
                "sample",
                "strategy",
            ]:
                view[field] = source[field]
            view.update(
                condition_id=identifier, measurement_id=canonical, is_alias=True
            )
            views.append(view)
        for a in attempts:
            a["latest"] = a["attempt"] == latest[a["condition_id"]][1].parent.name
        table(output / "tables/cells.tsv", rows)
        table(output / "tables/plot_views.tsv", views)
        table(output / "tables/attempts.tsv", attempts)
        # Collect frozen cache provenance and recipe definitions; no payload downloads.
        paths = set()
        for task in primary["tasks"]:
            paths.add(task["cache"]["canonical_directory"] + "/annotation.yaml")
            paths.add(task["cache"]["canonical_directory"] + "/params.snapshot.yaml")
            paths.update(
                a["path"]
                for a in task["annotation_assets"]
                if "provenance" in a["path"] and a["path"].endswith(".json")
            )
        supplemental = root / "cache_metadata"
        supplemental.mkdir()
        index = []
        resolved = controller_paths(
            sorted(p for p in paths if p.startswith("/results/"))
        )
        for number, path in enumerate(sorted(paths)):
            # All selected resources are on shared /results after canonical staging.
            if not path.startswith("/results/"):
                continue
            name = f"{number:03d}-" + Path(path).name
            copy_remote(resolved[path], supplemental / name)
            index.append(
                dict(
                    path=path,
                    archive_path="cache_metadata/" + name,
                    sha256=sha256sum(supplemental / name),
                )
            )
        write_json_atomic(supplemental / "INDEX.json", index)
        # Exact original strategy metadata carries cohort seed identities/build provenance.
        strategy_paths = set()
        for task in primary["tasks"]:
            for role in ["direct", "cached"]:
                path = task["baseline"].get(role, {}).get("path", "")
                if "/measured/" in path:
                    strategy_paths.add(
                        path.split("/measured/")[0] + "/manifests/strategies.json"
                    )
        for number, path in enumerate(sorted(strategy_paths)):
            copy_remote(path, supplemental / f"strategies-{number}.json")
        archive = output / "source/terminal-evidence.tar.gz"
        with tarfile.open(archive, "w:gz") as bundle:
            for path in sorted(root.iterdir()):
                bundle.add(path, arcname=path.name)
        # A compact, directly inspectable recipe/source inventory accompanies the archive.
        table(output / "tables/cache_metadata_index.tsv", index)
        summary = dict(
            created_utc=datetime.now(timezone.utc).isoformat(),
            remote_main=MAIN,
            remote_repair=REPAIR,
            manifest_sha256=sha256sum(main / "manifest.json"),
            repair_manifest_sha256=sha256sum(repair / "manifest.json"),
            conditions=len(rows),
            plot_views=len(views),
            aliases=len(aliases),
            attempts=len(attempts),
            states=dict(Counter(r["state"] for r in rows)),
            evidence_classes=dict(Counter(r["evidence_class"] for r in rows)),
            archive_sha256=sha256sum(archive),
            last_finished_unix=max(s["finished_unix"] for s, d in latest.values()),
            note="Historical direct timing reuse and synthetic composed-reference limitations retained; no failed/excluded row promoted to validated",
        )
        write_json_atomic(output / "source/SNAPSHOT.json", summary)
        print(json.dumps(summary, indent=2))


def main():
    """Create a NEW immutable local evidence snapshot and analysis tables."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--resume-preparation", action="store_true")
    args = parser.parse_args()
    collect(args.output, args.resume_preparation)


if __name__ == "__main__":
    main()
