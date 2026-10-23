"""Offline-safe collection: missing work is not success and overlaps are not summed."""

import csv
import hashlib
import json
import shutil
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

from benchmarks.engine_profile.run import latest, read, verify
from benchmarks.run_pilot import write_json_atomic


def union_seconds(spans):
    """Duration of an interval union, not the sum of parallel/nested workers."""
    total, end = 0, None
    for start, finish in sorted(spans):
        if finish < start:
            raise ValueError("negative span")
        total += max(0, finish - max(start, end if end is not None else start))
        end = max(finish, end if end is not None else finish)
    return total / 1e9


def table(path, rows, fields):
    """Write portable LF-separated TSV, including a header for empty data."""
    with path.open("w") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def job_states(root):
    """Snapshot live Slurm when available; retain unavailable/expired states honestly."""
    batches = (
        read(root / "submissions.json") if (root / "submissions.json").exists() else []
    )
    result = {}
    for batch in batches:
        for job in batch["jobs"]:
            identifier = str(job["job_id"])
            if shutil.which("scontrol"):
                check = subprocess.run(
                    ["scontrol", "show", "job", identifier, "-o"],
                    capture_output=True,
                    text=True,
                )
                fields = dict(
                    token.split("=", 1)
                    for token in check.stdout.split()
                    if "=" in token
                )
                result[identifier] = dict(
                    group=job["group"],
                    state=fields.get("JobState", "unavailable"),
                    reason=fields.get("Reason"),
                    exit_code=fields.get("ExitCode"),
                )
            else:
                result[identifier] = dict(group=job["group"], state="unavailable")
    return result


def collect(root, write=True):
    """Collect all attempts and produce a content-addressed, immutable review."""
    manifest = verify(root)
    rows, phases, signatures = [], [], []
    for group in manifest["groups"]:
        for cell in group["cells"]:
            files = sorted((root / "attempts").glob(f"*/{cell['id']}/state.json"))
            state = latest(root, cell["id"]) or {}
            metrics = state.get("metrics", {})
            if not metrics and files and (files[-1].parent / "resources.json").exists():
                metrics = read(files[-1].parent / "resources.json")
            row = dict(
                cell=cell["id"],
                group=group["id"],
                sample=cell["sample"],
                recipe=cell["recipe"],
                assay=cell["assay"],
                purpose=cell["purpose"],
                role=cell["role"],
                cpus=cell["cpus"],
                workers=cell["workers"],
                memory_mb=cell["memory_mb"],
                profile=cell["profile"],
                requested_engine=cell["engine"],
                effective_engine=state.get("effective_engine"),
                state=state.get("state", "missing"),
                reason=state.get("reason", state.get("error", "")),
                direct_pass=state.get("direct_pass"),
                instrumentation_equivalent=state.get("instrumentation_equivalent"),
                wall_seconds=metrics.get("wall_seconds"),
                cpu_seconds=metrics.get("cpu_seconds"),
                peak_bytes=(metrics.get("cgroup_after") or {}).get("memory_peak_bytes"),
                cgroup=(metrics.get("cgroup_after") or {}).get("path"),
                attempts=len(files),
            )
            signatures.append(dict(cell=cell, states=[read(p) for p in files]))
            before, after = (
                metrics.get("cgroup_before") or {},
                metrics.get("cgroup_after") or {},
            )
            row["aggregate_cpu_seconds"] = (
                (
                    after.get("cpu_stat", {}).get("usage_usec", 0)
                    - before.get("cpu_stat", {}).get("usage_usec", 0)
                )
                / 1e6
                if after
                else None
            )
            for counter in ("rbytes", "wbytes"):

                def io_total(snapshot, counter=counter):
                    return sum(
                        int(token.split("=")[1])
                        for line in snapshot.get("io_stat", [])
                        for token in line.split()
                        if token.startswith(counter + "=")
                    )

                row[counter] = io_total(after) - io_total(before) if after else None
            if row["wall_seconds"] and row["cpu_seconds"] is not None:
                row["average_cpus"] = row["cpu_seconds"] / row["wall_seconds"]
            else:
                row["average_cpus"] = None
            events = []
            broken = []
            for trace in (
                sorted(files[-1].parent.glob("trace/*.jsonl")) if files else []
            ):
                for line in trace.read_text().splitlines():
                    try:
                        event = json.loads(line)
                        if event["duration_ns"] < 0:
                            raise ValueError("negative duration")
                        events.append(event)
                    except (ValueError, KeyError) as error:
                        broken.append(f"{trace.name}: {error}")
            signatures.append(dict(events=events, broken=broken))
            spans = defaultdict(list)
            for event in events:
                start = event["start_unix_ns"]
                spans[event["name"]].append((start, start + event["duration_ns"]))
            for name, intervals in spans.items():
                phases.append(
                    dict(
                        cell=cell["id"],
                        recipe=cell["recipe"],
                        assay=cell["assay"],
                        phase=name,
                        union_seconds=union_seconds(intervals),
                        summed_seconds=sum(b - a for a, b in intervals) / 1e9,
                        occurrences=len(intervals),
                        requested_engine=cell["engine"],
                        effective_engine=row["effective_engine"],
                    )
                )
            # Container/wait envelopes would falsely make unexplained work zero.
            containers = {
                "cli.import_and_run",
                "native.attempt",
                "adaptive.attempt",
                "native.execute",
                "native.await_misses",
                "native.await_annotation",
                "native.lookup.worker",
            }
            intervals = [
                item
                for name, values in spans.items()
                if name not in containers
                for item in values
            ]
            row["covered_phase_seconds"] = union_seconds(intervals)
            row["unattributed_seconds"] = (
                max(0, row["wall_seconds"] - row["covered_phase_seconds"])
                if row["wall_seconds"] is not None
                else None
            )
            row["profile_errors"] = "; ".join(broken)
            rows.append(row)
    scheduler = job_states(root)
    # Preserve terminal scheduler evidence before controller retention expires.
    cached = (
        read(root / "scheduler-status.json")
        if (root / "scheduler-status.json").exists()
        else {}
    )
    for identifier, value in scheduler.items():
        if value["state"] != "unavailable" or identifier not in cached:
            cached[identifier] = value
    batches = (
        read(root / "submissions.json") if (root / "submissions.json").exists() else []
    )
    terminal = {
        "COMPLETED",
        "FAILED",
        "CANCELLED",
        "TIMEOUT",
        "OUT_OF_MEMORY",
        "NODE_FAIL",
        "PREEMPTED",
        "BOOT_FAIL",
        "DEADLINE",
    }
    for row in rows:
        jobs = [
            job
            for batch in batches
            for job in batch["jobs"]
            if row["cell"] in job["cells"]
        ]
        if jobs and row["state"] in {"running", "missing"}:
            job = cached.get(str(jobs[-1]["job_id"]), {})
            if job.get("state") in terminal:
                row["reason"] = f"No final cell record; scheduler state {job['state']}"
                row["state"] = "incomplete_" + job["state"].lower()
    counts = dict(Counter(r["state"] for r in rows))
    summary = dict(cells=len(rows), counts=counts, scheduler=cached)
    print(json.dumps(summary, indent=2))
    if not write:
        return summary
    write_json_atomic(root / "scheduler-status.json", cached)
    digest = hashlib.sha256(
        json.dumps(
            dict(signatures=signatures, scheduler=cached), sort_keys=True
        ).encode()
    ).hexdigest()[:16]
    destination = root / "reports" / digest
    if destination.exists() and (destination / "COMPLETE").exists():
        return summary
    destination.mkdir(parents=True, exist_ok=True)
    table(destination / "cells.tsv", rows, list(rows[0]))
    table(
        destination / "phases.tsv",
        phases,
        [
            "cell",
            "recipe",
            "assay",
            "phase",
            "union_seconds",
            "summed_seconds",
            "occurrences",
            "requested_engine",
            "effective_engine",
        ],
    )
    write_json_atomic(destination / "summary.json", summary)
    calibration = []
    pairs = defaultdict(dict)
    for row in rows:
        if row["purpose"].startswith("instrumentation-pair-") and row[
            "state"
        ].startswith("complete"):
            pairs[row["purpose"]][row["profile"]] = row
    for pair, values in pairs.items():
        if True in values and False in values:
            on, off = values[True], values[False]
            calibration.append(
                dict(
                    pair=pair,
                    on_seconds=on["wall_seconds"],
                    off_seconds=off["wall_seconds"],
                    relative_overhead=on["wall_seconds"] / off["wall_seconds"] - 1,
                )
            )
    table(
        destination / "instrumentation.tsv",
        calibration,
        ["pair", "on_seconds", "off_seconds", "relative_overhead"],
    )
    duplicates = [
        key
        for key, n in Counter(r["cgroup"] for r in rows if r["cgroup"]).items()
        if n > 1
    ]
    lines = [
        "# Engine overhead profile",
        "",
        "Engineering diagnostics, not new publication speedup claims.",
        "",
        f"Cells: {len(rows)}. States: {json.dumps(counts, sort_keys=True)}.",
        "",
        "Read `cells.tsv` for requested/effective backends, strict audits and resources; `phases.tsv` for inclusive stage timings.",
        "Intervals are unioned within each phase. Different phases can overlap; do not stack or sum them as disjoint wall time.",
        "Unattributed time excludes container/wait envelopes and includes process startup and uninstrumented boundaries. Compression and decoding CPU attribution requires a usable stack profile; coarse overlapping wall spans cannot establish it.",
        "",
        "All attempts, commands, partial traces and raw audit failures are retained. Missing or running cells are not successes.",
        "Cold storage is NOT claimed: staging and hashing warm files; first-use-after-copy is an annotation-access condition only.",
        "",
        "## Integrity gates",
        "",
        f"Repeated cgroup paths (condition peak-memory values invalid if nonempty): {duplicates}",
        f"Trace parse failures: {sum(bool(r['profile_errors']) for r in rows)}.",
        "Direct self-checks are reference sanity checks, not independent equivalence proof. `complete_with_direct_difference` remains a failed direct comparison; no fields are ignored.",
        "",
        "## Optimization screening",
        "",
        f"Completed instrumentation pairs: {len(calibration)}/3. See `instrumentation.tsv`; absolute paired shifts above 5% require timing-noise/instrumentation review before interpreting stage fractions.",
        "Rank measured cache-side phases only after checking instrumentation-on/off pairs and CPU/worker controls. A removable phase is an optimistic bound, not a promised speedup.",
        "| Cell | Phase | Union wall s | Fraction of end-to-end |",
        "| --- | --- | ---: | ---: |",
    ]
    candidates = []
    cells_by_id = {r["cell"]: r for r in rows}
    for phase in phases:
        row = cells_by_id[phase["cell"]]
        if (
            row["purpose"] == "primary"
            and row["wall_seconds"]
            and row["instrumentation_equivalent"] is True
            and phase["phase"]
            in {
                "native.lookup",
                "native.construct_encode_output",
                "native.prepare_decode_split_stream",
                "native.finalize_output",
                "pipeline.preprocess",
                "pipeline.route",
                "pipeline.assemble",
                "pipeline.index",
            }
        ):
            fraction = phase["union_seconds"] / row["wall_seconds"]
            if fraction >= 0.05:
                candidates.append((fraction, phase))
    for fraction, phase in sorted(candidates, key=lambda x: x[0], reverse=True)[:12]:
        lines.append(
            f"| {phase['cell']} | {phase['phase']} | {phase['union_seconds']:.3f} | {fraction:.1%} |"
        )
    lines += [
        "",
        "Candidates require output equivalence, reproducible end-to-end benefit and a separate implementation review. Native fallback waste and memory pathology remain relevant even below a 5% phase threshold. No optimization is implemented by this campaign.",
        "",
        "## Review days later",
        "",
        "Run `python -m benchmarks.engine_profile.run collect CAMPAIGN` again after transferring the campaign directory if desired. No worker scratch or Slurm accounting is required for recorded cells.",
        "Plot with `Rscript benchmarks/engine_profile/plot.R REPORT_DIRECTORY`; only R/ggplot2 is used.",
        "",
    ]
    (destination / "REPORT.md").write_text("\n".join(lines))
    script = Path(__file__).with_name("plot.R")
    if shutil.which("Rscript"):
        with (destination / "plot.log").open("w") as log:
            result = subprocess.run(
                ["Rscript", str(script), str(destination)], stdout=log, stderr=log
            )
        (destination / "plot-status.txt").write_text(
            f"Rscript exit {result.returncode}\n"
        )
    else:
        (destination / "plot-status.txt").write_text(
            "Rscript unavailable on collector; run plot.R on the R workstation after collection.\n"
        )
    (destination / "COMPLETE").write_text(digest + "\n")
    write_json_atomic(
        root / "latest-review.json", dict(report=str(destination.relative_to(root)))
    )
    return summary
