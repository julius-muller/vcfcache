"""Consolidate the bounded WES pilot without choosing fastest technical repeats."""

import argparse
import csv
import json
import re
from collections import Counter
from pathlib import Path

from benchmarks.performance_sprint.repair_headers import direct_binding
from benchmarks.run_pilot import sha256sum


def read(path):
    """Read one durable checkpoint."""
    return json.loads(path.read_text())


def table(path, rows):
    """Write complete scalar records; preserve unavailable values as missing."""
    if not rows:
        return
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w") as stream:
        writer = csv.DictWriter(stream, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        # Explicit missing values preserve rectangular TSVs without empty
        # trailing cells being mistaken for accidental whitespace by Git.
        writer.writerows(
            {
                key: "NA" if value is None or value == "" else value
                for key, value in row.items()
            }
            for row in rows
        )


def collect(root):
    """Keep all attempts and separate CLI, batch-warm and loading costs."""
    tables = root / "tables"
    tables.mkdir(exist_ok=True)
    rows = []
    attempts = []
    warm = []
    cumulative = []
    for name in [
        "20261103-screen",
        "20261103-coalesced",
        "20261103-coalesced-auditfix",
        "20261103-confirm",
    ]:
        campaign = root / "source" / name
        if not campaign.exists():
            continue
        manifest = read(campaign / "manifest.json")
        assert (
            sha256sum(campaign / "manifest.json")
            == (campaign / "READY").read_text().strip()
        )
        tasks = {t["id"]: t for t in manifest["tasks"]}
        for state_path in sorted((campaign / "cases").glob("*/state.json")):
            state = read(state_path)
            task = tasks[state_path.parent.name]
            attempts.append(
                dict(
                    campaign=name,
                    id=task["id"],
                    tool=task["tool"],
                    sample=task["sample"],
                    state=state["state"],
                    job_id=state["job_id"],
                    error=state.get("error"),
                )
            )
            # The initial coalesced build omitted its marker from the summary.
            # All six failed audits are retained, never promoted or best-picked.
            if state["state"] != "complete":
                continue
            if state.get("reused_direct"):
                direct_seconds = state["reused_direct"]["wall_seconds"]
                direct_origin = state["reused_direct"]["path"]
            else:
                direct_seconds = next(
                    r["wall_seconds"] for r in state["rows"] if r["mode"] == "direct"
                )
                direct_origin = str(state_path.parent)
            for result in state["rows"]:
                mode = result["mode"]
                metrics = read(state_path.parent / mode / "resources.json")
                command = read(state_path.parent / mode / "command.json")
                assert direct_binding(command["task"]) == direct_binding(task)
                assert command["binary"] == manifest["binary"]
                stats = read(state_path.parent / mode / "statistics.json")
                counts = stats["variant_counts"]
                assert (
                    counts["input_variants"]
                    == counts["total_output"]
                    == task["input"]["records"]
                )
                if mode in {"resident", "coalesced"}:
                    assert stats["effective_engine"] == "native"
                if mode == "coalesced":
                    assert stats["coalesced_chunks"] is True
                phases = {
                    p["phase"]: p["wall_seconds"]
                    for p in stats.get("phase_timings", [])
                }
                annotation = sum(
                    p["wall_seconds"]
                    for p in stats.get("stage_timings", [])
                    if p["stage"] == "annotate-misses"
                )
                assert metrics["exit_code"] == 0 and len(metrics["affinity"]) == 8
                if mode != "direct":
                    assert read(state_path.parent / mode / "against_direct.json")[
                        "semantic_pass"
                    ]
                rows.append(
                    dict(
                        id=task["id"],
                        sample=task["sample"],
                        cohort=task["cohort"],
                        tool=task["tool"],
                        mode=mode,
                        input_variants=counts["input_variants"],
                        missing_variants=counts["missing_variants"],
                        hit_fraction=(
                            (1 - counts["missing_variants"] / counts["input_variants"])
                            if mode != "direct"
                            else None
                        ),
                        cached_seconds=result["wall_seconds"],
                        direct_seconds=direct_seconds,
                        ratio=direct_seconds / result["wall_seconds"],
                        cpu_seconds=metrics["cpu_seconds"],
                        peak_gib=metrics["cgroup_after"]["memory_peak_bytes"] / 1024**3,
                        campaign=name,
                        job_id=state["job_id"],
                        direct_origin=direct_origin,
                        annotation_seconds=annotation if mode != "direct" else None,
                        lookup_seconds=phases.get("lookup"),
                        zero_non_annotation_bound=(
                            direct_seconds / annotation if annotation else None
                        ),
                        screen_order=(
                            "direct-first"
                            if name == "20261103-screen"
                            else (
                                "reused screen direct"
                                if state.get("reused_direct")
                                else "counterbalanced"
                            )
                        ),
                    )
                )
    # One accepted timing for each distinct CLI condition/configuration.
    keys = [(r["id"], r["mode"]) for r in rows]
    assert len(keys) == len(set(keys))
    batch = root / "source/20261103-batch"
    batch_status = []
    if batch.exists():
        for folder in sorted(batch.glob("batch-*")):
            if not (folder / "state.json").exists():
                continue
            state = read(folder / "state.json")
            batch_status.append(state)
            if state["state"] != "complete":
                continue
            direct = read(folder / "direct-batch.json")["rows"]
            resident = read(folder / "resident-batch.json")["rows"]
            assert len(direct) == len(resident) == 12
            dm = read(folder / "direct-batch/resources.json")
            cm = read(folder / "resident-batch/resources.json")
            assert dm["exit_code"] == cm["exit_code"] == 0
            load = read(folder / "load.json")
            overhead = cm["wall_seconds"] - sum(r["wall_seconds"] for r in resident)
            assert overhead >= 0
            total_d = 0
            total_c = overhead
            for index, (d, c) in enumerate(zip(direct, resident, strict=True), 1):
                assert (
                    d["id"] == c["id"]
                    and read(folder / (c["id"] + "-comparison.json"))["semantic_pass"]
                )
                warm.append(
                    dict(
                        tool=state["tool"],
                        id=c["id"],
                        sample=c["sample"],
                        direct_seconds=d["wall_seconds"],
                        warm_seconds=c["wall_seconds"],
                        warm_ratio=d["wall_seconds"] / c["wall_seconds"],
                        loader_seconds=load["controller_elapsed_seconds"],
                        session_overhead_seconds=overhead,
                        first_sample_with_session_overhead_seconds=(
                            overhead + c["wall_seconds"] if index == 1 else None
                        ),
                        session_peak_gib=cm["cgroup_after"]["memory_peak_bytes"]
                        / 1024**3,
                    )
                )
                total_d += d["wall_seconds"]
                total_c += c["wall_seconds"]
                cumulative.append(
                    dict(
                        tool=state["tool"],
                        samples=index,
                        direct_seconds=total_d,
                        cached_seconds_including_load=total_c,
                        ratio=total_d / total_c,
                    )
                )
    table(tables / "cli_conditions.tsv", rows)
    table(tables / "attempts.tsv", attempts)
    table(tables / "resident_warm.tsv", warm)
    table(tables / "resident_cumulative.tsv", cumulative)
    capacity = []
    for label, representation in [
        ("20261103-probes", "eager INFO"),
        ("20261103-lazy-probes", "lazy BCF"),
    ]:
        probe = root / "source" / label
        for path in sorted(probe.glob("*/resources.json")):
            metrics = read(path)
            text = (path.parent / "probe.out").read_text()
            match = re.search(r"after (\d+) records", text)
            capacity.append(
                dict(
                    tool=path.parent.name,
                    representation=representation,
                    exit_code=metrics["exit_code"],
                    elapsed_seconds=metrics["wall_seconds"],
                    cgroup_peak_gib=metrics["cgroup_after"]["memory_peak_bytes"]
                    / 1024**3,
                    child_peak_rss_gib=metrics["largest_child_rss_kib"] / 1024**2,
                    records_before_stop=int(match[1]) if match else None,
                    outcome=(
                        "20-GiB RSS guard"
                        if "20-GiB RSS budget" in text
                        else text[-200:].strip()
                    ),
                )
            )
    table(tables / "capacity_probes.tsv", capacity)
    tuning = []
    campaigns = sorted((root / "source").glob("20261103-workers-w*")) + sorted(
        (root / "source").glob("20261103-workers8-*")
    )
    for campaign in campaigns:
        manifest = read(campaign / "manifest.json")
        assert (
            sha256sum(campaign / "manifest.json")
            == (campaign / "READY").read_text().strip()
        )
        tasks = {t["id"]: t for t in manifest["tasks"]}
        for state_path in sorted((campaign / "cases").glob("*/state.json")):
            state = read(state_path)
            assert state["state"] == "complete"
            task = tasks[state_path.parent.name]
            assert len(state["rows"]) == 1
            row = state["rows"][0]
            mode = row["mode"]
            stats = read(state_path.parent / mode / "statistics.json")
            metrics = read(state_path.parent / mode / "resources.json")
            command = read(state_path.parent / mode / "command.json")
            assert direct_binding(command["task"]) == direct_binding(task)
            assert command["binary"] == manifest["binary"]
            assert metrics["exit_code"] == 0 and len(metrics["affinity"]) == 8
            counts = stats["variant_counts"]
            assert (
                counts["input_variants"]
                == counts["total_output"]
                == task["input"]["records"]
            )
            assert stats["effective_engine"] == "native"
            assert int(mode.removeprefix("resident-w")) == stats["workers"]
            assert read(state_path.parent / mode / "against_direct.json")[
                "semantic_pass"
            ]
            control = next(
                r
                for r in rows
                if r["id"] == state_path.parent.name and r["mode"] == "resident"
            )
            phases = {p["phase"]: p["wall_seconds"] for p in stats["phase_timings"]}
            annotation = sum(
                p["wall_seconds"]
                for p in stats.get("stage_timings", [])
                if p["stage"] == "annotate-misses"
            )
            tuning.append(
                dict(
                    id=state_path.parent.name,
                    tool=task["tool"],
                    sample=task["sample"],
                    cohort=task["cohort"],
                    campaign=campaign.name,
                    mode=mode,
                    input_variants=counts["input_variants"],
                    missing_variants=counts["missing_variants"],
                    hit_fraction=1
                    - counts["missing_variants"] / counts["input_variants"],
                    workers=stats["workers"],
                    cached_seconds=row["wall_seconds"],
                    direct_seconds=state["reused_direct"]["wall_seconds"],
                    ratio=state["reused_direct"]["wall_seconds"] / row["wall_seconds"],
                    gain_over_resident4=control["cached_seconds"] / row["wall_seconds"],
                    lookup_seconds=phases["lookup"],
                    control_lookup_seconds=control["lookup_seconds"],
                    annotation_seconds=annotation,
                    zero_non_annotation_bound=(
                        state["reused_direct"]["wall_seconds"] / annotation
                        if annotation
                        else None
                    ),
                    peak_gib=metrics["cgroup_after"]["memory_peak_bytes"] / 1024**3,
                    cpu_seconds=metrics["cpu_seconds"],
                    job_id=state["job_id"],
                )
            )
    assert len({(r["id"], r["workers"]) for r in tuning}) == len(tuning)
    table(tables / "worker_tuning.tsv", tuning)
    # Freeze small terminal evidence independently of disposable output BCFs.
    inventory = [
        dict(path=str(p.relative_to(root)), sha256=sha256sum(p), bytes=p.stat().st_size)
        for p in sorted((root / "source").rglob("*"))
        if p.is_file()
    ]
    (root / "EVIDENCE_INVENTORY.json").write_text(
        json.dumps(inventory, indent=2) + "\n"
    )
    (root / "COLLECTION.json").write_text(
        json.dumps(
            dict(
                cli_rows=len(rows),
                worker_allocation_rows=len(tuning),
                attempts=dict(Counter(a["state"] for a in attempts)),
                batch=batch_status,
            ),
            indent=2,
        )
        + "\n"
    )
    print(
        json.dumps(
            dict(
                cli_rows=len(rows),
                attempts=dict(Counter(a["state"] for a in attempts)),
                warm_rows=len(warm),
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    collect(parser.parse_args().root)
