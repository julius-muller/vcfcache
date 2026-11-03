"""Consolidate the bounded WES pilot without choosing fastest technical repeats."""

import argparse
import csv
import json
from collections import Counter
from pathlib import Path


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
        writer.writerows(rows)


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
                        cached_seconds=result["wall_seconds"],
                        direct_seconds=direct_seconds,
                        ratio=direct_seconds / result["wall_seconds"],
                        cpu_seconds=metrics["cpu_seconds"],
                        peak_gib=metrics["cgroup_after"]["memory_peak_bytes"] / 1024**3,
                        campaign=name,
                        job_id=state["job_id"],
                        direct_origin=direct_origin,
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
    (root / "COLLECTION.json").write_text(
        json.dumps(
            dict(
                cli_rows=len(rows),
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
