"""Audit finished pinned-node blocks; never promote partial repeats or best timings."""

import argparse
import csv
import json
import tarfile
from collections import defaultdict
from pathlib import Path
from statistics import median

from benchmarks.pass_fusion_pilot.run import MODES, manifest
from benchmarks.performance_sprint.repair_headers import direct_binding
from benchmarks.run_pilot import sha256sum, write_json_atomic


def read(path):
    """Read one durable checkpoint."""
    return json.loads(path.read_text())


def table(path, rows):
    """Keep missing fields explicit and TSVs rectangular."""
    if not rows:
        return
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w") as stream:
        writer = csv.DictWriter(stream, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {k: ("NA" if row.get(k) in {None, ""} else row[k]) for k in fields}
            )


def io_deltas(metrics):
    """Preserve devices separately: summing stacked block devices can double-count."""

    def parse(snapshot):
        result = {}
        for line in (snapshot or {}).get("io_stat", []):
            device, *items = line.split()
            result[device] = {k: int(v) for k, v in (item.split("=") for item in items)}
        return result

    before = parse(metrics.get("cgroup_before"))
    after = parse(metrics.get("cgroup_after"))
    return json.dumps(
        {
            device: {
                key: value - before.get(device, {}).get(key, 0)
                for key, value in values.items()
            }
            for device, values in after.items()
        },
        sort_keys=True,
    )


def summarize(root):
    """Export only complete, strictly validated blocks, preserving failure statuses."""
    source = root / "source" if (root / "source/manifest.json").exists() else root
    data = manifest(source)
    tasks = {t["sample"]: t for t in data["tasks"]}
    conditions = []
    attempts = []
    for case in data["cases"]:
        folder = source / "cases" / case["id"]
        state_path = folder / "state.json"
        state = read(state_path) if state_path.exists() else dict(state="not_started")
        recovery_path = folder / "AUDIT_RECOVERY.json"
        recovered = False
        if recovery_path.exists():
            recovery = read(recovery_path)
            assert (
                recovery["accepted"]
                and recovery["audit_only"]
                and recovery["no_timing_reruns"]
            )
            assert recovery["original_state_sha256"] == sha256sum(state_path)
            assert recovery["source_manifest_sha256"] == sha256sum(
                source / "manifest.json"
            )
            for name, digest in recovery["comparison_sha256"].items():
                assert sha256sum(folder / name) == digest
            for binding in recovery["timed_bindings"]:
                directory = folder / f"r{binding['repeat']:02d}" / binding["mode"]
                assert (
                    sha256sum(directory / "command.json") == binding["command_sha256"]
                )
                assert (
                    sha256sum(directory / "resources.json")
                    == binding["resources_sha256"]
                )
            state = dict(state, state="complete")
            recovered = True
        attempts.append(
            dict(
                case=case["id"],
                node=case["node"],
                sample=case["sample"],
                state=state["state"],
                job_id=state.get("job_id"),
                error=state.get("error"),
                audit_recovered=recovered,
            )
        )
        if state["state"] != "complete":
            continue
        assert (
            state["node"] == case["node"]
            and len(state["rows"]) == len(MODES) * data["repetitions"]
        )
        hardware = read(folder / "hardware.json")
        cpu = next(
            item["data"]
            for item in hardware["lscpu"]["lscpu"]
            if item["field"].rstrip(":") == "Model name"
        )
        task = tasks[case["sample"]]
        keys = set()
        for row in state["rows"]:
            repeat, mode = row["repeat"], row["mode"]
            key = (repeat, mode)
            assert key not in keys
            keys.add(key)
            directory = folder / f"r{repeat:02d}" / mode
            command = read(directory / "command.json")
            stats = read(directory / "statistics.json")
            metrics = read(directory / "resources.json")
            assert direct_binding(command["task"]) == direct_binding(task)
            assert (
                command["source_commit"] == data["commit"]
                and command["parameters"]["threads"] == 8
            )
            assert metrics["exit_code"] == 0 and len(metrics["affinity"]) == 8
            assert command["experimental_environment"] == {
                "VCFCACHE_EXPERIMENT_FUSED_PUBLISH": (
                    "1" if mode.startswith("fused") else "0"
                )
            }
            counts = stats["variant_counts"]
            assert (
                counts["input_variants"]
                == counts["total_output"]
                == task["input"]["records"]
            )
            if mode != "direct":
                comparison = read(directory / "against_direct.json")
                assert comparison["records_compared"] == task["input"]["records"]
                assert (
                    comparison["semantic_pass"] and comparison["ignored_fields"] == []
                )
                assert stats["engine"] == "adaptive" and stats[
                    "fused_publish"
                ] == mode.startswith("fused")
                if mode.startswith("fused"):
                    assert (
                        stats["shards"] * (2 + stats["publish_threads_per_shard"]) <= 8
                    )
            conditions.append(
                dict(
                    case=case["id"],
                    node=case["node"],
                    cpu_model=cpu,
                    sample=case["sample"],
                    repeat=repeat,
                    mode=mode,
                    wall_seconds=metrics["wall_seconds"],
                    cpu_seconds=metrics["cpu_seconds"],
                    actual_shards=stats.get("shards"),
                    input_variants=counts["input_variants"],
                    missing_variants=counts.get("missing_variants"),
                    hit_fraction=(
                        1 - counts["missing_variants"] / counts["input_variants"]
                        if mode != "direct"
                        else None
                    ),
                    cgroup_peak_gib=metrics["cgroup_after"]["memory_peak_bytes"]
                    / 1024**3,
                    local_block_io_deltas=io_deltas(metrics),
                    job_id=state["job_id"],
                )
            )
    by_repeat = defaultdict(dict)
    for row in conditions:
        by_repeat[(row["case"], row["repeat"])][row["mode"]] = row
    pairs = []
    for (case, repeat), values in by_repeat.items():
        assert set(values) == set(MODES)
        direct = values["direct"]["wall_seconds"]
        auto = values["control-auto"]["wall_seconds"]
        four = values["control4"]["wall_seconds"]
        misses = {
            v["missing_variants"] for mode, v in values.items() if mode != "direct"
        }
        assert len(misses) == 1
        for mode, value in values.items():
            if mode == "direct":
                continue
            pairs.append(
                dict(
                    case=case,
                    repeat=repeat,
                    node=value["node"],
                    cpu_model=value["cpu_model"],
                    sample=value["sample"],
                    mode=mode,
                    direct_seconds=direct,
                    cached_seconds=value["wall_seconds"],
                    speedup_vs_direct=direct / value["wall_seconds"],
                    gain_over_default=auto / value["wall_seconds"],
                    gain_over_control4=four / value["wall_seconds"],
                    seconds_saved_vs_default=auto - value["wall_seconds"],
                    actual_shards=value["actual_shards"],
                )
            )
    grouped = defaultdict(list)
    for row in pairs:
        grouped[(row["case"], row["mode"])].append(row)
    summary = []
    for (case, mode), values in grouped.items():
        assert len(values) == data["repetitions"]
        summary.append(
            dict(
                case=case,
                node=values[0]["node"],
                cpu_model=values[0]["cpu_model"],
                sample=values[0]["sample"],
                mode=mode,
                technical_repeats=len(values),
                median_cached_seconds=median(v["cached_seconds"] for v in values),
                median_speedup_vs_direct=median(v["speedup_vs_direct"] for v in values),
                median_gain_over_default=median(v["gain_over_default"] for v in values),
                median_gain_over_control4=median(
                    v["gain_over_control4"] for v in values
                ),
            )
        )
    table(root / "conditions.tsv", conditions)
    table(root / "paired_ratios.tsv", pairs)
    table(root / "summary_by_node.tsv", summary)
    table(root / "attempts.tsv", attempts)
    inventory = [
        dict(path=str(p.relative_to(root)), sha256=sha256sum(p))
        for p in sorted(source.rglob("*"))
        if p.is_file()
        and p.name not in {"EVIDENCE_INVENTORY.json", "ACCEPTANCE.json"}
        and p.suffix in {".json", ".log"}
    ]
    write_json_atomic(root / "EVIDENCE_INVENTORY.json", inventory)
    result = dict(
        accepted_blocks=sum(row["state"] == "complete" for row in attempts),
        conditions=len(conditions),
        planned_blocks=len(data["cases"]),
        attempts=attempts,
    )
    write_json_atomic(root / "ACCEPTANCE.json", result)
    # Freeze only terminal experiment evidence. Never overwrite an earlier
    # archive or include large disposable BCFs / recursively include ourselves.
    if all(row["state"] in {"complete", "failed", "skipped"} for row in attempts):
        archive_path = root / "terminal-evidence.tar.gz"
        if not archive_path.exists():
            with tarfile.open(archive_path, "x:gz") as archive:
                for name in [
                    "manifest.json",
                    "READY",
                    "submissions.json",
                    "review_submission.json",
                    "REVIEW.json",
                    "cases",
                    "logs",
                ]:
                    path = source / name
                    if path.exists():
                        archive.add(path, arcname="source/" + name)
        write_json_atomic(
            root / "SNAPSHOT.json",
            dict(
                archive_sha256=sha256sum(archive_path),
                source_manifest_sha256=sha256sum(source / "manifest.json"),
                source_commit=data["commit"],
                accepted_blocks=result["accepted_blocks"],
                planned_blocks=result["planned_blocks"],
            ),
        )
    print(json.dumps({k: v for k, v in result.items() if k != "attempts"}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    summarize(parser.parse_args().root)
