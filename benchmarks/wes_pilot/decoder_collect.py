"""Audit the fixed-worker decoder-pool experiment and retain every paired cell."""

import argparse
import json
from statistics import median
from pathlib import Path

from benchmarks.performance_sprint.repair_headers import direct_binding
from benchmarks.run_pilot import sha256sum
from benchmarks.wes_pilot.collect import read, table


def collect(root):
    """Require unchanged paired work, CPU allocation and strict output equivalence."""
    source = root / "source"
    manifest = read(source / "manifest.json")
    assert sha256sum(source / "manifest.json") == (source / "READY").read_text().strip()
    tasks = {task["id"]: task for task in manifest["tasks"]}
    conditions, pairs, attempts = [], [], []
    for identifier, task in tasks.items():
        folder = source / "cases" / identifier
        state = read(folder / "state.json")
        attempts.append(
            dict(id=identifier, state=state["state"], job_id=state["job_id"])
        )
        assert (
            state["state"] == "complete"
        ), f"unfinished or rejected pair: {identifier}"
        assert {row["mode"] for row in state["rows"]} == {
            "resident-w4-c0",
            "resident-w4-c4",
        }
        pair = {}
        for row in state["rows"]:
            mode = row["mode"]
            codec = int(mode.rsplit("-c", 1)[1])
            stats = read(folder / mode / "statistics.json")
            metrics = read(folder / mode / "resources.json")
            command = read(folder / mode / "command.json")
            comparison = read(folder / mode / "against_direct.json")
            assert direct_binding(command["task"]) == direct_binding(task)
            assert command["binary"] == manifest["binary"]
            assert command["experimental_environment"] == {
                "VCFCACHE_EXPERIMENT_CODEC_THREADS": str(codec)
            }
            assert stats["effective_engine"] == "native"
            assert stats["workers"] == 4 and stats["lookup_codec_threads"] == codec
            assert stats["cpu_budget"] == 8 and stats["experimental_window_bp"] == 16384
            assert stats["coalesced_chunks"] is False
            assert metrics["exit_code"] == 0 and len(metrics["affinity"]) == 8
            assert comparison["semantic_pass"] and comparison["ignored_fields"] == []
            counts = stats["variant_counts"]
            assert (
                counts["input_variants"]
                == counts["total_output"]
                == task["input"]["records"]
            )
            phases = {p["phase"]: p["wall_seconds"] for p in stats["phase_timings"]}
            annotation = sum(
                p["wall_seconds"]
                for p in stats["stage_timings"]
                if p["stage"] == "annotate-misses"
            )
            value = dict(
                id=identifier,
                tool=task["tool"],
                sample=task["sample"],
                cohort=task["cohort"],
                mode=mode,
                lookup_workers=4,
                codec_threads=codec,
                cpu_budget=8,
                cli_seconds=metrics["wall_seconds"],
                lookup_seconds=phases["lookup"],
                annotation_seconds=annotation,
                cpu_seconds=metrics["cpu_seconds"],
                input_variants=counts["input_variants"],
                missing_variants=counts["missing_variants"],
                candidate_records=stats["candidate_records"],
                lookup_queries=stats["lookup_queries"],
                direct_seconds=state["reused_direct"]["wall_seconds"],
                speedup_over_reused_direct=state["reused_direct"]["wall_seconds"]
                / metrics["wall_seconds"],
                semantic_sha256=comparison["cached_semantic_sha256"],
                peak_gib=metrics["cgroup_after"]["memory_peak_bytes"] / 1024**3,
                job_id=state["job_id"],
                first_mode=state["modes"][0],
            )
            conditions.append(value)
            pair[codec] = value
        control, candidate = pair[4], pair[0]
        for key in [
            "input_variants",
            "missing_variants",
            "candidate_records",
            "lookup_queries",
            "semantic_sha256",
        ]:
            assert control[key] == candidate[key], f"pair changed work/output: {key}"
        pairs.append(
            dict(
                id=identifier,
                tool=task["tool"],
                sample=task["sample"],
                job_id=state["job_id"],
                first_codec_threads=int(state["modes"][0].rsplit("-c", 1)[1]),
                pooled_lookup_seconds=control["lookup_seconds"],
                synchronous_lookup_seconds=candidate["lookup_seconds"],
                lookup_gain=control["lookup_seconds"] / candidate["lookup_seconds"],
                pooled_cli_seconds=control["cli_seconds"],
                synchronous_cli_seconds=candidate["cli_seconds"],
                cli_gain=control["cli_seconds"] / candidate["cli_seconds"],
                pooled_annotation_seconds=control["annotation_seconds"],
                synchronous_annotation_seconds=candidate["annotation_seconds"],
                lookup_seconds_saved=control["lookup_seconds"]
                - candidate["lookup_seconds"],
                cli_seconds_saved=control["cli_seconds"] - candidate["cli_seconds"],
            )
        )
    assert len(conditions) == 2 * len(tasks)
    table(root / "conditions.tsv", conditions)
    table(root / "pairs.tsv", pairs)
    summary = []
    for tool in sorted({p["tool"] for p in pairs}):
        group = [p for p in pairs if p["tool"] == tool]
        summary.append(
            dict(
                tool=tool,
                n=len(group),
                median_lookup_gain=median(p["lookup_gain"] for p in group),
                lookup_wins=sum(p["lookup_gain"] > 1 for p in group),
                median_lookup_seconds_saved=median(
                    p["lookup_seconds_saved"] for p in group
                ),
                median_cli_gain=median(p["cli_gain"] for p in group),
                cli_wins=sum(p["cli_gain"] > 1 for p in group),
                median_cli_seconds_saved=median(p["cli_seconds_saved"] for p in group),
            )
        )
    table(root / "summary.tsv", summary)
    table(root / "attempts.tsv", attempts)
    inventory = [
        dict(path=str(p.relative_to(root)), sha256=sha256sum(p))
        for p in sorted(source.rglob("*"))
        if p.is_file()
    ]
    (root / "EVIDENCE_INVENTORY.json").write_text(
        json.dumps(inventory, indent=2) + "\n"
    )
    print(
        json.dumps(
            dict(pairs=len(pairs), conditions=len(conditions), all_strict_pass=True),
            indent=2,
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    collect(parser.parse_args().root)
