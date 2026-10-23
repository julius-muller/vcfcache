#!/usr/bin/env python3
"""Recompute CLI-only medians and the experimental store gate from raw evidence."""

import argparse
import csv
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path

import yaml


def collect(root):
    """Prefer the in-step CLI timer; older outer timers include Slurm launch."""
    rows = []
    for measurements in sorted(root.rglob("measurements.json")):
        manifest_path = measurements.parent / "manifest.json"
        if not manifest_path.is_file():
            # Engine-only paired controls use a different schema and timer.
            continue
        manifest = json.loads(manifest_path.read_text())
        for value in json.loads(measurements.read_text()):
            sample = Path(value["input"]).stem
            condition = measurements.parent / f"{sample}-{value['arm']}-{value['rep']}"
            resources = json.loads((condition / "resources.json").read_text())
            stats_paths = list(condition.glob("*_vcstats/compare_stats.yaml"))
            stats = yaml.safe_load(stats_paths[0].read_text()) if stats_paths else {}
            annotation = sum(
                stage["wall_seconds"]
                for stage in stats.get("stage_timings", [])
                if stage["stage"] == "annotate-misses"
            )
            counts = stats.get("variant_counts", {})
            n, miss = counts.get("input_variants"), counts.get("missing_variants")
            peak = (resources.get("cgroup_after") or {}).get("memory.peak")
            rows.append(
                dict(
                    experiment=str(measurements.parent.relative_to(root)),
                    sample=sample,
                    cpus=manifest["cpus"],
                    arm=value["arm"],
                    rep=value["rep"],
                    cli_seconds=resources["wall_seconds"],
                    recipe_seconds=(
                        annotation
                        if value["arm"] not in {"legacy", "original-fused", "direct"}
                        else None
                    ),
                    cache_side_seconds=(
                        resources["wall_seconds"] - annotation
                        if value["arm"] not in {"legacy", "original-fused", "direct"}
                        else None
                    ),
                    cpu_seconds=resources["cpu_seconds"],
                    peak_step_memory_bytes=int(peak) if peak else None,
                    input_records=n,
                    cache_hit_rate=1 - miss / n if n and miss is not None else None,
                    semantic_pass=value["semantic_pass"],
                    source=str(measurements.relative_to(root)),
                )
            )
    return rows


def main():
    """Write a new report directory, refusing to overwrite prior summaries."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = collect(args.archive)
    args.output.mkdir(parents=True, exist_ok=False)
    with (args.output / "conditions.tsv").open("w") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    groups = defaultdict(list)
    for row in rows:
        groups[row["experiment"], row["sample"], row["cpus"], row["arm"]].append(row)
    medians = []
    for (experiment, sample, cpus, arm), values in sorted(groups.items()):
        medians.append(
            dict(
                experiment=experiment,
                sample=sample,
                cpus=cpus,
                arm=arm,
                runs=len(values),
                cli_seconds=statistics.median(r["cli_seconds"] for r in values),
                cache_side_seconds=(
                    statistics.median(r["cache_side_seconds"] for r in values)
                    if values[0]["cache_side_seconds"] is not None
                    else None
                ),
                semantic_pass=all(r["semantic_pass"] for r in values),
            )
        )
    (args.output / "medians.json").write_text(json.dumps(medians, indent=2) + "\n")
    gates = []
    for row in medians:
        if row["arm"] != "store":
            continue
        baseline = next(
            r
            for r in medians
            if r["sample"] == row["sample"]
            and r["experiment"] == row["experiment"]
            and r["cpus"] == row["cpus"]
            and r["arm"] == "auto"
        )
        runtime_gain = baseline["cli_seconds"] / row["cli_seconds"]
        cache_gain = baseline["cache_side_seconds"] / row["cache_side_seconds"]
        store_manifest = next(
            (args.archive / rows_item["source"]).parent.parent / "store_manifest.json"
            for rows_item in rows
            if rows_item["sample"] == row["sample"]
            and rows_item["arm"] == "store"
            and rows_item["experiment"] == row["experiment"]
        )
        build = json.loads(store_manifest.read_text())
        saving = baseline["cli_seconds"] - row["cli_seconds"]
        gates.append(
            dict(
                sample=row["sample"],
                cache_handling_speedup=cache_gain,
                end_to_end_speedup=runtime_gain,
                integration_gate_pass=cache_gain >= 2
                and runtime_gain >= 1.5
                and row["semantic_pass"],
                build_seconds=build["build_seconds"],
                store_bytes=build["payload_bytes"]
                + build["key_bytes"]
                + build["index_bytes"],
                conversion_payback_samples=(
                    math.ceil(build["build_seconds"] / saving) if saving > 0 else None
                ),
            )
        )
    (args.output / "store_gate.json").write_text(json.dumps(gates, indent=2) + "\n")
    print(
        json.dumps(
            dict(
                conditions=len(rows),
                all_semantic_pass=all(r["semantic_pass"] for r in rows),
                store_gate=gates,
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
