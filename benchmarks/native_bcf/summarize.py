#!/usr/bin/env python3
"""Recompute native-pilot tables, retaining failures and experiment boundaries."""

import argparse
import csv
import json
import runpy
import statistics
from collections import defaultdict
from pathlib import Path

collect = runpy.run_path(
    str(Path(__file__).parents[1] / "adaptive_routing/summarize.py")
)["collect"]


def main():
    """Write a fresh summary without modifying raw measurements or old figures."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = collect(args.archive)
    args.output.mkdir(parents=True, exist_ok=False)
    groups = defaultdict(list)
    for row in rows:
        experiment = row["experiment"].split("/")[0]
        row["role"] = (
            "diagnostic_failed_pilot"
            if experiment == "fastvep-1270"
            else (
                "calibration"
                if experiment.startswith("calibration-")
                else (
                    "control"
                    if experiment.startswith("controls-")
                    else (
                        "final_source_check"
                        if experiment.startswith("final-")
                        else "primary"
                    )
                )
            )
        )
        groups[row["experiment"], row["sample"], row["arm"]].append(row)
    with (args.output / "conditions.tsv").open("w") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    medians = []
    for (experiment, sample, arm), values in sorted(groups.items()):
        medians.append(
            dict(
                experiment=experiment,
                sample=sample,
                arm=arm,
                runs=len(values),
                cli_seconds=statistics.median(v["cli_seconds"] for v in values),
                cache_side_seconds=(
                    statistics.median(v["cache_side_seconds"] for v in values)
                    if values[0]["cache_side_seconds"] is not None
                    else None
                ),
                semantic_pass=all(v["semantic_pass"] for v in values),
            )
        )
    (args.output / "medians.json").write_text(json.dumps(medians, indent=2) + "\n")
    outcomes = dict(
        conditions=len(rows),
        failed_semantic_conditions=[
            {
                k: v
                for k, v in row.items()
                if k in {"experiment", "sample", "arm", "semantic_pass"}
            }
            for row in rows
            if not row["semantic_pass"]
        ],
        note="Clinical mismatch is retained, not excluded or converted into a speedup claim. Distinct snapshots/hardware/input scopes are not pooled.",
    )
    (args.output / "outcomes.json").write_text(json.dumps(outcomes, indent=2) + "\n")
    print(json.dumps(outcomes, indent=2))


if __name__ == "__main__":
    main()
