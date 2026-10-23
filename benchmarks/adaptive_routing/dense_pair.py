#!/usr/bin/env python3
"""Small paired throughput control against the original fused engine."""

import argparse
import json
import os
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path

from benchmarks.run_pilot import strict_semantic_compare


def main():
    """Alternate arm order on one existing synthetic dataset and fixed affinity."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--original-root", type=Path, required=True)
    parser.add_argument("--original-native", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cpus", type=int, default=14)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    root = Path(__file__).resolve().parents[2]
    affinity = sorted(os.sched_getaffinity(0))[: args.cpus]
    work = Path(tempfile.mkdtemp(prefix="vcfcache-dense-pair-"))
    measurements = []
    for rep in range(3):
        arms = ["original", "adaptive"] if rep % 2 == 0 else ["adaptive", "original"]
        outputs = {}
        for arm in arms:
            destination = work / f"{arm}-{rep}"
            environment = dict(
                os.environ,
                PYTHONPATH=(
                    f"{args.original_root}:{root}" if arm == "original" else str(root)
                ),
                VCFCACHE_ENGINE=str(
                    args.original_native
                    if arm == "original"
                    else root / "native/target/release/vcfcache-engine"
                ),
            )
            subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "benchmarks.fused_engine",
                    "--cell",
                    "--engine",
                    "fused",
                    "--data",
                    str(args.data),
                    "--cell-output",
                    str(destination),
                    "--affinity",
                    *map(str, affinity),
                ],
                cwd=work,
                env=environment,
                check=True,
            )
            value = json.loads((destination / "measurement.json").read_text())
            value.update(arm=arm, rep=rep)
            measurements.append(value)
            outputs[arm] = destination / "output.bcf"
        # Compare the first pair fully; subsequent controls use the same fixed
        # input, commands and binaries and quantify timing variability only.
        if rep == 0:
            audit = strict_semantic_compare(outputs["original"], outputs["adaptive"])
            (args.output / "comparison.json").write_text(
                json.dumps(audit, indent=2) + "\n"
            )
            if not audit["semantic_pass"]:
                raise RuntimeError("dense pair changed complete records")
    old = statistics.median(
        r["wall_seconds"] for r in measurements if r["arm"] == "original"
    )
    new = statistics.median(
        r["wall_seconds"] for r in measurements if r["arm"] == "adaptive"
    )
    summary = dict(
        original_median_seconds=old,
        adaptive_median_seconds=new,
        runtime_ratio=new / old,
        within_10_percent=new / old <= 1.1,
        affinity=affinity,
        work_directory=str(work),
        note="synthetic engine-only control, not CLI time or a real annotator; three paired repetitions",
    )
    (args.output / "measurements.json").write_text(
        json.dumps(measurements, indent=2) + "\n"
    )
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
