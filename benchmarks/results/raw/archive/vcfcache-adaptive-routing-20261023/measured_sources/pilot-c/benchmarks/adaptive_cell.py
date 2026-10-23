#!/usr/bin/env python3
"""Capture resource counters inside one independently allocated benchmark step."""

import argparse
import json
import resource
import shutil
import subprocess
import time
from pathlib import Path


def snapshot():
    """Read current cgroup counters; absence is explicit, never a fabricated zero."""
    memberships = Path("/proc/self/cgroup").read_text().splitlines()
    relative = next((line[3:] for line in memberships if line.startswith("0::")), None)
    if relative is None:
        return None
    root = Path("/sys/fs/cgroup") / relative.lstrip("/")
    return {
        name: (root / name).read_text().strip()
        for name in ("cpu.stat", "io.stat", "memory.peak")
        if (root / name).exists()
    }


def main():
    """Run a complete CLI and save wall/CPU/step-aggregate memory counters."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metrics", type=Path, required=True)
    parser.add_argument("--store", type=Path)
    parser.add_argument("--native", type=Path)
    parser.add_argument("--source", type=Path)
    parser.add_argument("--cache", type=Path)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    before = snapshot()
    start = time.perf_counter()
    selection_seconds = 0
    if args.store:
        # All conversion, candidate construction and temporary indexing for this
        # sample are timed. Store construction itself is a separate one-time cost.
        candidate_root = args.metrics.parent / "candidate"
        candidate = candidate_root / "cache" / args.cache.name
        candidate.mkdir(parents=True)
        (candidate_root / "blueprint").mkdir()
        (candidate_root / "workflow").mkdir()
        shutil.copyfile(
            args.cache.parents[1] / "workflow/init.yaml",
            candidate_root / "workflow/init.yaml",
        )
        for name in ("annotation.yaml", "params.snapshot.yaml"):
            shutil.copyfile(args.cache / name, candidate / name)
        normalized = candidate_root / "input.bcf"
        producer = subprocess.Popen(
            ["bcftools", "norm", "-m-", str(args.source), "-Ou"], stdout=subprocess.PIPE
        )
        try:
            subprocess.run(
                [
                    "bcftools",
                    "view",
                    "-e",
                    'ALT="*"',
                    "-Ob0",
                    "-W",
                    "-o",
                    str(normalized),
                ],
                stdin=producer.stdout,
                check=True,
            )
        finally:
            producer.stdout.close()
            rc = producer.wait()
        if rc:
            raise RuntimeError("store input normalization failed")
        producer = subprocess.Popen(
            [
                str(args.native),
                "store-query",
                "--store",
                str(args.store),
                "--input",
                str(normalized),
            ],
            stdout=subprocess.PIPE,
        )
        try:
            subprocess.run(
                [
                    "bcftools",
                    "view",
                    "-Ob0",
                    "-W",
                    "-o",
                    str(candidate / "vcfcache_annotated.bcf"),
                ],
                stdin=producer.stdout,
                check=True,
            )
        finally:
            producer.stdout.close()
            rc = producer.wait()
        if rc:
            raise RuntimeError("store candidate extraction failed")
        command[command.index("-a") + 1] = str(candidate)
        command[command.index("-i") + 1] = str(normalized)
        command += ["--skip-split-multiallelic"]
        selection_seconds = time.perf_counter() - start
    result = subprocess.run(command, check=False)
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    args.metrics.write_text(
        json.dumps(
            dict(
                wall_seconds=time.perf_counter() - start,
                cpu_seconds=usage.ru_utime + usage.ru_stime,
                largest_child_rss_kib=usage.ru_maxrss,
                cgroup_before=before,
                cgroup_after=snapshot(),
                exit_code=result.returncode,
                store_preparation_seconds=selection_seconds,
            ),
            indent=2,
        )
        + "\n"
    )
    raise SystemExit(result.returncode)


if __name__ == "__main__":
    main()
