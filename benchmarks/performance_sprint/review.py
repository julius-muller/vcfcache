"""Collect small durable sprint evidence, without annotation outputs or mutation."""

import argparse
import contextlib
import io
import json
import subprocess
from collections import Counter
from pathlib import Path

from benchmarks.engine_profile.deploy import HOST
from benchmarks.engine_profile.review import table
from benchmarks.engine_profile.run import read, verify
from benchmarks.performance_sprint.pilot import report
from benchmarks.run_pilot import REPO_ROOT, sha256sum, write_json_atomic


def snapshot(destination, receipts):
    """Copy immutable evidence to a NEW local snapshot; never modify campaigns."""
    destination.mkdir(parents=True, exist_ok=False)
    rows = []
    for path in receipts:
        receipt = read(path)
        campaign = receipt["campaign"]
        mapped = campaign.replace("/results/", "/mnt/data/slurm-results/", 1)
        local = destination / Path(campaign).name
        local.mkdir()
        # Explicit filters exclude annotation files and dense resource sampling;
        # state, command, comparison, phase and deployment identities remain.
        subprocess.run(
            [
                "rsync",
                "-a",
                "--rsync-path=sudo -u appuser rsync",
                "--exclude=assets/",
                "--exclude=references/",
                "--exclude=resources.jsonl",
                "--include=*/",
                "--include=*.json",
                "--include=*.jsonl",
                "--include=*.tsv",
                "--include=*.log",
                "--include=READY",
                "--exclude=*",
                f"{HOST}:{mapped}/",
                str(local) + "/",
            ],
            check=True,
        )
        verify(local)
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            report(local)
        cells = json.loads(buffer.getvalue())
        write_json_atomic(local / "screening-report.json", cells)
        for cell in cells:
            rows.append(dict(campaign=Path(campaign).name, **cell))
        write_json_atomic(local / "launch-receipt.json", receipt)
    table(destination / "cells.tsv", rows, list(rows[0]) if rows else ["campaign"])
    summary = dict(
        cells=len(rows),
        states=dict(Counter(r["state"] for r in rows)),
        note="Snapshot, not automatic promotion. Single-run reductions require confirmation. Running states may become stale.",
    )
    write_json_atomic(destination / "SUMMARY.json", summary)
    hashes = {
        str(p.relative_to(destination)): sha256sum(p)
        for p in sorted(destination.rglob("*"))
        if p.is_file()
    }
    write_json_atomic(destination / "SHA256.json", hashes)
    return summary


def main():
    """Collect only explicitly recorded launches, with no implicit retries."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    receipts = sorted(
        (REPO_ROOT / "benchmarks/performance_sprint").glob("launch*.json")
    )
    print(json.dumps(snapshot(args.output, receipts), indent=2))


if __name__ == "__main__":
    main()
