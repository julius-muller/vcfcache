#!/usr/bin/env python3
"""Prepare, submit, inspect, and collect VCFcache Slurm cohort tasks."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA_ROOT = Path("/mnt/data/vcfcache_benchmarks")
DEFAULT_QC = DEFAULT_DATA_ROOT / "qc/sample_qc.tsv"
DEFAULT_TASKS = DEFAULT_DATA_ROOT / "manifests/slurm_tasks.tsv"
DEFAULT_RESULTS = Path("/results")
PAIR_SCRIPT = REPO_ROOT / "benchmarks/slurm_pair.sh"
MODES = ("cached", "uncached")


@dataclass(frozen=True)
class CohortTask:
    """One paired sample/replicate benchmark scheduled as a Slurm array task."""

    task_id: int
    sample: str
    input_vcf: str
    replicate: int
    first_mode: str
    second_mode: str
    randomization_key: str


def mode_order(sample: str, replicate: int, seed: str) -> tuple[str, str, str]:
    """Return a deterministic balanced-looking order and its audit key."""
    key = hashlib.sha256(f"{seed}:{sample}:{replicate}".encode()).hexdigest()
    order = MODES if int(key, 16) % 2 == 0 else tuple(reversed(MODES))
    return order[0], order[1], key


def read_eligible_samples(
    qc_path: Path, selected_sample: str | None = None
) -> list[tuple[str, str]]:
    """Read publication-ready 1000 Genomes samples from the frozen QC table."""
    with qc_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    samples = [
        (row["sample"], row["path"])
        for row in rows
        if row["cohort"] == "1000g"
        and row["status"] == "PASS"
        and (selected_sample is None or row["sample"] == selected_sample)
    ]
    if selected_sample and not samples:
        raise ValueError(f"Sample is not publication-ready: {selected_sample}")
    if selected_sample is None and len(samples) != 50:
        raise ValueError(
            f"Expected 50 eligible 1000 Genomes samples, found {len(samples)}"
        )
    return sorted(samples)


def build_tasks(
    qc_path: Path,
    *,
    replicates: int,
    seed: str,
    selected_sample: str | None = None,
) -> list[CohortTask]:
    """Build the stable Slurm array manifest."""
    if replicates < 1:
        raise ValueError("replicates must be positive")
    tasks: list[CohortTask] = []
    for sample, input_vcf in read_eligible_samples(qc_path, selected_sample):
        for replicate in range(1, replicates + 1):
            first, second, key = mode_order(sample, replicate, seed)
            tasks.append(
                CohortTask(
                    task_id=len(tasks),
                    sample=sample,
                    input_vcf=input_vcf,
                    replicate=replicate,
                    first_mode=first,
                    second_mode=second,
                    randomization_key=key,
                )
            )
    return tasks


def write_tasks(path: Path, tasks: Sequence[CohortTask]) -> None:
    """Write the task manifest atomically."""
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_suffix(path.suffix + ".partial")
    with partial.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(asdict(tasks[0]).keys()),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(asdict(task) for task in tasks)
    partial.replace(path)


def prepare(args: argparse.Namespace) -> None:
    """Create and summarize a deterministic task manifest."""
    tasks = build_tasks(
        args.qc,
        replicates=args.replicates,
        seed=args.seed,
        selected_sample=args.sample,
    )
    write_tasks(args.tasks, tasks)
    summary = {
        "task_manifest": str(args.tasks),
        "tasks": len(tasks),
        "samples": len({task.sample for task in tasks}),
        "replicates": args.replicates,
        "seed": args.seed,
    }
    print(json.dumps(summary, indent=2, sort_keys=True))


def submit(args: argparse.Namespace) -> None:
    """Submit the prepared manifest as a bounded Slurm array."""
    if shutil.which("sbatch") is None:
        raise RuntimeError("sbatch is not available")
    with args.tasks.open() as handle:
        task_count = sum(1 for _ in handle) - 1
    if task_count < 1:
        raise RuntimeError("Task manifest is empty")
    if not args.results.is_dir():
        raise RuntimeError(f"Result mount is unavailable: {args.results}")
    command = [
        "sbatch",
        f"--array=0-{task_count - 1}%{args.concurrency}",
        f"--export=ALL,VCFCACHE_TASK_MANIFEST={args.tasks},"
        f"VCFCACHE_RESULTS_ROOT={args.results},VCFCACHE_REPO_ROOT={REPO_ROOT}",
        str(PAIR_SCRIPT),
    ]
    completed = subprocess.run(command, check=True, capture_output=True, text=True)
    print(completed.stdout.strip())


def status(_: argparse.Namespace) -> None:
    """Show active VCFcache benchmark jobs."""
    subprocess.run(
        [
            "squeue",
            "--name=vcfcache-pair",
            "--format=%.18i %.9P %.28j %.8T %.10M %.20R",
        ],
        check=True,
    )


def collect(args: argparse.Namespace) -> None:
    """Report completed task summaries without modifying results."""
    summaries = sorted(args.results.glob("tasks/task-*/**/summary_r*.json"))
    semantic_passes = 0
    for path in summaries:
        value = json.loads(path.read_text())
        if value.get("semantic_comparison", {}).get("semantic_pass") is True:
            semantic_passes += 1
    print(
        json.dumps(
            {
                "completed_tasks": len(summaries),
                "semantic_passes": semantic_passes,
                "results_root": str(args.results),
            },
            indent=2,
            sort_keys=True,
        )
    )


def parser() -> argparse.ArgumentParser:
    """Build the cohort command-line interface."""
    root = argparse.ArgumentParser(description=__doc__)
    subparsers = root.add_subparsers(dest="command", required=True)

    prepare_parser = subparsers.add_parser("prepare")
    prepare_parser.add_argument("--qc", type=Path, default=DEFAULT_QC)
    prepare_parser.add_argument("--tasks", type=Path, default=DEFAULT_TASKS)
    prepare_parser.add_argument("--replicates", type=int, default=3)
    prepare_parser.add_argument("--seed", default="vcfcache-paper-v1")
    prepare_parser.add_argument("--sample")
    prepare_parser.set_defaults(function=prepare)

    submit_parser = subparsers.add_parser("submit")
    submit_parser.add_argument("--tasks", type=Path, default=DEFAULT_TASKS)
    submit_parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    submit_parser.add_argument("--concurrency", type=int, default=6)
    submit_parser.set_defaults(function=submit)

    status_parser = subparsers.add_parser("status")
    status_parser.set_defaults(function=status)

    collect_parser = subparsers.add_parser("collect")
    collect_parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    collect_parser.set_defaults(function=collect)
    return root


def main() -> None:
    """Run the selected cohort action."""
    args = parser().parse_args()
    args.function(args)


if __name__ == "__main__":
    main()
