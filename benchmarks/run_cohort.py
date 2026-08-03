#!/usr/bin/env python3
"""Prepare, submit, inspect, and collect VCFcache Slurm campaigns."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Sequence

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA_ROOT = Path("/mnt/data/vcfcache_benchmarks")
DEFAULT_QC = DEFAULT_DATA_ROOT / "qc/sample_qc.tsv"
DEFAULT_CONTROLLER_RESULTS = Path("/mnt/data/slurm-results")
DEFAULT_WORKER_RESULTS = Path("/results")
PAIR_SCRIPT = REPO_ROOT / "benchmarks/slurm_pair.sh"
MODES = ("cached", "uncached")
PHASES = ("smoke", "warmup", "measured")
DEFAULT_SEED = "vcfcache-paper-primary-wgs-v1"
DEFAULT_MEASURED_REPLICATES = 1


@dataclass(frozen=True)
class CohortTask:
    """One paired sample/replicate benchmark scheduled as an array task."""

    task_id: int
    phase: str
    measured: str
    sample: str
    population: str
    superpopulation: str
    sex: str
    input_vcf: str
    input_records: int
    input_sha256: str
    replicate: int
    first_mode: str
    second_mode: str
    randomization_key: str


def git_output(*args: str) -> str:
    """Return one Git command's stripped stdout."""
    return subprocess.run(
        ["git", "-C", REPO_ROOT, *args],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def sha256sum(path: Path) -> str:
    """Return a streaming SHA-256 digest."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 << 20):
            digest.update(chunk)
    return digest.hexdigest()


def write_json_atomic(path: Path, value: object) -> None:
    """Write JSON without exposing a partial file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_suffix(path.suffix + ".partial")
    partial.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    partial.replace(path)


def read_eligible_samples(
    qc_path: Path, selected_sample: str | None = None
) -> list[dict[str, str]]:
    """Read publication-ready 1000 Genomes samples from frozen QC."""
    with qc_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    samples = [
        row
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
    return sorted(samples, key=lambda row: row["sample"])


def _sample_patterns(samples: Sequence[dict[str, str]], seed: str) -> dict[str, str]:
    """Assign exactly half of a full cohort to each alternating order."""
    ranked = sorted(
        samples,
        key=lambda row: hashlib.sha256(
            f"{seed}:order-block:{row['sample']}".encode()
        ).hexdigest(),
    )
    split = len(ranked) // 2
    return {
        row["sample"]: ("cached" if index < split else "uncached")
        for index, row in enumerate(ranked)
    }


def mode_order(
    sample: str,
    replicate: int,
    seed: str,
    *,
    initial_mode: str | None = None,
) -> tuple[str, str, str]:
    """Return a deterministic alternating order and its audit key."""
    key = hashlib.sha256(f"{seed}:{sample}:{replicate}".encode()).hexdigest()
    initial = initial_mode or MODES[int(key, 16) % 2]
    first = (
        initial
        if replicate % 2 == 1
        else next(mode for mode in MODES if mode != initial)
    )
    second = next(mode for mode in MODES if mode != first)
    return first, second, key


def build_tasks(
    qc_path: Path,
    *,
    phase: str = "measured",
    replicates: int = 3,
    seed: str = DEFAULT_SEED,
    selected_sample: str | None = None,
) -> list[CohortTask]:
    """Build one stable campaign-phase task manifest."""
    if phase not in PHASES:
        raise ValueError(f"Unknown phase: {phase}")
    if replicates < 1:
        raise ValueError("replicates must be positive")
    samples = read_eligible_samples(qc_path, selected_sample)
    patterns = _sample_patterns(samples, seed)
    tasks: list[CohortTask] = []
    for row in samples:
        for replicate in range(1, replicates + 1):
            first, second, key = mode_order(
                row["sample"],
                replicate,
                seed,
                initial_mode=patterns[row["sample"]],
            )
            tasks.append(
                CohortTask(
                    task_id=len(tasks),
                    phase=phase,
                    measured="true" if phase == "measured" else "false",
                    sample=row["sample"],
                    population=row["population"],
                    superpopulation=row["superpopulation"],
                    sex=row["sex"],
                    input_vcf=row["path"],
                    input_records=int(row["records"]),
                    input_sha256=row["sha256"],
                    replicate=replicate,
                    first_mode=first,
                    second_mode=second,
                    randomization_key=key,
                )
            )
    return tasks


def write_tasks(path: Path, tasks: Sequence[CohortTask]) -> None:
    """Write a nonempty task manifest atomically."""
    if not tasks:
        raise ValueError("Task manifest must not be empty")
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


def campaign_root(controller_results: Path, campaign_id: str) -> Path:
    """Return a validated controller-side campaign path."""
    if not campaign_id or any(
        char not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_."
        for char in campaign_id
    ):
        raise ValueError(
            "Campaign ID may contain only letters, digits, dot, dash, and underscore"
        )
    return controller_results / "campaigns" / campaign_id


def worker_path(
    controller_path: Path, controller_results: Path, worker_results: Path
) -> Path:
    """Translate a controller export path to its worker mount path."""
    return worker_results / controller_path.relative_to(controller_results)


def prepare_campaign(args: argparse.Namespace) -> None:
    """Freeze smoke, diagnostic warm-up, and measured manifests.

    The publication design uses one paired measurement per biological sample.
    The warm-up manifest is retained for diagnostics and compatibility with
    already archived campaigns, but the default submission chain does not run
    it. A full-cohort warm-up would execute the same scientific comparison a
    second time without adding an independent sample.
    """
    root = campaign_root(args.controller_results, args.campaign_id)
    if root.exists() and any(root.iterdir()):
        raise FileExistsError(f"Campaign already exists and is nonempty: {root}")
    manifests = root / "manifests"
    specs = {
        "smoke": (args.smoke_sample, 1),
        "warmup": (None, 1),
        "measured": (None, DEFAULT_MEASURED_REPLICATES),
    }
    phase_values: dict[str, object] = {}
    for phase, (sample, replicates) in specs.items():
        # Create shared parents once on the controller. Concurrent mkdir -p
        # calls from different workers can race on the NFSv4 export.
        (root / phase / "tasks").mkdir(parents=True, exist_ok=True)
        (root / phase / "attempts").mkdir(parents=True, exist_ok=True)
        tasks = build_tasks(
            args.qc,
            phase=phase,
            replicates=replicates,
            seed=args.seed,
            selected_sample=sample,
        )
        path = manifests / f"{phase}.tsv"
        write_tasks(path, tasks)
        phase_values[phase] = {
            "tasks": len(tasks),
            "samples": len({task.sample for task in tasks}),
            "replicates": replicates,
            "manifest": str(path),
            "manifest_sha256": sha256sum(path),
        }
    metadata = {
        "campaign_id": args.campaign_id,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "commit": git_output("rev-parse", "HEAD"),
        "commit_short": git_output("rev-parse", "--short=12", "HEAD"),
        "tracked_tree_clean": subprocess.run(
            ["git", "-C", REPO_ROOT, "diff", "--quiet"], check=False
        ).returncode
        == 0,
        "seed": args.seed,
        "qc": str(args.qc),
        "qc_sha256": sha256sum(args.qc),
        "controller_results": str(args.controller_results),
        "worker_results": str(args.worker_results),
        "phases": phase_values,
    }
    if not metadata["tracked_tree_clean"]:
        raise RuntimeError("Tracked worktree must be clean before campaign preparation")
    write_json_atomic(root / "campaign.json", metadata)
    print(json.dumps(metadata, indent=2, sort_keys=True))


def _task_count(path: Path) -> int:
    with path.open() as handle:
        count = sum(1 for _ in handle) - 1
    if count < 1:
        raise RuntimeError(f"Task manifest is empty: {path}")
    return count


def submit_phase(
    *,
    campaign_id: str,
    phase: str,
    controller_results: Path,
    worker_results: Path,
    concurrency: int,
    dependency: str | None = None,
    task_ids: str | None = None,
) -> tuple[str, list[str]]:
    """Submit one campaign phase and return its parsable Slurm job ID."""
    if shutil.which("sbatch") is None:
        raise RuntimeError("sbatch is not available")
    root = campaign_root(controller_results, campaign_id)
    manifest = root / "manifests" / f"{phase}.tsv"
    task_count = _task_count(manifest)
    (root / "logs").mkdir(parents=True, exist_ok=True)
    worker_manifest = worker_path(manifest, controller_results, worker_results)
    array = task_ids or f"0-{task_count - 1}"
    command = [
        "sbatch",
        "--parsable",
        f"--job-name=vcfcache-{phase}",
        f"--array={array}%{concurrency}",
        f"--chdir={REPO_ROOT}",
        f"--output={worker_results}/campaigns/{campaign_id}/logs/{phase}-%A_%a.out",
        f"--export=ALL,VCFCACHE_CAMPAIGN_ID={campaign_id},"
        f"VCFCACHE_PHASE={phase},VCFCACHE_TASK_MANIFEST={worker_manifest},"
        f"VCFCACHE_RESULTS_ROOT={worker_results},VCFCACHE_REPO_ROOT={REPO_ROOT}",
    ]
    if dependency:
        command.append(f"--dependency=afterok:{dependency}")
    command.append(str(PAIR_SCRIPT))
    completed = subprocess.run(command, check=True, capture_output=True, text=True)
    job_id = completed.stdout.strip().split(";", maxsplit=1)[0]
    if not job_id.isdigit():
        raise RuntimeError(f"Unexpected sbatch response: {completed.stdout!r}")
    return job_id, command


def submit_chain(args: argparse.Namespace) -> None:
    """Submit one validation smoke followed by one measurement per sample."""
    root = campaign_root(args.controller_results, args.campaign_id)
    campaign = json.loads((root / "campaign.json").read_text())
    current_commit = git_output("rev-parse", "HEAD")
    if campaign["commit"] != current_commit:
        raise RuntimeError("Campaign commit does not match the checked-out commit")
    submissions: dict[str, object] = {}
    smoke_id, smoke_command = submit_phase(
        campaign_id=args.campaign_id,
        phase="smoke",
        controller_results=args.controller_results,
        worker_results=args.worker_results,
        concurrency=1,
    )
    submissions["smoke"] = {"job_id": smoke_id, "command": smoke_command}
    measured_id, measured_command = submit_phase(
        campaign_id=args.campaign_id,
        phase="measured",
        controller_results=args.controller_results,
        worker_results=args.worker_results,
        concurrency=args.concurrency,
        dependency=smoke_id,
    )
    submissions["measured"] = {
        "job_id": measured_id,
        "dependency": smoke_id,
        "command": measured_command,
    }
    value = {
        "campaign_id": args.campaign_id,
        "submitted_at": datetime.now(timezone.utc).isoformat(),
        "commit": current_commit,
        "hostname": os.uname().nodename,
        "concurrency": args.concurrency,
        "submissions": submissions,
    }
    write_json_atomic(root / "submission.json", value)
    print(json.dumps(value, indent=2, sort_keys=True))


def submit(args: argparse.Namespace) -> None:
    """Submit or retry one phase, optionally with sparse array task IDs."""
    job_id, command = submit_phase(
        campaign_id=args.campaign_id,
        phase=args.phase,
        controller_results=args.controller_results,
        worker_results=args.worker_results,
        concurrency=args.concurrency,
        dependency=args.dependency,
        task_ids=args.task_ids,
    )
    print(json.dumps({"job_id": job_id, "command": command}, indent=2))


def status(args: argparse.Namespace) -> None:
    """Show active jobs for one campaign's recorded submission chain."""
    root = campaign_root(args.controller_results, args.campaign_id)
    submission = json.loads((root / "submission.json").read_text())
    job_ids = [str(value["job_id"]) for value in submission["submissions"].values()]
    subprocess.run(
        [
            "squeue",
            "--jobs",
            ",".join(job_ids),
            "--format=%.18i %.24j %.8T %.10M %.24R",
        ],
        check=True,
    )


def collect(args: argparse.Namespace) -> None:
    """Report archived task summaries without modifying results."""
    root = campaign_root(args.controller_results, args.campaign_id)
    phases: dict[str, object] = {}
    for phase in PHASES:
        summaries = sorted((root / phase / "tasks").glob("task-*/**/summary_r*.json"))
        passes = 0
        for path in summaries:
            value = json.loads(path.read_text())
            if value.get("semantic_comparison", {}).get("semantic_pass") is True:
                passes += 1
        phases[phase] = {
            "completed_tasks": len(summaries),
            "semantic_passes": passes,
        }
    print(
        json.dumps(
            {
                "campaign_id": args.campaign_id,
                "campaign_root": str(root),
                "phases": phases,
            },
            indent=2,
            sort_keys=True,
        )
    )


def _add_campaign_paths(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--campaign-id", required=True)
    parser.add_argument(
        "--controller-results", type=Path, default=DEFAULT_CONTROLLER_RESULTS
    )
    parser.add_argument("--worker-results", type=Path, default=DEFAULT_WORKER_RESULTS)


def parser() -> argparse.ArgumentParser:
    """Build the cohort command-line interface."""
    root = argparse.ArgumentParser(description=__doc__)
    subparsers = root.add_subparsers(dest="command", required=True)

    prepare_parser = subparsers.add_parser("prepare")
    _add_campaign_paths(prepare_parser)
    prepare_parser.add_argument("--qc", type=Path, default=DEFAULT_QC)
    prepare_parser.add_argument("--seed", default=DEFAULT_SEED)
    prepare_parser.add_argument("--smoke-sample", default="HG02079")
    prepare_parser.set_defaults(function=prepare_campaign)

    chain_parser = subparsers.add_parser("submit-chain")
    _add_campaign_paths(chain_parser)
    chain_parser.add_argument("--concurrency", type=int, default=6)
    chain_parser.set_defaults(function=submit_chain)

    submit_parser = subparsers.add_parser("submit")
    _add_campaign_paths(submit_parser)
    submit_parser.add_argument("--phase", choices=PHASES, required=True)
    submit_parser.add_argument("--concurrency", type=int, default=6)
    submit_parser.add_argument("--dependency")
    submit_parser.add_argument("--task-ids")
    submit_parser.set_defaults(function=submit)

    status_parser = subparsers.add_parser("status")
    _add_campaign_paths(status_parser)
    status_parser.set_defaults(function=status)

    collect_parser = subparsers.add_parser("collect")
    _add_campaign_paths(collect_parser)
    collect_parser.set_defaults(function=collect)
    return root


def main() -> None:
    """Run the selected cohort action."""
    args = parser().parse_args()
    if getattr(args, "concurrency", 1) < 1:
        raise ValueError("concurrency must be positive")
    args.function(args)


if __name__ == "__main__":
    main()
