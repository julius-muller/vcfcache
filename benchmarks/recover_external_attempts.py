#!/usr/bin/env python3
"""Recover fully computed external-WGS tasks rejected by an old comparator."""

from __future__ import annotations

import argparse
import json
import os
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from benchmarks.run_external_task import (
    load_strategies,
    read_task,
    summarize_completed_runs,
)
from benchmarks.run_pilot import git_output, write_json_atomic

CONDITIONS = ("uncached", "gnomad_af_0.1", "gnomad_af_0.01", "cohort_3_genomes")


def _run_directory(
    attempt_run: Path, condition: str, sample: str, replicate: int
) -> Path | None:
    mode = "uncached" if condition == "uncached" else "cached"
    candidates = list(
        (attempt_run / "runs" / condition / "pilot" / sample).glob(
            f"*/{mode}_r{replicate:02d}"
        )
    )
    complete = [
        path
        for path in candidates
        if all(
            required.exists()
            for required in (
                path / "metrics.json",
                path / "output.bcf",
                path / "output.bcf.csi",
            )
        )
    ]
    if len(complete) > 1:
        raise RuntimeError(f"Multiple complete {condition} runs under {attempt_run}")
    return complete[0] if complete else None


def completed_attempt(
    campaign_root: Path, phase: str, task_id: int, sample: str, replicate: int
) -> tuple[Path, dict[str, Path]]:
    """Find the newest failed attempt containing all four completed runs."""
    task_attempts = campaign_root / phase / "attempts" / f"task-{task_id}"
    failures = sorted(
        task_attempts.glob("job-*/failure.json"),
        key=lambda path: int(path.parent.name.removeprefix("job-")),
        reverse=True,
    )
    for failure in failures:
        attempt_run = failure.parent / "run"
        run_dirs = {
            condition: _run_directory(attempt_run, condition, sample, replicate)
            for condition in CONDITIONS
        }
        if all(run_dirs.values()):
            return attempt_run, {
                condition: path
                for condition, path in run_dirs.items()
                if path is not None
            }
    raise RuntimeError(
        f"No failed attempt with four complete runs for {phase} task {task_id}"
    )


def promote_attempt(
    *,
    attempt_run: Path,
    result_dir: Path,
    task_id: int,
    phase: str,
    source_failure: Path,
) -> None:
    """Atomically hard-link one validated attempt into the immutable task area."""
    job_id = os.environ.get("SLURM_JOB_ID", "manual")
    partial = result_dir.with_name(f"{result_dir.name}.partial-recovery-{job_id}")
    if result_dir.exists():
        return
    if partial.exists():
        raise FileExistsError(f"Recovery partial already exists: {partial}")
    result_dir.parent.mkdir(parents=True, exist_ok=True)
    try:
        shutil.copytree(attempt_run, partial, copy_function=os.link)
        write_json_atomic(
            partial / "slurm-task.json",
            {
                "campaign_id": result_dir.parents[2].name,
                "phase": phase,
                "task_id": task_id,
                "slurm_job_id": job_id,
                "slurm_array_job_id": os.environ.get("SLURM_ARRAY_JOB_ID", job_id),
                "slurm_node": os.environ.get("SLURMD_NODENAME", "manual"),
                "recovered": True,
            },
        )
        write_json_atomic(
            partial / "recovery.json",
            {
                "recovered_at": datetime.now(timezone.utc).isoformat(),
                "comparison_commit": git_output("rev-parse", "HEAD"),
                "source_attempt": str(attempt_run),
                "source_failure": str(source_failure),
                "storage": "hardlinks",
            },
        )
        partial.replace(result_dir)
    except BaseException:
        if partial.exists():
            shutil.rmtree(partial)
        raise


def recover(args: argparse.Namespace) -> dict[str, Any]:
    """Revalidate and promote one fully computed failed task."""
    campaign_root = args.campaign_root.resolve()
    result_dir = campaign_root / args.phase / "tasks" / f"task-{args.task_id}"
    if result_dir.exists():
        return {"task_id": args.task_id, "status": "already_complete"}
    manifest = args.task_manifest or campaign_root / "manifests" / f"{args.phase}.tsv"
    strategies_path = args.strategies or campaign_root / "manifests/strategies.json"
    task = read_task(manifest, args.task_id)
    document, strategies = load_strategies(
        strategies_path, task["cohort"], task["assembly"]
    )
    replicate = int(task["replicate"])
    attempt_run, run_dirs = completed_attempt(
        campaign_root, args.phase, args.task_id, task["sample"], replicate
    )
    summary = summarize_completed_runs(
        task=task,
        document=document,
        strategies=strategies,
        execution_order=task["strategy_order"].split(","),
        run_dirs=run_dirs,
        run_root=attempt_run,
        comparison_workers=args.comparison_workers,
    )
    failure = attempt_run.parent / "failure.json"
    promote_attempt(
        attempt_run=attempt_run,
        result_dir=result_dir,
        task_id=args.task_id,
        phase=args.phase,
        source_failure=failure,
    )
    return {
        "task_id": args.task_id,
        "status": "recovered",
        "result_dir": str(result_dir),
        "semantic_passes": sum(row["semantic_pass"] for row in summary["rows"]),
    }


def parser() -> argparse.ArgumentParser:
    """Build the recovery CLI."""
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--campaign-root", required=True, type=Path)
    result.add_argument("--phase", default="warmup")
    result.add_argument("--task-id", required=True, type=int)
    result.add_argument("--task-manifest", type=Path)
    result.add_argument("--strategies", type=Path)
    result.add_argument("--comparison-workers", type=int, default=3)
    return result


def main() -> None:
    """Recover one task and print its machine-readable outcome."""
    print(json.dumps(recover(parser().parse_args()), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
