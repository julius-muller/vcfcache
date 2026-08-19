#!/usr/bin/env python3
"""Prepare, submit, and collect the controlled runtime Slurm campaign."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence

from benchmarks.prepare_controlled_runtime import HIT_RATES, PIPELINES
from benchmarks.run_cohort import campaign_root, git_output, sha256sum, worker_path
from benchmarks.run_pilot import PUBLICATION_STATISTICS_MODE, write_json_atomic

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONTROLLED_ROOT = Path("/mnt/data/vcfcache_benchmarks/controlled_runtime")
DEFAULT_RESULTS = Path("/mnt/data/slurm-results")
DEFAULT_WORKER_RESULTS = Path("/results")
SLURM_SCRIPT = REPO_ROOT / "benchmarks/slurm_controlled_runtime.sh"
DEFAULT_HIT_RATES = (50, 90)


@dataclass(frozen=True)
class ControlledTask:
    """One immutable controlled-runtime task."""

    task_id: int
    phase: str
    pipeline: str
    mode: str
    target_hit_rate: float
    delay_us: int
    input_vcf: str
    input_sha256: str
    cache_dir: str
    params_file: str
    baseline_result: str


def write_tasks(path: Path, tasks: Sequence[ControlledTask]) -> None:
    """Write a stable controlled task manifest."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(asdict(tasks[0])),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(asdict(task) for task in tasks)


def prepare(args: argparse.Namespace) -> None:
    """Freeze baseline and cached task manifests."""
    selected_hit_rates = tuple(dict.fromkeys(args.hit_rates))
    if not selected_hit_rates:
        raise ValueError("At least one controlled hit rate is required")
    invalid_hit_rates = sorted(set(selected_hit_rates) - set(HIT_RATES))
    if invalid_hit_rates:
        raise ValueError(
            f"Controlled caches are unavailable for hit rates {invalid_hit_rates}; "
            f"choose from {list(HIT_RATES)}"
        )
    controlled = args.controlled_root.resolve()
    ready = json.loads((controlled / "READY.json").read_text())
    root = campaign_root(args.controller_results, args.campaign_id)
    if root.exists() and any(root.iterdir()):
        raise FileExistsError(f"Campaign is nonempty: {root}")
    (root / "logs").mkdir(parents=True, exist_ok=True)
    input_vcf = Path(ready["input_vcf"])
    input_sha = str(ready["input_sha256"])
    baseline_tasks = []
    cached_tasks = []
    baseline_worker_results = {}
    for pipeline, delay in PIPELINES.items():
        cache = controlled / "caches" / pipeline / "hit-050"
        _validate_cache(cache, pipeline, 0.5)
        task_id = len(baseline_tasks)
        baseline_worker_results[pipeline] = (
            args.worker_results
            / "campaigns"
            / args.campaign_id
            / "baseline"
            / "tasks"
            / f"task-{task_id}"
        )
        baseline_tasks.append(
            ControlledTask(
                task_id=task_id,
                phase="baseline",
                pipeline=pipeline,
                mode="uncached",
                target_hit_rate=0.5,
                delay_us=int(delay or 0),
                input_vcf=str(input_vcf),
                input_sha256=input_sha,
                cache_dir=str(cache),
                params_file=str(cache / "params.snapshot.yaml"),
                baseline_result="",
            )
        )
    for pipeline, delay in PIPELINES.items():
        for hit_rate in selected_hit_rates:
            cache = controlled / "caches" / pipeline / f"hit-{hit_rate:03d}"
            _validate_cache(cache, pipeline, hit_rate / 100)
            cached_tasks.append(
                ControlledTask(
                    task_id=len(cached_tasks),
                    phase="cached",
                    pipeline=pipeline,
                    mode="cached",
                    target_hit_rate=hit_rate / 100,
                    delay_us=int(delay or 0),
                    input_vcf=str(input_vcf),
                    input_sha256=input_sha,
                    cache_dir=str(cache),
                    params_file=str(cache / "params.snapshot.yaml"),
                    baseline_result=str(baseline_worker_results[pipeline]),
                )
            )
    for phase, tasks in (("baseline", baseline_tasks), ("cached", cached_tasks)):
        (root / phase / "tasks").mkdir(parents=True)
        (root / phase / "attempts").mkdir(parents=True)
        write_tasks(root / "manifests" / f"{phase}.tsv", tasks)
    metadata = {
        "campaign_id": args.campaign_id,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "commit": git_output("rev-parse", "HEAD"),
        "controlled_root": str(controlled),
        "ready_sha256": sha256sum(controlled / "READY.json"),
        "input_vcf": str(input_vcf),
        "input_sha256": input_sha,
        "pipelines": PIPELINES,
        "hit_rates": list(selected_hit_rates),
        "statistics_mode": PUBLICATION_STATISTICS_MODE,
        "source_files": {
            path.name: sha256sum(path)
            for path in (
                Path(__file__),
                REPO_ROOT / "benchmarks/run_controlled_task.py",
                REPO_ROOT / "benchmarks/run_pilot.py",
                SLURM_SCRIPT,
                REPO_ROOT / "benchmarks/vep_plugins/SyntheticDelay.pm",
            )
        },
        "phases": {"baseline": len(baseline_tasks), "cached": len(cached_tasks)},
    }
    write_json_atomic(root / "campaign.json", metadata)
    print(json.dumps(metadata, indent=2, sort_keys=True))


def _validate_cache(cache: Path, pipeline: str, hit_rate: float) -> None:
    """Reject incomplete or mismatched controlled caches before submission."""
    required = (
        "annotation.yaml",
        "params.snapshot.yaml",
        "vcfcache_annotated.bcf",
        "vcfcache_annotated.bcf.csi",
        "controlled_cache.json",
    )
    missing = [name for name in required if not (cache / name).is_file()]
    if missing:
        raise FileNotFoundError(f"Incomplete controlled cache {cache}: {missing}")
    metadata = json.loads((cache / "controlled_cache.json").read_text())
    if metadata.get("pipeline") != pipeline:
        raise RuntimeError(f"Controlled cache pipeline mismatch: {cache}")
    if float(metadata.get("target_hit_rate", -1)) != hit_rate:
        raise RuntimeError(f"Controlled cache hit-rate mismatch: {cache}")


def _submit(
    args: argparse.Namespace, phase: str, concurrency: int, dependency: str | None
) -> tuple[str, list[str]]:
    root = campaign_root(args.controller_results, args.campaign_id)
    manifest = root / "manifests" / f"{phase}.tsv"
    with manifest.open() as handle:
        count = sum(1 for _ in handle) - 1
    command = [
        "sbatch",
        "--parsable",
        f"--job-name=controlled-{phase}",
        f"--array=0-{count - 1}%{concurrency}",
        f"--chdir={REPO_ROOT}",
        f"--output={args.worker_results}/campaigns/{args.campaign_id}/logs/{phase}-%A_%a.out",
        "--export=ALL,"
        f"VCFCACHE_CAMPAIGN_ID={args.campaign_id},VCFCACHE_PHASE={phase},"
        f"VCFCACHE_TASK_MANIFEST={worker_path(manifest, args.controller_results, args.worker_results)},"
        f"VCFCACHE_RESULTS_ROOT={args.worker_results},VCFCACHE_REPO_ROOT={REPO_ROOT}",
    ]
    if dependency:
        command.append(f"--dependency=afterok:{dependency}")
    command.append(str(SLURM_SCRIPT))
    completed = subprocess.run(command, check=True, capture_output=True, text=True)
    job = completed.stdout.strip().split(";", maxsplit=1)[0]
    if not job.isdigit():
        raise RuntimeError(f"Unexpected sbatch response: {completed.stdout!r}")
    return job, command


def submit(args: argparse.Namespace) -> None:
    """Submit four baselines followed by the selected cached conditions."""
    baseline, baseline_command = _submit(
        args, "baseline", min(4, args.concurrency), args.start_after_job
    )
    cached, cached_command = _submit(args, "cached", args.concurrency, baseline)
    value = {
        "submitted_at": datetime.now(timezone.utc).isoformat(),
        "baseline": {"job_id": baseline, "command": baseline_command},
        "cached": {"job_id": cached, "dependency": baseline, "command": cached_command},
    }
    root = campaign_root(args.controller_results, args.campaign_id)
    write_json_atomic(root / "submission.json", value)
    print(json.dumps(value, indent=2, sort_keys=True))


def collect(args: argparse.Namespace) -> Path:
    """Collect controlled cached summaries into one TSV."""
    root = campaign_root(args.controller_results, args.campaign_id)
    rows: list[dict[str, Any]] = []
    for path in sorted((root / "cached/tasks").glob("task-*/controlled_summary.json")):
        value = json.loads(path.read_text())
        rows.append(
            {
                "pipeline": value["pipeline"],
                "delay_us": value["delay_us"],
                "target_hit_rate": value["target_hit_rate"],
                "observed_hit_rate": value["observed_hit_rate"],
                "cached_wall_seconds": value["cached_wall_seconds"],
                "uncached_wall_seconds": value["uncached_wall_seconds"],
                "relative_runtime": value["relative_runtime"],
                "speedup": value["speedup"],
                "semantic_pass": value["semantic_pass"],
                "statistics_mode": value["metrics"].get("statistics_mode"),
            }
        )
    expected = json.loads((root / "campaign.json").read_text())["phases"]["cached"]
    if len(rows) != expected:
        raise RuntimeError(f"Expected {expected} controlled rows, found {len(rows)}")
    output = root / "publication/controlled_runtime_metrics.tsv"
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(output)
    return output


def _paths(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--campaign-id", required=True)
    parser.add_argument("--controller-results", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--worker-results", type=Path, default=DEFAULT_WORKER_RESULTS)


def parser() -> argparse.ArgumentParser:
    """Build the controlled campaign CLI."""
    result = argparse.ArgumentParser(description=__doc__)
    commands = result.add_subparsers(dest="command", required=True)
    prepare_parser = commands.add_parser("prepare")
    _paths(prepare_parser)
    prepare_parser.add_argument(
        "--controlled-root", type=Path, default=DEFAULT_CONTROLLED_ROOT
    )
    prepare_parser.add_argument(
        "--hit-rates",
        type=int,
        nargs="+",
        choices=HIT_RATES,
        default=DEFAULT_HIT_RATES,
        help="Prepared deterministic hit-rate caches to benchmark (default: 50 90)",
    )
    prepare_parser.set_defaults(function=prepare)
    submit_parser = commands.add_parser("submit")
    _paths(submit_parser)
    submit_parser.add_argument("--concurrency", type=int, default=6)
    submit_parser.add_argument("--start-after-job")
    submit_parser.set_defaults(function=submit)
    collect_parser = commands.add_parser("collect")
    _paths(collect_parser)
    collect_parser.set_defaults(function=collect)
    return result


def main() -> None:
    """Dispatch controlled campaign commands."""
    args = parser().parse_args()
    args.function(args)


if __name__ == "__main__":
    main()
