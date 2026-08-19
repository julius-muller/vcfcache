#!/usr/bin/env python3
"""Prepare, submit, and collect the real-WGS virtual pipeline spectrum."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from benchmarks.run_cohort import campaign_root, git_output, sha256sum, worker_path
from benchmarks.run_pilot import PUBLICATION_STATISTICS_MODE, write_json_atomic

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SPECTRUM_ROOT = Path(
    "/mnt/data/vcfcache_benchmarks/wgs_pipeline_spectrum"
)
DEFAULT_RESULTS = Path("/mnt/data/slurm-results")
DEFAULT_WORKER_RESULTS = Path("/results")
SLURM_SCRIPT = REPO_ROOT / "benchmarks/slurm_wgs_pipeline_spectrum.sh"


@dataclass(frozen=True)
class SpectrumTask:
    """One paired real-WGS virtual pipeline load."""

    task_id: int
    pipeline: str
    delay_us: int
    input_vcf: str
    input_sha256: str
    cache_dir: str
    params_file: str
    first_mode: str
    second_mode: str


def write_tasks(path: Path, tasks: list[SpectrumTask]) -> None:
    """Write the stable task table."""
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
    """Freeze six paired virtual-load tasks for one real WGS."""
    spectrum = args.spectrum_root.resolve()
    ready = json.loads((spectrum / "READY.json").read_text())
    root = campaign_root(args.controller_results, args.campaign_id)
    if root.exists() and any(root.iterdir()):
        raise FileExistsError(f"Campaign is nonempty: {root}")
    (root / "logs").mkdir(parents=True)
    (root / "tasks").mkdir()
    (root / "attempts").mkdir()
    tasks = []
    for task_id, pipeline in enumerate(ready["pipelines"]):
        cache = Path(pipeline["cache_dir"])
        first, second = (
            ("cached", "uncached")
            if task_id % 2 == 0
            else ("uncached", "cached")
        )
        tasks.append(
            SpectrumTask(
                task_id=task_id,
                pipeline=pipeline["name"],
                delay_us=int(pipeline["delay_us"]),
                input_vcf=ready["input_vcf"],
                input_sha256=ready["input_sha256"],
                cache_dir=str(cache),
                params_file=str(cache / "params.snapshot.yaml"),
                first_mode=first,
                second_mode=second,
            )
        )
    write_tasks(root / "manifests/tasks.tsv", tasks)
    metadata = {
        "campaign_id": args.campaign_id,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "commit": git_output("rev-parse", "HEAD"),
        "statistics_mode": PUBLICATION_STATISTICS_MODE,
        "spectrum_root": str(spectrum),
        "ready_sha256": sha256sum(spectrum / "READY.json"),
        "input_vcf": ready["input_vcf"],
        "input_sha256": ready["input_sha256"],
        "input_records": ready["input_records"],
        "sample": ready["sample"],
        "cohort": ready["cohort"],
        "assembly": ready["assembly"],
        "source_alias": ready["source_alias"],
        "source_doi": ready["source_doi"],
        "prior_observed_hit_rate": ready["prior_observed_hit_rate"],
        "delays_us": ready["delays_us"],
        "task_count": len(tasks),
        "technical_repeats": 0,
        "source_files": {
            path.name: sha256sum(path)
            for path in (
                Path(__file__),
                REPO_ROOT / "benchmarks/run_wgs_pipeline_spectrum_task.py",
                REPO_ROOT / "benchmarks/run_pilot.py",
                REPO_ROOT / "benchmarks/vep_plugins/SyntheticDelay.pm",
                SLURM_SCRIPT,
            )
        },
    }
    write_json_atomic(root / "campaign.json", metadata)
    print(json.dumps(metadata, indent=2, sort_keys=True))


def submit(args: argparse.Namespace) -> None:
    """Submit all six virtual loads concurrently on exclusive nodes."""
    root = campaign_root(args.controller_results, args.campaign_id)
    manifest = root / "manifests/tasks.tsv"
    count = sum(1 for _ in manifest.open()) - 1
    command = [
        "sbatch",
        "--parsable",
        "--job-name=wgs-pipeline-spectrum",
        f"--array=0-{count - 1}%{count}",
        f"--chdir={REPO_ROOT}",
        f"--output={args.worker_results}/campaigns/{args.campaign_id}/logs/task-%A_%a.out",
        "--export=ALL,"
        f"VCFCACHE_CAMPAIGN_ID={args.campaign_id},"
        f"VCFCACHE_TASK_MANIFEST={worker_path(manifest, args.controller_results, args.worker_results)},"
        f"VCFCACHE_RESULTS_ROOT={args.worker_results},VCFCACHE_REPO_ROOT={REPO_ROOT}",
        str(SLURM_SCRIPT),
    ]
    completed = subprocess.run(command, check=True, capture_output=True, text=True)
    job_id = completed.stdout.strip().split(";", maxsplit=1)[0]
    if not job_id.isdigit():
        raise RuntimeError(f"Unexpected sbatch response: {completed.stdout!r}")
    value = {
        "submitted_at": datetime.now(timezone.utc).isoformat(),
        "job_id": job_id,
        "command": command,
    }
    write_json_atomic(root / "submission.json", value)
    print(json.dumps(value, indent=2, sort_keys=True))


def collect(args: argparse.Namespace) -> Path:
    """Collect the six validated paired summaries into one table."""
    root = campaign_root(args.controller_results, args.campaign_id)
    rows: list[dict[str, Any]] = []
    for path in sorted((root / "tasks").glob("task-*/wgs_pipeline_summary.json")):
        value = json.loads(path.read_text())
        rows.append(
            {
                "sample": value["sample"],
                "cohort": value["cohort"],
                "assembly": value["assembly"],
                "pipeline": value["pipeline"],
                "delay_us": value["delay_us"],
                "input_records": value["input_records"],
                "cache_hit_rate": value["cache_hit_rate"],
                "uncached_wall_seconds": value["uncached_wall_seconds"],
                "cached_wall_seconds": value["cached_wall_seconds"],
                "relative_runtime": value["relative_runtime"],
                "speedup": value["speedup"],
                "time_saved_seconds": value["time_saved_seconds"],
                "semantic_pass": value["semantic_pass"],
                "semantic_comparator": value["semantic_comparator"],
                "statistics_mode": value["statistics_mode"],
                "first_mode": value["first_mode"],
                "second_mode": value["second_mode"],
            }
        )
    expected = json.loads((root / "campaign.json").read_text())["task_count"]
    if len(rows) != expected:
        raise RuntimeError(f"Expected {expected} WGS spectrum rows, found {len(rows)}")
    output = root / "publication/wgs_pipeline_spectrum.tsv"
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(output)
    return output


def add_paths(parser: argparse.ArgumentParser) -> None:
    """Add shared campaign path options."""
    parser.add_argument("--campaign-id", required=True)
    parser.add_argument("--controller-results", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--worker-results", type=Path, default=DEFAULT_WORKER_RESULTS)


def main() -> None:
    """Dispatch spectrum campaign commands."""
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    prepare_parser = commands.add_parser("prepare")
    add_paths(prepare_parser)
    prepare_parser.add_argument(
        "--spectrum-root", type=Path, default=DEFAULT_SPECTRUM_ROOT
    )
    prepare_parser.set_defaults(function=prepare)
    submit_parser = commands.add_parser("submit")
    add_paths(submit_parser)
    submit_parser.set_defaults(function=submit)
    collect_parser = commands.add_parser("collect")
    add_paths(collect_parser)
    collect_parser.set_defaults(function=collect)
    args = parser.parse_args()
    args.function(args)


if __name__ == "__main__":
    main()
