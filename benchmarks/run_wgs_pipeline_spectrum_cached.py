#!/usr/bin/env python3
"""Run an additional public-cache spectrum without repeating direct cells."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from benchmarks.run_cohort import campaign_root, git_output, worker_path
from benchmarks.run_pilot import (
    PUBLICATION_STATISTICS_MODE,
    sha256sum,
    write_json_atomic,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RESULTS = Path("/mnt/data/slurm-results")
DEFAULT_WORKER_RESULTS = Path("/results")
SLURM_SCRIPT = REPO_ROOT / "benchmarks/slurm_wgs_pipeline_spectrum_cached.sh"
EXPECTED_SOURCE_ALIAS = "cache-gnomad-v4.1-GRCh38-joint-af001-vep115.2-e"
EXPECTED_SOURCE_DOI = "10.5281/zenodo.18190046"
EXPECTED_STRATEGY = "gnomad_af_0.01"


@dataclass(frozen=True)
class CachedSpectrumTask:
    """One new cached timing paired to one immutable direct baseline."""

    task_id: int
    sample: str
    cohort: str
    assembly: str
    pipeline: str
    delay_us: int
    input_records: int
    input_vcf: str
    input_sha256: str
    cache_dir: str
    params_file: str
    expected_hit_rate: float
    baseline_uncached_wall_seconds: float
    baseline_campaign_id: str
    baseline_summary_sha256: str
    source_semantic_campaign_id: str
    source_semantic_task: str
    source_semantic_pass: bool
    source_semantic_comparator: str


def write_tasks(path: Path, tasks: list[CachedSpectrumTask]) -> None:
    """Write a stable task table."""
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


def _source_reference(path: Path, ready: dict[str, Any]) -> dict[str, Any]:
    """Validate the prior direct semantic comparison for this cache and input."""
    summary = json.loads(path.read_text())
    if summary.get("statistics_mode") != PUBLICATION_STATISTICS_MODE:
        raise RuntimeError("Source semantic reference did not use light statistics")
    task = summary.get("task", {})
    if (
        task.get("sample") != ready["sample"]
        or task.get("input_sha256") != ready["input_sha256"]
    ):
        raise RuntimeError("Source semantic reference uses a different input")
    rows = [
        row for row in summary.get("rows", []) if row.get("strategy") == EXPECTED_STRATEGY
    ]
    if len(rows) != 1:
        raise RuntimeError("AF >=1% source semantic row is absent or ambiguous")
    row = rows[0]
    if (
        row.get("bundled_alias") != EXPECTED_SOURCE_ALIAS
        or row.get("zenodo_doi") != EXPECTED_SOURCE_DOI
        or row.get("semantic_pass") is not True
        or row.get("statistics_mode") != PUBLICATION_STATISTICS_MODE
    ):
        raise RuntimeError("AF >=1% source semantic reference is invalid")
    if abs(float(row["cache_hit_rate"]) - float(ready["prior_observed_hit_rate"])) > 1e-12:
        raise RuntimeError("Source semantic reference hit rate changed")
    return {
        "campaign_id": path.parents[3].name,
        "task_path": str(path.parent),
        "summary_sha256": sha256sum(path),
        "row": row,
    }


def prepare(args: argparse.Namespace) -> None:
    """Freeze six cached tasks and their direct/semantic references."""
    spectrum = args.spectrum_root.resolve()
    ready_path = spectrum / "READY.json"
    ready = json.loads(ready_path.read_text())
    if ready.get("source_alias") != EXPECTED_SOURCE_ALIAS:
        raise RuntimeError("Cached spectrum is not based on the bundled AF >=1% cache")
    if ready.get("source_doi") != EXPECTED_SOURCE_DOI:
        raise RuntimeError("Cached spectrum has the wrong Zenodo DOI")
    baseline = args.baseline_campaign_root.resolve()
    baseline_campaign = json.loads((baseline / "campaign.json").read_text())
    if (
        baseline_campaign.get("sample") != ready["sample"]
        or baseline_campaign.get("input_sha256") != ready["input_sha256"]
        or baseline_campaign.get("statistics_mode") != PUBLICATION_STATISTICS_MODE
    ):
        raise RuntimeError("Direct baseline campaign is not comparable")
    baseline_rows: dict[int, tuple[Path, dict[str, Any]]] = {}
    for path in sorted((baseline / "tasks").glob("task-*/wgs_pipeline_summary.json")):
        value = json.loads(path.read_text())
        if value.get("semantic_pass") is not True:
            raise RuntimeError(f"Baseline semantic failure: {path}")
        baseline_rows[int(value["delay_us"])] = (path, value)
    expected_delays = {int(value) for value in ready["delays_us"]}
    if set(baseline_rows) != expected_delays:
        raise RuntimeError("Direct baseline does not contain the same six loads")
    source = _source_reference(args.source_semantic_summary.resolve(), ready)

    root = campaign_root(args.controller_results, args.campaign_id)
    if root.exists() and any(root.iterdir()):
        raise FileExistsError(f"Campaign is nonempty: {root}")
    (root / "logs").mkdir(parents=True)
    (root / "tasks").mkdir()
    (root / "attempts").mkdir()
    tasks: list[CachedSpectrumTask] = []
    for task_id, pipeline in enumerate(ready["pipelines"]):
        delay = int(pipeline["delay_us"])
        baseline_path, baseline_row = baseline_rows[delay]
        cache = Path(pipeline["cache_dir"])
        tasks.append(
            CachedSpectrumTask(
                task_id=task_id,
                sample=ready["sample"],
                cohort=ready["cohort"],
                assembly=ready["assembly"],
                pipeline=pipeline["name"],
                delay_us=delay,
                input_records=int(ready["input_records"]),
                input_vcf=ready["input_vcf"],
                input_sha256=ready["input_sha256"],
                cache_dir=str(cache),
                params_file=str(cache / "params.snapshot.yaml"),
                expected_hit_rate=float(ready["prior_observed_hit_rate"]),
                baseline_uncached_wall_seconds=float(
                    baseline_row["uncached_wall_seconds"]
                ),
                baseline_campaign_id=baseline_campaign["campaign_id"],
                baseline_summary_sha256=sha256sum(baseline_path),
                source_semantic_campaign_id=source["campaign_id"],
                source_semantic_task=source["task_path"],
                source_semantic_pass=source["row"]["semantic_pass"],
                source_semantic_comparator=source["row"]["semantic_comparator"],
            )
        )
    write_tasks(root / "manifests/tasks.tsv", tasks)
    metadata = {
        "campaign_id": args.campaign_id,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "commit": git_output("rev-parse", "HEAD"),
        "statistics_mode": PUBLICATION_STATISTICS_MODE,
        "sample": ready["sample"],
        "cohort": ready["cohort"],
        "assembly": ready["assembly"],
        "input_vcf": ready["input_vcf"],
        "input_sha256": ready["input_sha256"],
        "input_records": ready["input_records"],
        "source_alias": ready["source_alias"],
        "source_doi": ready["source_doi"],
        "prior_observed_hit_rate": ready["prior_observed_hit_rate"],
        "delays_us": ready["delays_us"],
        "task_count": len(tasks),
        "new_timed_cells": len(tasks),
        "reused_direct_cells": len(tasks),
        "technical_repeats": 0,
        "baseline_campaign": {
            "campaign_id": baseline_campaign["campaign_id"],
            "commit": baseline_campaign["commit"],
            "root": str(baseline),
        },
        "source_semantic_reference": source,
        "semantic_design": {
            "kind": "composed_reference_plus_cross_delay_fingerprint",
            "source_cache_direct_comparison": True,
            "plugin_emits_annotations": False,
            "new_outputs_require_identical_canonical_fingerprints": True,
        },
        "source_files": {
            path.name: sha256sum(path)
            for path in (
                Path(__file__),
                REPO_ROOT / "benchmarks/run_wgs_pipeline_spectrum_task.py",
                REPO_ROOT / "benchmarks/run_pilot.py",
                REPO_ROOT / "benchmarks/vep_plugins/SyntheticDelay.pm",
                SLURM_SCRIPT,
                ready_path,
            )
        },
    }
    write_json_atomic(root / "campaign.json", metadata)
    print(json.dumps(metadata, indent=2, sort_keys=True))


def submit(args: argparse.Namespace) -> None:
    """Submit the six new cached timings concurrently."""
    root = campaign_root(args.controller_results, args.campaign_id)
    manifest = root / "manifests/tasks.tsv"
    count = sum(1 for _ in manifest.open()) - 1
    command = [
        "sbatch",
        "--parsable",
        "--job-name=wgs-spectrum-af001",
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
    """Validate fingerprints and combine new timings with frozen baselines."""
    root = campaign_root(args.controller_results, args.campaign_id)
    campaign = json.loads((root / "campaign.json").read_text())
    summaries = []
    for path in sorted((root / "tasks").glob("task-*/wgs_pipeline_cached_summary.json")):
        summaries.append(json.loads(path.read_text()))
    if len(summaries) != campaign["task_count"]:
        raise RuntimeError(
            f"Expected {campaign['task_count']} cached rows, found {len(summaries)}"
        )
    digests = {
        value["semantic_fingerprint"]["semantic_sha256"] for value in summaries
    }
    headers = {
        value["semantic_fingerprint"]["csq_header_sha256"] for value in summaries
    }
    fingerprint_records = {
        int(value["semantic_fingerprint"]["records"]) for value in summaries
    }
    semantic_pass = (
        len(digests) == 1
        and len(headers) == 1
        and fingerprint_records == {int(campaign["input_records"])}
        and all(value["source_semantic_pass"] is True for value in summaries)
    )
    validation = {
        "semantic_pass": semantic_pass,
        "comparator": "composed_reference_plus_cross_delay_fingerprint_v1",
        "canonical_semantic_sha256": next(iter(digests)) if len(digests) == 1 else None,
        "csq_header_sha256": next(iter(headers)) if len(headers) == 1 else None,
        "records": next(iter(fingerprint_records)) if len(fingerprint_records) == 1 else None,
        "new_cached_outputs": len(summaries),
        "source_semantic_reference": campaign["source_semantic_reference"],
        "baseline_campaign": campaign["baseline_campaign"],
    }
    write_json_atomic(root / "publication/semantic_validation.json", validation)
    if not semantic_pass:
        raise RuntimeError(f"Composed semantic validation failed: {validation}")
    rows = []
    for value in sorted(summaries, key=lambda row: int(row["delay_us"])):
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
                "semantic_pass": True,
                "semantic_comparator": validation["comparator"],
                "statistics_mode": value["statistics_mode"],
                "baseline_reused": True,
                "baseline_campaign_id": value["baseline_campaign_id"],
                "source_alias": campaign["source_alias"],
                "source_doi": campaign["source_doi"],
            }
        )
    output = root / "publication/wgs_pipeline_spectrum_af001.tsv"
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(output)
    return output


def add_paths(parser: argparse.ArgumentParser) -> None:
    """Add paths shared by all cached-extension commands."""
    parser.add_argument("--campaign-id", required=True)
    parser.add_argument("--controller-results", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--worker-results", type=Path, default=DEFAULT_WORKER_RESULTS)


def main() -> None:
    """Dispatch cached-extension campaign commands."""
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    prepare_parser = commands.add_parser("prepare")
    add_paths(prepare_parser)
    prepare_parser.add_argument("--spectrum-root", type=Path, required=True)
    prepare_parser.add_argument("--baseline-campaign-root", type=Path, required=True)
    prepare_parser.add_argument("--source-semantic-summary", type=Path, required=True)
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
