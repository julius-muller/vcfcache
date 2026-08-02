#!/usr/bin/env python3
"""Prepare, submit, inspect, and collect external-WGS Slurm campaigns."""

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
from typing import Any, Sequence

from benchmarks.run_cohort import campaign_root, git_output, sha256sum, worker_path
from benchmarks.run_pilot import write_json_atomic
from benchmarks.run_strategy_comparison import (
    BUNDLED_CACHE_SPECS,
    VEP_CACHE_NAME,
    public_strategies,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA_ROOT = Path("/mnt/data/vcfcache_benchmarks")
DEFAULT_EXTERNAL_ROOT = DEFAULT_DATA_ROOT / "external_wgs"
DEFAULT_QC = DEFAULT_EXTERNAL_ROOT / "qc/external_wgs_qc.tsv"
DEFAULT_CACHE_ROOT = DEFAULT_DATA_ROOT / "bundled_zenodo_caches"
DEFAULT_CONTROLLER_RESULTS = Path("/mnt/data/slurm-results")
DEFAULT_WORKER_RESULTS = Path("/results")
SLURM_SCRIPT = REPO_ROOT / "benchmarks/slurm_external_multi.sh"
PHASES = ("smoke", "warmup", "measured")
STRATEGIES = ("uncached", "gnomad_af_0.1", "gnomad_af_0.01", "cohort_3_genomes")
DEFAULT_SEED = "vcfcache-paper-external-wgs-v1"
EXPECTED_EVALUATION = {"kpgp": 20, "sgdp": 20, "pgp": 12}


@dataclass(frozen=True)
class ExternalTask:
    """One sample/replicate task containing all four benchmark conditions."""

    task_id: int
    phase: str
    measured: str
    cohort: str
    assembly: str
    sample: str
    population: str
    superpopulation: str
    sex: str
    provider: str
    input_vcf: str
    input_records: int
    input_sha256: str
    replicate: int
    strategy_order: str
    randomization_key: str


def read_evaluation_samples(
    qc_path: Path, *, require_complete: bool = True
) -> list[dict[str, str]]:
    """Read passing, unrelated, held-out external genomes."""
    with qc_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    eligible = [
        row
        for row in rows
        if row["role"] == "evaluation"
        and row["status"] == "PASS"
        and row["relatedness_status"] == "PASS"
        and row["documented_overlap"] == "none"
    ]
    if require_complete:
        counts = {
            cohort: sum(row["cohort"] == cohort for row in eligible)
            for cohort in EXPECTED_EVALUATION
        }
        present = {row["cohort"] for row in rows}
        expected = {
            cohort: count
            for cohort, count in EXPECTED_EVALUATION.items()
            if cohort in present
        }
        differences = {
            cohort: (counts[cohort], count)
            for cohort, count in expected.items()
            if counts[cohort] != count
        }
        if differences:
            raise RuntimeError(
                f"External evaluation cohort is incomplete: {differences}"
            )
    if not eligible:
        raise RuntimeError("No unrelated external evaluation genomes are ready")
    return sorted(eligible, key=lambda row: (row["cohort"], row["sample"]))


def condition_order(sample_rank: int, replicate: int) -> tuple[list[str], str]:
    """Return balanced cyclic condition order and an audit hash."""
    rotation = (sample_rank + replicate - 1) % len(STRATEGIES)
    order = [*STRATEGIES[rotation:], *STRATEGIES[:rotation]]
    key = hashlib.sha256(
        f"{DEFAULT_SEED}:{sample_rank}:{replicate}:{','.join(order)}".encode()
    ).hexdigest()
    return order, key


def build_tasks(
    samples: Sequence[dict[str, str]], *, phase: str, replicates: int
) -> list[ExternalTask]:
    """Build a deterministic task matrix with balanced first conditions."""
    if phase not in PHASES or replicates < 1:
        raise ValueError("Invalid phase or replicate count")
    ranked = sorted(
        samples,
        key=lambda row: hashlib.sha256(
            f"{DEFAULT_SEED}:sample:{row['cohort']}:{row['sample']}".encode()
        ).hexdigest(),
    )
    tasks: list[ExternalTask] = []
    for sample_rank, row in enumerate(ranked):
        for replicate in range(1, replicates + 1):
            order, key = condition_order(sample_rank, replicate)
            tasks.append(
                ExternalTask(
                    task_id=len(tasks),
                    phase=phase,
                    measured="true" if phase == "measured" else "false",
                    cohort=row["cohort"],
                    assembly=row["assembly"],
                    sample=row["sample"],
                    population=row["population"],
                    superpopulation=row["superpopulation"],
                    sex=row["sex"],
                    provider=row["provider"],
                    input_vcf=row["path"],
                    input_records=int(row["records"]),
                    input_sha256=row["sha256"],
                    replicate=replicate,
                    strategy_order=",".join(order),
                    randomization_key=key,
                )
            )
    return tasks


def write_tasks(path: Path, tasks: Sequence[ExternalTask]) -> None:
    """Atomically write a nonempty external task manifest."""
    if not tasks:
        raise ValueError("External task manifest must not be empty")
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_suffix(path.suffix + ".partial")
    with partial.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(asdict(tasks[0])),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(asdict(task) for task in tasks)
    partial.replace(path)


def _bundled_manifest(cache_root: Path, assembly: str) -> list[dict[str, Any]]:
    strategies = public_strategies(cache_root, assembly)
    specs = {spec.name: spec for spec in BUNDLED_CACHE_SPECS if spec.genome == assembly}
    values = []
    for strategy in strategies:
        spec = specs[strategy.name]
        cache = strategy.cache_dir
        provenance = cache.parents[1] / "zenodo_provenance.json"
        values.append(
            {
                "name": strategy.name,
                "kind": "bundled_zenodo",
                "assembly": assembly,
                "alias": spec.alias,
                "doi": spec.doi,
                "cache_dir": str(cache),
                "annotation_yaml_sha256": sha256sum(cache / "annotation.yaml"),
                "provenance_path": str(provenance),
                "provenance_expected": {
                    "alias": spec.alias,
                    "doi": spec.doi,
                    "archive_md5": spec.archive_md5,
                    "archive_md5_verified": True,
                    "source": "zenodo_production",
                },
            }
        )
    return values


def _custom_manifest(external_root: Path, cohort: str, assembly: str) -> dict[str, Any]:
    root = external_root / "cohort_caches" / cohort
    cache = root / "cohort3_database/cache" / VEP_CACHE_NAME
    provenance_path = root / "cohort3_cache_provenance.json"
    provenance = json.loads(provenance_path.read_text())
    expected = {
        "cohort": cohort,
        "assembly": assembly,
        "kind": "custom_cohort_three_disjoint_genomes",
        "complete": True,
    }
    for key, value in expected.items():
        if provenance.get(key) != value:
            raise RuntimeError(f"Invalid custom cache provenance for {cohort}: {key}")
    return {
        "name": "cohort_3_genomes",
        "kind": "custom_cohort",
        "assembly": assembly,
        "cache_dir": str(cache),
        "annotation_yaml_sha256": sha256sum(cache / "annotation.yaml"),
        "provenance_path": str(provenance_path),
        "provenance_expected": expected,
        "training_samples": provenance["training_samples"],
        "build_wall_seconds": provenance["wall_seconds"],
    }


def prepare_campaign(args: argparse.Namespace) -> None:
    """Freeze task, strategy, and provenance manifests."""
    root = campaign_root(args.controller_results, args.campaign_id)
    if root.exists() and any(root.iterdir()):
        raise FileExistsError(f"Campaign already exists and is nonempty: {root}")
    if subprocess.run(["git", "diff", "--quiet"], check=False).returncode:
        raise RuntimeError("Tracked worktree must be clean before campaign preparation")
    samples = read_evaluation_samples(args.qc)
    smoke = min(
        samples,
        key=lambda row: hashlib.sha256(
            f"{args.seed}:smoke:{row['cohort']}:{row['sample']}".encode()
        ).hexdigest(),
    )
    phase_specs = {
        "smoke": ([smoke], 1),
        "warmup": (samples, 1),
        "measured": (samples, 3),
    }
    phase_values = {}
    for phase, (phase_samples, replicates) in phase_specs.items():
        (root / phase / "tasks").mkdir(parents=True, exist_ok=True)
        (root / phase / "attempts").mkdir(parents=True, exist_ok=True)
        tasks = build_tasks(phase_samples, phase=phase, replicates=replicates)
        manifest = root / "manifests" / f"{phase}.tsv"
        write_tasks(manifest, tasks)
        phase_values[phase] = {
            "tasks": len(tasks),
            "samples": len({task.sample for task in tasks}),
            "replicates": replicates,
            "manifest": str(manifest),
            "manifest_sha256": sha256sum(manifest),
        }
    cohorts = sorted({row["cohort"] for row in samples})
    cohort_assemblies = {
        cohort: next(row["assembly"] for row in samples if row["cohort"] == cohort)
        for cohort in cohorts
    }
    if any(
        {row["assembly"] for row in samples if row["cohort"] == cohort}
        != {cohort_assemblies[cohort]}
        for cohort in cohorts
    ):
        raise RuntimeError("Each external cohort must use exactly one assembly")
    bundled: dict[str, list[dict[str, Any]]] = {
        assembly: _bundled_manifest(args.cache_root, assembly)
        for assembly in sorted(set(cohort_assemblies.values()))
    }
    cohort_strategies = {
        cohort: _custom_manifest(args.external_root, cohort, cohort_assemblies[cohort])
        for cohort in cohorts
    }
    strategies = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "commit": git_output("rev-parse", "HEAD"),
        "bundled_strategies_by_assembly": bundled,
        "cohort_assemblies": cohort_assemblies,
        "cohort_strategies": cohort_strategies,
    }
    strategy_path = root / "manifests/strategies.json"
    write_json_atomic(strategy_path, strategies)
    for cohort, assembly in cohort_assemblies.items():
        annotation_hashes = {
            item["annotation_yaml_sha256"]
            for item in [
                *bundled[assembly],
                cohort_strategies[cohort],
            ]
        }
        if len(annotation_hashes) != 1:
            raise RuntimeError(
                f"All four {cohort} conditions must use an identical recipe"
            )
    metadata = {
        "campaign_id": args.campaign_id,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "commit": git_output("rev-parse", "HEAD"),
        "commit_short": git_output("rev-parse", "--short=12", "HEAD"),
        "seed": args.seed,
        "qc": str(args.qc),
        "qc_sha256": sha256sum(args.qc),
        "strategies": str(strategy_path),
        "strategies_sha256": sha256sum(strategy_path),
        "conditions": list(STRATEGIES),
        "cohort_assemblies": cohort_assemblies,
        "evaluation_counts": {
            cohort: sum(row["cohort"] == cohort for row in samples)
            for cohort in cohorts
        },
        "phases": phase_values,
    }
    write_json_atomic(root / "campaign.json", metadata)
    print(json.dumps(metadata, indent=2, sort_keys=True))


def _task_count(path: Path) -> int:
    with path.open() as handle:
        count = sum(1 for _line in handle) - 1
    if count < 1:
        raise RuntimeError(f"Empty task manifest: {path}")
    return count


def submit_phase(
    args: argparse.Namespace, phase: str, dependency: str | None = None
) -> tuple[str, list[str]]:
    """Submit one external phase and return the Slurm job ID and command."""
    if shutil.which("sbatch") is None:
        raise RuntimeError("sbatch is unavailable")
    root = campaign_root(args.controller_results, args.campaign_id)
    manifest = root / "manifests" / f"{phase}.tsv"
    strategies = root / "manifests/strategies.json"
    (root / "logs").mkdir(parents=True, exist_ok=True)
    count = _task_count(manifest)
    command = [
        "sbatch",
        "--parsable",
        f"--job-name=external-{phase}",
        f"--array=0-{count - 1}%{args.concurrency}",
        f"--chdir={REPO_ROOT}",
        f"--output={args.worker_results}/campaigns/{args.campaign_id}/logs/{phase}-%A_%a.out",
        "--export=ALL,"
        f"VCFCACHE_CAMPAIGN_ID={args.campaign_id},"
        f"VCFCACHE_PHASE={phase},"
        f"VCFCACHE_TASK_MANIFEST={worker_path(manifest, args.controller_results, args.worker_results)},"
        f"VCFCACHE_STRATEGIES={worker_path(strategies, args.controller_results, args.worker_results)},"
        f"VCFCACHE_RESULTS_ROOT={args.worker_results},"
        f"VCFCACHE_REPO_ROOT={REPO_ROOT}",
    ]
    if dependency:
        command.append(f"--dependency=afterok:{dependency}")
    command.append(str(SLURM_SCRIPT))
    completed = subprocess.run(command, check=True, capture_output=True, text=True)
    job_id = completed.stdout.strip().split(";", maxsplit=1)[0]
    if not job_id.isdigit():
        raise RuntimeError(f"Unexpected sbatch response: {completed.stdout!r}")
    return job_id, command


def submit_chain(args: argparse.Namespace) -> None:
    """Submit smoke, warmup, and measured arrays with fail-closed dependencies."""
    root = campaign_root(args.controller_results, args.campaign_id)
    campaign = json.loads((root / "campaign.json").read_text())
    if campaign["commit"] != git_output("rev-parse", "HEAD"):
        raise RuntimeError("Prepared campaign commit differs from checkout")
    dependency = args.start_after_job
    submissions = {}
    for phase in PHASES:
        job_id, command = submit_phase(args, phase, dependency)
        submissions[phase] = {
            "job_id": job_id,
            "dependency": dependency,
            "command": command,
        }
        dependency = job_id
    value = {
        "campaign_id": args.campaign_id,
        "submitted_at": datetime.now(timezone.utc).isoformat(),
        "commit": campaign["commit"],
        "hostname": os.uname().nodename,
        "concurrency": args.concurrency,
        "submissions": submissions,
    }
    write_json_atomic(root / "submission.json", value)
    print(json.dumps(value, indent=2, sort_keys=True))


def status(args: argparse.Namespace) -> None:
    """Show queue state for the external campaign."""
    root = campaign_root(args.controller_results, args.campaign_id)
    submission = json.loads((root / "submission.json").read_text())
    jobs = ",".join(
        str(value["job_id"]) for value in submission["submissions"].values()
    )
    subprocess.run(
        ["squeue", "--jobs", jobs, "--format=%.18i %.24j %.8T %.10M %.24R"],
        check=True,
    )


def collect(args: argparse.Namespace) -> Path:
    """Collect measured summaries into a stable publication-source TSV."""
    root = campaign_root(args.controller_results, args.campaign_id)
    summaries = sorted((root / "measured/tasks").glob("task-*/external_summary.json"))
    rows = []
    for path in summaries:
        document = json.loads(path.read_text())
        task = document["task"]
        for value in document["rows"]:
            rows.append(
                {
                    **value,
                    "population": task["population"],
                    "superpopulation": task["superpopulation"],
                    "provider": task["provider"],
                    "input_sha256": task["input_sha256"],
                    "randomization_key": task["randomization_key"],
                }
            )
    campaign = json.loads((root / "campaign.json").read_text())
    expected = campaign["phases"]["measured"]["tasks"] * 3
    if len(rows) != expected:
        raise RuntimeError(f"Expected {expected} strategy rows, found {len(rows)}")
    output = root / "publication/external_wgs_metrics.tsv"
    output.parent.mkdir(parents=True, exist_ok=True)
    partial = output.with_suffix(".tsv.partial")
    with partial.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    partial.replace(output)
    print(f"Collected {len(rows)} measured strategy rows -> {output}")
    return output


def _paths(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--campaign-id", required=True)
    parser.add_argument(
        "--controller-results", type=Path, default=DEFAULT_CONTROLLER_RESULTS
    )
    parser.add_argument("--worker-results", type=Path, default=DEFAULT_WORKER_RESULTS)


def parser() -> argparse.ArgumentParser:
    """Build the external campaign CLI."""
    root = argparse.ArgumentParser(description=__doc__)
    commands = root.add_subparsers(dest="command", required=True)
    prepare = commands.add_parser("prepare")
    _paths(prepare)
    prepare.add_argument("--qc", type=Path, default=DEFAULT_QC)
    prepare.add_argument("--external-root", type=Path, default=DEFAULT_EXTERNAL_ROOT)
    prepare.add_argument("--cache-root", type=Path, default=DEFAULT_CACHE_ROOT)
    prepare.add_argument("--seed", default=DEFAULT_SEED)
    prepare.set_defaults(function=prepare_campaign)
    submit = commands.add_parser("submit-chain")
    _paths(submit)
    submit.add_argument("--concurrency", type=int, default=6)
    submit.add_argument("--start-after-job")
    submit.set_defaults(function=submit_chain)
    inspect = commands.add_parser("status")
    _paths(inspect)
    inspect.set_defaults(function=status)
    collection = commands.add_parser("collect")
    _paths(collection)
    collection.set_defaults(function=collect)
    return root


def main() -> None:
    """Run one external campaign action."""
    args = parser().parse_args()
    for name in (
        "controller_results",
        "worker_results",
        "qc",
        "external_root",
        "cache_root",
    ):
        if hasattr(args, name):
            setattr(args, name, getattr(args, name).expanduser().resolve())
    if getattr(args, "concurrency", 1) < 1:
        raise ValueError("Concurrency must be positive")
    args.function(args)


if __name__ == "__main__":
    main()
