#!/usr/bin/env python3
"""Prepare, submit, inspect, and collect external-WGS Slurm campaigns."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
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
PHASES = ("measured",)
LEGACY_COLLECTION_PHASES = ("warmup", "measured")
STRATEGIES = ("uncached", "gnomad_af_0.1", "gnomad_af_0.01", "cohort_3_genomes")
DEFAULT_SEED = "vcfcache-paper-external-wgs-v1"
EXPECTED_EVALUATION = {"kpgp": 20, "sgdp": 20, "pgp": 12}
DEFAULT_VEP_BUFFER = 100_000
DEFAULT_MEASURED_REPLICATES = 1


def calibration_subset(
    samples: Sequence[dict[str, str]], per_cohort: int, seed: str
) -> list[dict[str, str]]:
    """Select a deterministic small light-statistics calibration by cohort."""
    if per_cohort < 1:
        return list(samples)
    selected: list[dict[str, str]] = []
    for cohort in sorted({row["cohort"] for row in samples}):
        cohort_rows = [row for row in samples if row["cohort"] == cohort]
        ranked = sorted(
            cohort_rows,
            key=lambda row: hashlib.sha256(
                f"{seed}:light-calibration:{cohort}:{row['sample']}".encode()
            ).hexdigest(),
        )
        if len(ranked) < per_cohort:
            raise RuntimeError(
                f"Only {len(ranked)} {cohort} samples for calibration of {per_cohort}"
            )
        selected.extend(ranked[:per_cohort])
    return selected


@dataclass(frozen=True)
class ExternalTask:
    """One sample/tool task containing all four benchmark conditions."""

    task_id: int
    tool: str
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


def condition_order(
    sample_rank: int, replicate: int = 1, seed: str = DEFAULT_SEED
) -> tuple[list[str], str]:
    """Return balanced cyclic condition order and an audit hash."""
    rotation = (sample_rank + replicate - 1) % len(STRATEGIES)
    order = [*STRATEGIES[rotation:], *STRATEGIES[:rotation]]
    key = hashlib.sha256(
        f"{seed}:{sample_rank}:{replicate}:{','.join(order)}".encode()
    ).hexdigest()
    return order, key


def build_tasks(
    samples: Sequence[dict[str, str]],
    *,
    phase: str = "measured",
    replicates: int = 1,
    tool: str = "vep",
    seed: str = DEFAULT_SEED,
) -> list[ExternalTask]:
    """Build a deterministic task matrix with balanced first conditions."""
    if phase not in PHASES:
        raise ValueError(f"Unknown publication phase: {phase}")
    if replicates != 1:
        raise ValueError(
            "Publication campaigns permit exactly one run per sample/tool/condition"
        )
    if tool not in {"vep", "fastvep"}:
        raise ValueError(f"Unsupported annotator: {tool}")
    ranked = sorted(
        samples,
        key=lambda row: hashlib.sha256(
            f"{DEFAULT_SEED}:sample:{row['cohort']}:{row['sample']}".encode()
        ).hexdigest(),
    )
    tasks: list[ExternalTask] = []
    for sample_rank, row in enumerate(ranked):
        order, key = condition_order(sample_rank, 1, seed)
        tasks.append(
            ExternalTask(
                task_id=len(tasks),
                tool=tool,
                phase=phase,
                measured="true",
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
                replicate=1,
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


def _runtime_params(
    root: Path,
    assembly: str,
    source: Path,
    vep_buffer: int,
    published_path: Path | None = None,
) -> dict[str, str | int]:
    """Freeze lower-memory runtime params without modifying a bundled cache."""
    updated, replacements = re.subn(
        r"(?m)^vep_buffer:\s*\d+\s*$",
        f"vep_buffer: {vep_buffer}",
        source.read_text(),
    )
    if replacements != 1:
        raise RuntimeError(f"Expected one vep_buffer entry in {source}")
    destination = root / "manifests" / f"runtime_params_{assembly}.yaml"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(updated)
    return {
        "path": str(published_path or destination),
        "controller_path": str(destination),
        "sha256": sha256sum(destination),
        "source": str(source),
        "source_sha256": sha256sum(source),
        "vep_buffer": vep_buffer,
    }


def validate_strategy_document(
    document: dict[str, Any],
    *,
    tool: str,
    cohort_assemblies: dict[str, str],
) -> None:
    """Fail closed on an annotator-specific three-cache strategy bundle."""
    if document.get("tool", "vep") != tool:
        raise RuntimeError(
            f"Strategy bundle is for {document.get('tool', 'vep')}, not {tool}"
        )
    if document.get("cohort_assemblies") != cohort_assemblies:
        raise RuntimeError(
            "Strategy bundle does not cover the exact evaluation cohorts"
        )
    bundled = document.get("bundled_strategies_by_assembly", {})
    custom = document.get("cohort_strategies", {})
    runtime = document.get("runtime_params_by_assembly", {})
    for cohort, assembly in cohort_assemblies.items():
        values = [*bundled.get(assembly, []), custom.get(cohort)]
        if any(value is None for value in values):
            raise RuntimeError(f"Strategy bundle is incomplete for {cohort}")
        names = [value["name"] for value in values]
        if names != ["gnomad_af_0.1", "gnomad_af_0.01", "cohort_3_genomes"]:
            raise RuntimeError(f"Unexpected strategy set for {cohort}: {names}")
        recipe_hashes = set()
        for value in values:
            if value.get("assembly") != assembly:
                raise RuntimeError(f"Strategy assembly mismatch for {value['name']}")
            cache = Path(value.get("controller_cache_dir", value["cache_dir"]))
            required = (
                cache / "annotation.yaml",
                cache / "params.snapshot.yaml",
                cache / "vcfcache_annotated.bcf",
                cache / "vcfcache_annotated.bcf.csi",
            )
            missing = [str(path) for path in required if not path.exists()]
            if missing:
                raise FileNotFoundError(f"Incomplete {tool} cache: {missing}")
            observed = sha256sum(cache / "annotation.yaml")
            if observed != value.get("annotation_yaml_sha256"):
                raise RuntimeError(f"Recipe checksum mismatch for {value['name']}")
            recipe_hashes.add(observed)
        if len(recipe_hashes) != 1:
            raise RuntimeError(
                f"All four {tool}/{cohort} conditions must use one recipe"
            )
        params = runtime.get(assembly)
        if not params:
            raise RuntimeError(f"Runtime params are absent for {assembly}")
        params_path = Path(params.get("controller_path", params["path"]))
        if not params_path.exists() or sha256sum(params_path) != params.get("sha256"):
            raise RuntimeError(f"Runtime params are invalid for {assembly}")


def _external_strategy_bundle(
    args: argparse.Namespace,
    root: Path,
    cohort_assemblies: dict[str, str],
) -> dict[str, Any]:
    """Create VEP strategies or freeze a prepared fastVEP strategy bundle."""
    tool = getattr(args, "tool", "vep")
    if tool == "fastvep":
        source = getattr(args, "strategy_manifest", None)
        if source is None:
            raise ValueError("--strategy-manifest is required for fastvep")
        document = json.loads(source.read_text())
        validate_strategy_document(
            document, tool="fastvep", cohort_assemblies=cohort_assemblies
        )
        return {
            **document,
            "commit": git_output("rev-parse", "HEAD"),
            "source_manifest": str(source),
            "source_manifest_sha256": sha256sum(source),
        }

    if getattr(args, "strategy_manifest", None) is not None:
        raise ValueError("--strategy-manifest is only used with --tool fastvep")
    bundled: dict[str, list[dict[str, Any]]] = {
        assembly: _bundled_manifest(args.cache_root, assembly)
        for assembly in sorted(set(cohort_assemblies.values()))
    }
    cohort_strategies = {
        cohort: _custom_manifest(args.external_root, cohort, cohort_assemblies[cohort])
        for cohort in cohort_assemblies
    }
    runtime_params_by_assembly = {
        assembly: _runtime_params(
            root,
            assembly,
            Path(bundled[assembly][0]["cache_dir"]) / "params.snapshot.yaml",
            args.vep_buffer,
            worker_path(
                root / "manifests" / f"runtime_params_{assembly}.yaml",
                args.controller_results,
                args.worker_results,
            ),
        )
        for assembly in bundled
    }
    document = {
        "tool": "vep",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "commit": git_output("rev-parse", "HEAD"),
        "bundled_strategies_by_assembly": bundled,
        "cohort_assemblies": cohort_assemblies,
        "cohort_strategies": cohort_strategies,
        "runtime_params_by_assembly": runtime_params_by_assembly,
    }
    validate_strategy_document(
        document, tool="vep", cohort_assemblies=cohort_assemblies
    )
    return document


def prepare_campaign(args: argparse.Namespace) -> None:
    """Freeze task, strategy, and provenance manifests."""
    root = campaign_root(args.controller_results, args.campaign_id)
    if root.exists() and any(root.iterdir()):
        raise FileExistsError(f"Campaign already exists and is nonempty: {root}")
    if subprocess.run(
        ["git", "-C", str(REPO_ROOT), "diff", "--quiet"], check=False
    ).returncode:
        raise RuntimeError("Tracked worktree must be clean before campaign preparation")
    samples = read_evaluation_samples(args.qc)
    calibration_per_cohort = getattr(args, "calibration_per_cohort", 0)
    if calibration_per_cohort:
        if getattr(args, "tool", "vep") != "vep":
            raise ValueError("--calibration-per-cohort is reserved for VEP calibration")
        samples = calibration_subset(samples, calibration_per_cohort, args.seed)
    phase_specs = {"measured": (samples, DEFAULT_MEASURED_REPLICATES)}
    phase_values = {}
    for phase, (phase_samples, replicates) in phase_specs.items():
        (root / phase / "tasks").mkdir(parents=True, exist_ok=True)
        (root / phase / "attempts").mkdir(parents=True, exist_ok=True)
        tasks = build_tasks(
            phase_samples,
            phase=phase,
            replicates=replicates,
            tool=getattr(args, "tool", "vep"),
            seed=args.seed,
        )
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
    strategies = _external_strategy_bundle(args, root, cohort_assemblies)
    strategy_path = root / "manifests/strategies.json"
    write_json_atomic(strategy_path, strategies)
    metadata = {
        "campaign_id": args.campaign_id,
        "tool": getattr(args, "tool", "vep"),
        "created_at": datetime.now(timezone.utc).isoformat(),
        "commit": git_output("rev-parse", "HEAD"),
        "commit_short": git_output("rev-parse", "--short=12", "HEAD"),
        "seed": args.seed,
        "qc": str(args.qc),
        "qc_sha256": sha256sum(args.qc),
        "strategies": str(strategy_path),
        "strategies_sha256": sha256sum(strategy_path),
        "conditions": list(STRATEGIES),
        "statistics_mode": "light",
        "calibration_per_cohort": calibration_per_cohort,
        "cohort_assemblies": cohort_assemblies,
        "runtime_params_by_assembly": strategies["runtime_params_by_assembly"],
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
    args: argparse.Namespace,
    phase: str,
    dependency: str | None = None,
    *,
    array_spec: str | None = None,
    job_suffix: str | None = None,
) -> tuple[str, list[str]]:
    """Submit one external phase and return the Slurm job ID and command."""
    if shutil.which("sbatch") is None:
        raise RuntimeError("sbatch is unavailable")
    root = campaign_root(args.controller_results, args.campaign_id)
    manifest = root / "manifests" / f"{phase}.tsv"
    strategies = root / "manifests/strategies.json"
    tool = json.loads((root / "campaign.json").read_text()).get("tool", "vep")
    (root / "logs").mkdir(parents=True, exist_ok=True)
    count = _task_count(manifest)
    task_array = array_spec or f"0-{count - 1}%{args.concurrency}"
    command = [
        "sbatch",
        "--parsable",
        f"--job-name=external-{tool}-{job_suffix or phase}",
        f"--array={task_array}",
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
    """Submit one four-condition measurement per sample."""
    root = campaign_root(args.controller_results, args.campaign_id)
    campaign = json.loads((root / "campaign.json").read_text())
    if campaign["commit"] != git_output("rev-parse", "HEAD"):
        raise RuntimeError("Prepared campaign commit differs from checkout")
    dependency = args.start_after_job
    submissions = {}
    if getattr(args, "smoke_first", False):
        phase = "measured"
        smoke_id, smoke_command = submit_phase(
            args,
            phase,
            dependency,
            array_spec="0-0",
            job_suffix="smoke",
        )
        submissions["smoke"] = {
            "job_id": smoke_id,
            "dependency": dependency,
            "command": smoke_command,
        }
        count = _task_count(root / "manifests" / f"{phase}.tsv")
        if count > 1:
            dependency = smoke_id
            measured_id, measured_command = submit_phase(
                args,
                phase,
                dependency,
                array_spec=f"1-{count - 1}%{args.concurrency}",
            )
            submissions[phase] = {
                "job_id": measured_id,
                "dependency": dependency,
                "command": measured_command,
            }
    else:
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
    """Collect one complete source phase into a publication-source TSV.

    ``--phase warmup`` allows a complete, validated first pass from a legacy
    campaign to be promoted without rerunning identical work. The source phase
    remains recorded in every output row for auditability.
    """
    root = campaign_root(args.controller_results, args.campaign_id)
    summaries = sorted(
        (root / args.phase / "tasks").glob("task-*/external_summary.json")
    )
    rows = []
    for path in summaries:
        document = json.loads(path.read_text())
        task = document["task"]
        for value in document["rows"]:
            rows.append(
                {
                    **value,
                    "statistics_mode": value.get(
                        "statistics_mode",
                        document.get("statistics_mode", "legacy_full_rescan"),
                    ),
                    "semantic_comparator": value.get(
                        "semantic_comparator",
                        document.get("semantic_comparator", "legacy_vep_semantic"),
                    ),
                    "population": task["population"],
                    "superpopulation": task["superpopulation"],
                    "provider": task["provider"],
                    "input_sha256": task["input_sha256"],
                    "randomization_key": task["randomization_key"],
                }
            )
    campaign = json.loads((root / "campaign.json").read_text())
    expected = campaign["phases"][args.phase]["tasks"] * 3
    if len(rows) != expected:
        raise RuntimeError(f"Expected {expected} strategy rows, found {len(rows)}")
    output = root / "publication" / args.output_name
    output.parent.mkdir(parents=True, exist_ok=True)
    partial = output.with_suffix(".tsv.partial")
    with partial.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    partial.replace(output)
    print(f"Collected {len(rows)} {args.phase} strategy rows -> {output}")
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
    prepare.add_argument("--tool", choices=("vep", "fastvep"), default="vep")
    prepare.add_argument(
        "--strategy-manifest",
        type=Path,
        help="Prepared annotator-specific cache bundle (required for fastvep)",
    )
    prepare.add_argument("--vep-buffer", type=int, default=DEFAULT_VEP_BUFFER)
    prepare.add_argument(
        "--calibration-per-cohort",
        type=int,
        default=0,
        help="Prepare a deterministic VEP light-mode calibration subset",
    )
    prepare.set_defaults(function=prepare_campaign)
    submit = commands.add_parser("submit-chain")
    _paths(submit)
    submit.add_argument("--concurrency", type=int, default=6)
    submit.add_argument("--start-after-job")
    submit.add_argument(
        "--smoke-first",
        action="store_true",
        help="Run task 0 as a gate before submitting the remaining array",
    )
    submit.set_defaults(function=submit_chain)
    inspect = commands.add_parser("status")
    _paths(inspect)
    inspect.set_defaults(function=status)
    collection = commands.add_parser("collect")
    _paths(collection)
    collection.add_argument(
        "--phase", choices=LEGACY_COLLECTION_PHASES, default="measured"
    )
    collection.add_argument(
        "--output-name",
        default="external_wgs_metrics.tsv",
        help="Output filename below publication/; use a new name for amendments",
    )
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
        "strategy_manifest",
    ):
        if hasattr(args, name) and getattr(args, name) is not None:
            setattr(args, name, getattr(args, name).expanduser().resolve())
    if getattr(args, "concurrency", 1) < 1:
        raise ValueError("Concurrency must be positive")
    if getattr(args, "vep_buffer", 1) < 1:
        raise ValueError("VEP buffer must be positive")
    if getattr(args, "calibration_per_cohort", 0) < 0:
        raise ValueError("Calibration sample count cannot be negative")
    if hasattr(args, "output_name") and (
        Path(args.output_name).name != args.output_name
        or not args.output_name.endswith(".tsv")
    ):
        raise ValueError("--output-name must be a simple .tsv filename")
    args.function(args)


if __name__ == "__main__":
    main()
