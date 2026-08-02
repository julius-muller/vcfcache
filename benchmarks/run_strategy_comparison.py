#!/usr/bin/env python3
"""Run and plot the small publication cache-strategy comparison."""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib
import json
import shutil
import statistics
import subprocess
import tempfile
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence

from benchmarks.run_pilot import (
    PilotConfig,
    preflight,
    run_one,
    semantic_compare,
    write_json_atomic,
)

DEFAULT_DATA_ROOT = Path("/mnt/data/vcfcache_benchmarks")
DEFAULT_CACHE_ROOT = DEFAULT_DATA_ROOT / "bundled_zenodo_caches"
DEFAULT_BLUEPRINT_ROOT = DEFAULT_DATA_ROOT / "bundled_zenodo_blueprints"
DEFAULT_ROOT = DEFAULT_DATA_ROOT / "strategy_comparison_zenodo_v1"
DEFAULT_SELECTION = (
    Path(__file__).resolve().parent / "manifests/selected_1000g_samples.tsv"
)
SELECTION_SEED = "vcfcache-paper-cache-strategy-v1"
VEP_CACHE_NAME = "vep115.2_everything"


@dataclass(frozen=True)
class BundledCacheSpec:
    """One published cache bundle discoverable through VCFcache and Zenodo."""

    name: str
    genome: str
    alias: str
    doi: str
    root_name: str
    archive_name: str
    archive_md5: str
    af_threshold: float


BUNDLED_CACHE_SPECS = (
    BundledCacheSpec(
        name="gnomad_af_0.1",
        genome="GRCh38",
        alias="cache-gnomad-v4.1-GRCh38-joint-af01-vep115.2-e",
        doi="10.5281/zenodo.18189447",
        root_name="gnomad_v4.1_GRCh38_joint_af010",
        archive_name="cache_gnomad_v4.1_GRCh38_joint_af010.tar.gz",
        archive_md5="088cf426461a51b77bdfcd5dcd2233f4",
        af_threshold=0.1,
    ),
    BundledCacheSpec(
        name="gnomad_af_0.01",
        genome="GRCh38",
        alias="cache-gnomad-v4.1-GRCh38-joint-af001-vep115.2-e",
        doi="10.5281/zenodo.18190046",
        root_name="gnomad_v4.1_GRCh38_joint_af001",
        archive_name="cache_gnomad_v4.1_GRCh38_joint_af001.tar.gz",
        archive_md5="3ac438461eac0cf42c75717156d7b2d4",
        af_threshold=0.01,
    ),
    BundledCacheSpec(
        name="gnomad_af_0.1",
        genome="GRCh37",
        alias="cache-gnomad-v4.1-GRCh37-joint-af01-vep115.2-e",
        doi="10.5281/zenodo.18189051",
        root_name="gnomad_v4.1_GRCh37_joint_af010",
        archive_name="cache_gnomad_v4.1_GRCh37_joint_af010.tar.gz",
        archive_md5="96bb1edd0e067d9c933256bd112e4589",
        af_threshold=0.1,
    ),
    BundledCacheSpec(
        name="gnomad_af_0.01",
        genome="GRCh37",
        alias="cache-gnomad-v4.1-GRCh37-joint-af001-vep115.2-e",
        doi="10.5281/zenodo.18189348",
        root_name="gnomad_v4.1_GRCh37_joint_af001",
        archive_name="cache_gnomad_v4.1_GRCh37_joint_af001.tar.gz",
        archive_md5="f7d246a7adf44b778d6dc1383153eff2",
        af_threshold=0.01,
    ),
)


def bundled_cache_specs(genome: str = "GRCh38") -> tuple[BundledCacheSpec, ...]:
    """Return the frozen two-cache strategy for one reference assembly."""
    values = tuple(spec for spec in BUNDLED_CACHE_SPECS if spec.genome == genome)
    if len(values) != 2:
        raise RuntimeError(f"Expected two bundled cache strategies for {genome}")
    return values


@dataclass(frozen=True)
class BundledBlueprintSpec:
    """One official Zenodo blueprint allowed for local scenario cache builds."""

    name: str
    alias: str
    doi: str
    root_name: str
    archive_name: str
    archive_bytes: int
    archive_md5: str
    af_threshold: float


BUNDLED_BLUEPRINT_SPECS = (
    BundledBlueprintSpec(
        name="gnomad_af_0.1",
        alias="bp-gnomad-v4.1-GRCh38-joint-af01",
        doi="10.5281/zenodo.18190424",
        root_name="gnomad_v4.1_GRCh38_joint_af010",
        archive_name="bp_gnomad_v4.1_GRCh38_joint_af010.tar.gz",
        archive_bytes=45_756_242,
        archive_md5="c3d1ea67acd62b3fd9f30ea132d98a41",
        af_threshold=0.1,
    ),
    BundledBlueprintSpec(
        name="gnomad_af_0.01",
        alias="bp-gnomad-v4.1-GRCh38-joint-af001",
        doi="10.5281/zenodo.18190436",
        root_name="gnomad_v4.1_GRCh38_joint_af001",
        archive_name="bp_gnomad_v4.1_GRCh38_joint_af001.tar.gz",
        archive_bytes=103_598_148,
        archive_md5="6b7403ff03815500ba49c52ad285746c",
        af_threshold=0.01,
    ),
    BundledBlueprintSpec(
        name="gnomad_af_0.001",
        alias="bp-gnomad-v4.1-GRCh38-joint-af0001",
        doi="10.5281/zenodo.18190666",
        root_name="gnomad_v4.1_GRCh38_joint_af0001",
        archive_name="bp_gnomad_v4.1_GRCh38_joint_af0001.tar.gz",
        archive_bytes=496_477_799,
        archive_md5="1e44e7c08c8fb6aec6913eb2914ffabc",
        af_threshold=0.001,
    ),
)


@dataclass(frozen=True)
class Strategy:
    """One cache strategy evaluated against the common uncached baseline."""

    name: str
    kind: str
    cache_dir: Path
    blueprint: Path
    af_threshold: float | None = None
    alias: str | None = None
    doi: str | None = None


def run_checked(args: Sequence[str | Path], *, log: Path | None = None) -> None:
    """Run one preparation command, optionally capturing a persistent log."""
    command = [str(value) for value in args]
    if log is None:
        subprocess.run(command, check=True)
        return
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("w") as handle:
        subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, check=True)


def stable_rank(value: str, *, namespace: str) -> str:
    """Return an auditable deterministic selection key."""
    return hashlib.sha256(f"{SELECTION_SEED}:{namespace}:{value}".encode()).hexdigest()


def select_samples(selection: Path) -> dict[str, Any]:
    """Select three training and three held-out samples without cherry-picking."""
    with selection.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    populations = sorted(
        {row["superpopulation"] for row in rows},
        key=lambda value: stable_rank(value, namespace="population"),
    )[:3]
    training: list[dict[str, str]] = []
    evaluation: list[dict[str, str]] = []
    for population in populations:
        candidates = sorted(
            (row for row in rows if row["superpopulation"] == population),
            key=lambda row: stable_rank(row["sample"], namespace="sample"),
        )
        training.append(candidates[0])
        evaluation.append(candidates[1])
    return {
        "selection_seed": SELECTION_SEED,
        "superpopulations": populations,
        "training": training,
        "evaluation": evaluation,
    }


def sample_path(data_root: Path, sample: str, superpopulation: str) -> Path:
    """Resolve one frozen primary-cohort VCF."""
    return (
        data_root
        / "samples/GRCh38/1000g"
        / superpopulation
        / f"{sample}.GRCh38.small_variants.vcf.gz"
    )


def provenance_path(cache_root: Path, spec: BundledCacheSpec) -> Path:
    """Return the retained provenance record for a downloaded bundle."""
    return cache_root / spec.root_name / "zenodo_provenance.json"


def verify_bundled_cache(cache_root: Path, spec: BundledCacheSpec) -> None:
    """Fail closed unless a cache has verified Zenodo-download provenance."""
    root = cache_root / spec.root_name
    provenance = provenance_path(cache_root, spec)
    required = (
        root / "blueprint/vcfcache.bcf",
        root / f"cache/{VEP_CACHE_NAME}/vcfcache_annotated.bcf",
        root / f"cache/{VEP_CACHE_NAME}/params.snapshot.yaml",
        provenance,
    )
    require_paths(required)
    value = json.loads(provenance.read_text())
    expected = {
        "alias": spec.alias,
        "doi": spec.doi,
        "archive_name": spec.archive_name,
        "archive_md5": spec.archive_md5,
        "archive_md5_verified": True,
        "source": "zenodo_production",
    }
    differences = {
        key: (value.get(key), wanted)
        for key, wanted in expected.items()
        if value.get(key) != wanted
    }
    if differences:
        raise RuntimeError(
            f"Bundled cache lacks matching Zenodo provenance: {provenance}: "
            f"{differences}"
        )


def fetch_bundled(args: argparse.Namespace) -> None:
    """Download and checksum-verify two published assembly-matched bundles."""
    zenodo = importlib.import_module("vcfcache.integrations.zenodo")
    archive_utils = importlib.import_module("vcfcache.utils.archive")
    args.cache_root.mkdir(parents=True, exist_ok=True)
    archive_root = args.cache_root / "archives"
    archive_root.mkdir(parents=True, exist_ok=True)
    for spec in bundled_cache_specs(args.genome):
        try:
            verify_bundled_cache(args.cache_root, spec)
            print(f"Verified existing Zenodo bundle: {spec.alias}")
            continue
        except (FileNotFoundError, RuntimeError, json.JSONDecodeError):
            pass

        root = args.cache_root / spec.root_name
        if root.exists():
            raise RuntimeError(
                f"Refusing to replace unverified cache directory: {root}. "
                "Move it aside or select a clean --cache-root."
            )
        archive = archive_root / spec.archive_name
        if not archive.exists():
            partial = archive.with_suffix(archive.suffix + ".partial")
            if partial.exists():
                partial.unlink()
            zenodo.download_doi(spec.doi, partial, sandbox=False)
            partial_md5 = archive_utils.file_md5(partial)
            if partial_md5 != spec.archive_md5:
                raise RuntimeError(
                    f"Zenodo archive checksum mismatch for {spec.doi}: "
                    f"{partial_md5} != {spec.archive_md5}"
                )
            partial.replace(archive)
        observed_md5 = archive_utils.file_md5(archive)
        if observed_md5 != spec.archive_md5:
            raise RuntimeError(
                f"Zenodo archive checksum mismatch for {spec.doi}: "
                f"{observed_md5} != {spec.archive_md5}"
            )
        temporary = Path(tempfile.mkdtemp(prefix=f".{spec.name}-", dir=args.cache_root))
        try:
            extracted = archive_utils.extract_cache(archive, temporary)
            if extracted.name != spec.root_name:
                raise RuntimeError(
                    f"Unexpected root in Zenodo bundle {spec.doi}: {extracted.name}"
                )
            extracted.rename(root)
        finally:
            shutil.rmtree(temporary)
        write_json_atomic(
            provenance_path(args.cache_root, spec),
            {
                "alias": spec.alias,
                "doi": spec.doi,
                "source": "zenodo_production",
                "downloaded_at": datetime.now(timezone.utc).isoformat(),
                "archive_name": spec.archive_name,
                "archive_md5": observed_md5,
                "archive_md5_verified": True,
            },
        )
        verify_bundled_cache(args.cache_root, spec)
        print(f"Downloaded and verified Zenodo bundle: {spec.alias}")


def blueprint_provenance_path(blueprint_root: Path, spec: BundledBlueprintSpec) -> Path:
    """Return the retained provenance record for a downloaded blueprint."""
    return blueprint_root / spec.root_name / "zenodo_blueprint_provenance.json"


def verify_bundled_blueprint(blueprint_root: Path, spec: BundledBlueprintSpec) -> None:
    """Fail closed unless a blueprint came from its frozen Zenodo archive."""
    root = blueprint_root / spec.root_name
    provenance = blueprint_provenance_path(blueprint_root, spec)
    require_paths(
        (
            root / "blueprint/vcfcache.bcf",
            root / "blueprint/vcfcache.bcf.csi",
            provenance,
        )
    )
    value = json.loads(provenance.read_text())
    expected = {
        "alias": spec.alias,
        "doi": spec.doi,
        "archive_name": spec.archive_name,
        "archive_bytes": spec.archive_bytes,
        "archive_md5": spec.archive_md5,
        "archive_md5_verified": True,
        "source": "zenodo_production",
        "artifact_role": "blueprint_source_for_local_annotation_scenarios",
    }
    differences = {
        key: (value.get(key), wanted)
        for key, wanted in expected.items()
        if value.get(key) != wanted
    }
    if differences:
        raise RuntimeError(
            f"Blueprint lacks matching Zenodo provenance: {provenance}: {differences}"
        )


def fetch_blueprints(args: argparse.Namespace) -> None:
    """Download and checksum-verify official GRCh38 blueprint inputs."""
    zenodo = importlib.import_module("vcfcache.integrations.zenodo")
    archive_utils = importlib.import_module("vcfcache.utils.archive")
    args.blueprint_root.mkdir(parents=True, exist_ok=True)
    archive_root = args.blueprint_root / "archives"
    archive_root.mkdir(parents=True, exist_ok=True)
    for spec in BUNDLED_BLUEPRINT_SPECS:
        try:
            verify_bundled_blueprint(args.blueprint_root, spec)
            print(f"Verified existing Zenodo blueprint: {spec.alias}")
            continue
        except (FileNotFoundError, RuntimeError, json.JSONDecodeError):
            pass

        root = args.blueprint_root / spec.root_name
        if root.exists():
            raise RuntimeError(
                f"Refusing to replace unverified blueprint directory: {root}. "
                "Move it aside or select a clean --blueprint-root."
            )
        archive = archive_root / spec.archive_name
        if not archive.exists():
            partial = archive.with_suffix(archive.suffix + ".partial")
            if partial.exists():
                partial.unlink()
            zenodo.download_doi(spec.doi, partial, sandbox=False)
            partial_md5 = archive_utils.file_md5(partial)
            if partial_md5 != spec.archive_md5:
                raise RuntimeError(
                    f"Zenodo blueprint checksum mismatch for {spec.doi}: "
                    f"{partial_md5} != {spec.archive_md5}"
                )
            if partial.stat().st_size != spec.archive_bytes:
                raise RuntimeError(
                    f"Zenodo blueprint size mismatch for {spec.doi}: "
                    f"{partial.stat().st_size} != {spec.archive_bytes}"
                )
            partial.replace(archive)
        observed_md5 = archive_utils.file_md5(archive)
        if (
            observed_md5 != spec.archive_md5
            or archive.stat().st_size != spec.archive_bytes
        ):
            raise RuntimeError(
                f"Zenodo blueprint archive mismatch for {spec.doi}: "
                f"bytes={archive.stat().st_size}, md5={observed_md5}"
            )
        temporary = Path(
            tempfile.mkdtemp(prefix=f".{spec.name}-", dir=args.blueprint_root)
        )
        try:
            extracted = archive_utils.extract_cache(archive, temporary)
            if extracted.name != spec.root_name:
                raise RuntimeError(
                    f"Unexpected root in Zenodo blueprint {spec.doi}: {extracted.name}"
                )
            extracted.rename(root)
        finally:
            shutil.rmtree(temporary)
        write_json_atomic(
            blueprint_provenance_path(args.blueprint_root, spec),
            {
                "alias": spec.alias,
                "doi": spec.doi,
                "source": "zenodo_production",
                "artifact_role": "blueprint_source_for_local_annotation_scenarios",
                "downloaded_at": datetime.now(timezone.utc).isoformat(),
                "archive_name": spec.archive_name,
                "archive_bytes": spec.archive_bytes,
                "archive_md5": observed_md5,
                "archive_md5_verified": True,
            },
        )
        verify_bundled_blueprint(args.blueprint_root, spec)
        print(f"Downloaded and verified Zenodo blueprint: {spec.alias}")


def public_strategies(cache_root: Path, genome: str = "GRCh38") -> list[Strategy]:
    """Resolve only published cache bundles downloaded from Zenodo."""
    values: list[Strategy] = []
    for spec in bundled_cache_specs(genome):
        verify_bundled_cache(cache_root, spec)
        root = cache_root / spec.root_name
        values.append(
            Strategy(
                name=spec.name,
                kind="bundled_zenodo",
                cache_dir=root / "cache" / VEP_CACHE_NAME,
                blueprint=root / "blueprint/vcfcache.bcf",
                af_threshold=spec.af_threshold,
                alias=spec.alias,
                doi=spec.doi,
            )
        )
    return values


def cohort_strategy(root: Path) -> Strategy:
    """Resolve the custom cache built from the three training genomes."""
    database = root / "cohort3_database"
    return Strategy(
        name="cohort_3_genomes",
        kind="custom_cohort",
        cache_dir=database / "cache" / VEP_CACHE_NAME,
        blueprint=database / "blueprint/vcfcache.bcf",
    )


def require_paths(paths: Sequence[Path]) -> None:
    """Reject an incomplete strategy setup early."""
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing strategy prerequisites: {missing}")


def prepare(args: argparse.Namespace) -> None:
    """Freeze the design, merge training genomes, and build the custom cache."""
    args.root.mkdir(parents=True, exist_ok=True)
    design_path = args.root / "design.json"
    design = select_samples(args.selection)
    design.update(
        {
            "created_at": datetime.now(timezone.utc).isoformat(),
            "data_root": str(args.data_root),
            "cache_root": str(args.cache_root),
            "evaluation_is_disjoint": True,
            "bundled_caches": [
                {
                    "name": spec.name,
                    "alias": spec.alias,
                    "doi": spec.doi,
                    "archive_md5": spec.archive_md5,
                }
                for spec in bundled_cache_specs(args.genome)
            ],
        }
    )
    if design_path.exists():
        existing = json.loads(design_path.read_text())
        for key in (
            "selection_seed",
            "superpopulations",
            "training",
            "evaluation",
            "bundled_caches",
        ):
            if existing[key] != design[key]:
                raise RuntimeError(f"Frozen design differs at {key}")
        design = existing
    else:
        write_json_atomic(design_path, design)

    training_paths = [
        sample_path(args.data_root, row["sample"], row["superpopulation"])
        for row in design["training"]
    ]
    public = public_strategies(args.cache_root, args.genome)
    source_cache = next(
        strategy.cache_dir for strategy in public if strategy.name == "gnomad_af_0.01"
    )
    require_paths(
        [
            *training_paths,
            *(Path(f"{path}.tbi") for path in training_paths),
            *(strategy.cache_dir / "vcfcache_annotated.bcf" for strategy in public),
            *(strategy.blueprint for strategy in public),
            source_cache / "annotation.yaml",
            source_cache / "params.snapshot.yaml",
        ]
    )

    union = args.root / "cohort3_training_union.bcf"
    if not union.exists():
        run_checked(
            [
                "bcftools",
                "merge",
                "--merge",
                "none",
                "--threads",
                str(args.threads),
                "-Ob",
                "-o",
                union,
                *training_paths,
            ],
            log=args.root / "logs/merge_training.log",
        )
        run_checked(["bcftools", "index", "--threads", str(args.threads), union])
    require_paths([union, Path(f"{union}.csi")])

    custom = cohort_strategy(args.root)
    database = args.root / "cohort3_database"
    if not custom.blueprint.exists():
        run_checked(
            [
                Path(__file__).resolve().parents[1] / ".venv/bin/vcfcache",
                "blueprint-init",
                "--vcf",
                union,
                "--output",
                database,
                "--normalize",
            ],
            log=args.root / "logs/build_blueprint.log",
        )

    cache_bcf = custom.cache_dir / "vcfcache_annotated.bcf"
    if not cache_bcf.exists():
        started_at = datetime.now(timezone.utc)
        started = time.monotonic()
        try:
            run_checked(
                [
                    Path(__file__).resolve().parents[1] / ".venv/bin/vcfcache",
                    "cache-build",
                    "--db",
                    database,
                    "--name",
                    VEP_CACHE_NAME,
                    "--anno-config",
                    source_cache / "annotation.yaml",
                    "--params",
                    source_cache / "params.snapshot.yaml",
                ],
                log=args.root / "logs/build_custom_cache.log",
            )
        finally:
            write_json_atomic(
                args.root / "custom_cache_build.json",
                {
                    "started_at": started_at.isoformat(),
                    "completed_at": datetime.now(timezone.utc).isoformat(),
                    "wall_seconds": time.monotonic() - started,
                    "complete": cache_bcf.exists(),
                },
            )
    require_paths([cache_bcf, Path(f"{cache_bcf}.csi")])
    print(f"Strategy preparation complete: {args.root}")


def load_design(root: Path) -> dict[str, Any]:
    """Load the immutable train/test allocation."""
    path = root / "design.json"
    if not path.exists():
        raise FileNotFoundError(f"Run prepare first: {path}")
    return json.loads(path.read_text())


def config_for(
    root: Path,
    sample: Path,
    strategy: str,
    cache: Path,
    params: Path,
) -> PilotConfig:
    """Create an isolated pilot configuration for one matrix cell."""
    return PilotConfig(
        data_root=root / "runs" / sample.name.split(".", maxsplit=1)[0] / strategy,
        input_vcf=sample,
        cache_dir=cache,
        params_file=params,
        replicate=1,
    )


def completed_metrics(config: PilotConfig, mode: str) -> dict[str, Any] | None:
    """Return metrics for an already completed resumable matrix cell."""
    path = config.run_dir(mode) / "metrics.json"
    return json.loads(path.read_text()) if path.exists() else None


def execute(args: argparse.Namespace) -> None:
    """Run one uncached baseline and three cache strategies per held-out sample."""
    design = load_design(args.root)
    strategies = [
        *public_strategies(args.cache_root, args.genome),
        cohort_strategy(args.root),
    ]
    require_paths(
        [
            *(strategy.cache_dir / "vcfcache_annotated.bcf" for strategy in strategies),
            *(strategy.cache_dir / "params.snapshot.yaml" for strategy in strategies),
        ]
    )
    baseline_cache = strategies[1].cache_dir
    for row in design["evaluation"]:
        sample = sample_path(args.data_root, row["sample"], row["superpopulation"])
        baseline = config_for(
            args.root,
            sample,
            "uncached_baseline",
            baseline_cache,
            baseline_cache / "params.snapshot.yaml",
        )
        if completed_metrics(baseline, "uncached") is None:
            preflight(baseline)
            run_one(baseline, "uncached")
        baseline_bcf = baseline.run_dir("uncached") / "output.bcf"

        for strategy in strategies:
            config = config_for(
                args.root,
                sample,
                strategy.name,
                strategy.cache_dir,
                strategy.cache_dir / "params.snapshot.yaml",
            )
            if completed_metrics(config, "cached") is None:
                preflight(config)
                run_one(config, "cached")
            comparison_path = config.pilot_root / "comparison_to_baseline.json"
            if comparison_path.exists():
                comparison = json.loads(comparison_path.read_text())
            else:
                comparison = semantic_compare(
                    config.run_dir("cached") / "output.bcf", baseline_bcf
                )
                write_json_atomic(comparison_path, comparison)
            if comparison.get("semantic_pass") is not True:
                raise RuntimeError(
                    f"Strategy output failed semantic comparison: {comparison_path}"
                )
    print(f"Strategy execution complete: {args.root}")


def cache_records(path: Path) -> int:
    """Return an indexed blueprint's variant count."""
    output = subprocess.run(
        ["bcftools", "index", "--nrecords", str(path)],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    return int(output.strip())


def collect_rows(args: argparse.Namespace) -> list[dict[str, Any]]:
    """Collect complete matrix cells into publication source rows."""
    design = load_design(args.root)
    strategies = [
        *public_strategies(args.cache_root, args.genome),
        cohort_strategy(args.root),
    ]
    training = ",".join(row["sample"] for row in design["training"])
    rows: list[dict[str, Any]] = []
    for evaluation in design["evaluation"]:
        sample = sample_path(
            args.data_root, evaluation["sample"], evaluation["superpopulation"]
        )
        baseline = config_for(
            args.root,
            sample,
            "uncached_baseline",
            strategies[1].cache_dir,
            strategies[1].cache_dir / "params.snapshot.yaml",
        )
        baseline_metrics = completed_metrics(baseline, "uncached")
        if baseline_metrics is None:
            continue
        baseline_seconds = float(baseline_metrics["wall_seconds"])
        for strategy in strategies:
            config = config_for(
                args.root,
                sample,
                strategy.name,
                strategy.cache_dir,
                strategy.cache_dir / "params.snapshot.yaml",
            )
            metrics = completed_metrics(config, "cached")
            comparison_path = config.pilot_root / "comparison_to_baseline.json"
            if metrics is None or not comparison_path.exists():
                continue
            comparison = json.loads(comparison_path.read_text())
            cached_seconds = float(metrics["wall_seconds"])
            rows.append(
                {
                    "strategy": strategy.name,
                    "strategy_kind": strategy.kind,
                    "bundled_alias": strategy.alias or "",
                    "zenodo_doi": strategy.doi or "",
                    "af_threshold": strategy.af_threshold or "",
                    "cache_records": cache_records(strategy.blueprint),
                    "cache_bytes": (strategy.cache_dir / "vcfcache_annotated.bcf")
                    .stat()
                    .st_size,
                    "training_samples": (
                        training if strategy.kind == "custom_cohort" else ""
                    ),
                    "evaluation_sample": evaluation["sample"],
                    "evaluation_superpopulation": evaluation["superpopulation"],
                    "input_records": metrics["output_records"],
                    "cache_hit_rate": metrics["cache_hit_rate"],
                    "cached_wall_seconds": cached_seconds,
                    "uncached_wall_seconds": baseline_seconds,
                    "speedup": baseline_seconds / cached_seconds,
                    "fraction_wall_time_saved": 1 - cached_seconds / baseline_seconds,
                    "semantic_pass": comparison["semantic_pass"],
                    "raw_annotation_mismatches": comparison.get(
                        "raw_annotation_mismatches", 0
                    ),
                    "ignored_annotation_mismatches": comparison.get(
                        "ignored_annotation_mismatches", 0
                    ),
                }
            )
    return rows


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    """Write figure source data with a stable schema."""
    if not rows:
        raise RuntimeError("No completed strategy rows to collect")
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_suffix(path.suffix + ".partial")
    with partial.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    partial.replace(path)


def xml_escape(value: str) -> str:
    """Escape text inserted into the generated SVG."""
    return value.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def make_svg(path: Path, rows: list[dict[str, Any]]) -> None:
    """Create a dependency-free two-panel SVG with medians and sample points."""
    labels = {
        "gnomad_af_0.1": "Bundled gnomAD\nany-stratum AF ≥ 10%",
        "gnomad_af_0.01": "Bundled gnomAD\nany-stratum AF ≥ 1%",
        "cohort_3_genomes": "3-genome\ncohort cache",
    }
    order = list(labels)
    colors = {
        "gnomad_af_0.1": "#7AA6C2",
        "gnomad_af_0.01": "#4C78A8",
        "cohort_3_genomes": "#E07A5F",
    }
    grouped = {name: [row for row in rows if row["strategy"] == name] for name in order}
    if any(not grouped[name] for name in order):
        raise RuntimeError("Figure requires complete results for all three strategies")
    width, height = 1200, 560
    panel_width = 480
    lefts = (85, 685)
    top, bottom = 65, 420
    panel_height = bottom - top
    values = {
        "hit": {
            name: [100 * float(row["cache_hit_rate"]) for row in grouped[name]]
            for name in order
        },
        "speedup": {
            name: [float(row["speedup"]) for row in grouped[name]] for name in order
        },
    }
    maxima = {
        "hit": 100.0,
        "speedup": max(1.0, max(map(max, values["speedup"].values()))) * 1.12,
    }
    titles = {
        "hit": "A  Variants served from cache",
        "speedup": "B  End-to-end speedup",
    }
    ylabels = {"hit": "Cache hit rate (%)", "speedup": "Speedup vs uncached (×)"}
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        "<style>text{font-family:Arial,Helvetica,sans-serif;fill:#222}.title{font-size:20px;font-weight:700}.axis{font-size:14px}.label{font-size:13px}.value{font-size:13px;font-weight:700}</style>",
    ]
    offsets = (-12, 0, 12)
    for panel, key in enumerate(("hit", "speedup")):
        left = lefts[panel]
        maximum = maxima[key]
        svg.append(f'<text class="title" x="{left}" y="30">{titles[key]}</text>')
        for tick in range(6):
            value = maximum * tick / 5
            y = bottom - panel_height * value / maximum
            svg.append(
                f'<line x1="{left}" y1="{y:.1f}" x2="{left + panel_width}" y2="{y:.1f}" stroke="#E5E5E5"/>'
            )
            shown = f"{value:.0f}" if key == "hit" else f"{value:.1f}"
            svg.append(
                f'<text class="axis" text-anchor="end" x="{left - 10}" y="{y + 5:.1f}">{shown}</text>'
            )
        svg.append(
            f'<line x1="{left}" y1="{bottom}" x2="{left + panel_width}" y2="{bottom}" stroke="#222"/>'
        )
        svg.append(
            f'<line x1="{left}" y1="{top}" x2="{left}" y2="{bottom}" stroke="#222"/>'
        )
        svg.append(
            f'<text class="axis" transform="translate({left - 58},{(top + bottom) / 2}) rotate(-90)" text-anchor="middle">{ylabels[key]}</text>'
        )
        if key == "speedup":
            y = bottom - panel_height / maximum
            svg.append(
                f'<line x1="{left}" y1="{y:.1f}" x2="{left + panel_width}" y2="{y:.1f}" stroke="#555" stroke-dasharray="5,5"/>'
            )
        slot = panel_width / len(order)
        for index, name in enumerate(order):
            x = left + slot * (index + 0.5)
            sample_values = values[key][name]
            median = statistics.median(sample_values)
            y = bottom - panel_height * median / maximum
            bar_width = 66
            svg.append(
                f'<rect x="{x - bar_width / 2:.1f}" y="{y:.1f}" width="{bar_width}" height="{bottom - y:.1f}" fill="{colors[name]}" opacity="0.82"/>'
            )
            svg.append(
                f'<text class="value" text-anchor="middle" x="{x:.1f}" y="{max(top + 13, y - 8):.1f}">{median:.1f}</text>'
            )
            for point_index, point in enumerate(sample_values):
                point_y = bottom - panel_height * point / maximum
                offset = offsets[point_index] if point_index < len(offsets) else 0
                svg.append(
                    f'<circle cx="{x + offset:.1f}" cy="{point_y:.1f}" r="4.5" fill="#111" stroke="white" stroke-width="1"/>'
                )
            lines = labels[name].split("\n")
            for line_index, line in enumerate(lines):
                svg.append(
                    f'<text class="label" text-anchor="middle" x="{x:.1f}" y="{447 + 17 * line_index}">{xml_escape(line)}</text>'
                )
    svg.append(
        '<text class="label" x="600" y="525" text-anchor="middle">Points: three held-out genomes. Bundled caches are verified Zenodo downloads; the custom cache uses disjoint training genomes.</text>'
    )
    svg.append("</svg>")
    path.write_text("\n".join(svg) + "\n")


def collect(args: argparse.Namespace) -> None:
    """Write the strategy TSV and render its publication figure."""
    rows = collect_rows(args)
    expected = 3 * 3
    if len(rows) != expected:
        raise RuntimeError(f"Expected {expected} complete rows, found {len(rows)}")
    tsv = args.root / "figure/cache_strategy.tsv"
    svg = args.root / "figure/cache_strategy.svg"
    write_tsv(tsv, rows)
    make_svg(svg, rows)
    print(f"Wrote {tsv}")
    print(f"Wrote {svg}")


def build_parser() -> argparse.ArgumentParser:
    """Build the strategy comparison CLI."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "command",
        choices=(
            "fetch-bundled",
            "fetch-blueprints",
            "prepare",
            "execute",
            "collect",
            "all",
        ),
    )
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    parser.add_argument("--cache-root", type=Path, default=DEFAULT_CACHE_ROOT)
    parser.add_argument("--blueprint-root", type=Path, default=DEFAULT_BLUEPRINT_ROOT)
    parser.add_argument("--selection", type=Path, default=DEFAULT_SELECTION)
    parser.add_argument("--genome", choices=("GRCh37", "GRCh38"), default="GRCh38")
    parser.add_argument("--threads", type=int, default=8)
    return parser


def main() -> None:
    """Run the selected resumable strategy step."""
    args = build_parser().parse_args()
    for name in (
        "root",
        "data_root",
        "cache_root",
        "blueprint_root",
        "selection",
    ):
        setattr(args, name, getattr(args, name).expanduser().resolve())
    if args.threads < 1:
        raise ValueError("--threads must be positive")
    if args.command in ("fetch-bundled", "all"):
        fetch_bundled(args)
    if args.command in ("fetch-blueprints", "all"):
        fetch_blueprints(args)
    if args.command in ("prepare", "all"):
        prepare(args)
    if args.command in ("execute", "all"):
        execute(args)
    if args.command in ("collect", "all"):
        collect(args)


if __name__ == "__main__":
    main()
