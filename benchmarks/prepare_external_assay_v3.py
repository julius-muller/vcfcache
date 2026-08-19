#!/usr/bin/env python3
"""Prepare frozen manifests for the independent Panel/WES annotator extension."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_METRICS = (
    REPO_ROOT
    / "benchmarks/figures/source_data/final"
    / "2026-08-09T223718Z-light-matched-final"
    / "annotator_external_wgs_metrics.tsv"
)
DEFAULT_CAMPAIGN = "external-assay-v3-independent-af001-warm"
DEFAULT_CONTROLLER_ROOT = Path("/mnt/data/slurm-results/campaigns")
DEFAULT_WORKER_ROOT = Path("/results/campaigns")
SEED = "vcfcache-paper-external-assay-v3"
COHORTS = ("kpgp", "sgdp")
WARMUP_SAMPLE = "KPGP-00319"
ASSAYS = {
    "panel": Path(
        "/mnt/data/vcfcache_benchmarks/regions/GRCh38/"
        "acmg_sf_v3.3.ensembl115.mane_select_cds_pad20.bed"
    ),
    "wes": Path(
        "/mnt/data/vcfcache_benchmarks/regions/GRCh38/"
        "twist_human_core_exome_hg38.chr1-22-X.pad125.merged.bed"
    ),
}
TOOL_CONFIG = {
    "vep": {
        "cache_dir": (
            "/mnt/data/vcfcache_benchmarks/bundled_zenodo_caches/"
            "gnomad_v4.1_GRCh38_joint_af001/cache/vep115.2_everything"
        ),
        "params_file": (
            "/results/campaigns/external-vep-light-full52-bbef41d7db8e/"
            "manifests/runtime_params_GRCh38.yaml"
        ),
        "cache_kind": "bundled_zenodo",
        "cache_alias": "cache-gnomad-v4.1-GRCh38-joint-af001-vep115.2-e",
    },
    "fastvep": {
        "cache_dir": (
            "/mnt/data/vcfcache_benchmarks/fastvep_external_caches/databases/"
            "GRCh38/bundled-GRCh38-gnomad_af_0.01/cache/"
            "cache-fastvep-0.3.0-publication"
        ),
        "params_file": (
            "/mnt/data/vcfcache_benchmarks/fastvep_external_caches/"
            "config/runtime_params_GRCh38.yaml"
        ),
        "cache_kind": "locally_built_fastvep_from_bundled_blueprint",
        "cache_alias": "",
    },
}


def sha256sum(path: Path) -> str:
    """Return a streaming SHA-256 digest."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def selected_samples(metrics: Path, per_cohort: int) -> list[dict[str, str]]:
    """Select deterministic independent GRCh38 genomes from the frozen campaign."""
    with metrics.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    unique: dict[tuple[str, str], dict[str, str]] = {}
    for row in rows:
        if (
            row["tool"] == "vep"
            and row["strategy"] == "gnomad_af_0.01"
            and row["cohort"] in COHORTS
            and row["assembly"] == "GRCh38"
        ):
            unique[(row["cohort"], row["sample"])] = row
    selected: list[dict[str, str]] = []
    for cohort in COHORTS:
        candidates = [row for (group, _sample), row in unique.items() if group == cohort]
        ranked = sorted(
            candidates,
            key=lambda row: hashlib.sha256(
                f"{SEED}:{cohort}:{row['sample']}".encode()
            ).hexdigest(),
        )
        if len(ranked) < per_cohort:
            raise RuntimeError(f"Only {len(ranked)} independent {cohort} genomes")
        selected.extend(ranked[:per_cohort])
    return selected


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    """Write a nonempty tab-separated manifest atomically."""
    if not rows:
        raise ValueError(f"Cannot write an empty manifest: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_suffix(path.suffix + ".partial")
    with partial.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    partial.replace(path)


def prepare(args: argparse.Namespace) -> None:
    """Freeze source-sample, input-preparation, and annotation task manifests."""
    samples = selected_samples(args.metrics, args.per_cohort)
    controller = args.controller_root / args.campaign_id
    worker = args.worker_root / args.campaign_id
    manifests = controller / "manifests"
    source_rows: list[dict[str, object]] = []
    prep_rows: list[dict[str, object]] = []
    benchmark_rows: list[dict[str, object]] = []
    warmup_inputs = {
        assay: worker / "inputs" / "untimed_warmup" / f"{WARMUP_SAMPLE}.{assay}.vcf.gz"
        for assay in ASSAYS
    }

    for sample_row in samples:
        cohort = sample_row["cohort"]
        sample = sample_row["sample"]
        source_vcf = Path(
            f"/mnt/data/vcfcache_benchmarks/external_wgs/samples/GRCh38/"
            f"{cohort}/{sample}.GRCh38.small_variants.vcf.gz"
        )
        source_rows.append(
            {
                "cohort": cohort,
                "sample": sample,
                "assembly": "GRCh38",
                "source_vcf": source_vcf,
                "wgs_input_records": sample_row["input_records"],
                "wgs_cache_hit_rate": sample_row["cache_hit_rate"],
            }
        )
        for assay, bed in ASSAYS.items():
            prep_id = len(prep_rows)
            input_vcf = worker / "inputs" / assay / f"{sample}.GRCh38.{assay}.vcf.gz"
            worker_bed = worker / "regions" / bed.name
            prep_rows.append(
                {
                    "task_id": prep_id,
                    "cohort": cohort,
                    "sample": sample,
                    "assay": assay,
                    "source_vcf": source_vcf,
                    "bed": worker_bed,
                    "output_vcf": input_vcf,
                }
            )
            for tool, config in TOOL_CONFIG.items():
                task_id = len(benchmark_rows)
                order_key = hashlib.sha256(
                    f"{SEED}:{tool}:{assay}:{sample}".encode()
                ).digest()[0]
                condition_order = (
                    "cached,uncached" if order_key % 2 else "uncached,cached"
                )
                benchmark_rows.append(
                    {
                        "task_id": task_id,
                        "tool": tool,
                        "assay": assay,
                        "cohort": cohort,
                        "sample": sample,
                        "assembly": "GRCh38",
                        "input_vcf": input_vcf,
                        "warmup_input_vcf": warmup_inputs[assay],
                        "cache_strategy": "gnomad_af_0.01",
                        "cache_kind": config["cache_kind"],
                        "cache_alias": config["cache_alias"],
                        "cache_dir": config["cache_dir"],
                        "params_file": config["params_file"],
                        "condition_order": condition_order,
                        "replicate": 1,
                    }
                )

    warmup_source = Path(
        "/mnt/data/vcfcache_benchmarks/external_wgs/samples/GRCh38/"
        f"kpgp/{WARMUP_SAMPLE}.GRCh38.small_variants.vcf.gz"
    )
    for assay, bed in ASSAYS.items():
        prep_rows.append(
            {
                "task_id": len(prep_rows),
                "cohort": "untimed_warmup",
                "sample": WARMUP_SAMPLE,
                "assay": assay,
                "source_vcf": warmup_source,
                "bed": worker / "regions" / bed.name,
                "output_vcf": warmup_inputs[assay],
            }
        )

    write_tsv(manifests / "source_samples.tsv", source_rows)
    write_tsv(manifests / "prepare_tasks.tsv", prep_rows)
    write_tsv(manifests / "benchmark_tasks.tsv", benchmark_rows)
    campaign = {
        "campaign_id": args.campaign_id,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "design": "matched independent Panel/WES extension to frozen WGS campaign",
        "seed": SEED,
        "source_metrics": str(args.metrics),
        "source_metrics_sha256": sha256sum(args.metrics),
        "cohorts": list(COHORTS),
        "samples_per_cohort": args.per_cohort,
        "samples": len(samples),
        "assays": list(ASSAYS),
        "source_interval_beds": {key: str(value) for key, value in ASSAYS.items()},
        "worker_interval_beds": {
            key: str(worker / "regions" / value.name)
            for key, value in ASSAYS.items()
        },
        "tools": list(TOOL_CONFIG),
        "cache_strategy": "gnomad_af_0.01",
        "technical_repeats": 0,
        "untimed_warmup_sample": WARMUP_SAMPLE,
        "untimed_warmup_before_each_task": True,
        "warmup_in_reported_wall_time": False,
        "prep_tasks": len(prep_rows),
        "benchmark_tasks": len(benchmark_rows),
        "worker_root": str(worker),
        "tool_config": TOOL_CONFIG,
        "implementation_sha256": {
            relative: sha256sum(REPO_ROOT / relative)
            for relative in (
                "benchmarks/prepare_external_assay_v3.py",
                "benchmarks/run_external_assay_task.py",
                "benchmarks/slurm_prepare_external_assay_v3.sh",
                "benchmarks/slurm_external_assay_v3.sh",
                "benchmarks/collect_external_assay_v3.py",
            )
        },
    }
    (controller / "campaign.json").write_text(
        json.dumps(campaign, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(campaign, indent=2, sort_keys=True))


def parser() -> argparse.ArgumentParser:
    """Build the command-line parser."""
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--metrics", type=Path, default=DEFAULT_METRICS)
    result.add_argument("--campaign-id", default=DEFAULT_CAMPAIGN)
    result.add_argument("--controller-root", type=Path, default=DEFAULT_CONTROLLER_ROOT)
    result.add_argument("--worker-root", type=Path, default=DEFAULT_WORKER_ROOT)
    result.add_argument("--per-cohort", type=int, default=6)
    return result


if __name__ == "__main__":
    prepare(parser().parse_args())
