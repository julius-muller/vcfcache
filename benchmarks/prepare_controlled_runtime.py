#!/usr/bin/env python3
"""Prepare deterministic self-caches for the controlled runtime experiment."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from benchmarks.run_pilot import (
    PilotConfig,
    preflight,
    run_one,
    sha256sum,
    write_json_atomic,
)

DEFAULT_DATA_ROOT = Path("/mnt/data/vcfcache_benchmarks")
DEFAULT_OUTPUT = DEFAULT_DATA_ROOT / "controlled_runtime"
DEFAULT_INPUT = (
    DEFAULT_DATA_ROOT
    / "samples/GRCh38/wes_twist_core/HG02374.GRCh38.wes_twist_core.vcf.gz"
)
DEFAULT_BUNDLED_CACHE = (
    DEFAULT_DATA_ROOT
    / "bundled_zenodo_caches/gnomad_v4.1_GRCh38_joint_af001/cache/vep115.2_everything"
)
HIT_RATES = (50, 80, 90, 95, 100)
PIPELINES = {
    "vanilla": None,
    "everything": None,
    "delay_5ms": 5_000,
    "delay_20ms": 20_000,
}


def annotation_yaml(source: str, delay_us: int | None = None) -> str:
    """Return a vanilla or synthetic-delay recipe derived from `--everything`."""
    lines = [
        line
        for line in source.splitlines()
        if not line.strip().startswith(("--hgvsg", "--everything"))
    ]
    if delay_us is not None:
        marker = next(
            index for index, line in enumerate(lines) if "--stats_file" in line
        )
        lines[marker:marker] = [
            "    --dir_plugins /mnt/data/vcfcache_benchmarks/controlled_runtime/plugins \\",
            f"    --plugin SyntheticDelay,delay_us={delay_us} \\",
        ]
    return "\n".join(lines) + "\n"


def _run(args: list[str | Path]) -> None:
    subprocess.run([str(value) for value in args], check=True)


def _records(path: Path) -> int:
    return int(
        subprocess.run(
            ["bcftools", "index", "--nrecords", str(path)],
            check=True,
            capture_output=True,
            text=True,
        ).stdout
    )


def _placeholder_cache(
    root: Path, annotation: str, params: Path, bundled_cache: Path
) -> Path:
    cache = root / "cache" / "placeholder"
    cache.mkdir(parents=True)
    (root / "blueprint").mkdir()
    workflow = root / "workflow"
    workflow.mkdir()
    shutil.copy2(
        bundled_cache.parents[1] / "workflow/init.yaml", workflow / "init.yaml"
    )
    (cache / "annotation.yaml").write_text(annotation)
    shutil.copy2(params, cache / "params.snapshot.yaml")
    os.symlink(
        bundled_cache / "vcfcache_annotated.bcf",
        cache / "vcfcache_annotated.bcf",
    )
    os.symlink(
        bundled_cache / "vcfcache_annotated.bcf.csi",
        cache / "vcfcache_annotated.bcf.csi",
    )
    return cache


def _uncached_reference(work: Path, input_vcf: Path, cache: Path, params: Path) -> Path:
    config = PilotConfig(
        data_root=work,
        input_vcf=input_vcf,
        cache_dir=cache,
        params_file=params,
        replicate=1,
    )
    preflight(config)
    run_one(config, "uncached")
    return config.run_dir("uncached") / "output.bcf"


def _sanitize_cache(source: Path, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    view = subprocess.Popen(
        ["bcftools", "view", "-G", "-Ou", str(source)], stdout=subprocess.PIPE
    )
    assert view.stdout is not None
    annotate = subprocess.run(
        [
            "bcftools",
            "annotate",
            "-x",
            "^INFO/CSQ",
            "-Ob",
            "-o",
            str(output),
        ],
        stdin=view.stdout,
        check=True,
    )
    view.stdout.close()
    if view.wait() or annotate.returncode:
        raise RuntimeError(f"Failed to sanitize controlled cache source: {source}")
    _run(["bcftools", "index", "--force", output])


def _subset(source: Path, output: Path, hit_rate: int) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    if hit_rate == 100:
        os.link(source, output)
    else:
        _run(
            [
                "bcftools",
                "view",
                "-i",
                f"POS%100<{hit_rate}",
                "-Ob",
                "-o",
                output,
                source,
            ]
        )
    _run(["bcftools", "index", "--force", output])


def prepare(args: argparse.Namespace) -> None:
    """Build full self-annotations and nested deterministic hit-rate caches."""
    output = args.output.resolve()
    partial = output.with_name(output.name + ".partial")
    if output.exists():
        ready = output / "READY.json"
        if ready.exists():
            print(f"Controlled runtime caches already ready: {output}")
            return
        raise FileExistsError(f"Incomplete controlled runtime root exists: {output}")
    if partial.exists():
        raise FileExistsError(f"Partial controlled runtime root exists: {partial}")
    partial.mkdir(parents=True)
    input_vcf = args.input.resolve()
    bundled = args.bundled_cache.resolve()
    params = bundled / "params.snapshot.yaml"
    source_annotation = (bundled / "annotation.yaml").read_text()
    recipes = {
        "everything": source_annotation,
        "vanilla": annotation_yaml(source_annotation),
        "delay_5ms": annotation_yaml(source_annotation, PIPELINES["delay_5ms"]),
        "delay_20ms": annotation_yaml(source_annotation, PIPELINES["delay_20ms"]),
    }
    plugins = partial / "plugins"
    plugins.mkdir()
    shutil.copy2(
        Path(__file__).with_name("vep_plugins") / "SyntheticDelay.pm",
        plugins / "SyntheticDelay.pm",
    )
    work = partial / "work"
    everything_output = _uncached_reference(
        work / "everything", input_vcf, bundled, params
    )
    vanilla_placeholder = _placeholder_cache(
        work / "vanilla-cache", recipes["vanilla"], params, bundled
    )
    vanilla_output = _uncached_reference(
        work / "vanilla", input_vcf, vanilla_placeholder, params
    )
    full_sources = {}
    for name, source in (
        ("everything", everything_output),
        ("vanilla", vanilla_output),
    ):
        destination = partial / "full" / name / "vcfcache_annotated.bcf"
        _sanitize_cache(source, destination)
        full_sources[name] = destination

    input_records = _records(input_vcf)
    cache_rows = []
    for pipeline, delay_us in PIPELINES.items():
        annotation_source = "everything" if pipeline == "everything" else "vanilla"
        for hit_rate in HIT_RATES:
            cache = partial / "caches" / pipeline / f"hit-{hit_rate:03d}"
            cache.mkdir(parents=True)
            (cache / "annotation.yaml").write_text(recipes[pipeline])
            shutil.copy2(params, cache / "params.snapshot.yaml")
            bcf = cache / "vcfcache_annotated.bcf"
            source_bcf = full_sources[annotation_source]
            if pipeline.startswith("delay_"):
                vanilla_bcf = (
                    partial
                    / "caches"
                    / "vanilla"
                    / f"hit-{hit_rate:03d}"
                    / "vcfcache_annotated.bcf"
                )
                os.link(vanilla_bcf, bcf)
                os.link(Path(f"{vanilla_bcf}.csi"), Path(f"{bcf}.csi"))
            else:
                _subset(source_bcf, bcf, hit_rate)
            records = _records(bcf)
            metadata = {
                "pipeline": pipeline,
                "delay_us_per_transcript_consequence": delay_us,
                "target_hit_rate": hit_rate / 100,
                "cache_records": records,
                "input_records": input_records,
                "construction_rule": (
                    "all records" if hit_rate == 100 else f"POS modulo 100 < {hit_rate}"
                ),
                "annotation_yaml_sha256": sha256sum(cache / "annotation.yaml"),
                "cache_bcf_sha256": sha256sum(bcf),
            }
            write_json_atomic(cache / "controlled_cache.json", metadata)
            cache_rows.append(metadata)
    shutil.rmtree(work)
    ready = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "input_vcf": str(input_vcf),
        "input_sha256": sha256sum(input_vcf),
        "input_records": input_records,
        "hit_rates": list(HIT_RATES),
        "pipelines": PIPELINES,
        "cache_count": len(cache_rows),
        "complete": True,
    }
    write_json_atomic(partial / "READY.json", ready)
    partial.replace(output)
    print(json.dumps(ready, indent=2, sort_keys=True))


def parser() -> argparse.ArgumentParser:
    """Build the controlled-cache preparation CLI."""
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    result.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    result.add_argument("--bundled-cache", type=Path, default=DEFAULT_BUNDLED_CACHE)
    result.set_defaults(function=prepare)
    return result


def main() -> None:
    """Prepare controlled caches."""
    args = parser().parse_args()
    args.function(args)


if __name__ == "__main__":
    main()
