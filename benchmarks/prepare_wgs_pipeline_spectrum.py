#!/usr/bin/env python3
"""Prepare Zenodo-cache wrappers for a real-WGS pipeline-cost spectrum."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import yaml

from benchmarks.run_pilot import sha256sum, write_json_atomic

DEFAULT_ROOT = Path("/mnt/data/vcfcache_benchmarks/wgs_pipeline_spectrum")
DEFAULT_INPUT = Path(
    "/mnt/data/vcfcache_benchmarks/external_wgs/samples/GRCh38/kpgp/"
    "KPGP-00319.GRCh38.small_variants.vcf.gz"
)
DEFAULT_SOURCE_CACHE = Path(
    "/mnt/data/vcfcache_benchmarks/bundled_zenodo_caches/"
    "gnomad_v4.1_GRCh38_joint_af010/cache/vep115.2_everything"
)
DELAYS_US = (500, 1_000, 2_000, 4_000, 7_000, 10_000)
DEFAULT_EXPECTED_HIT_RATE = 0.8023469328603546
DEFAULT_SOURCE_ALIAS = "cache-gnomad-v4.1-GRCh38-joint-af01-vep115.2-e"
DEFAULT_SOURCE_DOI = "10.5281/zenodo.18189447"


def delayed_annotation(source: str, delay_us: int, plugin_dir: Path) -> str:
    """Add the no-output delay plugin while preserving VEP --everything."""
    lines = source.splitlines()
    marker = next(
        index for index, line in enumerate(lines) if "--stats_file" in line
    )
    lines[marker:marker] = [
        f"    --dir_plugins {plugin_dir} \\",
        f"    --plugin SyntheticDelay,delay_us={delay_us} \\",
    ]
    return "\n".join(lines) + "\n"


def record_count(path: Path) -> int:
    """Return the number of indexed input records."""
    return int(
        subprocess.run(
            ["bcftools", "index", "--nrecords", str(path)],
            check=True,
            capture_output=True,
            text=True,
        ).stdout
    )


def prepare(args: argparse.Namespace) -> None:
    """Create lightweight cache wrappers and immutable provenance."""
    root = args.output.resolve()
    if root.exists():
        if (root / "READY.json").is_file():
            print(f"WGS pipeline spectrum already ready: {root}")
            return
        raise FileExistsError(f"Incomplete WGS pipeline spectrum root: {root}")
    partial = root.with_name(root.name + ".partial")
    if partial.exists():
        raise FileExistsError(partial)
    source_cache = args.source_cache.resolve()
    input_vcf = args.input.resolve()
    required = (
        input_vcf,
        Path(f"{input_vcf}.tbi"),
        source_cache / "annotation.yaml",
        source_cache / "params.snapshot.yaml",
        source_cache / "vcfcache_annotated.bcf",
        source_cache / "vcfcache_annotated.bcf.csi",
        source_cache.parents[1] / "workflow/init.yaml",
    )
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(missing)

    (partial / "blueprint").mkdir(parents=True)
    (partial / "workflow").mkdir()
    (partial / "cache").mkdir()
    (partial / "plugins").mkdir()
    shutil.copy2(
        source_cache.parents[1] / "workflow/init.yaml",
        partial / "workflow/init.yaml",
    )
    shutil.copy2(
        Path(__file__).with_name("vep_plugins") / "SyntheticDelay.pm",
        partial / "plugins/SyntheticDelay.pm",
    )
    source_annotation = (source_cache / "annotation.yaml").read_text()
    params = yaml.safe_load((source_cache / "params.snapshot.yaml").read_text())
    params["vep_buffer"] = 100_000
    rows = []
    for delay_us in DELAYS_US:
        name = f"everything_delay_{delay_us:05d}us"
        cache = partial / "cache" / name
        cache.mkdir()
        annotation = delayed_annotation(
            source_annotation,
            delay_us,
            root / "plugins",
        )
        (cache / "annotation.yaml").write_text(annotation)
        (cache / "params.snapshot.yaml").write_text(
            yaml.safe_dump(params, sort_keys=False)
        )
        os.symlink(
            source_cache / "vcfcache_annotated.bcf",
            cache / "vcfcache_annotated.bcf",
        )
        os.symlink(
            source_cache / "vcfcache_annotated.bcf.csi",
            cache / "vcfcache_annotated.bcf.csi",
        )
        provenance = {
            "kind": "virtual_pipeline_load_from_bundled_zenodo_cache",
            "source_alias": args.source_alias,
            "source_doi": args.source_doi,
            "source_cache": str(source_cache),
            "source_annotation_yaml_sha256": sha256sum(
                source_cache / "annotation.yaml"
            ),
            "delay_us_per_transcript_consequence": delay_us,
            "plugin_emits_annotations": False,
            "expected_hit_rate_from_prior_light_run": args.expected_hit_rate,
        }
        write_json_atomic(cache / "spectrum_provenance.json", provenance)
        rows.append(
            {
                "name": name,
                "delay_us": delay_us,
                "cache_dir": str(cache).replace(str(partial), str(root), 1),
                "annotation_yaml_sha256": sha256sum(cache / "annotation.yaml"),
                "params_yaml_sha256": sha256sum(cache / "params.snapshot.yaml"),
            }
        )
    ready = {
        "complete": True,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "input_vcf": str(input_vcf),
        "input_sha256": sha256sum(input_vcf),
        "input_records": record_count(input_vcf),
        "sample": "KPGP-00319",
        "cohort": "KPGP",
        "assembly": "GRCh38",
        "source_alias": args.source_alias,
        "source_doi": args.source_doi,
        "source_cache": str(source_cache),
        "prior_observed_hit_rate": args.expected_hit_rate,
        "delays_us": list(DELAYS_US),
        "pipelines": rows,
    }
    write_json_atomic(partial / "READY.json", ready)
    partial.replace(root)
    print(json.dumps(ready, indent=2, sort_keys=True))


def main() -> None:
    """Parse CLI arguments and prepare the spectrum wrappers."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--source-cache", type=Path, default=DEFAULT_SOURCE_CACHE)
    parser.add_argument("--source-alias", default=DEFAULT_SOURCE_ALIAS)
    parser.add_argument("--source-doi", default=DEFAULT_SOURCE_DOI)
    parser.add_argument(
        "--expected-hit-rate", type=float, default=DEFAULT_EXPECTED_HIT_RATE
    )
    prepare(parser.parse_args())


if __name__ == "__main__":
    main()
