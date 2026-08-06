#!/usr/bin/env python3
"""Run a dense supplementary-annotation stress pilot for fastVEP/VCFcache."""

# This is an operational benchmark runner, not an installed public API.
# ruff: noqa: D101,D103,E501

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shlex
import shutil
import subprocess
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from itertools import zip_longest
from pathlib import Path
from typing import Sequence

import yaml  # type: ignore[import-untyped]

CUSTOM_SOURCES = {
    "dense_af_01": "AF",
    "dense_af_02": "AF",
    "dense_counts_01": "AC,AN",
    "dense_counts_02": "AC,AN",
    "dense_population_01": "AF,AC,AN",
    "dense_population_02": "AF,AC,AN",
}
PROJECTED_INFO = ("CSQ", "FV_CLINVAR", "FV_GNOMAD", "FV_1KG", "FV_TOPMED")
HIT_RATES = (90, 100)


@dataclass(frozen=True)
class Config:
    pilot: Path
    repo: Path
    source: Path
    clinvar: Path
    threads: int

    @property
    def root(self) -> Path:
        return self.pilot / "heavy"

    @property
    def fastvep(self) -> Path:
        return self.pilot / "tools/fastVEP/target/release/fastvep"

    @property
    def environment(self) -> dict[str, object]:
        return json.loads((self.pilot / "environment.json").read_text())

    @property
    def input_vcf(self) -> Path:
        return Path(str(self.environment["input_vcf"]))

    @property
    def transcript_cache(self) -> Path:
        return Path(str(self.environment["transcript_cache"]))

    @property
    def fasta(self) -> Path:
        return Path(str(self.environment["fasta"]))

    @property
    def sa_dir(self) -> Path:
        return self.root / "sa"

    @property
    def recipe(self) -> Path:
        return self.root / "config/annotation.yaml"

    @property
    def params(self) -> Path:
        return self.root / "config/params.yaml"

    @property
    def python(self) -> Path:
        return self.repo / ".venv/bin/python"


def run(
    command: Sequence[str | Path],
    *,
    env: dict[str, str] | None = None,
    stdout: int | None = subprocess.PIPE,
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [str(value) for value in command],
        env=env,
        stdout=stdout,
        stderr=subprocess.PIPE,
        check=True,
        text=True,
    )


def vcfcache_env(config: Config) -> dict[str, str]:
    env = os.environ.copy()
    env.update(
        {
            "PYTHONPATH": str(config.source),
            "RAYON_NUM_THREADS": str(config.threads),
            "LC_ALL": "C",
            "LANG": "C",
            "TMPDIR": str(config.root / "tmp"),
        }
    )
    return env


def write_configs(config: Config) -> None:
    config.recipe.parent.mkdir(parents=True, exist_ok=True)
    annotation = """genome_build: "GRCh38"
annotation_cmd: |
  ${params.bcftools_cmd} view \\${INPUT_BCF} -Ov | \\
  env RAYON_NUM_THREADS=${params.rayon_threads} ${params.annotation_tool_cmd} annotate \\
    --input - \\
    --output - \\
    --output-format vcf \\
    --transcript-cache ${params.transcript_cache} \\
    --fasta ${params.fasta} \\
    --sa-dir ${params.sa_dir} \\
    --hgvs --symbol --canonical --no-progress \\
    2> \\${AUXILIARY_DIR}/fastvep_stderr.txt | \\
  ${params.bcftools_cmd} view -Ob -o \\${OUTPUT_BCF} -W --threads ${params.threads}

must_contain_info_tag: CSQ
required_tool_version: "0.3.0"
optional_checks: {}
"""
    params = {
        "genome_build": "GRCh38",
        "annotation_tool_cmd": str(config.fastvep),
        "tool_version_command": f"{config.fastvep} --version",
        "bcftools_cmd": "bcftools",
        "temp_dir": str(config.root / "tmp"),
        "threads": config.threads,
        "rayon_threads": config.threads,
        "transcript_cache": str(config.transcript_cache),
        "fasta": str(config.fasta),
        "sa_dir": str(config.sa_dir),
        "optional_checks": {},
    }
    config.recipe.write_text(annotation)
    config.params.write_text(yaml.safe_dump(params, sort_keys=False))


def build_sa(
    config: Config,
    source: str,
    input_path: Path,
    output: Path,
    *,
    name: str | None = None,
    info_fields: str | None = None,
) -> None:
    if output.with_suffix(".osa").exists() or output.with_suffix(".osa2").exists():
        return
    command: list[str | Path] = [
        config.fastvep,
        "sa-build",
        "--source",
        source,
        "--input",
        input_path,
        "--output",
        output,
        "--format",
        "auto",
        "--no-progress",
    ]
    if name:
        command.extend(["--name", name])
    if info_fields:
        command.extend(["--info-fields", info_fields])
    started = time.monotonic()
    completed = run(command)
    with (config.root / "setup.log").open("a") as handle:
        handle.write(
            f"{shlex.join(map(str, command))}\n"
            f"wall_seconds={time.monotonic() - started:.6f}\n"
            f"{completed.stdout}{completed.stderr}\n"
        )


def fastvep_command(
    config: Config, input_path: Path, output: Path, fmt: str
) -> list[str | Path]:
    return [
        config.fastvep,
        "annotate",
        "--input",
        input_path,
        "--output",
        output,
        "--output-format",
        fmt,
        "--transcript-cache",
        config.transcript_cache,
        "--fasta",
        config.fasta,
        "--sa-dir",
        config.sa_dir,
        "--hgvs",
        "--symbol",
        "--canonical",
        "--no-progress",
    ]


def setup(config: Config) -> None:
    for directory in (
        config.root,
        config.sa_dir,
        config.root / "tmp",
        config.root / "smoke",
        config.root / "blueprint",
        config.root / "cache",
        config.root / "workflow",
    ):
        directory.mkdir(parents=True, exist_ok=True)
    if not config.source.joinpath("vcfcache/__init__.py").exists():
        raise FileNotFoundError(
            f"Benchmark source snapshot missing: {config.source}; sync the checkout first"
        )
    write_configs(config)
    shutil.copy2(config.params, config.root / "workflow/init.yaml")
    for name, fields in CUSTOM_SOURCES.items():
        build_sa(
            config,
            "custom_vcf",
            config.input_vcf,
            config.sa_dir / name,
            name=name,
            info_fields=fields,
        )
    for source in ("gnomad", "onekg", "topmed"):
        build_sa(config, source, config.input_vcf, config.sa_dir / f"dense_{source}")
    build_sa(config, "clinvar", config.clinvar, config.sa_dir / "clinvar_real")

    smoke_input = config.pilot / "smoke/input.chr22_20_25mb.v3.vcf.gz"
    smoke_vcf = config.root / "smoke/stress.vcf"
    smoke_json = config.root / "smoke/stress.json"
    if not smoke_vcf.exists():
        run(fastvep_command(config, smoke_input, smoke_vcf, "vcf"))
    if not smoke_json.exists():
        run(fastvep_command(config, smoke_input, smoke_json, "json"))
    header = run(["bcftools", "view", "-h", smoke_vcf]).stdout
    json_text = smoke_json.read_text()
    projected_info_present = {
        tag: f"##INFO=<ID={tag}," in header for tag in PROJECTED_INFO
    }
    custom_sources_present = {name: f'"{name}"' in json_text for name in CUSTOM_SOURCES}
    gate = {
        "projected_info_present": projected_info_present,
        "custom_sources_present_in_json": custom_sources_present,
        "note": (
            "fastVEP 0.3.0 computes custom_vcf sources and exposes them in JSON, "
            "but does not project arbitrary custom sources into VCF; official sources "
            "exercise VCF output construction."
        ),
    }
    gate["pass"] = all(projected_info_present.values()) and all(
        custom_sources_present.values()
    )
    (config.root / "smoke/gate.json").write_text(
        json.dumps(gate, indent=2, sort_keys=True) + "\n"
    )
    if not gate["pass"]:
        raise RuntimeError("Heavy fastVEP smoke gate failed")
    manifest = {
        "input": str(config.input_vcf),
        "clinvar": str(config.clinvar),
        "threads": config.threads,
        "custom_dense_sources": CUSTOM_SOURCES,
        "synthetic_projected_sources": ["gnomad", "onekg", "topmed"],
        "real_projected_sources": ["clinvar"],
        "purpose": (
            "A deliberately dense upper-bound stress configuration. The personal WGS "
            "is reused to make lookup coverage deterministic; it is not a biological "
            "population-reference benchmark."
        ),
    }
    (config.root / "source_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )


def placeholder_cache(config: Config) -> Path:
    cache = config.root / "cache/placeholder"
    output = cache / "vcfcache_annotated.bcf"
    if Path(f"{output}.csi").exists():
        return cache
    cache.mkdir(parents=True, exist_ok=True)
    shutil.copy2(config.recipe, cache / "annotation.yaml")
    shutil.copy2(config.params, cache / "params.snapshot.yaml")
    smoke = config.root / "smoke/stress.vcf"
    keep = ",".join(f"INFO/{tag}" for tag in PROJECTED_INFO)
    shell = (
        f"bcftools view -G -Ou {shlex.quote(str(smoke))} | "
        f"bcftools annotate -x {shlex.quote('^' + keep)} -Ob -o {shlex.quote(str(output))}"
    )
    subprocess.run(["bash", "-euo", "pipefail", "-c", shell], check=True)
    run(["bcftools", "index", "--force", output])
    return cache


def vcfcache_command(
    config: Config, cache: Path, run_dir: Path, *, uncached: bool
) -> list[str | Path]:
    command: list[str | Path] = [
        config.python,
        "-m",
        "vcfcache",
        "annotate",
        "-a",
        cache,
        "-i",
        config.input_vcf,
        "-o",
        run_dir / "output.bcf",
        "--stats-dir",
        run_dir / "stats",
        "--statistics",
        "light",
        "-y",
        config.params,
        "--skip-split-multiallelic",
        "--force",
    ]
    if uncached:
        command.append("--uncached")
    return command


def timed_cell(
    config: Config, name: str, cache: Path, *, uncached: bool
) -> dict[str, object]:
    run_dir = config.root / "runs" / name
    metrics_path = run_dir / "metrics.json"
    if metrics_path.exists():
        previous = json.loads(metrics_path.read_text())
        if previous.get("returncode") == 0:
            return previous
        failed_root = config.root / "runs_failed"
        failed_root.mkdir(parents=True, exist_ok=True)
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        run_dir.rename(failed_root / f"{name}_{stamp}")
    if run_dir.exists():
        raise FileExistsError(f"Incomplete run exists: {run_dir}")
    run_dir.mkdir(parents=True)
    command = vcfcache_command(config, cache, run_dir, uncached=uncached)
    started = time.monotonic()
    with (run_dir / "command.log").open("w") as log:
        completed = subprocess.run(
            [
                "/usr/bin/time",
                "--verbose",
                "--output",
                run_dir / "resource_usage.txt",
                *map(str, command),
            ],
            env=vcfcache_env(config),
            cwd=config.root,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
            text=True,
        )
    metrics = {
        "command": shlex.join(map(str, command)),
        "returncode": completed.returncode,
        "wall_seconds": time.monotonic() - started,
        "output": str(run_dir / "output.bcf"),
    }
    metrics_path.write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n")
    if completed.returncode:
        raise RuntimeError(f"Heavy cell {name} failed; see {run_dir / 'command.log'}")
    return metrics


def build_hit_caches(config: Config, direct: Path) -> None:
    full = config.root / "cache/full/vcfcache_annotated.bcf"
    if not Path(f"{full}.csi").exists():
        full.parent.mkdir(parents=True, exist_ok=True)
        keep = ",".join(f"INFO/{tag}" for tag in PROJECTED_INFO)
        shell = (
            f"bcftools view -G -Ou {shlex.quote(str(direct))} | "
            f"bcftools annotate -x {shlex.quote('^' + keep)} -Ob -o {shlex.quote(str(full))}"
        )
        subprocess.run(["bash", "-euo", "pipefail", "-c", shell], check=True)
        run(["bcftools", "index", "--force", full])
    for rate in HIT_RATES:
        cache = config.root / f"cache/hit-{rate:03d}"
        output = cache / "vcfcache_annotated.bcf"
        if Path(f"{output}.csi").exists():
            continue
        cache.mkdir(parents=True, exist_ok=True)
        shutil.copy2(config.recipe, cache / "annotation.yaml")
        shutil.copy2(config.params, cache / "params.snapshot.yaml")
        if rate == 100:
            os.link(full, output)
        else:
            run(
                ["bcftools", "view", "-i", f"POS%100<{rate}", "-Ob", "-o", output, full]
            )
        run(["bcftools", "index", "--force", output])


def canonical_info(value: str) -> str:
    fields = []
    for item in value.split(";") if value not in ("", ".") else []:
        if item.startswith("CSQ="):
            item = "CSQ=" + ",".join(sorted(item[4:].split(",")))
        fields.append(item)
    return ";".join(sorted(fields)) if fields else value


def compare(direct: Path, cached: Path) -> dict[str, object]:
    processes = [
        subprocess.Popen(
            ["bcftools", "view", "-H", path], stdout=subprocess.PIPE, text=True
        )
        for path in (direct, cached)
    ]
    assert processes[0].stdout is not None and processes[1].stdout is not None
    digests = [hashlib.sha256(), hashlib.sha256()]
    mismatches = 0
    records = 0
    for pair in zip_longest(processes[0].stdout, processes[1].stdout):
        records += 1
        normalized = []
        for line in pair:
            if line is None:
                normalized.append(None)
                continue
            columns = line.rstrip("\n").split("\t")
            columns[7] = canonical_info(columns[7])
            normalized.append("\t".join(columns) + "\n")
        for digest, line in zip(digests, normalized, strict=True):
            if line is not None:
                digest.update(line.encode())
        mismatches += normalized[0] != normalized[1]
    returncodes = [process.wait() for process in processes]
    result = {
        "records": records,
        "mismatches": mismatches,
        "sha256": [digest.hexdigest() for digest in digests],
        "returncodes": returncodes,
    }
    result["pass"] = mismatches == 0 and result["sha256"][0] == result["sha256"][1]  # type: ignore[index]
    return result


def benchmark(config: Config) -> None:
    placeholder = placeholder_cache(config)
    direct = timed_cell(config, "direct", placeholder, uncached=True)
    direct_output = Path(str(direct["output"]))
    build_hit_caches(config, direct_output)
    rows = [
        {"condition": "direct", "wall_seconds": direct["wall_seconds"], "speedup": 1.0}
    ]
    for rate in HIT_RATES:
        cell = timed_cell(
            config,
            f"hit-{rate:03d}",
            config.root / f"cache/hit-{rate:03d}",
            uncached=False,
        )
        equality = compare(direct_output, Path(str(cell["output"])))
        equality_path = config.root / f"runs/hit-{rate:03d}/equality.json"
        equality_path.write_text(json.dumps(equality, indent=2, sort_keys=True) + "\n")
        if not equality["pass"]:
            raise RuntimeError(f"Heavy {rate}% output equality failed")
        rows.append(
            {
                "condition": f"hit-{rate:03d}",
                "wall_seconds": cell["wall_seconds"],
                "speedup": float(direct["wall_seconds"]) / float(cell["wall_seconds"]),
                "semantic_pass": True,
            }
        )
    summary = {"exploratory": True, "rows": rows}
    (config.root / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    with (config.root / "summary.tsv").open("w") as handle:
        handle.write("condition\twall_seconds\tspeedup\tsemantic_pass\n")
        for row in rows:
            handle.write(
                f"{row['condition']}\t{row['wall_seconds']}\t{row['speedup']}\t"
                f"{row.get('semantic_pass', True)}\n"
            )
    print(json.dumps(summary, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("setup", "benchmark", "all", "status"))
    parser.add_argument(
        "--pilot-root",
        type=Path,
        default=Path("/mnt/data/vcfcache_benchmarks/fastvep_pilot"),
    )
    parser.add_argument(
        "--repo", type=Path, default=Path("/mnt/data/home/appuser-projects/vcfcache")
    )
    parser.add_argument("--source", type=Path)
    parser.add_argument(
        "--clinvar",
        type=Path,
        default=Path("/mnt/data/resources/clinvar/vcf_GRCh38/clinvar_20241111.vcf.gz"),
    )
    parser.add_argument("--threads", type=int, default=16)
    args = parser.parse_args()
    source = args.source or args.pilot_root / "heavy/source"
    config = Config(
        args.pilot_root.resolve(),
        args.repo.resolve(),
        source.resolve(),
        args.clinvar.resolve(),
        args.threads,
    )
    if args.command in ("setup", "all"):
        setup(config)
    if args.command in ("benchmark", "all"):
        benchmark(config)
    if args.command == "status":
        print(
            json.dumps(
                {
                    "setup": (config.root / "source_manifest.json").exists(),
                    "direct": (config.root / "runs/direct/metrics.json").exists(),
                    "hit_090": (config.root / "runs/hit-090/metrics.json").exists(),
                    "hit_100": (config.root / "runs/hit-100/metrics.json").exists(),
                    "summary": (config.root / "summary.json").exists(),
                },
                indent=2,
            )
        )


if __name__ == "__main__":
    main()
