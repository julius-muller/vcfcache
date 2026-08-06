#!/usr/bin/env python3
"""Prepare, run, and summarize the exploratory fastVEP/VCFcache pilot."""

# Ruff's public-docstring rules target the installed package API. This file is
# an operational script whose command phases are documented in README.md.
# ruff: noqa: D101, D103

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import platform
import re
import shlex
import shutil
import socket
import subprocess
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from itertools import zip_longest
from pathlib import Path
from typing import Iterable, Sequence

import yaml  # type: ignore[import-untyped]

FASTVEP_REPOSITORY = "https://github.com/Huang-lab/fastVEP.git"
FASTVEP_COMMIT = "e47216cebe3abcd8dff722b7fb0ab1b19d4fcc80"
FASTVEP_VERSION = "0.3.0"
GFF_URL = (
    "https://ftp.ensembl.org/pub/release-115/gff3/homo_sapiens/"
    "Homo_sapiens.GRCh38.115.gff3.gz"
)
HIT_RATES = (100, 90, 80)
CORE_FIELDS = (
    "Allele",
    "Consequence",
    "IMPACT",
    "SYMBOL",
    "Gene",
    "Feature_type",
    "Feature",
    "BIOTYPE",
    "EXON",
    "INTRON",
    "HGVSc",
    "HGVSp",
    "cDNA_position",
    "CDS_position",
    "Protein_position",
    "Amino_acids",
    "Codons",
    "CANONICAL",
)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_suffix(path.suffix + ".partial")
    partial.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    partial.replace(path)


def sha256sum(path: Path, chunk_size: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def command_text(args: Sequence[str | Path]) -> str:
    return shlex.join([str(value) for value in args])


def run_checked(
    args: Sequence[str | Path],
    *,
    cwd: Path | None = None,
    env: dict[str, str] | None = None,
    stdout_path: Path | None = None,
    stderr_path: Path | None = None,
) -> subprocess.CompletedProcess[str]:
    stdout_handle = stdout_path.open("w") if stdout_path else subprocess.PIPE
    stderr_handle = stderr_path.open("w") if stderr_path else subprocess.PIPE
    try:
        return subprocess.run(
            [str(value) for value in args],
            cwd=cwd,
            env=env,
            stdout=stdout_handle,
            stderr=stderr_handle,
            check=True,
            text=True,
        )
    finally:
        if stdout_path:
            stdout_handle.close()  # type: ignore[union-attr]
        if stderr_path:
            stderr_handle.close()  # type: ignore[union-attr]


def run_pipeline(script: str, *, cwd: Path, log: Path) -> None:
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("w") as handle:
        subprocess.run(
            ["/bin/bash", "-euo", "pipefail", "-c", script],
            cwd=cwd,
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=True,
            text=True,
        )


@dataclass(frozen=True)
class Config:
    root: Path
    repo: Path
    input_vcf: Path
    fasta_gz: Path
    vep_sif: Path
    threads: int

    @property
    def tools(self) -> Path:
        return self.root / "tools"

    @property
    def fastvep_source(self) -> Path:
        return self.tools / "fastVEP"

    @property
    def fastvep(self) -> Path:
        return self.fastvep_source / "target/release/fastvep"

    @property
    def cargo_home(self) -> Path:
        return self.tools / "cargo"

    @property
    def rustup_home(self) -> Path:
        return self.tools / "rustup"

    @property
    def cargo(self) -> Path:
        return self.cargo_home / "bin/cargo"

    @property
    def data(self) -> Path:
        return self.root / "data"

    @property
    def gff_gz(self) -> Path:
        return self.data / "Homo_sapiens.GRCh38.115.gff3.gz"

    @property
    def gff(self) -> Path:
        return self.data / "Homo_sapiens.GRCh38.115.gff3"

    @property
    def gff_sorted_bgz(self) -> Path:
        return self.data / "Homo_sapiens.GRCh38.115.sorted.gff3.bgz"

    @property
    def fasta(self) -> Path:
        return self.data / "Homo_sapiens.GRCh38.dna.toplevel.fa"

    @property
    def transcript_cache(self) -> Path:
        return self.data / "Homo_sapiens.GRCh38.115.fastvep.cache"

    @property
    def recipe(self) -> Path:
        return self.root / "config/annotation.yaml"

    @property
    def params(self) -> Path:
        return self.root / "config/params.yaml"

    @property
    def vcfcache(self) -> Path:
        return self.repo / ".venv/bin/vcfcache"

    @property
    def smoke_input(self) -> Path:
        return self.root / "smoke/input.chr22_20_25mb.v3.vcf.gz"


def require_commands(names: Iterable[str]) -> None:
    missing = [name for name in names if shutil.which(name) is None]
    if missing:
        raise RuntimeError(f"Missing required commands: {', '.join(missing)}")


def decompress_once(source: Path, destination: Path) -> None:
    if destination.exists():
        return
    partial = destination.with_suffix(destination.suffix + ".partial")
    with gzip.open(source, "rb") as reader, partial.open("wb") as writer:
        shutil.copyfileobj(reader, writer, length=8 << 20)
    partial.replace(destination)


def prepare_tabix_gff(source: Path, destination: Path) -> None:
    if destination.exists() and Path(f"{destination}.tbi").exists():
        return
    partial = Path(f"{destination}.partial")
    script = (
        f"(grep '^#' {shlex.quote(str(source))}; "
        f"grep -v '^#' {shlex.quote(str(source))} | "
        "LC_ALL=C sort -k1,1V -k4,4n -k5,5n) | "
        f"bgzip -c > {shlex.quote(str(partial))}"
    )
    run_pipeline(
        script, cwd=destination.parent, log=destination.parent / "sort_gff.log"
    )
    partial.replace(destination)
    run_checked(["tabix", "-f", "-p", "gff", destination])


def write_configs(config: Config) -> None:
    config.recipe.parent.mkdir(parents=True, exist_ok=True)
    (config.root / "blueprint").mkdir(parents=True, exist_ok=True)
    (config.root / "cache").mkdir(parents=True, exist_ok=True)
    (config.root / "workflow").mkdir(parents=True, exist_ok=True)
    annotation = f"""genome_build: "GRCh38"
annotation_cmd: |
  ${{params.bcftools_cmd}} view \\${{INPUT_BCF}} -Ov | \\
  env RAYON_NUM_THREADS=${{params.rayon_threads}} ${{params.annotation_tool_cmd}} annotate \\
    --input - \\
    --output - \\
    --output-format vcf \\
    --transcript-cache ${{params.transcript_cache}} \\
    --fasta ${{params.fasta}} \\
    --hgvs --symbol --canonical --no-progress \\
    2> \\${{AUXILIARY_DIR}}/fastvep_stderr.txt | \\
  ${{params.bcftools_cmd}} view -Ob -o \\${{OUTPUT_BCF}} -W --threads ${{params.threads}}

must_contain_info_tag: CSQ
required_tool_version: "{FASTVEP_VERSION}"
optional_checks: {{}}
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
        "optional_checks": {},
    }
    config.recipe.write_text(annotation)
    config.params.write_text(yaml.safe_dump(params, sort_keys=False))
    shutil.copy2(config.params, config.root / "workflow/init.yaml")


def rust_environment(config: Config) -> dict[str, str]:
    environment = os.environ.copy()
    environment.update(
        {
            "CARGO_HOME": str(config.cargo_home),
            "RUSTUP_HOME": str(config.rustup_home),
            "PATH": f"{config.cargo_home / 'bin'}:{environment['PATH']}",
        }
    )
    return environment


def ensure_rust_toolchain(config: Config) -> None:
    if config.cargo.exists():
        return
    rustup_init = config.tools / "rustup-init"
    if not rustup_init.exists():
        run_checked(
            [
                "curl",
                "-L",
                "--fail",
                "--retry",
                "3",
                "-o",
                rustup_init,
                "https://static.rust-lang.org/rustup/dist/x86_64-unknown-linux-gnu/rustup-init",
            ]
        )
        rustup_init.chmod(0o755)
    run_checked(
        [
            rustup_init,
            "-y",
            "--profile",
            "minimal",
            "--default-toolchain",
            "stable",
            "--no-modify-path",
        ],
        env=rust_environment(config),
    )


def setup(config: Config) -> None:
    require_commands(("git", "curl", "bcftools", "samtools", "bgzip", "tabix"))
    for directory in (config.root, config.tools, config.data, config.root / "tmp"):
        directory.mkdir(parents=True, exist_ok=True)

    print("[setup] Preparing pinned fastVEP source", flush=True)
    if not config.fastvep_source.exists():
        run_checked(["git", "clone", FASTVEP_REPOSITORY, config.fastvep_source])
    run_checked(["git", "fetch", "origin", FASTVEP_COMMIT], cwd=config.fastvep_source)
    run_checked(
        ["git", "checkout", "--detach", FASTVEP_COMMIT], cwd=config.fastvep_source
    )
    actual_commit = run_checked(
        ["git", "rev-parse", "HEAD"], cwd=config.fastvep_source
    ).stdout.strip()
    if actual_commit != FASTVEP_COMMIT:
        raise RuntimeError(f"Unexpected fastVEP commit: {actual_commit}")

    print("[setup] Preparing Ensembl 115 GFF3 and FASTA", flush=True)
    if not config.gff_gz.exists():
        partial = Path(f"{config.gff_gz}.partial")
        run_checked(["curl", "-L", "--fail", "--retry", "3", "-o", partial, GFF_URL])
        partial.replace(config.gff_gz)
    decompress_once(config.gff_gz, config.gff)
    prepare_tabix_gff(config.gff, config.gff_sorted_bgz)
    decompress_once(config.fasta_gz, config.fasta)
    if not Path(f"{config.fasta}.fai").exists():
        run_checked(["samtools", "faidx", config.fasta])

    print("[setup] Preparing a current user-local Rust toolchain", flush=True)
    ensure_rust_toolchain(config)
    print("[setup] Building fastVEP release binary", flush=True)
    if not config.fastvep.exists():
        run_checked(
            [config.cargo, "build", "--release", "--bin", "fastvep"],
            cwd=config.fastvep_source,
            env=rust_environment(config),
        )
    version = run_checked([config.fastvep, "--version"]).stdout.strip()
    if FASTVEP_VERSION not in version:
        raise RuntimeError(f"Unexpected fastVEP version: {version}")
    print("[setup] Building fastVEP transcript cache", flush=True)
    if not config.transcript_cache.exists():
        run_checked(
            [
                config.fastvep,
                "cache",
                "--gff3",
                config.gff,
                "--fasta",
                config.fasta,
                "--output",
                config.transcript_cache,
                "--no-progress",
            ]
        )
    write_configs(config)
    manifest = {
        "created_at": utc_now(),
        "hostname": socket.gethostname(),
        "platform": platform.platform(),
        "cpu_count": os.cpu_count(),
        "threads": config.threads,
        "fastvep_repository": FASTVEP_REPOSITORY,
        "fastvep_commit": actual_commit,
        "fastvep_version": version,
        "cargo_version": run_checked(
            [config.cargo, "--version"], env=rust_environment(config)
        ).stdout.strip(),
        "rustc_version": run_checked(
            [config.cargo_home / "bin/rustc", "--version"],
            env=rust_environment(config),
        ).stdout.strip(),
        "input_vcf": str(config.input_vcf),
        "input_records": int(
            run_checked(["bcftools", "index", "-n", config.input_vcf]).stdout
        ),
        "gff": str(config.gff),
        "fasta": str(config.fasta),
        "transcript_cache": str(config.transcript_cache),
        "bcftools_version": run_checked(["bcftools", "--version"]).stdout.splitlines()[
            0
        ],
    }
    write_json(config.root / "environment.json", manifest)
    print(f"Setup ready: {config.root}")


def fastvep_command(
    config: Config, input_vcf: Path, output_vcf: Path, flags: Sequence[str]
) -> list[str | Path]:
    return [
        config.fastvep,
        "annotate",
        "--input",
        input_vcf,
        "--output",
        output_vcf,
        "--output-format",
        "vcf",
        "--transcript-cache",
        config.transcript_cache,
        "--fasta",
        config.fasta,
        "--no-progress",
        *flags,
    ]


def csq_header(path: Path) -> tuple[str, tuple[str, ...]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("##INFO=<ID=CSQ"):
                match = re.search(r"Format: ([^\"]+)", line)
                return line.rstrip(), tuple(match.group(1).split("|")) if match else ()
            if line.startswith("#CHROM"):
                break
    return "", ()


def vcf_body_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        for line in handle:
            if not line.startswith("#"):
                digest.update(line.encode())
    return digest.hexdigest()


def parse_csq(
    path: Path,
) -> tuple[
    set[tuple[str, ...]], dict[tuple[str, ...], dict[str, str]], tuple[str, ...]
]:
    _, fields = csq_header(path)
    if not fields:
        raise RuntimeError(f"No parseable CSQ header in {path}")
    field_index = {name: index for index, name in enumerate(fields)}
    feature_index = field_index.get("Feature")
    allele_index = field_index.get("Allele", 0)
    if feature_index is None:
        raise RuntimeError(f"No Feature field in CSQ header: {path}")
    variants: set[tuple[str, ...]] = set()
    pairs: dict[tuple[str, ...], dict[str, str]] = {}
    with path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            columns = line.rstrip().split("\t")
            key = (columns[0], columns[1], columns[3], columns[4])
            variants.add(key)
            info = columns[7].split(";")
            csq = next((item[4:] for item in info if item.startswith("CSQ=")), "")
            for entry in csq.split(",") if csq else ():
                values = entry.split("|")
                values.extend([""] * (len(fields) - len(values)))
                pair = key + (values[allele_index], values[feature_index])
                pairs[pair] = dict(zip(fields, values, strict=True))
    return variants, pairs, fields


def compare_annotators(fastvep_vcf: Path, vep_vcf: Path) -> dict[str, object]:
    fast_variants, fast_pairs, fast_fields = parse_csq(fastvep_vcf)
    vep_variants, vep_pairs, vep_fields = parse_csq(vep_vcf)
    shared_pairs = sorted(fast_pairs.keys() & vep_pairs.keys())
    common_core = [
        field for field in CORE_FIELDS if field in fast_fields and field in vep_fields
    ]
    field_comparison: dict[str, dict[str, int | float]] = {}
    for field in common_core:
        matches = sum(
            fast_pairs[pair].get(field, "") == vep_pairs[pair].get(field, "")
            for pair in shared_pairs
        )
        total = len(shared_pairs)
        field_comparison[field] = {
            "matches": matches,
            "total": total,
            "fraction": matches / total if total else 0.0,
        }
    return {
        "fastvep_variants": len(fast_variants),
        "vep_variants": len(vep_variants),
        "shared_variants": len(fast_variants & vep_variants),
        "fastvep_only_variants": len(fast_variants - vep_variants),
        "vep_only_variants": len(vep_variants - fast_variants),
        "fastvep_allele_transcript_pairs": len(fast_pairs),
        "vep_allele_transcript_pairs": len(vep_pairs),
        "shared_allele_transcript_pairs": len(shared_pairs),
        "fastvep_only_allele_transcript_pairs": len(
            fast_pairs.keys() - vep_pairs.keys()
        ),
        "vep_only_allele_transcript_pairs": len(vep_pairs.keys() - fast_pairs.keys()),
        "common_core_fields": common_core,
        "field_comparison": field_comparison,
    }


def smoke(config: Config) -> None:
    if not (config.root / "environment.json").exists():
        raise RuntimeError("Run setup before smoke")
    smoke_root = config.root / "smoke"
    smoke_root.mkdir(parents=True, exist_ok=True)
    print("[smoke] Preparing chr22 test input", flush=True)
    if not config.smoke_input.exists():
        run_checked(
            [
                "bcftools",
                "view",
                "--regions",
                "chr22:20000000-25000000",
                "-Oz",
                "-o",
                config.smoke_input,
                config.input_vcf,
            ]
        )
        run_checked(["bcftools", "index", "--tbi", config.smoke_input])
    smoke_records = int(
        run_checked(["bcftools", "index", "-n", config.smoke_input]).stdout
    )
    if smoke_records == 0:
        raise RuntimeError(f"Smoke input contains no records: {config.smoke_input}")
    base = smoke_root / "fastvep_hgvs.chr22_20_25mb.v3.vcf"
    flags = smoke_root / "fastvep_flags.chr22_20_25mb.v3.vcf"
    print("[smoke] Running fastVEP flag comparison", flush=True)
    if not base.exists():
        run_checked(fastvep_command(config, config.smoke_input, base, ["--hgvs"]))
    if not flags.exists():
        run_checked(
            fastvep_command(
                config,
                config.smoke_input,
                flags,
                ["--hgvs", "--symbol", "--canonical", "--everything"],
            )
        )
    base_header, _ = csq_header(base)
    flags_header, _ = csq_header(flags)
    flag_report = {
        "base_command": command_text(
            fastvep_command(config, config.smoke_input, base, ["--hgvs"])
        ),
        "flags_command": command_text(
            fastvep_command(
                config,
                config.smoke_input,
                flags,
                ["--hgvs", "--symbol", "--canonical", "--everything"],
            )
        ),
        "csq_headers_equal": base_header == flags_header,
        "bodies_equal": vcf_body_sha256(base) == vcf_body_sha256(flags),
        "base_body_sha256": vcf_body_sha256(base),
        "flags_body_sha256": vcf_body_sha256(flags),
    }
    write_json(smoke_root / "flag_behavior.json", flag_report)

    vep_output = smoke_root / "vep_gff.chr22_20_25mb.v3.vcf"
    print("[smoke] Running matched Ensembl VEP GFF comparison", flush=True)
    if not vep_output.exists():
        command: list[str | Path] = [
            "apptainer",
            "exec",
            "-B",
            "/mnt/data:/mnt/data",
            config.vep_sif,
            "vep",
            "--gff",
            config.gff_sorted_bgz,
            "--fasta",
            config.fasta,
            "--species",
            "homo_sapiens",
            "--assembly",
            "GRCh38",
            "--format",
            "vcf",
            "--vcf",
            "--hgvs",
            "--symbol",
            "--canonical",
            "--force_overwrite",
            "--no_stats",
            "--input_file",
            config.smoke_input,
            "--output_file",
            vep_output,
        ]
        run_checked(
            command,
            stdout_path=smoke_root / "vep_stdout.log",
            stderr_path=smoke_root / "vep_stderr.log",
        )
        (smoke_root / "vep_command.txt").write_text(command_text(command) + "\n")
    write_json(
        smoke_root / "annotator_comparison.json", compare_annotators(flags, vep_output)
    )
    make_placeholder_cache(config, flags)
    print(f"Smoke complete: {smoke_root}")


def make_placeholder_cache(config: Config, source_vcf: Path) -> Path:
    cache = config.root / "cache/placeholder"
    if (cache / "vcfcache_annotated.bcf.csi").exists():
        return cache
    cache.mkdir(parents=True, exist_ok=True)
    shutil.copy2(config.recipe, cache / "annotation.yaml")
    shutil.copy2(config.params, cache / "params.snapshot.yaml")
    output = cache / "vcfcache_annotated.bcf"
    script = (
        f"bcftools view -G -Ou {shlex.quote(str(source_vcf))} | "
        f"bcftools annotate -x '^INFO/CSQ' -Ob -o {shlex.quote(str(output))}"
    )
    run_pipeline(script, cwd=cache, log=cache / "prepare.log")
    run_checked(["bcftools", "index", "--force", output])
    return cache


def parse_gnu_time(path: Path) -> dict[str, float | int]:
    values: dict[str, str] = {}
    for line in path.read_text().splitlines():
        if ": " in line:
            key, value = line.strip().rsplit(": ", maxsplit=1)
            values[key] = value
    elapsed_key = next(key for key in values if key.startswith("Elapsed (wall clock)"))
    parts = values[elapsed_key].split(":")
    if len(parts) == 3:
        elapsed = int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
    elif len(parts) == 2:
        elapsed = int(parts[0]) * 60 + float(parts[1])
    else:
        elapsed = float(parts[0])
    return {
        "wall_seconds_gnu": elapsed,
        "user_seconds": float(values["User time (seconds)"]),
        "system_seconds": float(values["System time (seconds)"]),
        "cpu_percent": int(values["Percent of CPU this job got"].rstrip("%")),
        "max_rss_kib": int(values["Maximum resident set size (kbytes)"]),
        "filesystem_inputs": int(values["File system inputs"]),
        "filesystem_outputs": int(values["File system outputs"]),
    }


def collect_stage_timings(run_dir: Path) -> list[dict[str, object]]:
    values: list[dict[str, object]] = []
    for path in sorted(run_dir.rglob("timing.txt")):
        for line_number, line in enumerate(path.read_text().splitlines(), start=1):
            fields = line.split("\t")
            if len(fields) == 2:
                values.append(
                    {
                        "file": str(path.relative_to(run_dir)),
                        "line": line_number,
                        "command": fields[0],
                        "seconds": float(fields[1]),
                    }
                )
    return values


def parse_command_log_timings(path: Path) -> dict[str, object]:
    commands: list[dict[str, object]] = []
    workflow_seconds: float | None = None
    for line in path.read_text().splitlines():
        command_match = re.search(r"Command completed in ([0-9.]+)s: (.+)$", line)
        if command_match:
            commands.append(
                {
                    "command": command_match.group(2),
                    "seconds": float(command_match.group(1)),
                }
            )
        workflow_match = re.search(r"Annotation completed in ([0-9.]+)s", line)
        if workflow_match:
            workflow_seconds = float(workflow_match.group(1))
    return {"commands": commands, "workflow_seconds": workflow_seconds}


def vcfcache_command(
    config: Config, cache: Path, run_dir: Path, *, uncached: bool
) -> list[str | Path]:
    command: list[str | Path] = [
        config.vcfcache,
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
        "full",
        "-y",
        config.params,
        "--skip-split-multiallelic",
        "--force",
    ]
    if uncached:
        command.append("--uncached")
    return command


def run_wgs_cell(
    config: Config, name: str, cache: Path, *, uncached: bool = False
) -> dict[str, object]:
    run_dir = config.root / "runs" / name
    metrics_path = run_dir / "metrics.json"
    if metrics_path.exists():
        print(f"Already complete: {name}")
        return json.loads(metrics_path.read_text())
    if run_dir.exists():
        raise FileExistsError(f"Incomplete run exists: {run_dir}")
    run_dir.mkdir(parents=True)
    command = vcfcache_command(config, cache, run_dir, uncached=uncached)
    write_json(
        run_dir / "command.json",
        {"created_at": utc_now(), "command": command_text(command)},
    )
    time_file = run_dir / "resource_usage.txt"
    log_file = run_dir / "command.log"
    env = os.environ.copy()
    env.update(
        {
            "LC_ALL": "C",
            "LANG": "C",
            "RAYON_NUM_THREADS": str(config.threads),
            "TMPDIR": str(config.root / "tmp"),
        }
    )
    started = time.monotonic()
    with log_file.open("w") as log:
        completed = subprocess.run(
            ["/usr/bin/time", "--verbose", "--output", time_file, *map(str, command)],
            stdout=log,
            stderr=subprocess.STDOUT,
            env=env,
            check=False,
            text=True,
        )
    wall = time.monotonic() - started
    if completed.returncode:
        write_json(
            run_dir / "failure.json",
            {"returncode": completed.returncode, "wall_seconds": wall},
        )
        raise RuntimeError(f"WGS cell failed: {name}; see {log_file}")
    output = run_dir / "output.bcf"
    records = int(run_checked(["bcftools", "index", "-n", output]).stdout)
    metrics: dict[str, object] = {
        "name": name,
        "uncached": uncached,
        "completed_at": utc_now(),
        "wall_seconds": wall,
        "output": str(output),
        "output_records": records,
        "output_bytes": output.stat().st_size,
        "stage_timings": collect_stage_timings(run_dir),
        "log_timings": parse_command_log_timings(log_file),
        **parse_gnu_time(time_file),
    }
    write_json(metrics_path, metrics)
    print(f"Completed {name}: {wall:.1f}s")
    return metrics


def build_hit_caches(config: Config, direct_output: Path) -> None:
    full = config.root / "cache/full/vcfcache_annotated.bcf"
    if not Path(f"{full}.csi").exists():
        full.parent.mkdir(parents=True, exist_ok=True)
        script = (
            f"bcftools view -G -Ou {shlex.quote(str(direct_output))} | "
            f"bcftools annotate -x '^INFO/CSQ' -Ob -o {shlex.quote(str(full))}"
        )
        run_pipeline(script, cwd=full.parent, log=full.parent / "prepare.log")
        run_checked(["bcftools", "index", "--force", full])
    for hit_rate in HIT_RATES:
        cache = config.root / "cache" / f"hit-{hit_rate:03d}"
        output = cache / "vcfcache_annotated.bcf"
        if Path(f"{output}.csi").exists():
            continue
        cache.mkdir(parents=True, exist_ok=True)
        shutil.copy2(config.recipe, cache / "annotation.yaml")
        shutil.copy2(config.params, cache / "params.snapshot.yaml")
        if hit_rate == 100:
            os.link(full, output)
        else:
            run_checked(
                [
                    "bcftools",
                    "view",
                    "-i",
                    f"POS%100<{hit_rate}",
                    "-Ob",
                    "-o",
                    output,
                    full,
                ]
            )
        run_checked(["bcftools", "index", "--force", output])
        write_json(
            cache / "cache.json",
            {
                "target_hit_rate": hit_rate / 100,
                "construction": (
                    "all records" if hit_rate == 100 else f"POS modulo 100 < {hit_rate}"
                ),
                "records": int(run_checked(["bcftools", "index", "-n", output]).stdout),
            },
        )


def canonical_info(value: str) -> str:
    if value in ("", "."):
        return value
    fields: list[str] = []
    for item in value.split(";"):
        if item.startswith("CSQ="):
            item = "CSQ=" + ",".join(sorted(item[4:].split(",")))
        fields.append(item)
    return ";".join(sorted(fields))


def canonical_vcf_stream(path: Path) -> subprocess.Popen[str]:
    return subprocess.Popen(
        ["bcftools", "view", "-H", str(path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )


def compare_bcf(
    direct: Path, cached: Path, mismatch_limit: int = 20
) -> dict[str, object]:
    direct_process = canonical_vcf_stream(direct)
    cached_process = canonical_vcf_stream(cached)
    assert direct_process.stdout is not None and cached_process.stdout is not None
    direct_digest = hashlib.sha256()
    cached_digest = hashlib.sha256()
    mismatches = 0
    records = 0
    examples: list[dict[str, object]] = []
    sentinel = object()
    for number, pair in enumerate(
        zip_longest(direct_process.stdout, cached_process.stdout, fillvalue=sentinel),
        start=1,
    ):
        direct_line, cached_line = pair
        records += 1
        canonical: list[str | object] = []
        for line in (direct_line, cached_line):
            if line is sentinel:
                canonical.append(line)
                continue
            columns = str(line).rstrip("\n").split("\t")
            columns[7] = canonical_info(columns[7])
            canonical.append("\t".join(columns) + "\n")
        if canonical[0] is not sentinel:
            direct_digest.update(str(canonical[0]).encode())
        if canonical[1] is not sentinel:
            cached_digest.update(str(canonical[1]).encode())
        if canonical[0] != canonical[1]:
            mismatches += 1
            if len(examples) < mismatch_limit:
                examples.append(
                    {
                        "record": number,
                        "direct": (
                            None
                            if canonical[0] is sentinel
                            else str(canonical[0])[:1000]
                        ),
                        "cached": (
                            None
                            if canonical[1] is sentinel
                            else str(canonical[1])[:1000]
                        ),
                    }
                )
    direct_stderr = direct_process.communicate()[1]
    cached_stderr = cached_process.communicate()[1]
    if direct_process.returncode or cached_process.returncode:
        raise RuntimeError(
            f"bcftools comparison failed: direct={direct_stderr}; cached={cached_stderr}"
        )
    direct_header, _ = csq_header_from_bcf(direct)
    cached_header, _ = csq_header_from_bcf(cached)
    return {
        "semantic_pass": (
            mismatches == 0
            and direct_digest.hexdigest() == cached_digest.hexdigest()
            and direct_header == cached_header
        ),
        "records_compared": records,
        "record_mismatches": mismatches,
        "csq_headers_equal": direct_header == cached_header,
        "direct_canonical_sha256": direct_digest.hexdigest(),
        "cached_canonical_sha256": cached_digest.hexdigest(),
        "examples": examples,
    }


def csq_header_from_bcf(path: Path) -> tuple[str, tuple[str, ...]]:
    header = run_checked(["bcftools", "view", "-h", path]).stdout
    line = next(
        (value for value in header.splitlines() if value.startswith("##INFO=<ID=CSQ")),
        "",
    )
    match = re.search(r"Format: ([^\"]+)", line)
    return line, tuple(match.group(1).split("|")) if match else ()


def observed_hit_rate(run_dir: Path, output_records: int) -> float | None:
    candidates = list(run_dir.rglob("compare_stats.yaml"))
    if len(candidates) != 1:
        return None
    value = yaml.safe_load(candidates[0].read_text()) or {}
    counts = value.get("variant_counts", {}) or {}
    total = counts.get("input_variants") or counts.get("total_output") or output_records
    misses = counts.get("tool_annotated")
    if isinstance(total, int) and isinstance(misses, int) and total:
        return 1 - misses / total
    return None


def wgs(config: Config) -> None:
    placeholder = config.root / "cache/placeholder"
    if not (placeholder / "vcfcache_annotated.bcf.csi").exists():
        raise RuntimeError("Run smoke before wgs")
    print("[wgs] Running matched direct fastVEP baseline", flush=True)
    direct = run_wgs_cell(config, "direct", placeholder, uncached=True)
    direct_output = Path(str(direct["output"]))
    build_hit_caches(config, direct_output)
    for hit_rate in HIT_RATES:
        print(f"[wgs] Running cached cell at target {hit_rate}% hits", flush=True)
        name = f"hit-{hit_rate:03d}"
        cache = config.root / "cache" / name
        metrics = run_wgs_cell(config, name, cache)
        comparison = compare_bcf(direct_output, Path(str(metrics["output"])))
        write_json(config.root / "runs" / name / "equality.json", comparison)
        if not comparison["semantic_pass"]:
            raise RuntimeError(f"Output equality failed at {hit_rate}% hits")


def collect(config: Config) -> dict[str, object]:
    direct_path = config.root / "runs/direct/metrics.json"
    if not direct_path.exists():
        raise RuntimeError("No completed direct run")
    direct = json.loads(direct_path.read_text())
    if not direct.get("log_timings"):
        direct["log_timings"] = parse_command_log_timings(
            config.root / "runs/direct/command.log"
        )
        write_json(direct_path, direct)
    direct_wall = float(direct["wall_seconds"])
    direct_workflow = float(direct["log_timings"]["workflow_seconds"])
    rows: list[dict[str, object]] = [
        {
            "condition": "direct",
            "target_hit_rate": 0.0,
            "observed_hit_rate": 0.0,
            "wall_seconds": direct_wall,
            "workflow_seconds": direct_workflow,
            "speedup": 1.0,
            "workflow_speedup": 1.0,
            "semantic_pass": True,
        }
    ]
    for hit_rate in sorted(HIT_RATES):
        name = f"hit-{hit_rate:03d}"
        run_dir = config.root / "runs" / name
        metrics = json.loads((run_dir / "metrics.json").read_text())
        if not metrics.get("log_timings"):
            metrics["log_timings"] = parse_command_log_timings(run_dir / "command.log")
            write_json(run_dir / "metrics.json", metrics)
        equality = json.loads((run_dir / "equality.json").read_text())
        wall = float(metrics["wall_seconds"])
        workflow = float(metrics["log_timings"]["workflow_seconds"])
        rows.append(
            {
                "condition": name,
                "target_hit_rate": hit_rate / 100,
                "observed_hit_rate": observed_hit_rate(
                    run_dir, int(metrics["output_records"])
                ),
                "wall_seconds": wall,
                "workflow_seconds": workflow,
                "speedup": direct_wall / wall,
                "workflow_speedup": direct_workflow / workflow,
                "semantic_pass": bool(equality["semantic_pass"]),
            }
        )
    by_hit = {int(round(float(row["target_hit_rate"]) * 100)): row for row in rows}
    gate = {
        "speedup_80_at_least_2x": float(by_hit[80]["speedup"]) >= 2,
        "speedup_90_at_least_3x": float(by_hit[90]["speedup"]) >= 3,
        "all_outputs_equal": all(bool(row["semantic_pass"]) for row in rows),
    }
    gate["pass"] = all(gate.values())
    summary = {
        "created_at": utc_now(),
        "exploratory": True,
        "rows": rows,
        "gate": gate,
        "smoke": {
            "flag_behavior": json.loads(
                (config.root / "smoke/flag_behavior.json").read_text()
            ),
            "annotator_comparison": json.loads(
                (config.root / "smoke/annotator_comparison.json").read_text()
            ),
        },
    }
    write_json(config.root / "summary.json", summary)
    columns = (
        "condition",
        "target_hit_rate",
        "observed_hit_rate",
        "wall_seconds",
        "workflow_seconds",
        "speedup",
        "workflow_speedup",
        "semantic_pass",
    )
    with (config.root / "summary.tsv").open("w") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write(
                "\t".join(str(row.get(column, "")) for column in columns) + "\n"
            )
    comparison = summary["smoke"]["annotator_comparison"]  # type: ignore[index]
    result_lines = [
        "# Exploratory fastVEP pilot results",
        "",
        f"Generated: {summary['created_at']}",
        "",
        "## VCFcache performance",
        "",
        "| Condition | Observed hit rate | Wall time | Speedup | Equal output |",
        "|---|---:|---:|---:|:---:|",
    ]
    for row in rows:
        observed = row["observed_hit_rate"]
        observed_text = "n/a" if observed is None else f"{100 * float(observed):.2f}%"
        result_lines.append(
            f"| {row['condition']} | {observed_text} | {float(row['wall_seconds']):.1f} s | "
            f"{float(row['speedup']):.2f}x | {'yes' if row['semantic_pass'] else 'no'} |"
        )
    result_lines.extend(
        [
            "",
            f"Performance gate: **{'PASS' if gate['pass'] else 'FAIL'}**.",
            "",
            "The default statistics scan is included above. Internal workflow-only "
            "speedups were "
            + ", ".join(
                f"{int(float(row['target_hit_rate']) * 100)}%: "
                f"{float(row['workflow_speedup']):.2f}x"
                for row in rows
                if row["condition"] != "direct"
            )
            + ". These are diagnostic trace ratios, not separately timed `--no-stats` runs.",
            "",
            "## Small-input VEP comparison",
            "",
            f"- Shared variants: {comparison['shared_variants']}.",  # type: ignore[index]
            f"- fastVEP-only variants: {comparison['fastvep_only_variants']}.",  # type: ignore[index]
            f"- VEP-only variants: {comparison['vep_only_variants']}.",  # type: ignore[index]
            f"- Shared allele-transcript pairs: {comparison['shared_allele_transcript_pairs']}.",  # type: ignore[index]
            f"- fastVEP-only pairs: {comparison['fastvep_only_allele_transcript_pairs']}.",  # type: ignore[index]
            f"- VEP-only pairs: {comparison['vep_only_allele_transcript_pairs']}.",  # type: ignore[index]
            "",
            "This is an exploratory engineering result and is not a publication-grade fastVEP validation.",
            "",
        ]
    )
    (config.root / "FASTVEP_PILOT_RESULTS.md").write_text("\n".join(result_lines))
    print(json.dumps(gate, indent=2, sort_keys=True))
    return summary


def status(config: Config) -> None:
    values = {
        "setup": (config.root / "environment.json").exists(),
        "smoke": (config.root / "smoke/annotator_comparison.json").exists(),
        "direct": (config.root / "runs/direct/metrics.json").exists(),
        **{
            f"hit_{hit_rate}": (
                config.root / f"runs/hit-{hit_rate:03d}/metrics.json"
            ).exists()
            for hit_rate in HIT_RATES
        },
        "summary": (config.root / "summary.json").exists(),
    }
    print(json.dumps(values, indent=2, sort_keys=True))


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument(
        "command", choices=("setup", "smoke", "wgs", "collect", "status", "all")
    )
    result.add_argument(
        "--root", type=Path, default=Path("/mnt/data/vcfcache_benchmarks/fastvep_pilot")
    )
    result.add_argument(
        "--repo", type=Path, default=Path("/mnt/data/home/appuser-projects/vcfcache")
    )
    result.add_argument(
        "--input",
        type=Path,
        default=Path(
            "/mnt/data/vcfcache_benchmarks/samples/GRCh38/1000g/EAS/"
            "HG02079.GRCh38.small_variants.vcf.gz"
        ),
    )
    result.add_argument(
        "--fasta-gz",
        type=Path,
        default=Path(
            "/mnt/data/apps/ensembl-vep/115/cachedir/homo_sapiens/115_GRCh38/"
            "Homo_sapiens.GRCh38.dna.toplevel.fa.gz"
        ),
    )
    result.add_argument(
        "--vep-sif", type=Path, default=Path("/mnt/data/apps/ensembl-vep/115/vep.sif")
    )
    result.add_argument("--threads", type=int, default=16)
    return result


def main(argv: Sequence[str] | None = None) -> int:
    args = parser().parse_args(argv)
    config = Config(
        root=args.root.resolve(),
        repo=args.repo.resolve(),
        input_vcf=args.input.resolve(),
        fasta_gz=args.fasta_gz.resolve(),
        vep_sif=args.vep_sif.resolve(),
        threads=args.threads,
    )
    if args.command in ("setup", "all"):
        setup(config)
    if args.command in ("smoke", "all"):
        smoke(config)
    if args.command in ("wgs", "all"):
        wgs(config)
    if args.command in ("collect", "all"):
        collect(config)
    if args.command == "status":
        status(config)
    return 0


if __name__ == "__main__":
    sys.exit(main())
