#!/usr/bin/env python3
"""Prepare one Slurm worker for the matched external fastVEP campaign."""

from __future__ import annotations

import argparse
import gzip
import os
import shutil
import subprocess
from pathlib import Path
from types import SimpleNamespace

import yaml

from benchmarks.prepare_external_fastvep import prepare as prepare_caches
from benchmarks.run_cohort import sha256sum
from benchmarks.run_pilot import write_json_atomic

FASTVEP_VERSION = "0.3.0"
GFF_URLS = {
    "GRCh37": (
        "https://ftp.ensembl.org/pub/grch37/release-115/gff3/homo_sapiens/"
        "Homo_sapiens.GRCh37.87.gff3.gz"
    ),
    "GRCh38": (
        "https://ftp.ensembl.org/pub/release-115/gff3/homo_sapiens/"
        "Homo_sapiens.GRCh38.115.gff3.gz"
    ),
}
SOURCE_FASTA = {
    "GRCh37": Path(
        "/mnt/data/apps/ensembl-vep/115/cachedir/homo_sapiens/115_GRCh37/"
        "Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz"
    ),
    "GRCh38": Path(
        "/mnt/data/apps/ensembl-vep/115/cachedir/homo_sapiens/115_GRCh38/"
        "Homo_sapiens.GRCh38.dna.toplevel.fa.gz"
    ),
}


def run(command: list[str | Path]) -> None:
    """Run one setup command."""
    print("$", " ".join(map(str, command)), flush=True)
    subprocess.run(list(map(str, command)), check=True)


def decompress(source: Path, destination: Path) -> None:
    """Decompress once with an atomic destination."""
    if destination.exists():
        return
    partial = destination.with_suffix(destination.suffix + ".partial")
    with gzip.open(source, "rb") as reader, partial.open("wb") as writer:
        shutil.copyfileobj(reader, writer, length=8 << 20)
    partial.replace(destination)


def download(source: str, destination: Path) -> None:
    """Download once without accepting an interrupted file as complete."""
    if destination.exists():
        return
    partial = destination.with_suffix(destination.suffix + ".partial")
    run(["curl", "-L", "--fail", "--retry", "3", "-o", partial, source])
    partial.replace(destination)


def write_fasta_index(fasta: Path, destination: Path) -> None:
    """Write a samtools-compatible FAI for an uncompressed FASTA."""
    if destination.exists():
        return
    records: list[tuple[str, int, int, int, int]] = []
    name: str | None = None
    length = 0
    sequence_offset = 0
    line_bases = 0
    line_width = 0

    def finish_record() -> None:
        if name is None:
            return
        if length < 1 or line_bases < 1 or line_width < line_bases:
            raise RuntimeError(f"Invalid FASTA sequence {name!r} in {fasta}")
        records.append((name, length, sequence_offset, line_bases, line_width))

    with fasta.open("rb") as handle:
        while line := handle.readline():
            if line.startswith(b">"):
                finish_record()
                name = line[1:].split(maxsplit=1)[0].decode()
                if not name:
                    raise RuntimeError(f"Empty FASTA sequence name in {fasta}")
                length = 0
                sequence_offset = handle.tell()
                line_bases = 0
                line_width = 0
                continue
            if name is None:
                raise RuntimeError(f"FASTA sequence occurs before a header in {fasta}")
            bases = len(line.rstrip(b"\r\n"))
            if not bases:
                continue
            if line_bases == 0:
                line_bases = bases
                line_width = len(line)
            length += bases
    finish_record()
    if not records:
        raise RuntimeError(f"No FASTA sequences found in {fasta}")

    partial = destination.with_suffix(destination.suffix + ".partial")
    partial.write_text(
        "".join(
            f"{record_name}\t{record_length}\t{offset}\t{bases}\t{width}\n"
            for record_name, record_length, offset, bases, width in records
        )
    )
    partial.replace(destination)


def write_configs(
    config_root: Path,
    assembly: str,
    fastvep: Path,
    transcript_cache: Path,
    fasta: Path,
    threads: int,
) -> tuple[Path, Path]:
    """Write the immutable command recipe for one assembly."""
    config_root.mkdir(parents=True, exist_ok=True)
    recipe = config_root / f"annotation_{assembly}.yaml"
    params = config_root / f"params_{assembly}.yaml"
    recipe.write_text(f"""genome_build: "{assembly}"
annotation_cmd: |
  ${{params.bcftools_cmd}} view \\${{INPUT_BCF}} -Ov | \\
  env RAYON_NUM_THREADS=${{params.rayon_threads}} ${{params.annotation_tool_cmd}} annotate \\
    --input - \\
    --output - \\
    --output-format vcf \\
    --transcript-cache ${{params.transcript_cache}} \\
    --fasta ${{params.fasta}} \\
    --hgvs --no-progress \\
    2> \\${{AUXILIARY_DIR}}/fastvep_stderr.txt | \\
  ${{params.bcftools_cmd}} view -Ob -o \\${{OUTPUT_BCF}} -W --threads ${{params.threads}}

must_contain_info_tag: CSQ
required_tool_version: "{FASTVEP_VERSION}"
optional_checks: {{}}
""")
    params.write_text(
        yaml.safe_dump(
            {
                "genome_build": assembly,
                "annotation_tool_cmd": str(fastvep),
                "tool_version_command": f"{fastvep} --version",
                "bcftools_cmd": "bcftools",
                "temp_dir": "/mnt/data/tmp/fastvep-publication",
                "threads": threads,
                "rayon_threads": threads,
                "transcript_cache": str(transcript_cache),
                "fasta": str(fasta),
                "optional_checks": {},
            },
            sort_keys=False,
        )
    )
    return recipe, params


def prepare(args: argparse.Namespace) -> None:
    """Install assets, build transcript caches, and rebuild all blueprints."""
    args.root.mkdir(parents=True, exist_ok=True)
    binary = args.root / "bin/fastvep"
    binary.parent.mkdir(parents=True, exist_ok=True)
    if not binary.exists():
        shutil.copy2(args.fastvep_binary, binary)
        binary.chmod(0o755)
    version = subprocess.run(
        [binary, "--version"], check=True, capture_output=True, text=True
    ).stdout.strip()
    if FASTVEP_VERSION not in version:
        raise RuntimeError(f"Unexpected fastVEP version: {version}")

    configs: dict[str, tuple[Path, Path]] = {}
    assets: dict[str, object] = {}
    for assembly in ("GRCh37", "GRCh38"):
        data = args.root / "data" / assembly
        data.mkdir(parents=True, exist_ok=True)
        gff_gz = data / Path(GFF_URLS[assembly]).name
        download(GFF_URLS[assembly], gff_gz)
        gff = gff_gz.with_suffix("")
        decompress(gff_gz, gff)
        source_fasta = SOURCE_FASTA[assembly]
        if not source_fasta.exists():
            raise FileNotFoundError(source_fasta)
        fasta = data / source_fasta.name.removesuffix(".gz")
        decompress(source_fasta, fasta)
        write_fasta_index(fasta, Path(f"{fasta}.fai"))
        transcript_cache = data / f"Homo_sapiens.{assembly}.fastvep.cache"
        if not transcript_cache.exists():
            run(
                [
                    binary,
                    "cache",
                    "--gff3",
                    gff,
                    "--fasta",
                    fasta,
                    "--output",
                    transcript_cache,
                    "--no-progress",
                ]
            )
        configs[assembly] = write_configs(
            args.root / "config",
            assembly,
            binary,
            transcript_cache,
            fasta,
            args.threads,
        )
        assets[assembly] = {
            "gff": str(gff),
            "gff_sha256": sha256sum(gff),
            "fasta": str(fasta),
            "fasta_sha256": sha256sum(fasta),
            "transcript_cache": str(transcript_cache),
            "transcript_cache_sha256": sha256sum(transcript_cache),
        }

    cache_root = args.cache_root
    prepare_caches(
        SimpleNamespace(
            vep_strategies=args.vep_strategies,
            output_root=cache_root,
            published_root=cache_root,
            recipe_grch37=configs["GRCh37"][0],
            params_grch37=configs["GRCh37"][1],
            recipe_grch38=configs["GRCh38"][0],
            params_grch38=configs["GRCh38"][1],
            vcfcache=args.vcfcache,
            cache_name="cache-fastvep-0.3.0-publication",
        )
    )
    write_json_atomic(
        args.root / "node_ready.json",
        {
            "hostname": os.uname().nodename,
            "fastvep_version": version,
            "fastvep_sha256": sha256sum(binary),
            "assets": assets,
            "strategy_manifest": str(cache_root / "fastvep_strategies.json"),
            "complete": True,
        },
    )


def main() -> None:
    """Run node preparation."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--cache-root", type=Path, required=True)
    parser.add_argument("--fastvep-binary", type=Path, required=True)
    parser.add_argument("--vep-strategies", type=Path, required=True)
    parser.add_argument("--vcfcache", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=8)
    prepare(parser.parse_args())


if __name__ == "__main__":
    main()
