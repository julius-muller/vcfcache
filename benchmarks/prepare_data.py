#!/usr/bin/env python3
"""Prepare reproducible public VCFs for the VCFcache publication benchmarks."""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import re
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import requests  # type: ignore[import-untyped]

DEFAULT_ROOT = Path("/mnt/data/vcfcache_benchmarks")
SELECTION_SEED = "vcfcache-paper-v1"
SUPERPOPS = ("AFR", "AMR", "EAS", "EUR", "SAS")
AUTOSOMES = tuple(str(i) for i in range(1, 23))
GIAB_EXCLUSIONS = {
    "NA12878",
    "NA24385",
    "NA24149",
    "NA24143",
    "NA24631",
    "NA24694",
    "NA24695",
}

THOUSAND_GENOMES_BASE = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
    "1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV"
)
THOUSAND_GENOMES_MANIFEST_URL = f"{THOUSAND_GENOMES_BASE}/20220804_manifest.txt"
THOUSAND_GENOMES_PANEL_URL = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
    "integrated_call_samples_v3.20130502.ALL.panel"
)


@dataclass(frozen=True)
class SelectedSample:
    """Frozen metadata for a selected 1000 Genomes sample."""

    sample: str
    population: str
    superpopulation: str
    sex: str


@dataclass(frozen=True)
class SourceFile:
    """Download source with an optional upstream MD5."""

    cohort: str
    sample: str
    url: str
    md5: str = ""

    @property
    def filename(self) -> str:
        return self.url.rsplit("/", 1)[-1]


GIAB_SOURCES = (
    SourceFile(
        "giab",
        "HG001",
        "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/"
        "NA12878_HG001/latest/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
        "ec61d76ad725d808e3a9e1c8acec402a",
    ),
    SourceFile(
        "giab",
        "HG001",
        "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/"
        "NA12878_HG001/latest/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi",
        "476166561e6a486aee67d727dbd102ef",
    ),
    SourceFile(
        "giab",
        "HG002",
        "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/"
        "AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/"
        "HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
    ),
    SourceFile(
        "giab",
        "HG002",
        "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/"
        "AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/"
        "HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi",
    ),
    SourceFile(
        "giab",
        "HG003",
        "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/"
        "AshkenazimTrio/HG003_NA24149_father/latest/GRCh38/"
        "HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
    ),
    SourceFile(
        "giab",
        "HG003",
        "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/"
        "AshkenazimTrio/HG003_NA24149_father/latest/GRCh38/"
        "HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi",
    ),
    SourceFile(
        "giab",
        "HG004",
        "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/"
        "AshkenazimTrio/HG004_NA24143_mother/latest/GRCh38/"
        "HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
    ),
    SourceFile(
        "giab",
        "HG004",
        "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/"
        "AshkenazimTrio/HG004_NA24143_mother/latest/GRCh38/"
        "HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi",
    ),
    SourceFile(
        "giab",
        "HG005",
        "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/"
        "ChineseTrio/HG005_NA24631_son/latest/GRCh38/"
        "HG005_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
        "35fe853b7e8d979daa91033ed6863e98",
    ),
    SourceFile(
        "giab",
        "HG005",
        "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/"
        "ChineseTrio/HG005_NA24631_son/latest/GRCh38/"
        "HG005_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi",
        "ffd21b44e0bfe4133b4f34fba88ea57e",
    ),
    SourceFile(
        "giab",
        "HG006",
        "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/"
        "ChineseTrio/HG006_NA24694_father/latest/GRCh38/"
        "HG006_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
        "d02b15294b232f96506dea2dd8d6d82d",
    ),
    SourceFile(
        "giab",
        "HG006",
        "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/"
        "ChineseTrio/HG006_NA24694_father/latest/GRCh38/"
        "HG006_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi",
        "1d8dae44b3d23713f09693e970ced858",
    ),
    SourceFile(
        "giab",
        "HG007",
        "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/"
        "ChineseTrio/HG007_NA24695_mother/latest/GRCh38/"
        "HG007_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
        "336f8cf6806521618779b6c5f4eadeb7",
    ),
    SourceFile(
        "giab",
        "HG007",
        "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/"
        "ChineseTrio/HG007_NA24695_mother/latest/GRCh38/"
        "HG007_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi",
        "6001c661abf3b4c8eec1aa96fbb82fe9",
    ),
)


def fetch_text(url: str, timeout: int = 60) -> str:
    """Fetch a small UTF-8 source manifest."""
    response = requests.get(url, timeout=timeout)
    response.raise_for_status()
    return response.text


def parse_1000g_manifest(text: str) -> dict[str, str]:
    """Parse the upstream whitespace-delimited filename/MD5 pairs."""
    tokens = text.split()
    if len(tokens) % 2:
        raise ValueError("1000 Genomes manifest has an odd number of fields")
    parsed: dict[str, str] = {}
    for filename, checksum in zip(tokens[0::2], tokens[1::2], strict=True):
        if not re.fullmatch(r"[0-9a-fA-F]{32}", checksum):
            raise ValueError(f"Invalid MD5 for {filename}: {checksum}")
        parsed[filename] = checksum.lower()
    return parsed


def parse_population_panel(text: str) -> list[SelectedSample]:
    """Parse the phase-three unrelated sample panel."""
    rows = csv.DictReader(text.splitlines(), delimiter="\t")
    samples = []
    for row in rows:
        samples.append(
            SelectedSample(
                sample=row["sample"],
                population=row["pop"],
                superpopulation=row["super_pop"],
                sex=row["gender"],
            )
        )
    return samples


def deterministic_sample_selection(
    candidates: Sequence[SelectedSample],
    *,
    per_superpopulation: int = 10,
    seed: str = SELECTION_SEED,
    exclusions: set[str] | None = None,
) -> list[SelectedSample]:
    """Select a deterministic sex-balanced cohort."""
    exclusions = GIAB_EXCLUSIONS if exclusions is None else exclusions
    if per_superpopulation % 2:
        raise ValueError("per_superpopulation must be even for sex balance")
    per_sex = per_superpopulation // 2
    selected: list[SelectedSample] = []
    for superpopulation in SUPERPOPS:
        for sex in ("female", "male"):
            eligible = [
                sample
                for sample in candidates
                if sample.superpopulation == superpopulation
                and sample.sex == sex
                and sample.sample not in exclusions
            ]
            eligible.sort(
                key=lambda item: hashlib.sha256(
                    f"{seed}:{item.sample}".encode()
                ).hexdigest()
            )
            if len(eligible) < per_sex:
                raise ValueError(
                    f"Not enough {sex} samples in {superpopulation}: {len(eligible)}"
                )
            selected.extend(eligible[:per_sex])
    return sorted(selected, key=lambda item: (item.superpopulation, item.sample))


def md5sum(path: Path, chunk_size: int = 8 << 20) -> str:
    """Return a streaming MD5 for upstream verification."""
    digest = hashlib.md5()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def sha256sum(path: Path, chunk_size: int = 8 << 20) -> str:
    """Return a streaming SHA-256 for prepared-data provenance."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def run_command(
    args: Sequence[str],
    *,
    cwd: Path | None = None,
    stdout=None,
    stderr=None,
) -> subprocess.CompletedProcess:
    """Run and check an external command."""
    return subprocess.run(
        list(args),
        check=True,
        cwd=cwd,
        stdout=stdout,
        stderr=stderr,
        text=stdout == subprocess.PIPE or stderr == subprocess.PIPE,
    )


def ensure_layout(root: Path) -> None:
    """Create the external benchmark-data directory layout."""
    paths = (
        "sources/1000g/vcf",
        "sources/1000g/index",
        "sources/1000g/metadata",
        "sources/giab",
        "samples/GRCh38/1000g",
        "samples/GRCh38/giab",
        "manifests",
        "qc/bcftools_stats",
        "qc/failures",
        "work/partial_downloads",
        "work/chromosome_shards",
        "work/sort",
        "work/tmp",
        "logs",
    )
    for relative in paths:
        (root / relative).mkdir(parents=True, exist_ok=True)


def preflight(root: Path, minimum_free_gib: int = 100) -> None:
    """Validate tools, storage, and root placement."""
    ensure_layout(root)
    if root.resolve().is_relative_to(
        Path("/").resolve()
    ) and not root.resolve().is_relative_to(Path("/mnt/data").resolve()):
        raise RuntimeError("Benchmark root must be placed under /mnt/data")
    for command in ("bcftools", "bgzip", "tabix", "curl"):
        if shutil.which(command) is None:
            raise RuntimeError(f"Required command not found: {command}")
    free = shutil.disk_usage(root).free
    required = minimum_free_gib * (1 << 30)
    if free < required:
        raise RuntimeError(
            f"Insufficient free space: {free / (1 << 30):.1f} GiB; "
            f"need at least {minimum_free_gib} GiB"
        )
    version = run_command(
        ["bcftools", "--version"], stdout=subprocess.PIPE
    ).stdout.splitlines()[0]
    print(f"Preflight OK: root={root}, free={free / (1 << 30):.1f} GiB, {version}")


def write_selected_samples(root: Path, selected: Sequence[SelectedSample]) -> Path:
    """Write the frozen selected-sample manifest."""
    output = root / "manifests/selected_samples.tsv"
    temporary = output.with_suffix(".tsv.partial")
    with temporary.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("sample", "population", "superpopulation", "sex", "seed"))
        for item in selected:
            writer.writerow(
                (
                    item.sample,
                    item.population,
                    item.superpopulation,
                    item.sex,
                    SELECTION_SEED,
                )
            )
    temporary.replace(output)
    sample_list = root / "manifests/selected_sample_ids.txt"
    sample_list.write_text("".join(f"{item.sample}\n" for item in selected))
    print(f"Selected {len(selected)} samples -> {output}")
    return output


def select_samples(root: Path) -> list[SelectedSample]:
    """Fetch population metadata and freeze the selected cohort."""
    ensure_layout(root)
    text = fetch_text(THOUSAND_GENOMES_PANEL_URL)
    metadata_path = root / "sources/1000g/metadata/phase3_unrelated.panel"
    metadata_path.write_text(text)
    selected = deterministic_sample_selection(parse_population_panel(text))
    write_selected_samples(root, selected)
    return selected


def read_selected_samples(root: Path) -> list[SelectedSample]:
    """Read the frozen selected cohort."""
    path = root / "manifests/selected_samples.tsv"
    if not path.exists():
        return select_samples(root)
    with path.open() as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        return [
            SelectedSample(
                row["sample"], row["population"], row["superpopulation"], row["sex"]
            )
            for row in rows
        ]


def download_file(source: SourceFile, destination: Path) -> Path:
    """Download a source atomically with curl resume and checksum verification."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        if not source.md5 or md5sum(destination) == source.md5:
            print(f"Verified existing: {destination}")
            return destination
        raise RuntimeError(f"Existing file has wrong checksum: {destination}")
    partial = destination.with_name(destination.name + ".partial")
    run_command(
        [
            "curl",
            "--location",
            "--fail",
            "--show-error",
            "--silent",
            "--retry",
            "8",
            "--retry-all-errors",
            "--continue-at",
            "-",
            "--output",
            str(partial),
            source.url,
        ]
    )
    if source.md5:
        observed = md5sum(partial)
        if observed != source.md5:
            raise RuntimeError(
                f"MD5 mismatch for {source.url}: expected {source.md5}, got {observed}"
            )
    partial.replace(destination)
    print(f"Downloaded: {destination}")
    return destination


def download_many(
    jobs: Iterable[tuple[SourceFile, Path]], *, workers: int = 4
) -> list[Path]:
    """Download multiple independent files concurrently."""
    results = []
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(download_file, source, destination): destination
            for source, destination in jobs
        }
        for future in as_completed(futures):
            results.append(future.result())
    return results


def thousand_genomes_sources(root: Path) -> list[tuple[SourceFile, Path]]:
    """Resolve autosomal VCF and TBI sources from the upstream manifest."""
    text = fetch_text(THOUSAND_GENOMES_MANIFEST_URL)
    manifest_path = root / "sources/1000g/metadata/20220804_manifest.txt"
    manifest_path.write_text(text)
    checksums = parse_1000g_manifest(text)
    jobs = []
    for chromosome in AUTOSOMES:
        stem = (
            f"1kGP_high_coverage_Illumina.chr{chromosome}."
            "filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
        )
        for filename in (stem, f"{stem}.tbi"):
            source = SourceFile(
                "1000g",
                f"chr{chromosome}",
                f"{THOUSAND_GENOMES_BASE}/{filename}",
                checksums[filename],
            )
            # htslib discovers local indexes beside the data file by appending
            # .tbi, so retain each VCF/index pair in the same directory.
            jobs.append((source, root / f"sources/1000g/vcf/{filename}"))
    return jobs


def download_1000g(root: Path, workers: int = 4) -> list[Path]:
    """Download and verify all autosomal 1000 Genomes sources."""
    ensure_layout(root)
    return download_many(thousand_genomes_sources(root), workers=workers)


def download_giab(root: Path, workers: int = 4) -> list[Path]:
    """Download and verify the seven GIAB v4.2.1 sources."""
    ensure_layout(root)
    jobs = [
        (source, root / "sources/giab" / source.sample / source.filename)
        for source in GIAB_SOURCES
    ]
    return download_many(jobs, workers=workers)


def _chromosome_source(root: Path, chromosome: str) -> Path:
    name = (
        f"1kGP_high_coverage_Illumina.chr{chromosome}."
        "filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    )
    return root / "sources/1000g/vcf" / name


def prepare_1000g_chromosome(root: Path, chromosome: str) -> Path:
    """Split one chromosome into compact single-sample VCF shards."""
    ensure_layout(root)
    source = _chromosome_source(root, chromosome)
    if not source.exists():
        raise FileNotFoundError(source)
    output_dir = root / "work/chromosome_shards" / f"chr{chromosome}"
    output_dir.mkdir(parents=True, exist_ok=True)
    marker = output_dir / ".complete"
    if marker.exists():
        return output_dir
    sample_list = root / "manifests/selected_sample_ids.txt"
    log = root / "logs" / f"split_chr{chromosome}.log"
    view = subprocess.Popen(
        [
            "bcftools",
            "view",
            "--samples-file",
            str(sample_list),
            "--no-update",
            "--types",
            "snps,indels",
            "--exclude",
            'ALT~"[<>]" || INFO/SVTYPE!="."',
            "--output-type",
            "u",
            str(source),
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    annotate = subprocess.Popen(
        [
            "bcftools",
            "annotate",
            "--remove",
            "^INFO/AF,INFO/AC,INFO/AN,FORMAT",
            "--output-type",
            "u",
        ],
        stdin=view.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert view.stdout is not None
    view.stdout.close()
    split = subprocess.run(
        [
            "bcftools",
            "+split",
            "--output-type",
            "z",
            "--output",
            str(output_dir),
            "--samples-file",
            str(sample_list),
            "--include",
            'GT="alt"',
            "--keep-tags",
            "INFO,FMT/GT",
            "--write-index=tbi",
        ],
        stdin=annotate.stdout,
        capture_output=True,
    )
    assert annotate.stdout is not None
    annotate.stdout.close()
    annotate_stderr = annotate.communicate()[1]
    view_stderr = view.communicate()[1]
    log.write_bytes(
        b"[view]\n"
        + view_stderr
        + b"\n[annotate]\n"
        + annotate_stderr
        + b"\n[split]\n"
        + split.stderr
    )
    failures = {
        "view": view.returncode,
        "annotate": annotate.returncode,
        "split": split.returncode,
    }
    failures = {name: code for name, code in failures.items() if code}
    if failures:
        raise RuntimeError(
            f"Chromosome {chromosome} split failed: {failures}; see {log}"
        )
    selected = read_selected_samples(root)
    missing = [
        item.sample
        for item in selected
        if not (output_dir / f"{item.sample}.vcf.gz").exists()
    ]
    if missing:
        raise RuntimeError(f"Missing chr{chromosome} shards: {missing}")
    marker.write_text("complete\n")
    print(f"Prepared chromosome {chromosome}")
    return output_dir


def finalize_1000g_sample(root: Path, sample: SelectedSample) -> Path:
    """Concatenate, sort, compress, and index one 1000 Genomes sample."""
    destination_dir = root / f"samples/GRCh38/1000g/{sample.superpopulation}"
    destination_dir.mkdir(parents=True, exist_ok=True)
    destination = destination_dir / f"{sample.sample}.GRCh38.small_variants.vcf.gz"
    if destination.exists() and Path(f"{destination}.tbi").exists():
        return destination
    shards = [
        root / f"work/chromosome_shards/chr{chromosome}/{sample.sample}.vcf.gz"
        for chromosome in AUTOSOMES
    ]
    missing = [str(path) for path in shards if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing shards for {sample.sample}: {missing[:3]}")
    partial = destination.with_name(destination.name + ".partial")
    sort_dir = root / "work/sort" / sample.sample
    sort_dir.mkdir(parents=True, exist_ok=True)
    concat = subprocess.Popen(
        ["bcftools", "concat", "--allow-overlaps", "--output-type", "u"]
        + [str(path) for path in shards],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    sort = subprocess.run(
        [
            "bcftools",
            "sort",
            "--temp-dir",
            str(sort_dir),
            "--output-type",
            "z",
            "--output",
            str(partial),
        ],
        stdin=concat.stdout,
        capture_output=True,
    )
    assert concat.stdout is not None
    concat.stdout.close()
    concat_stderr = concat.communicate()[1]
    if concat.returncode or sort.returncode:
        raise RuntimeError(
            f"Finalize failed for {sample.sample}: concat={concat.returncode}, "
            f"sort={sort.returncode}\n{concat_stderr.decode()}\n{sort.stderr.decode()}"
        )
    partial.replace(destination)
    run_command(["tabix", "--preset", "vcf", str(destination)])
    print(f"Finalized: {destination}")
    return destination


def prepare_1000g(
    root: Path,
    *,
    chromosome_workers: int = 4,
    sample_workers: int = 8,
) -> list[Path]:
    """Prepare all 50 single-sample 1000 Genomes VCFs."""
    selected = read_selected_samples(root)
    with ThreadPoolExecutor(max_workers=chromosome_workers) as executor:
        futures = [
            executor.submit(prepare_1000g_chromosome, root, chromosome)
            for chromosome in AUTOSOMES
        ]
        for future in as_completed(futures):
            future.result()
    results = []
    with ThreadPoolExecutor(max_workers=sample_workers) as executor:
        sample_futures = {
            executor.submit(finalize_1000g_sample, root, sample): sample
            for sample in selected
        }
        for future in as_completed(sample_futures):
            results.append(future.result())
    return results


def prepare_one_giab(root: Path, sample: str) -> Path:
    """Sort and index one official GIAB VCF."""
    sources = [
        source
        for source in GIAB_SOURCES
        if source.sample == sample and source.filename.endswith(".vcf.gz")
    ]
    if len(sources) != 1:
        raise RuntimeError(f"Expected one GIAB VCF for {sample}")
    source = root / "sources/giab" / sample / sources[0].filename
    destination = root / f"samples/GRCh38/giab/{sample}.GRCh38.v4.2.1.vcf.gz"
    if destination.exists() and Path(f"{destination}.tbi").exists():
        return destination
    partial = destination.with_name(destination.name + ".partial")
    sort_dir = root / "work/sort" / sample
    sort_dir.mkdir(parents=True, exist_ok=True)
    run_command(
        [
            "bcftools",
            "sort",
            "--temp-dir",
            str(sort_dir),
            "--output-type",
            "z",
            "--output",
            str(partial),
            str(source),
        ]
    )
    partial.replace(destination)
    run_command(["tabix", "--preset", "vcf", str(destination)])
    print(f"Prepared GIAB: {destination}")
    return destination


def prepare_giab(root: Path, workers: int = 4) -> list[Path]:
    """Prepare all seven official GIAB files."""
    samples = [f"HG{i:03d}" for i in range(1, 8)]
    with ThreadPoolExecutor(max_workers=workers) as executor:
        return list(
            executor.map(lambda sample: prepare_one_giab(root, sample), samples)
        )


def _stats_summary(vcf: Path) -> dict[str, int]:
    result = run_command(["bcftools", "stats", str(vcf)], stdout=subprocess.PIPE)
    summary: dict[str, int] = {}
    mapping = {
        "number of records": "records",
        "number of SNPs": "snps",
        "number of indels": "indels",
        "number of MNPs": "mnps",
        "number of others": "others",
        "number of multiallelic sites": "multiallelic",
    }
    for line in result.stdout.splitlines():
        if not line.startswith("SN\t"):
            continue
        fields = line.split("\t")
        label = fields[2].rstrip(":")
        if label in mapping:
            summary[mapping[label]] = int(fields[3])
    return summary


def _record_contigs(vcf: Path) -> set[str]:
    result = run_command(
        ["bcftools", "index", "--stats", str(vcf)], stdout=subprocess.PIPE
    )
    return {
        line.split("\t", maxsplit=1)[0] for line in result.stdout.splitlines() if line
    }


def _header_tags(vcf: Path, kind: str) -> set[str]:
    header = run_command(
        ["bcftools", "view", "--header-only", str(vcf)], stdout=subprocess.PIPE
    ).stdout
    pattern = re.compile(rf"^##{kind}=<ID=([^,>]+)", re.MULTILINE)
    return set(pattern.findall(header))


def _has_non_alt_genotype(vcf: Path) -> bool:
    process = subprocess.Popen(
        ["bcftools", "view", "--exclude", 'GT="alt"', "--no-header", str(vcf)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert process.stdout is not None
    first = process.stdout.readline()
    process.terminate()
    process.communicate()
    return bool(first)


def validate_prepared_vcf(vcf: Path, *, cohort: str) -> dict[str, object]:
    """Comprehensively validate one prepared single-sample VCF."""
    errors: list[str] = []
    try:
        run_command(["bgzip", "--test", str(vcf)])
    except subprocess.CalledProcessError:
        errors.append("bgzip_test")
    index = Path(f"{vcf}.tbi")
    if not index.exists():
        errors.append("missing_tbi")
    else:
        try:
            run_command(
                ["bcftools", "index", "--nrecords", str(vcf)],
                stdout=subprocess.PIPE,
            )
        except subprocess.CalledProcessError:
            errors.append("invalid_index")
    samples = run_command(
        ["bcftools", "query", "--list-samples", str(vcf)], stdout=subprocess.PIPE
    ).stdout.splitlines()
    if len(samples) != 1:
        errors.append(f"sample_count={len(samples)}")
    stats = _stats_summary(vcf)
    contigs = _record_contigs(vcf)
    allowed_contigs = set(AUTOSOMES) | {f"chr{chromosome}" for chromosome in AUTOSOMES}
    if not contigs <= allowed_contigs:
        errors.append(f"unexpected_contigs={sorted(contigs - allowed_contigs)}")
    if cohort == "1000g":
        if stats.get("mnps", 0) or stats.get("others", 0):
            errors.append("non_small_variant_records")
        info_tags = _header_tags(vcf, "INFO")
        if info_tags != {"AF", "AC", "AN"}:
            errors.append(f"info_tags={sorted(info_tags)}")
        format_tags = _header_tags(vcf, "FORMAT")
        if format_tags != {"GT"}:
            errors.append(f"format_tags={sorted(format_tags)}")
        if _has_non_alt_genotype(vcf):
            errors.append("homref_or_missing_record")
    return {
        "status": "PASS" if not errors else "FAIL",
        "errors": ",".join(errors),
        "sample": samples[0] if len(samples) == 1 else "",
        "records": stats.get("records", 0),
        "snps": stats.get("snps", 0),
        "indels": stats.get("indels", 0),
        "multiallelic": stats.get("multiallelic", 0),
        "contigs": ",".join(sorted(contigs)),
        "bytes": vcf.stat().st_size,
        "sha256": sha256sum(vcf),
    }


def build_flat_sample_view(root: Path) -> Path:
    """Expose all final VCF/index pairs in one non-duplicating directory."""
    vcfs = sorted((root / "samples/GRCh38/1000g").glob("*/*.vcf.gz"))
    vcfs += sorted((root / "samples/GRCh38/giab").glob("*.vcf.gz"))
    if len(vcfs) != 57:
        raise RuntimeError(f"Expected 57 prepared VCFs, found {len(vcfs)}")
    output = root / "samples/GRCh38/all"
    output.mkdir(parents=True, exist_ok=True)
    for vcf in vcfs:
        for source in (vcf, Path(f"{vcf}.tbi")):
            if not source.exists():
                raise FileNotFoundError(source)
            destination = output / source.name
            relative_source = Path(os.path.relpath(source, output))
            if destination.is_symlink():
                if destination.readlink() != relative_source:
                    raise RuntimeError(f"Unexpected symlink target: {destination}")
                continue
            if destination.exists():
                raise RuntimeError(f"Refusing to replace existing path: {destination}")
            destination.symlink_to(relative_source)
    print(f"Flat sample view: {output}")
    return output


def qc_all(root: Path, *, workers: int = 8) -> Path:
    """Validate all final VCFs and write a machine-readable QC table."""
    selected = {item.sample: item for item in read_selected_samples(root)}
    vcfs = sorted((root / "samples/GRCh38/1000g").glob("*/*.vcf.gz"))
    vcfs += sorted((root / "samples/GRCh38/giab").glob("*.vcf.gz"))
    if len(vcfs) != 57:
        raise RuntimeError(f"Expected 57 prepared VCFs, found {len(vcfs)}")
    output = root / "qc/sample_qc.tsv"
    temporary = output.with_suffix(".tsv.partial")
    fields = (
        "cohort",
        "sample",
        "population",
        "superpopulation",
        "sex",
        "path",
        "records",
        "snps",
        "indels",
        "multiallelic",
        "contigs",
        "bytes",
        "sha256",
        "status",
        "errors",
    )
    work = [(vcf, "giab" if "/giab/" in str(vcf) else "1000g") for vcf in vcfs]
    with ThreadPoolExecutor(max_workers=workers) as executor:
        validated = list(
            executor.map(
                lambda item: validate_prepared_vcf(item[0], cohort=item[1]),
                work,
            )
        )
    failures = []
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for (vcf, cohort), result in zip(work, validated, strict=True):
            sample_id = str(result["sample"])
            metadata = selected.get(sample_id)
            row = {
                "cohort": cohort,
                "sample": sample_id,
                "population": metadata.population if metadata else "",
                "superpopulation": metadata.superpopulation if metadata else "",
                "sex": metadata.sex if metadata else "",
                "path": str(vcf),
                **result,
            }
            writer.writerow({field: row.get(field, "") for field in fields})
            if result["status"] != "PASS":
                failures.append((vcf, result["errors"]))
            print(f"QC {result['status']}: {vcf.name} ({result['records']:,} records)")
    temporary.replace(output)
    if failures:
        raise RuntimeError(f"{len(failures)} VCFs failed QC; see {output}")
    build_flat_sample_view(root)
    print(f"All 57 VCFs passed QC -> {output}")
    return output


def smoke_test(root: Path) -> None:
    """Exercise official metadata and a small remote chr22 region."""
    preflight(root)
    selected = read_selected_samples(root)
    test_samples = [
        next(item for item in selected if item.superpopulation == "AFR"),
        next(item for item in selected if item.superpopulation == "EUR"),
    ]
    source_url = (
        f"{THOUSAND_GENOMES_BASE}/"
        "1kGP_high_coverage_Illumina.chr22.filtered."
        "SNV_INDEL_SV_phased_panel.vcf.gz"
    )
    smoke_dir = root / "work/tmp/network_smoke"
    smoke_dir.mkdir(parents=True, exist_ok=True)
    output = smoke_dir / "smoke.vcf.gz"
    run_command(
        [
            "bcftools",
            "view",
            "--regions",
            "chr22:16050000-17000000",
            "--samples",
            ",".join(item.sample for item in test_samples),
            "--types",
            "snps,indels",
            "--output-type",
            "z",
            "--output",
            str(output),
            source_url,
        ],
        cwd=smoke_dir,
    )
    run_command(["tabix", "--force", "--preset", "vcf", str(output)])
    observed = run_command(
        ["bcftools", "query", "--list-samples", str(output)],
        stdout=subprocess.PIPE,
    ).stdout.splitlines()
    expected = [item.sample for item in test_samples]
    if observed != expected:
        raise RuntimeError(
            f"Smoke sample mismatch: expected {expected}, got {observed}"
        )
    print(f"Network smoke test passed: {output}")


def write_source_manifest(root: Path) -> Path:
    """Record the immutable source URLs and checksums."""
    output = root / "manifests/source_files.tsv"
    rows: list[SourceFile] = [source for source, _ in thousand_genomes_sources(root)]
    rows.extend(GIAB_SOURCES)
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("cohort", "sample_or_chromosome", "url", "upstream_md5"))
        for source in rows:
            writer.writerow((source.cohort, source.sample, source.url, source.md5))
    return output


def run_all(args: argparse.Namespace) -> None:
    """Run all resumable data preparation stages."""
    preflight(args.root)
    read_selected_samples(args.root)
    write_source_manifest(args.root)
    download_1000g(args.root, workers=args.download_workers)
    download_giab(args.root, workers=args.download_workers)
    prepare_1000g(
        args.root,
        chromosome_workers=args.chromosome_workers,
        sample_workers=args.sample_workers,
    )
    prepare_giab(args.root, workers=min(args.sample_workers, 7))
    qc_all(args.root, workers=args.sample_workers)


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root", type=Path, default=DEFAULT_ROOT, help="External benchmark-data root"
    )
    parser.add_argument("--download-workers", type=int, default=4)
    parser.add_argument("--chromosome-workers", type=int, default=4)
    parser.add_argument("--sample-workers", type=int, default=8)
    parser.add_argument(
        "command",
        choices=(
            "preflight",
            "select",
            "smoke",
            "download-1000g",
            "download-giab",
            "prepare-1000g",
            "prepare-giab",
            "qc",
            "all",
        ),
    )
    return parser


def main() -> None:
    """Command-line entry point."""
    args = build_parser().parse_args()
    args.root = args.root.expanduser().resolve()
    os.environ["TMPDIR"] = str(args.root / "work/tmp")
    commands = {
        "preflight": lambda: preflight(args.root),
        "select": lambda: select_samples(args.root),
        "smoke": lambda: smoke_test(args.root),
        "download-1000g": lambda: download_1000g(
            args.root, workers=args.download_workers
        ),
        "download-giab": lambda: download_giab(
            args.root, workers=args.download_workers
        ),
        "prepare-1000g": lambda: prepare_1000g(
            args.root,
            chromosome_workers=args.chromosome_workers,
            sample_workers=args.sample_workers,
        ),
        "prepare-giab": lambda: prepare_giab(
            args.root, workers=min(args.sample_workers, 7)
        ),
        "qc": lambda: qc_all(args.root, workers=args.sample_workers),
        "all": lambda: run_all(args),
    }
    commands[args.command]()


if __name__ == "__main__":
    main()
