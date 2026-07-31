#!/usr/bin/env python3
"""Prepare HPRC R2, matched WES, and matched small-panel benchmark VCFs."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import re
import shutil
import subprocess
import tempfile
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence, TextIO

from benchmarks.prepare_data import (
    GIAB_EXCLUSIONS,
    THOUSAND_GENOMES_BASE,
    THOUSAND_GENOMES_PANEL_URL,
    SourceFile,
    download_file,
    fetch_text,
    parse_population_panel,
    run_command,
    sha256sum,
)

DEFAULT_ROOT = Path("/mnt/data/vcfcache_benchmarks")
REFERENCE_FASTA = Path("/mnt/data/resources/reference/ucsc/hg38.fa.gz")
ACMG_GENE_FILE = Path(__file__).parent / "config/acmg_sf_v3.3_genes.txt"
HPRC_SELECTION_SEED = "vcfcache-paper-hprc-r2-v1"
HPRC_QUOTAS = {"AFR": 5, "AMR": 3, "EAS": 4, "EUR": 4, "SAS": 4}
AUTOSOMES = tuple(f"chr{number}" for number in range(1, 23))
ASSAY_CONTIGS = (*AUTOSOMES, "chrX")

HPRC_R2_URL = (
    "https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/"
    "scratch/2025_02_28_minigraph_cactus/hprc-v2.0-mc-grch38/"
    "hprc-v2.0-mc-grch38.wave.vcf.gz"
)
TWIST_HG38_BED_URL = (
    "https://www.twistbioscience.com/content/dam/twistbioscience/resources/"
    "2024-09/Twist_Exome_Core_Covered_Targets_hg38.bed"
)
ENSEMBL_115_GTF_URL = (
    "https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/"
    "Homo_sapiens.GRCh38.115.gtf.gz"
)
THOUSAND_GENOMES_X_NAME = (
    "1kGP_high_coverage_Illumina.chrX." "filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz"
)
THOUSAND_GENOMES_X_MD5 = "d2c80aa7b3bcb8f895f98fd5779fb448"
THOUSAND_GENOMES_X_TBI_MD5 = "6d4abec5fd7119e6f3152a640d758b13"


@dataclass(frozen=True)
class AssaySample:
    """Population metadata for a selected benchmark sample."""

    sample: str
    population: str
    superpopulation: str
    sex: str


def ensure_assay_layout(root: Path) -> None:
    """Create external directories without changing existing benchmark cohorts."""
    for relative in (
        "sources/hprc_r2",
        "sources/twist",
        "sources/ensembl",
        "sources/1000g/vcf",
        "samples/GRCh38/hprc_r2",
        "samples/GRCh38/wes_twist_core",
        "samples/GRCh38/panel_acmg_sf_v3.3",
        "regions/GRCh38",
        "manifests",
        "qc",
        "logs",
        "work/hprc_r2_shards",
        "work/chromosome_shards/chrX",
        "work/assays",
    ):
        (root / relative).mkdir(parents=True, exist_ok=True)


def preflight(root: Path, minimum_free_gib: int = 15) -> None:
    """Check storage and commands needed by the assay preparation."""
    ensure_assay_layout(root)
    if not root.resolve().is_relative_to(Path("/mnt/data").resolve()):
        raise RuntimeError("Benchmark root must be below /mnt/data")
    for command in ("bcftools", "bgzip", "tabix", "curl"):
        if shutil.which(command) is None:
            raise RuntimeError(f"Required command not found: {command}")
    free = shutil.disk_usage(root).free
    if free < minimum_free_gib * (1 << 30):
        raise RuntimeError(
            f"Only {free / (1 << 30):.1f} GiB free; "
            f"need at least {minimum_free_gib} GiB"
        )
    if not REFERENCE_FASTA.exists() or not Path(f"{REFERENCE_FASTA}.fai").exists():
        raise FileNotFoundError(
            f"Indexed GRCh38 reference is missing: {REFERENCE_FASTA}"
        )
    version = run_command(
        ["bcftools", "--version"], stdout=subprocess.PIPE
    ).stdout.splitlines()[0]
    print(
        f"Assay preflight OK: root={root}, free={free / (1 << 30):.1f} GiB, "
        f"{version}"
    )


def load_acmg_genes(path: Path = ACMG_GENE_FILE) -> list[str]:
    """Load and validate the frozen ACMG SF v3.3 symbols."""
    genes = [
        line.strip()
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ]
    invalid = [gene for gene in genes if not re.fullmatch(r"[A-Z0-9-]+", gene)]
    if invalid:
        raise ValueError(f"Invalid ACMG gene symbols: {invalid}")
    if len(genes) != 84 or len(set(genes)) != 84:
        raise ValueError(
            f"Expected 84 unique ACMG SF v3.3 genes, got {len(genes)} rows and "
            f"{len(set(genes))} unique symbols"
        )
    return genes


def _rank(sample: AssaySample, seed: str) -> str:
    return hashlib.sha256(f"{seed}:{sample.sample}".encode()).hexdigest()


def deterministic_hprc_selection(
    candidates: Sequence[AssaySample],
    *,
    exclusions: set[str],
    seed: str = HPRC_SELECTION_SEED,
) -> list[AssaySample]:
    """Select 20 HPRC samples across all superpopulations and both sexes."""
    selected: list[AssaySample] = []
    for superpopulation, quota in HPRC_QUOTAS.items():
        eligible = [
            sample
            for sample in candidates
            if sample.superpopulation == superpopulation
            and sample.sample not in exclusions
        ]
        by_sex = {
            sex: sorted(
                (sample for sample in eligible if sample.sex == sex),
                key=lambda sample: _rank(sample, seed),
            )
            for sex in ("female", "male")
        }
        per_sex = quota // 2
        if len(eligible) < quota or any(
            len(by_sex[sex]) < per_sex for sex in ("female", "male")
        ):
            raise ValueError(
                f"Insufficient HPRC {superpopulation} candidates for quota {quota}: "
                f"female={len(by_sex['female'])}, male={len(by_sex['male'])}"
            )
        chosen = by_sex["female"][:per_sex] + by_sex["male"][:per_sex]
        remaining = sorted(
            (sample for sample in eligible if sample not in chosen),
            key=lambda sample: _rank(sample, seed),
        )
        chosen.extend(remaining[: quota - len(chosen)])
        selected.extend(chosen)
    return sorted(selected, key=lambda sample: (sample.superpopulation, sample.sample))


def _hprc_source(root: Path) -> Path:
    return root / "sources/hprc_r2/hprc-v2.0-mc-grch38.wave.vcf.gz"


def _twist_source(root: Path) -> Path:
    return root / "sources/twist/Twist_Exome_Core_Covered_Targets_hg38.bed"


def _ensembl_source(root: Path) -> Path:
    return root / "sources/ensembl/Homo_sapiens.GRCh38.115.gtf.gz"


def _chromosome_x_source(root: Path) -> Path:
    return root / "sources/1000g/vcf" / THOUSAND_GENOMES_X_NAME


def assay_source_jobs(root: Path) -> list[tuple[SourceFile, Path]]:
    """Return immutable source downloads and any upstream checksums."""
    return [
        (SourceFile("hprc_r2", "joint", HPRC_R2_URL), _hprc_source(root)),
        (
            SourceFile("hprc_r2", "joint", f"{HPRC_R2_URL}.tbi"),
            Path(f"{_hprc_source(root)}.tbi"),
        ),
        (SourceFile("twist", "hg38", TWIST_HG38_BED_URL), _twist_source(root)),
        (
            SourceFile("ensembl", "GRCh38", ENSEMBL_115_GTF_URL),
            _ensembl_source(root),
        ),
        (
            SourceFile(
                "1000g_chrX",
                "chrX",
                f"{THOUSAND_GENOMES_BASE}/{THOUSAND_GENOMES_X_NAME}",
                THOUSAND_GENOMES_X_MD5,
            ),
            _chromosome_x_source(root),
        ),
        (
            SourceFile(
                "1000g_chrX",
                "chrX",
                f"{THOUSAND_GENOMES_BASE}/{THOUSAND_GENOMES_X_NAME}.tbi",
                THOUSAND_GENOMES_X_TBI_MD5,
            ),
            Path(f"{_chromosome_x_source(root)}.tbi"),
        ),
    ]


def write_assay_source_manifest(
    root: Path, jobs: Iterable[tuple[SourceFile, Path]]
) -> Path:
    """Freeze URLs, sizes, and checksums after successful downloads."""
    output = root / "manifests/assay_sources.tsv"
    partial = output.with_suffix(".tsv.partial")
    with partial.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            ("cohort", "sample", "url", "path", "bytes", "upstream_md5", "sha256")
        )
        for source, path in jobs:
            writer.writerow(
                (
                    source.cohort,
                    source.sample,
                    source.url,
                    path,
                    path.stat().st_size,
                    source.md5,
                    sha256sum(path),
                )
            )
    partial.replace(output)
    return output


def download_sources(root: Path, workers: int = 3) -> list[Path]:
    """Download all sources concurrently with resume and atomic rename."""
    ensure_assay_layout(root)
    jobs = assay_source_jobs(root)
    results: list[Path] = []
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(download_file, source, path): path for source, path in jobs
        }
        for future in as_completed(futures):
            results.append(future.result())
    write_assay_source_manifest(root, jobs)
    return results


def _existing_sample_exclusions(root: Path) -> set[str]:
    exclusions = set(GIAB_EXCLUSIONS)
    exclusions.update({f"HG{number:03d}" for number in range(1, 8)})
    selected = root / "manifests/selected_samples.tsv"
    if selected.exists():
        with selected.open() as handle:
            exclusions.update(
                row["sample"] for row in csv.DictReader(handle, delimiter="\t")
            )
    return exclusions


def _hprc_sample_ids(root: Path) -> set[str]:
    source = _hprc_source(root)
    target = str(source if source.exists() else HPRC_R2_URL)
    result = run_command(
        ["bcftools", "query", "--list-samples", target], stdout=subprocess.PIPE
    )
    return set(result.stdout.splitlines())


def select_hprc_samples(root: Path) -> list[AssaySample]:
    """Intersect HPRC R2 with official population metadata and freeze 20 IDs."""
    ensure_assay_layout(root)
    metadata = root / "sources/1000g/metadata/phase3_unrelated.panel"
    if metadata.exists():
        panel_text = metadata.read_text()
    else:
        panel_text = fetch_text(THOUSAND_GENOMES_PANEL_URL)
        metadata.parent.mkdir(parents=True, exist_ok=True)
        metadata.write_text(panel_text)
    hprc_ids = _hprc_sample_ids(root)
    candidates = [
        AssaySample(row.sample, row.population, row.superpopulation, row.sex)
        for row in parse_population_panel(panel_text)
        if row.sample in hprc_ids
    ]
    selected = deterministic_hprc_selection(
        candidates, exclusions=_existing_sample_exclusions(root)
    )
    output = root / "manifests/selected_hprc_r2_samples.tsv"
    partial = output.with_suffix(".tsv.partial")
    with partial.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("sample", "population", "superpopulation", "sex", "seed"))
        for sample in selected:
            writer.writerow(
                (
                    sample.sample,
                    sample.population,
                    sample.superpopulation,
                    sample.sex,
                    HPRC_SELECTION_SEED,
                )
            )
    partial.replace(output)
    ids = root / "manifests/selected_hprc_r2_sample_ids.txt"
    ids.write_text("".join(f"{sample.sample}\n" for sample in selected))
    print(f"Selected {len(selected)} HPRC R2 samples -> {output}")
    return selected


def read_hprc_samples(root: Path) -> list[AssaySample]:
    """Read or create the frozen HPRC R2 selection."""
    path = root / "manifests/selected_hprc_r2_samples.tsv"
    if not path.exists():
        return select_hprc_samples(root)
    with path.open() as handle:
        return [
            AssaySample(
                row["sample"], row["population"], row["superpopulation"], row["sex"]
            )
            for row in csv.DictReader(handle, delimiter="\t")
        ]


def parse_gtf_attributes(text: str) -> dict[str, list[str]]:
    """Parse repeated key/value attributes from one GTF record."""
    parsed: dict[str, list[str]] = defaultdict(list)
    for match in re.finditer(r'(\S+)\s+"([^"]*)"\s*;', text):
        parsed[match.group(1)].append(match.group(2))
    return dict(parsed)


def merge_intervals(
    intervals: Iterable[tuple[str, int, int]],
) -> list[tuple[str, int, int]]:
    """Sort and merge overlapping or bookended BED intervals."""
    contig_order = {contig: number for number, contig in enumerate(ASSAY_CONTIGS)}
    ordered = sorted(
        intervals,
        key=lambda item: (contig_order.get(item[0], 10_000), item[1], item[2]),
    )
    merged: list[tuple[str, int, int]] = []
    for contig, start, end in ordered:
        if not merged or merged[-1][0] != contig or start > merged[-1][2]:
            merged.append((contig, start, end))
        else:
            previous = merged[-1]
            merged[-1] = (contig, previous[1], max(previous[2], end))
    return merged


def _open_gtf(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def build_mane_panel_intervals(
    gtf: Path, genes: set[str], *, padding: int = 20
) -> tuple[list[tuple[str, int, int]], set[str]]:
    """Build merged MANE Select CDS intervals padded in both directions."""
    intervals = []
    observed = set()
    allowed = set(ASSAY_CONTIGS)
    with _open_gtf(gtf) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "CDS":
                continue
            attributes = parse_gtf_attributes(fields[8])
            gene = attributes.get("gene_name", [""])[0]
            if gene not in genes or "MANE_Select" not in attributes.get("tag", []):
                continue
            contig = fields[0]
            if not contig.startswith("chr"):
                contig = f"chr{contig}"
            if contig not in allowed:
                continue
            start = max(0, int(fields[3]) - 1 - padding)
            end = int(fields[4]) + padding
            intervals.append((contig, start, end))
            observed.add(gene)
    return merge_intervals(intervals), observed


def _write_bed(path: Path, intervals: Sequence[tuple[str, int, int]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_suffix(".bed.partial")
    with partial.open("w") as handle:
        for contig, start, end in intervals:
            handle.write(f"{contig}\t{start}\t{end}\n")
    partial.replace(path)
    return path


def normalize_twist_bed(source: Path, destination: Path) -> Path:
    """Freeze the vendor targets as merged chr1-22/X GRCh38 intervals."""
    allowed = set(ASSAY_CONTIGS)
    intervals = []
    with source.open() as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"Malformed Twist BED line {line_number}")
            contig, start, end = fields[0], int(fields[1]), int(fields[2])
            if contig in allowed:
                intervals.append((contig, start, end))
    if not intervals:
        raise ValueError("Twist BED had no chr1-22/X intervals")
    return _write_bed(destination, merge_intervals(intervals))


def prepare_regions(root: Path) -> tuple[Path, Path]:
    """Create frozen WES and ACMG panel interval files."""
    twist_source = _twist_source(root)
    ensembl_source = _ensembl_source(root)
    if not twist_source.exists() or not ensembl_source.exists():
        raise FileNotFoundError("Download the Twist BED and Ensembl GTF first")
    twist = normalize_twist_bed(
        twist_source,
        root / "regions/GRCh38/twist_human_core_exome_hg38.chr1-22-X.merged.bed",
    )
    genes = set(load_acmg_genes())
    intervals, observed = build_mane_panel_intervals(ensembl_source, genes, padding=20)
    missing = sorted(genes - observed)
    if missing:
        raise RuntimeError(f"No MANE Select CDS intervals for ACMG genes: {missing}")
    panel = _write_bed(
        root / "regions/GRCh38/acmg_sf_v3.3.ensembl115.mane_select_cds_pad20.bed",
        intervals,
    )
    output = root / "manifests/assay_regions.tsv"
    with output.with_suffix(".tsv.partial").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            ("assay", "bed", "regions", "target_bases", "sha256", "definition")
        )
        for assay, bed, definition in (
            (
                "wes_twist_core",
                twist,
                "Twist Human Core Exome hg38 targets; merged; chr1-22/X",
            ),
            (
                "panel_acmg_sf_v3.3",
                panel,
                "84 ACMG SF v3.3 genes; Ensembl 115 MANE Select CDS; +/-20 bp",
            ),
        ):
            rows = [line.split("\t") for line in bed.read_text().splitlines()]
            writer.writerow(
                (
                    assay,
                    bed,
                    len(rows),
                    sum(int(row[2]) - int(row[1]) for row in rows),
                    sha256sum(bed),
                    definition,
                )
            )
    output.with_suffix(".tsv.partial").replace(output)
    print(f"Prepared interval definitions: {twist}, {panel}")
    return twist, panel


def _run_split_pipeline(
    *,
    source: Path,
    sample_ids: Path,
    selected: Sequence[str],
    destination_dir: Path,
    log: Path,
    regions: str,
    max_allele_length: int | None,
) -> None:
    """Split a cohort VCF into compact non-reference single-sample shards."""
    work_parent = destination_dir.parent
    temporary = Path(tempfile.mkdtemp(prefix="split.partial.", dir=work_parent))
    commands: list[list[str]] = [
        [
            "bcftools",
            "view",
            "--samples-file",
            str(sample_ids),
            "--no-update",
            # Stream in source order. HPRC's contig order is lexical rather than
            # natural, so --regions would seek backwards through the 2.2-GB file.
            "--targets",
            regions,
            "--output-type",
            "u",
            str(source),
        ]
    ]
    if max_allele_length is not None:
        commands.extend(
            [
                ["bcftools", "norm", "--multiallelics", "-any", "--output-type", "u"],
                [
                    "bcftools",
                    "view",
                    "--types",
                    "snps,indels",
                    "--exclude",
                    f"strlen(REF)>{max_allele_length} || "
                    f'strlen(ALT)>{max_allele_length} || ALT="*"',
                    "--output-type",
                    "u",
                ],
            ]
        )
    else:
        commands.append(
            [
                "bcftools",
                "view",
                "--types",
                "snps,indels",
                "--exclude",
                'ALT~"[<>]" || INFO/SVTYPE!="."',
                "--output-type",
                "u",
            ]
        )
    commands.extend(
        [
            [
                "bcftools",
                "annotate",
                "--remove",
                "^INFO/AF,INFO/AC,INFO/AN,FORMAT",
                "--output-type",
                "u",
            ],
            [
                "bcftools",
                "+split",
                "--output-type",
                "z",
                "--output",
                str(temporary),
                "--samples-file",
                str(sample_ids),
                "--include",
                'GT="alt"',
                "--keep-tags",
                "INFO,FMT/GT",
                "--write-index=tbi",
            ],
        ]
    )
    log.parent.mkdir(parents=True, exist_ok=True)
    processes = []
    with log.open("wb") as log_handle:
        previous = None
        for command in commands:
            process = subprocess.Popen(
                command,
                stdin=previous,
                stdout=subprocess.PIPE if command is not commands[-1] else None,
                stderr=log_handle,
            )
            if previous is not None:
                previous.close()
            previous = process.stdout
            processes.append(process)
        returncodes = [process.wait() for process in processes]
    if any(returncodes):
        shutil.rmtree(temporary)
        raise RuntimeError(f"Split pipeline failed {returncodes}; see {log}")
    missing = [
        sample for sample in selected if not (temporary / f"{sample}.vcf.gz").exists()
    ]
    if missing:
        shutil.rmtree(temporary)
        raise RuntimeError(f"Split pipeline omitted samples: {missing}")
    destination_dir.mkdir(parents=True, exist_ok=True)
    for sample in selected:
        for suffix in (".vcf.gz", ".vcf.gz.tbi"):
            source_path = temporary / f"{sample}{suffix}"
            source_path.replace(destination_dir / source_path.name)
    shutil.rmtree(temporary)


def prepare_chromosome_x(root: Path) -> list[Path]:
    """Prepare chrX shards for the 50 matched WES/panel samples."""
    source = _chromosome_x_source(root)
    if not source.exists() or not Path(f"{source}.tbi").exists():
        raise FileNotFoundError("Download the 1000 Genomes chrX VCF and index first")
    sample_manifest = root / "manifests/selected_samples.tsv"
    sample_ids = root / "manifests/selected_sample_ids.txt"
    with sample_manifest.open() as handle:
        selected = [row["sample"] for row in csv.DictReader(handle, delimiter="\t")]
    destination = root / "work/chromosome_shards/chrX"
    marker = destination / ".complete"
    if not marker.exists():
        _run_split_pipeline(
            source=source,
            sample_ids=sample_ids,
            selected=selected,
            destination_dir=destination,
            log=root / "logs/split_chrX.log",
            regions="chrX",
            max_allele_length=None,
        )
        marker.write_text("complete\n")
    return [destination / f"{sample}.vcf.gz" for sample in selected]


def prepare_hprc(root: Path) -> list[Path]:
    """Prepare the 20 HPRC R2 robustness VCFs in a single source pass."""
    source = _hprc_source(root)
    if not source.exists() or not Path(f"{source}.tbi").exists():
        raise FileNotFoundError("Download the HPRC R2 VCF and index first")
    selected = read_hprc_samples(root)
    sample_ids = root / "manifests/selected_hprc_r2_sample_ids.txt"
    shards = root / "work/hprc_r2_shards"
    marker = shards / ".complete"
    if not marker.exists():
        _run_split_pipeline(
            source=source,
            sample_ids=sample_ids,
            selected=[sample.sample for sample in selected],
            destination_dir=shards,
            log=root / "logs/split_hprc_r2.log",
            regions=",".join(AUTOSOMES),
            max_allele_length=49,
        )
        marker.write_text("complete\n")
    outputs = []
    destination_dir = root / "samples/GRCh38/hprc_r2"
    for sample in selected:
        destination = (
            destination_dir / f"{sample.sample}.GRCh38.hprc_r2.small_variants.vcf.gz"
        )
        if not destination.exists():
            shutil.copy2(shards / f"{sample.sample}.vcf.gz", destination)
            shutil.copy2(
                shards / f"{sample.sample}.vcf.gz.tbi", Path(f"{destination}.tbi")
            )
        outputs.append(destination)
    return outputs


def _subset_vcf(source: Path, bed: Path, destination: Path) -> None:
    run_command(
        [
            "bcftools",
            "view",
            # Compact per-sample inputs are faster to stream once than to seek
            # through ~192k exome intervals. Overlap mode 1 matches -R BED.
            "--targets-file",
            str(bed),
            "--targets-overlap",
            "1",
            "--output-type",
            "z",
            "--output",
            str(destination),
            str(source),
        ]
    )
    run_command(["tabix", "--preset", "vcf", str(destination)])


def prepare_one_interval_sample(
    *,
    sample: str,
    autosomal_vcf: Path,
    chromosome_x_vcf: Path,
    bed: Path,
    destination: Path,
    work_dir: Path,
) -> Path:
    """Subset, combine, sort, BGZF-compress, and index one assay sample."""
    if destination.exists() and Path(f"{destination}.tbi").exists():
        return destination
    destination.parent.mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=f"{sample}.", dir=work_dir))
    autosomal_subset = temporary / "autosomal.vcf.gz"
    chromosome_x_subset = temporary / "chrX.vcf.gz"
    partial = destination.with_name(destination.name + ".partial")
    try:
        _subset_vcf(autosomal_vcf, bed, autosomal_subset)
        _subset_vcf(chromosome_x_vcf, bed, chromosome_x_subset)
        concat = subprocess.Popen(
            [
                "bcftools",
                "concat",
                "--allow-overlaps",
                "--output-type",
                "v",
                str(autosomal_subset),
                str(chromosome_x_subset),
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        sort = subprocess.run(
            [
                "bcftools",
                "sort",
                "--temp-dir",
                str(temporary / "sort.XXXXXX"),
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
                f"Assay preparation failed for {sample}: concat={concat.returncode}, "
                f"sort={sort.returncode}\n{concat_stderr.decode()}\n{sort.stderr.decode()}"
            )
        partial.replace(destination)
        run_command(["tabix", "--preset", "vcf", str(destination)])
    finally:
        if partial.exists():
            partial.unlink()
        shutil.rmtree(temporary)
    return destination


def _read_primary_samples(root: Path) -> list[AssaySample]:
    path = root / "manifests/selected_samples.tsv"
    with path.open() as handle:
        return [
            AssaySample(
                row["sample"], row["population"], row["superpopulation"], row["sex"]
            )
            for row in csv.DictReader(handle, delimiter="\t")
        ]


def _primary_autosomal_vcf(root: Path, sample: AssaySample) -> Path:
    matches = list(
        (root / f"samples/GRCh38/1000g/{sample.superpopulation}").glob(
            f"{sample.sample}.*.vcf.gz"
        )
    )
    if len(matches) != 1:
        raise RuntimeError(
            f"Expected one autosomal source for {sample.sample}: {matches}"
        )
    return matches[0]


def prepare_interval_cohort(
    root: Path, *, cohort: str, bed: Path, workers: int = 4
) -> list[Path]:
    """Prepare one matched 50-sample interval cohort."""
    samples = _read_primary_samples(root)
    output_dir = root / f"samples/GRCh38/{cohort}"
    work_dir = root / f"work/assays/{cohort}"
    futures = {}
    results = []
    with ThreadPoolExecutor(max_workers=workers) as executor:
        for sample in samples:
            destination = output_dir / f"{sample.sample}.GRCh38.{cohort}.vcf.gz"
            future = executor.submit(
                prepare_one_interval_sample,
                sample=sample.sample,
                autosomal_vcf=_primary_autosomal_vcf(root, sample),
                chromosome_x_vcf=root
                / f"work/chromosome_shards/chrX/{sample.sample}.vcf.gz",
                bed=bed,
                destination=destination,
                work_dir=work_dir,
            )
            futures[future] = sample.sample
        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            print(f"Prepared {cohort}: {result.name}")
    return sorted(results)


def _first_invalid_alt(vcf: Path) -> str:
    """Return the first symbolic or multiallelic ALT in one streaming scan."""
    process = subprocess.Popen(
        [
            "bcftools",
            "query",
            "--include",
            'ALT~"[<>*]" || N_ALT>1',
            "--format",
            "%ALT\n",
            str(vcf),
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert process.stdout is not None
    first = process.stdout.readline().decode().strip()
    process.terminate()
    process.communicate()
    return first


def validate_interval_vcf(vcf: Path, *, allowed_contigs: set[str]) -> dict[str, object]:
    """Validate a final sorted, indexed, one-sample interval VCF."""
    errors = []
    try:
        run_command(["bgzip", "--test", str(vcf)])
    except subprocess.CalledProcessError:
        errors.append("bgzip_test")
    index = Path(f"{vcf}.tbi")
    if not index.exists():
        errors.append("missing_index")
    samples = run_command(
        ["bcftools", "query", "--list-samples", str(vcf)], stdout=subprocess.PIPE
    ).stdout.splitlines()
    if len(samples) != 1:
        errors.append(f"sample_count={len(samples)}")
    records = int(
        run_command(
            ["bcftools", "index", "--nrecords", str(vcf)], stdout=subprocess.PIPE
        ).stdout
        or 0
    )
    stats = run_command(
        ["bcftools", "index", "--stats", str(vcf)], stdout=subprocess.PIPE
    ).stdout.splitlines()
    contigs = {line.split("\t", 1)[0] for line in stats if line}
    unexpected = sorted(contigs - allowed_contigs)
    if unexpected:
        errors.append(f"unexpected_contigs={','.join(unexpected)}")
    invalid_alt = _first_invalid_alt(vcf)
    if re.search(r"[<>*]", invalid_alt):
        errors.append("symbolic_allele")
    if "," in invalid_alt:
        errors.append("multiallelic")
    return {
        "status": "PASS" if not errors else "FAIL",
        "errors": ";".join(errors),
        "sample": samples[0] if len(samples) == 1 else "",
        "records": records,
        "contigs": ",".join(sorted(contigs)),
        "bytes": vcf.stat().st_size,
        "sha256": sha256sum(vcf),
    }


def qc_assay_cohorts(root: Path) -> Path:
    """Validate all 120 new VCFs and write one provenance-friendly QC table."""
    expected = {
        "hprc_r2": (20, set(AUTOSOMES)),
        "wes_twist_core": (50, set(ASSAY_CONTIGS)),
        "panel_acmg_sf_v3.3": (50, set(ASSAY_CONTIGS)),
    }
    output = root / "qc/assay_sample_qc.tsv"
    partial = output.with_suffix(".tsv.partial")
    failures = []
    with partial.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            (
                "cohort",
                "sample",
                "records",
                "bytes",
                "contigs",
                "sha256",
                "status",
                "errors",
                "vcf",
            )
        )
        for cohort, (count, allowed) in expected.items():
            vcfs = sorted((root / f"samples/GRCh38/{cohort}").glob("*.vcf.gz"))
            if len(vcfs) != count:
                raise RuntimeError(f"Expected {count} {cohort} VCFs, found {len(vcfs)}")
            for vcf in vcfs:
                result = validate_interval_vcf(vcf, allowed_contigs=allowed)
                writer.writerow(
                    (
                        cohort,
                        result["sample"],
                        result["records"],
                        result["bytes"],
                        result["contigs"],
                        result["sha256"],
                        result["status"],
                        result["errors"],
                        vcf,
                    )
                )
                if result["status"] != "PASS" or int(result["records"]) == 0:
                    failures.append(str(vcf))
    partial.replace(output)
    if failures:
        raise RuntimeError(f"Assay QC failures: {failures}")
    print(f"Assay QC PASS -> {output}")
    return output


def run_all(root: Path, workers: int) -> None:
    """Run the complete resumable assay preparation."""
    preflight(root)
    download_sources(root, workers=min(workers, 3))
    select_hprc_samples(root)
    prepare_regions(root)
    prepare_chromosome_x(root)
    prepare_hprc(root)
    prepare_interval_cohort(
        root,
        cohort="wes_twist_core",
        bed=root / "regions/GRCh38/twist_human_core_exome_hg38.chr1-22-X.merged.bed",
        workers=workers,
    )
    prepare_interval_cohort(
        root,
        cohort="panel_acmg_sf_v3.3",
        bed=root / "regions/GRCh38/acmg_sf_v3.3.ensembl115.mane_select_cds_pad20.bed",
        workers=workers,
    )
    qc_assay_cohorts(root)


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "command",
        choices=(
            "preflight",
            "download",
            "select-hprc",
            "regions",
            "prepare-x",
            "prepare-hprc",
            "prepare-wes",
            "prepare-panel",
            "qc",
            "all",
        ),
    )
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--workers", type=int, default=4)
    return parser


def main() -> None:
    """Run the selected preparation stage."""
    args = build_parser().parse_args()
    root: Path = args.root
    command = args.command
    if command == "preflight":
        preflight(root)
    elif command == "download":
        download_sources(root, workers=min(args.workers, 3))
    elif command == "select-hprc":
        select_hprc_samples(root)
    elif command == "regions":
        prepare_regions(root)
    elif command == "prepare-x":
        prepare_chromosome_x(root)
    elif command == "prepare-hprc":
        prepare_hprc(root)
    elif command == "prepare-wes":
        prepare_interval_cohort(
            root,
            cohort="wes_twist_core",
            bed=root
            / "regions/GRCh38/twist_human_core_exome_hg38.chr1-22-X.merged.bed",
            workers=args.workers,
        )
    elif command == "prepare-panel":
        prepare_interval_cohort(
            root,
            cohort="panel_acmg_sf_v3.3",
            bed=root
            / "regions/GRCh38/acmg_sf_v3.3.ensembl115.mane_select_cds_pad20.bed",
            workers=args.workers,
        )
    elif command == "qc":
        qc_assay_cohorts(root)
    else:
        run_all(root, args.workers)


if __name__ == "__main__":
    main()
