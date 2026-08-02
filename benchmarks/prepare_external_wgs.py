#!/usr/bin/env python3
"""Acquire and prepare independent real-world WGS benchmark cohorts."""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import re
import shutil
import subprocess
import tempfile
import time
import zipfile
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Sequence
from urllib.parse import urljoin
from xml.etree import ElementTree

import requests  # type: ignore[import-untyped]

from benchmarks.prepare_data import (
    download_file,
    run_command,
    sha256sum,
    validate_prepared_vcf,
)
from benchmarks.run_pilot import write_json_atomic
from benchmarks.run_strategy_comparison import VEP_CACHE_NAME, public_strategies

DEFAULT_ROOT = Path("/mnt/data/vcfcache_benchmarks/external_wgs")
DEFAULT_REFERENCE = Path("/mnt/data/resources/reference/ucsc/hg38.fa.gz")
DEFAULT_CACHE_ROOT = Path("/mnt/data/vcfcache_benchmarks/bundled_zenodo_caches")
DEFAULT_PLINK2 = DEFAULT_ROOT / "tools/plink2"
SELECTION_SEED = "vcfcache-paper-external-wgs-v1"
AUTOSOMES = tuple(f"chr{number}" for number in range(1, 23))
COHORT_COUNTS = {"kpgp": (3, 20), "sgdp": (3, 20), "pgp": (3, 12)}
DDBJ_BASE = "https://ddbj.nig.ac.jp/public/public-human-genomes/GRCh38"
DDBJ_COHORTS = {"kpgp": "KPGP", "sgdp": "SGDP"}
PGP_CATALOG_URL = "https://my.pgp-hms.org/public_genetic_data"
SGDP_METADATA_URL = (
    "https://discovery.ucl.ac.uk/id/eprint/1516771/12/"
    "Balloux_5_18_2016_SGDP_supplementary_data_table_1_"
    "sample_information_table.xlsx"
)


@dataclass(frozen=True)
class Candidate:
    """One remotely available source file considered for the benchmark."""

    cohort: str
    sample: str
    provider: str
    population: str
    region: str
    sex: str
    assembly: str
    source_kind: str
    url: str
    index_url: str
    source_name: str
    source_bytes: int
    upstream_md5: str = ""
    documented_overlap: str = "none"
    eligibility: str = "candidate"
    exclusion_reason: str = ""
    landing_url: str = ""


@dataclass(frozen=True)
class Selected:
    """One frozen training or held-out sample allocation."""

    cohort: str
    sample: str
    role: str
    provider: str
    population: str
    region: str
    sex: str
    assembly: str
    source_kind: str
    url: str
    index_url: str
    source_name: str
    source_bytes: int
    upstream_md5: str
    documented_overlap: str
    landing_url: str
    selection_seed: str
    selection_key: str


def stable_key(namespace: str, value: str, seed: str = SELECTION_SEED) -> str:
    """Return an auditable deterministic ranking key."""
    return hashlib.sha256(f"{seed}:{namespace}:{value}".encode()).hexdigest()


def _links(text: str) -> list[str]:
    return [html.unescape(value) for value in re.findall(r'href="([^"]+)"', text)]


def _fetch(url: str) -> requests.Response:
    response = requests.get(url, timeout=90)
    response.raise_for_status()
    return response


def _head_size(url: str) -> int:
    response = requests.head(url, allow_redirects=True, timeout=60)
    response.raise_for_status()
    return int(response.headers.get("content-length", 0))


def catalog_ddbj(cohort: str) -> list[Candidate]:
    """Discover predictable autosomal GRCh38 gVCFs in a DDBJ directory."""
    upstream = DDBJ_COHORTS[cohort]
    base = f"{DDBJ_BASE}/{upstream}/GVCF/"
    directories = sorted(
        value.rstrip("/")
        for value in _links(_fetch(base).text)
        if value.endswith("/") and not value.startswith(("/", "?"))
    )
    candidates = []
    for directory in directories:
        stem = directory
        name = f"{stem}.BQSR.autosome.g.vcf.gz"
        url = urljoin(base, f"{directory}/{name}")
        candidates.append(
            Candidate(
                cohort=cohort,
                sample=stem.removesuffix(".hs38DH"),
                provider="DDBJ_NIG",
                population="Korean" if cohort == "kpgp" else "",
                region="EastAsia" if cohort == "kpgp" else "",
                sex="",
                assembly="GRCh38",
                source_kind="gVCF",
                url=url,
                index_url=f"{url}.tbi",
                source_name=name,
                source_bytes=0,
            )
        )
    return candidates


def _xlsx_rows(content: bytes) -> list[dict[str, str]]:
    """Read the first XLSX worksheet using only the Python standard library."""
    namespace = {"x": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    with tempfile.NamedTemporaryFile(suffix=".xlsx") as handle:
        handle.write(content)
        handle.flush()
        with zipfile.ZipFile(handle.name) as archive:
            shared_root = ElementTree.fromstring(archive.read("xl/sharedStrings.xml"))
            shared = [
                "".join(node.text or "" for node in item.findall(".//x:t", namespace))
                for item in shared_root.findall("x:si", namespace)
            ]
            sheet = ElementTree.fromstring(archive.read("xl/worksheets/sheet1.xml"))

    def column(reference: str) -> int:
        letters = re.match(r"[A-Z]+", reference)
        if letters is None:
            raise ValueError(f"Invalid XLSX cell reference: {reference}")
        result = 0
        for character in letters.group():
            result = result * 26 + ord(character) - ord("A") + 1
        return result - 1

    values: list[list[str]] = []
    for row in sheet.findall(".//x:sheetData/x:row", namespace):
        cells: dict[int, str] = {}
        for cell in row.findall("x:c", namespace):
            value = cell.find("x:v", namespace)
            if value is None or value.text is None:
                continue
            text = value.text
            if cell.get("t") == "s":
                text = shared[int(text)]
            cells[column(cell.get("r", ""))] = text
        if cells:
            width = max(cells) + 1
            values.append([cells.get(index, "") for index in range(width)])
    if not values:
        raise ValueError("SGDP metadata worksheet is empty")
    headers = values[0]
    return [
        {
            headers[index]: value
            for index, value in enumerate(row)
            if index < len(headers)
        }
        for row in values[1:]
    ]


def enrich_sgdp(candidates: Sequence[Candidate]) -> list[Candidate]:
    """Join the primary SGDP sample table and reject documented overlaps."""
    rows = _xlsx_rows(_fetch(SGDP_METADATA_URL).content)
    metadata = {row.get("Sample ID (Illumina)", ""): row for row in rows}
    enriched = []
    for item in candidates:
        row = metadata.get(item.sample)
        if row is None:
            enriched.append(
                Candidate(
                    **{
                        **asdict(item),
                        "eligibility": "excluded",
                        "exclusion_reason": "missing_primary_metadata",
                    }
                )
            )
            continue
        collaborator = row.get("Sample ID (Collaborator)", "")
        public = row.get("Embargo level (X=Fully Public, Y=Signed Letter)", "")
        overlap = bool(re.match(r"(?:HGDP|NA\d+|HG\d+)", collaborator))
        reason = ""
        if public != "X":
            reason = "not_fully_public"
        elif overlap:
            reason = "documented_HGDP_or_1000G_identity"
        enriched.append(
            Candidate(
                **{
                    **asdict(item),
                    "population": row.get("Population ID", ""),
                    "region": row.get("Region", ""),
                    "sex": row.get("Genetic sex assignment", ""),
                    "documented_overlap": collaborator if overlap else "none",
                    "eligibility": "excluded" if reason else "candidate",
                    "exclusion_reason": reason,
                }
            )
        )
    return enriched


def _strip_markup(value: str) -> str:
    return " ".join(html.unescape(re.sub(r"<[^>]+>", " ", value)).split())


def _provider(text: str) -> str:
    for name in ("Sequencing.com", "Nebula", "Gencove", "Dante", "Veritas"):
        if name.lower() in text.lower():
            return name
    return "Participant"


def catalog_pgp() -> list[Candidate]:
    """Discover plausible participant-contributed WGS VCFs from Harvard PGP."""
    page = _fetch(PGP_CATALOG_URL).text
    rows = re.findall(r"<tr data-file-row.*?</tr>", page, re.I | re.S)
    candidates: list[Candidate] = []
    seen: set[str] = set()
    for row in rows:
        opening = row.split(">", maxsplit=1)[0]
        attrs = dict(re.findall(r'data-([\w-]+)="([^"]*)"', opening))
        participant = re.search(r'href="/profile/(hu[A-Fa-f0-9]+)"', row)
        hrefs = re.findall(r'<a[^>]+href="([^"]+)"[^>]*>Download</a>', row, re.I)
        text = _strip_markup(row)
        if participant is None or not hrefs or participant.group(1) in seen:
            continue
        if not re.search(r"\b(?:VCF|gVCF)\b", text, re.I):
            continue
        if not re.search(
            r"whole genome|\bWGS\b|Nebula|Dante|Veritas|Sequencing\.com",
            text,
            re.I,
        ):
            continue
        size = int(attrs.get("file-size", "0") or 0)
        if size < 100_000_000 or re.search(r"\b(?:SV|CNV)\b", text, re.I):
            continue
        sample = participant.group(1)
        seen.add(sample)
        assembly = "GRCh38" if re.search(r"hg38|GRCh38", text, re.I) else "unknown"
        landing_url = urljoin(PGP_CATALOG_URL, html.unescape(hrefs[0]))
        candidates.append(
            Candidate(
                cohort="pgp",
                sample=sample,
                provider=_provider(text),
                population="",
                region="",
                sex="",
                assembly=assembly,
                source_kind="participant_WGS_VCF",
                url=landing_url,
                index_url="",
                source_name=f"{sample}.source.vcf.gz",
                source_bytes=size,
                eligibility="header_probe_required",
                landing_url=landing_url,
            )
        )
    return candidates


def write_tsv(path: Path, rows: Sequence[object]) -> None:
    """Atomically write dataclass rows to a TSV."""
    if not rows:
        raise ValueError(f"Refusing to write empty manifest: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_suffix(path.suffix + ".partial")
    values = [asdict(row) for row in rows]
    with partial.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(values[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(values)
    partial.replace(path)


def read_candidates(path: Path) -> list[Candidate]:
    """Read a candidate manifest."""
    with path.open(newline="") as handle:
        return [
            Candidate(
                **{
                    **row,
                    "source_bytes": int(row["source_bytes"]),
                }
            )
            for row in csv.DictReader(handle, delimiter="\t")
        ]


def catalog(args: argparse.Namespace) -> Path:
    """Fetch source catalogues without downloading genome files."""
    args.root.mkdir(parents=True, exist_ok=True)
    path = args.root / "manifests/external_wgs_candidates.tsv"
    previous = (
        {row.sample: row for row in read_candidates(path) if row.cohort == "pgp"}
        if path.exists()
        else {}
    )
    candidates = catalog_ddbj("kpgp")
    candidates.extend(enrich_sgdp(catalog_ddbj("sgdp")))
    for item in catalog_pgp():
        old = previous.get(item.sample)
        preserve = old and (old.url == item.url or old.landing_url == item.url)
        retryable = old and old.exclusion_reason.startswith("header_probe_failed")
        verified = (
            old
            and old.eligibility == "candidate"
            and bool(old.landing_url)
            and old.url != old.landing_url
        )
        terminal_exclusion = old and old.eligibility == "excluded" and not retryable
        if preserve and (verified or terminal_exclusion):
            candidates.append(old)
        else:
            candidates.append(item)
    write_tsv(path, candidates)
    write_json_atomic(
        args.root / "manifests/catalog_provenance.json",
        {
            "accessed_at": datetime.now(timezone.utc).isoformat(),
            "selection_seed": args.seed,
            "sources": {
                "ddbj": DDBJ_BASE,
                "pgp": PGP_CATALOG_URL,
                "sgdp_metadata": SGDP_METADATA_URL,
            },
            "counts": {
                cohort: sum(item.cohort == cohort for item in candidates)
                for cohort in COHORT_COUNTS
            },
        },
    )
    print(f"Catalogued {len(candidates)} external WGS candidates -> {path}")
    return path


def _ranked(rows: Iterable[Candidate], namespace: str) -> list[Candidate]:
    return sorted(rows, key=lambda item: stable_key(namespace, item.sample))


def _select_sgdp(rows: Sequence[Candidate]) -> list[tuple[Candidate, str]]:
    groups = {
        "Africa": (1, 5),
        "WestEurasia": (1, 5),
        "EastAsia": (0, 5),
        "SouthCentralAsia": (1, 5),
    }
    selected: list[tuple[Candidate, str]] = []
    for group, (train_count, eval_count) in groups.items():
        eligible = [
            row
            for row in rows
            if row.eligibility == "candidate"
            and (
                row.region == group
                or (
                    group == "SouthCentralAsia"
                    and row.region in {"SouthAsia", "CentralAsiaSiberia"}
                )
            )
        ]
        eligible = _ranked(eligible, f"sgdp:{group}")
        if len(eligible) < train_count + eval_count:
            raise RuntimeError(f"Not enough eligible SGDP samples in {group}")
        selected.extend((row, "training") for row in eligible[:train_count])
        selected.extend(
            (row, "evaluation")
            for row in eligible[train_count : train_count + eval_count]
        )
    return selected


def _select_pgp(rows: Sequence[Candidate]) -> list[tuple[Candidate, str]]:
    """Select provider-diverse PGP training and held-out genomes."""
    eligible = _ranked((row for row in rows if row.eligibility == "candidate"), "pgp")
    training: list[Candidate] = []
    used_providers: set[str] = set()
    for row in eligible:
        if row.provider in used_providers:
            continue
        training.append(row)
        used_providers.add(row.provider)
        if len(training) == 3:
            break
    if len(training) < 3:
        raise RuntimeError("PGP needs three eligible sequencing providers for training")
    training_ids = {row.sample for row in training}
    evaluation: list[Candidate] = []
    provider_counts: dict[str, int] = {}
    for row in eligible:
        if row.sample in training_ids or provider_counts.get(row.provider, 0) >= 4:
            continue
        evaluation.append(row)
        provider_counts[row.provider] = provider_counts.get(row.provider, 0) + 1
        if len(evaluation) == 12:
            break
    if len(evaluation) < 12:
        raise RuntimeError(
            f"PGP needs 12 held-out genomes under the provider cap; found {len(evaluation)}"
        )
    return [
        *((row, "training") for row in training),
        *((row, "evaluation") for row in evaluation),
    ]


def select(args: argparse.Namespace) -> Path:
    """Freeze disjoint training and evaluation allocations."""
    source = args.root / "manifests/external_wgs_candidates.tsv"
    if not source.exists():
        catalog(args)
    candidates = read_candidates(source)
    allocations: list[tuple[Candidate, str]] = []
    for cohort in ("kpgp",):
        training, evaluation = COHORT_COUNTS[cohort]
        eligible = [
            row
            for row in candidates
            if row.cohort == cohort and row.eligibility == "candidate"
        ]
        eligible = _ranked(eligible, cohort)
        required = training + evaluation
        if len(eligible) < required:
            raise RuntimeError(f"Not enough eligible {cohort} samples: {len(eligible)}")
        allocations.extend((row, "training") for row in eligible[:training])
        allocations.extend((row, "evaluation") for row in eligible[training:required])
    allocations.extend(
        _select_sgdp([row for row in candidates if row.cohort == "sgdp"])
    )
    try:
        allocations.extend(
            _select_pgp([row for row in candidates if row.cohort == "pgp"])
        )
    except RuntimeError as error:
        if not args.allow_missing_pgp:
            raise RuntimeError(
                f"PGP gate is incomplete: {error}; run probe-pgp before selection"
            ) from error
        print(f"PGP gate pending: {error}; writing only DDBJ selections")
    rows = [
        Selected(
            **{
                **{
                    key: value
                    for key, value in asdict(item).items()
                    if key not in {"eligibility", "exclusion_reason"}
                },
                "role": role,
                "selection_seed": args.seed,
                "selection_key": stable_key(f"{item.cohort}:{role}", item.sample),
            }
        )
        for item, role in allocations
    ]
    output = args.root / "manifests/external_wgs_selection.tsv"
    write_tsv(output, sorted(rows, key=lambda row: (row.cohort, row.role, row.sample)))
    print(f"Selected {len(rows)} samples -> {output}")
    return output


def probe_pgp(args: argparse.Namespace) -> Path:
    """Download ranked PGP candidates until the GRCh38/provider gate is met."""
    path = args.root / "manifests/external_wgs_candidates.tsv"
    candidates = read_candidates(path)
    pending = _ranked(
        (
            row
            for row in candidates
            if row.cohort == "pgp" and row.eligibility == "header_probe_required"
        ),
        "pgp-probe",
    )
    from benchmarks.prepare_data import SourceFile

    updated = {row.sample: row for row in candidates if row.cohort == "pgp"}
    for row in pending:
        destination = args.root / "sources" / "pgp" / row.sample / row.source_name
        try:
            landing_url = row.landing_url or row.url
            download_file(SourceFile("pgp", row.sample, landing_url), destination)
            with destination.open("rb") as handle:
                magic = handle.read(32).lstrip()
            resolved_url = row.url
            if magic.startswith(b"<!DOCTYPE HTML"):
                landing = destination.with_suffix(destination.suffix + ".landing.html")
                destination.replace(landing)
                landing_text = html.unescape(landing.read_text(errors="replace"))
                base = re.search(
                    r"wget[^']*'(https://[^']+/_/)'", landing_text, re.I | re.S
                )
                files = re.findall(
                    r'<a class="item" href="\./([^"]+)">', landing_text, re.I
                )
                usable = [
                    name
                    for name in files
                    if re.search(r"g?vcf(?:\.gz|\.bgz)?$", name, re.I)
                    and not re.search(r"(?:^|[._-])(?:sv|cnv)(?:[._-]|$)", name, re.I)
                ]
                if base is None or len(usable) != 1:
                    raise RuntimeError("PGP landing page lacks one usable WGS VCF")
                resolved_url = urljoin(base.group(1), usable[0])
                download_file(SourceFile("pgp", row.sample, resolved_url), destination)
            _assert_grch38_header(destination)
            samples = run_command(
                ["bcftools", "query", "--list-samples", str(destination)],
                stdout=subprocess.PIPE,
            ).stdout.splitlines()
            if len(samples) != 1:
                raise RuntimeError(f"expected one VCF sample, found {len(samples)}")
            updated[row.sample] = Candidate(
                **{
                    **asdict(row),
                    "assembly": "GRCh38",
                    "url": resolved_url,
                    "landing_url": landing_url,
                    "eligibility": "candidate",
                    "exclusion_reason": "",
                }
            )
        except (RuntimeError, subprocess.CalledProcessError, OSError) as error:
            updated[row.sample] = Candidate(
                **{
                    **asdict(row),
                    "eligibility": "excluded",
                    "exclusion_reason": f"source_ineligible:{type(error).__name__}",
                }
            )
        non_pgp = [item for item in candidates if item.cohort != "pgp"]
        current = [*non_pgp, *updated.values()]
        write_tsv(path, current)
        try:
            _select_pgp(list(updated.values()))
            print("PGP GRCh38/provider gate satisfied")
            return path
        except RuntimeError:
            continue
    raise RuntimeError("PGP catalogue exhausted before the GRCh38/provider gate passed")


def _selected(path: Path) -> list[Selected]:
    with path.open(newline="") as handle:
        return [
            Selected(**{**row, "source_bytes": int(row["source_bytes"])})
            for row in csv.DictReader(handle, delimiter="\t")
        ]


def _source_path(root: Path, row: Selected) -> Path:
    return root / "sources" / row.cohort / row.sample / row.source_name


def download(args: argparse.Namespace) -> None:
    """Download and retain the frozen selected source files."""
    selection = args.root / "manifests/external_wgs_selection.tsv"
    rows = _selected(selection)
    jobs = []
    from benchmarks.prepare_data import SourceFile

    for row in rows:
        jobs.append(
            (
                SourceFile(row.cohort, row.sample, row.url, row.upstream_md5),
                _source_path(args.root, row),
            )
        )
        if row.index_url:
            jobs.append(
                (
                    SourceFile(row.cohort, row.sample, row.index_url),
                    Path(f"{_source_path(args.root, row)}.tbi"),
                )
            )
    with ThreadPoolExecutor(max_workers=args.download_workers) as executor:
        list(executor.map(lambda pair: download_file(*pair), jobs))


def _prepared_path(root: Path, row: Selected) -> Path:
    return (
        root
        / "samples/GRCh38"
        / row.cohort
        / f"{row.sample}.GRCh38.small_variants.vcf.gz"
    )


def _assert_grch38_header(path: Path) -> None:
    header = run_command(
        ["bcftools", "view", "--header-only", str(path)], stdout=subprocess.PIPE
    ).stdout
    if not re.search(r"GRCh38|hg38|assembly=GRCh38", header, re.I):
        lengths = dict(re.findall(r"##contig=<ID=(?:chr)?(\d+),length=(\d+)", header))
        if lengths.get("1") != "248956422":
            raise RuntimeError(f"Source does not establish GRCh38 coordinates: {path}")
    if not re.search(r"^##contig=<ID=chr1,", header, re.MULTILINE):
        raise RuntimeError(
            f"Source needs chr-prefixed GRCh38 contigs for the frozen cache: {path}"
        )


def prepare_one(root: Path, row: Selected, reference: Path) -> Path:
    """Normalize one source into a compact carried-autosomal small-variant VCF."""
    source = _source_path(root, row)
    destination = _prepared_path(root, row)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and Path(f"{destination}.tbi").exists():
        return destination
    _assert_grch38_header(source)
    partial = destination.with_name(destination.name + ".partial")
    first = subprocess.Popen(
        [
            "bcftools",
            "view",
            "--apply-filters",
            "PASS,.",
            "--include",
            'CHROM~"^chr([1-9]|1[0-9]|2[0-2])$"',
            "--output-type",
            "u",
            str(source),
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    normalize = subprocess.Popen(
        [
            "bcftools",
            "norm",
            "--fasta-ref",
            str(reference),
            "--multiallelics",
            "-any",
            "--output-type",
            "u",
        ],
        stdin=first.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert first.stdout is not None
    first.stdout.close()
    filtered = subprocess.Popen(
        [
            "bcftools",
            "view",
            "--include",
            'GT="alt" && ALT!="<NON_REF>" && ALT!="*" && TYPE="snp" || '
            '(GT="alt" && TYPE="indel" && strlen(REF)<=50 && strlen(ALT)<=50)',
            "--output-type",
            "u",
        ],
        stdin=normalize.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert normalize.stdout is not None
    normalize.stdout.close()
    compact = subprocess.run(
        [
            "bcftools",
            "annotate",
            "--remove",
            "^FORMAT/GT",
            "--output-type",
            "z",
            "--output",
            str(partial),
        ],
        stdin=filtered.stdout,
        capture_output=True,
    )
    assert filtered.stdout is not None
    filtered.stdout.close()
    errors = {
        "view": (first.wait(), first.stderr.read() if first.stderr else b""),
        "norm": (
            normalize.wait(),
            normalize.stderr.read() if normalize.stderr else b"",
        ),
        "filter": (filtered.wait(), filtered.stderr.read() if filtered.stderr else b""),
        "compact": (compact.returncode, compact.stderr),
    }
    failed = {name: value for name, value in errors.items() if value[0]}
    if failed:
        raise RuntimeError(f"Preparation failed for {row.sample}: {failed}")
    partial.replace(destination)
    run_command(["tabix", "--preset", "vcf", str(destination)])
    return destination


def prepare(args: argparse.Namespace) -> None:
    """Prepare all downloaded selected sources."""
    rows = _selected(args.root / "manifests/external_wgs_selection.tsv")
    if not args.reference.exists() or not Path(f"{args.reference}.fai").exists():
        raise FileNotFoundError(
            f"Indexed GRCh38 reference is missing: {args.reference}"
        )
    with ThreadPoolExecutor(max_workers=args.prepare_workers) as executor:
        list(
            executor.map(lambda row: prepare_one(args.root, row, args.reference), rows)
        )


def qc(args: argparse.Namespace) -> Path:
    """Validate prepared files and write the campaign input table."""
    rows = _selected(args.root / "manifests/external_wgs_selection.tsv")
    fields = (
        "cohort",
        "sample",
        "role",
        "population",
        "superpopulation",
        "sex",
        "provider",
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
        "relatedness_status",
        "documented_overlap",
    )
    results = []
    for row in rows:
        path = _prepared_path(args.root, row)
        result = validate_prepared_vcf(path, cohort="external")
        count = int(result["records"])
        errors = str(result["errors"])
        if not 2_500_000 <= count <= 6_500_000:
            errors = ",".join(filter(None, (errors, f"implausible_records={count}")))
            result["status"] = "FAIL"
        results.append(
            {
                "cohort": row.cohort,
                "sample": row.sample,
                "role": row.role,
                "population": row.population,
                "superpopulation": row.region,
                "sex": row.sex,
                "provider": row.provider,
                "path": str(path),
                **result,
                "errors": errors,
                "relatedness_status": "PENDING",
                "documented_overlap": row.documented_overlap,
            }
        )
    output = args.root / "qc/external_wgs_qc.tsv"
    output.parent.mkdir(parents=True, exist_ok=True)
    partial = output.with_suffix(".tsv.partial")
    with partial.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(
            {field: row.get(field, "") for field in fields} for row in results
        )
    partial.replace(output)
    failures = [row for row in results if row["status"] != "PASS"]
    if failures:
        raise RuntimeError(f"{len(failures)} external genomes failed QC; see {output}")
    print(f"External WGS QC passed for {len(results)} samples -> {output}")
    return output


def screen_relatedness(args: argparse.Namespace) -> Path:
    """Run PLINK2 KING screening and approve only a relationship-free design."""
    plink2 = args.plink2 if args.plink2.exists() else shutil.which("plink2")
    if plink2 is None:
        raise RuntimeError("plink2 is required for the publication relatedness gate")
    qc_path = args.root / "qc/external_wgs_qc.tsv"
    with qc_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    passing = [row for row in rows if row["status"] == "PASS"]
    if not passing:
        raise RuntimeError("No passing external genomes are available for screening")
    work = args.root / "work/relatedness"
    work.mkdir(parents=True, exist_ok=True)
    merged = work / "external_wgs_all.bcf"
    if not merged.exists():
        run_command(
            [
                "bcftools",
                "merge",
                "--missing-to-ref",
                "--force-samples",
                "--threads",
                str(args.threads),
                "--output-type",
                "b",
                "--output",
                str(merged),
                *(row["path"] for row in passing),
            ]
        )
        run_command(["bcftools", "index", str(merged)])
    prefix = work / "external_wgs"
    run_command(
        [
            plink2,
            "--bcf",
            str(merged),
            "--allow-extra-chr",
            "--snps-only",
            "just-acgt",
            "--max-alleles",
            "2",
            "--maf",
            "0.05",
            "--geno",
            "0.02",
            "--indep-pairwise",
            "200",
            "50",
            "0.2",
            "--out",
            str(prefix),
        ]
    )
    run_command(
        [
            plink2,
            "--bcf",
            str(merged),
            "--allow-extra-chr",
            "--extract",
            f"{prefix}.prune.in",
            "--make-king-table",
            "--king-table-filter",
            "0.0884",
            "--out",
            str(prefix),
        ]
    )
    kinship_path = Path(f"{prefix}.kin0")
    pairs: list[dict[str, str]] = []
    if kinship_path.exists() and kinship_path.stat().st_size:
        with kinship_path.open(newline="") as handle:
            pairs = list(csv.DictReader(handle, delimiter="\t"))
    report = work / "relatedness_report.json"
    write_json_atomic(
        report,
        {
            "created_at": datetime.now(timezone.utc).isoformat(),
            "method": "PLINK2 KING robust kinship on LD-pruned autosomal SNPs",
            "second_degree_threshold": 0.0884,
            "samples": len(passing),
            "related_pairs": pairs,
        },
    )
    if pairs:
        raise RuntimeError(
            f"Related external samples require deterministic replacement; see {report}"
        )
    fields = list(rows[0])
    partial = qc_path.with_suffix(".tsv.partial")
    with partial.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            if row["status"] == "PASS":
                row["relatedness_status"] = "PASS"
            writer.writerow(row)
    partial.replace(qc_path)
    print(f"Relatedness gate passed for {len(passing)} samples -> {report}")
    return report


def build_cache(args: argparse.Namespace) -> None:
    """Build one timed three-genome cache using the bundled AF>=1% recipe."""
    with (args.root / "qc/external_wgs_qc.tsv").open(newline="") as handle:
        rows = [
            row
            for row in csv.DictReader(handle, delimiter="\t")
            if row["cohort"] == args.cohort
            and row["role"] == "training"
            and row["status"] == "PASS"
        ]
    if len(rows) != 3:
        raise RuntimeError(
            f"Expected exactly three passing {args.cohort} training samples"
        )
    public = public_strategies(args.cache_root)
    source_cache = next(item for item in public if item.name == "gnomad_af_0.01")
    database = args.root / "cohort_caches" / args.cohort / "cohort3_database"
    database.parent.mkdir(parents=True, exist_ok=True)
    union = database.parent / "cohort3_training_union.bcf"
    started = time.monotonic()
    started_at = datetime.now(timezone.utc)
    if not union.exists():
        run_command(
            [
                "bcftools",
                "merge",
                "--merge",
                "none",
                "--threads",
                str(args.threads),
                "-Ob",
                "-o",
                str(union),
                *(row["path"] for row in rows),
            ]
        )
        run_command(["bcftools", "index", str(union)])
    vcfcache = Path(__file__).resolve().parents[1] / ".venv/bin/vcfcache"
    if not (database / "blueprint/vcfcache.bcf").exists():
        run_command(
            [
                vcfcache,
                "blueprint-init",
                "--vcf",
                union,
                "--output",
                database,
                "--normalize",
            ]
        )
    cache = database / "cache" / VEP_CACHE_NAME
    if not (cache / "vcfcache_annotated.bcf").exists():
        run_command(
            [
                vcfcache,
                "cache-build",
                "--db",
                database,
                "--name",
                VEP_CACHE_NAME,
                "--anno-config",
                source_cache.cache_dir / "annotation.yaml",
                "--params",
                source_cache.cache_dir / "params.snapshot.yaml",
            ]
        )
    provenance = {
        "cohort": args.cohort,
        "kind": "custom_cohort_three_disjoint_genomes",
        "training_samples": sorted(row["sample"] for row in rows),
        "training_input_sha256": {row["sample"]: row["sha256"] for row in rows},
        "annotation_yaml_sha256": sha256sum(cache / "annotation.yaml"),
        "params_yaml_sha256": sha256sum(cache / "params.snapshot.yaml"),
        "blueprint_sha256": sha256sum(database / "blueprint/vcfcache.bcf"),
        "cache_sha256": sha256sum(cache / "vcfcache_annotated.bcf"),
        "started_at": started_at.isoformat(),
        "completed_at": datetime.now(timezone.utc).isoformat(),
        "wall_seconds": time.monotonic() - started,
        "complete": True,
    }
    write_json_atomic(database.parent / "cohort3_cache_provenance.json", provenance)
    print(json.dumps(provenance, indent=2, sort_keys=True))


def preflight(args: argparse.Namespace) -> None:
    """Check data-root placement, tools, caches, reference, and free space."""
    args.root.mkdir(parents=True, exist_ok=True)
    if not args.root.resolve().is_relative_to(Path("/mnt/data")):
        raise RuntimeError("External WGS root must be under /mnt/data")
    for command in ("bcftools", "bgzip", "tabix", "curl"):
        if shutil.which(command) is None:
            raise RuntimeError(f"Required command is missing: {command}")
    free = shutil.disk_usage(args.root).free
    if free < args.minimum_free_gib * (1 << 30):
        raise RuntimeError(f"Only {free / (1 << 30):.1f} GiB free")
    public_strategies(args.cache_root)
    print(f"External WGS preflight OK: {free / (1 << 30):.1f} GiB free")


def parser() -> argparse.ArgumentParser:
    """Build the external-cohort command line."""
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument(
        "command",
        choices=(
            "preflight",
            "catalog",
            "select",
            "download",
            "probe-pgp",
            "prepare",
            "qc",
            "screen-relatedness",
            "build-cache",
            "all",
        ),
    )
    result.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    result.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE)
    result.add_argument("--cache-root", type=Path, default=DEFAULT_CACHE_ROOT)
    result.add_argument("--plink2", type=Path, default=DEFAULT_PLINK2)
    result.add_argument("--seed", default=SELECTION_SEED)
    result.add_argument("--download-workers", type=int, default=3)
    result.add_argument("--prepare-workers", type=int, default=3)
    result.add_argument("--threads", type=int, default=8)
    result.add_argument("--minimum-free-gib", type=int, default=200)
    result.add_argument("--allow-missing-pgp", action="store_true")
    result.add_argument("--cohort", choices=tuple(COHORT_COUNTS))
    return result


def main() -> None:
    """Run a resumable external WGS preparation stage."""
    args = parser().parse_args()
    for name in ("root", "reference", "cache_root", "plink2"):
        setattr(args, name, getattr(args, name).expanduser().resolve())
    if min(args.download_workers, args.prepare_workers, args.threads) < 1:
        raise ValueError("Worker and thread counts must be positive")
    if args.command == "preflight":
        preflight(args)
    elif args.command == "catalog":
        catalog(args)
    elif args.command == "select":
        select(args)
    elif args.command == "download":
        download(args)
    elif args.command == "probe-pgp":
        probe_pgp(args)
    elif args.command == "prepare":
        prepare(args)
    elif args.command == "qc":
        qc(args)
    elif args.command == "screen-relatedness":
        screen_relatedness(args)
    elif args.command == "build-cache":
        if not args.cohort:
            raise ValueError("--cohort is required for build-cache")
        build_cache(args)
    else:
        preflight(args)
        catalog(args)
        probe_pgp(args)
        select(args)
        download(args)
        prepare(args)
        qc(args)


if __name__ == "__main__":
    main()
