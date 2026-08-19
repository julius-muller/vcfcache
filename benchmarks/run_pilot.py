#!/usr/bin/env python3
"""Run and validate a paired cached/uncached VCFcache publication pilot."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import re
import shutil
import socket
import subprocess
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from itertools import zip_longest
from pathlib import Path
from typing import Iterable, Sequence

import yaml  # type: ignore[import-untyped]

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA_ROOT = Path("/mnt/data/vcfcache_benchmarks")
DEFAULT_INPUT = (
    DEFAULT_DATA_ROOT / "samples/GRCh38/1000g/EAS/HG02079.GRCh38.small_variants.vcf.gz"
)
DEFAULT_CACHE = Path(
    "/mnt/data/vcfcache_benchmarks/bundled_zenodo_caches/"
    "gnomad_v4.1_GRCh38_joint_af001/"
    "cache/vep115.2_everything"
)
DEFAULT_PARAMS = DEFAULT_CACHE / "params.snapshot.yaml"
DEFAULT_CACHE_PROVENANCE = DEFAULT_CACHE.parents[1] / "zenodo_provenance.json"
DEFAULT_CACHE_PROVENANCE_EXPECTED = {
    "alias": "cache-gnomad-v4.1-GRCh38-joint-af001-vep115.2-e",
    "doi": "10.5281/zenodo.18190046",
    "archive_md5": "3ac438461eac0cf42c75717156d7b2d4",
    "archive_md5_verified": True,
    "source": "zenodo_production",
}
VCFCACHE_CMD = REPO_ROOT / ".venv/bin/vcfcache"
TIME_CMD = Path("/usr/bin/time")
MODES = ("uncached", "cached")
PUBLICATION_STATISTICS_MODE = "light"
KNOWN_VEP_IGNORED_CSQ_FIELDS = ("HGNC_ID",)
KNOWN_VEP_UNORDERED_CSQ_FIELDS = ("DOMAINS",)
KNOWN_VEP_ISSUE_URL = "https://github.com/Ensembl/ensembl-vep/issues/1959"


@dataclass(frozen=True)
class PilotConfig:
    """Immutable paths and identifiers for one pilot."""

    data_root: Path
    input_vcf: Path
    cache_dir: Path
    params_file: Path
    replicate: int

    @property
    def sample(self) -> str:
        return self.input_vcf.name.split(".", maxsplit=1)[0]

    @property
    def commit(self) -> str:
        return git_output("rev-parse", "--short=12", "HEAD")

    @property
    def pilot_root(self) -> Path:
        return self.data_root / "pilot" / self.sample / self.commit

    def run_dir(self, mode: str) -> Path:
        return self.pilot_root / f"{mode}_r{self.replicate:02d}"

    @property
    def comparison_path(self) -> Path:
        return self.pilot_root / f"semantic_comparison_r{self.replicate:02d}.json"

    @property
    def summary_path(self) -> Path:
        return self.pilot_root / f"summary_r{self.replicate:02d}.json"


def run_checked(
    args: Sequence[str | Path], *, capture_output: bool = True
) -> subprocess.CompletedProcess[str]:
    """Run a command with stable text decoding."""
    return subprocess.run(
        [str(arg) for arg in args],
        check=True,
        capture_output=capture_output,
        text=True,
    )


def git_output(*args: str) -> str:
    """Return one Git command's stripped stdout."""
    return run_checked(["git", "-C", REPO_ROOT, *args]).stdout.strip()


def sha256sum(path: Path, chunk_size: int = 8 << 20) -> str:
    """Calculate a streaming SHA-256 digest."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def _key_value_file(path: Path) -> dict[str, int]:
    """Parse a cgroup key/value counter file."""
    values: dict[str, int] = {}
    for line in path.read_text().splitlines():
        fields = line.split()
        if len(fields) == 2:
            try:
                values[fields[0]] = int(fields[1])
            except ValueError:
                continue
    return values


def cgroup_v2_snapshot(
    *,
    proc_cgroup: Path = Path("/proc/self/cgroup"),
    cgroup_root: Path = Path("/sys/fs/cgroup"),
) -> dict[str, object] | None:
    """Read aggregate counters for the current Slurm cgroup-v2 step."""
    if not proc_cgroup.exists():
        return None
    relative = next(
        (
            line.removeprefix("0::")
            for line in proc_cgroup.read_text().splitlines()
            if line.startswith("0::")
        ),
        None,
    )
    if relative is None:
        return None
    group = cgroup_root / relative.lstrip("/")
    memory_peak = group / "memory.peak"
    if not memory_peak.exists():
        return None
    snapshot: dict[str, object] = {
        "path": relative,
        "memory_peak_bytes": int(memory_peak.read_text().strip()),
        "memory_current_bytes": int((group / "memory.current").read_text().strip()),
    }
    for name in ("memory.events", "cpu.stat"):
        path = group / name
        if path.exists():
            snapshot[name.replace(".", "_")] = _key_value_file(path)
    io_stat = group / "io.stat"
    if io_stat.exists():
        snapshot["io_stat"] = io_stat.read_text().splitlines()
    return snapshot


def tracked_tree_is_clean() -> bool:
    """Return whether tracked and staged files match HEAD."""
    unstaged = subprocess.run(
        ["git", "-C", REPO_ROOT, "diff", "--quiet"], check=False
    ).returncode
    staged = subprocess.run(
        ["git", "-C", REPO_ROOT, "diff", "--cached", "--quiet"], check=False
    ).returncode
    return unstaged == 0 and staged == 0


def _version(args: Sequence[str | Path]) -> str:
    result = run_checked(args)
    return (result.stdout or result.stderr).splitlines()[0]


def validate_default_cache_provenance(path: Path) -> dict[str, object]:
    """Validate the frozen Zenodo identity of the default publication cache."""
    value = json.loads(path.read_text())
    differences = {
        key: (value.get(key), wanted)
        for key, wanted in DEFAULT_CACHE_PROVENANCE_EXPECTED.items()
        if value.get(key) != wanted
    }
    if differences:
        raise RuntimeError(f"Invalid bundled-cache provenance at {path}: {differences}")
    return value


def validate_requirements_output(
    requirements: str, *, expected_tool_version: str
) -> None:
    """Require successful bcftools and annotation-tool version checks."""
    if not expected_tool_version.strip():
        raise RuntimeError("Cache annotation recipe has no required_tool_version")

    annotation_lines = [
        line.strip()
        for line in requirements.splitlines()
        if line.strip().startswith("annotation tool:")
    ]
    annotation_status = annotation_lines[0] if len(annotation_lines) == 1 else ""
    if "✓" not in annotation_status or expected_tool_version not in annotation_status:
        observed = annotation_status or "missing annotation-tool status"
        raise RuntimeError(
            "Annotation tool requirement check did not pass: "
            f"expected version {expected_tool_version!r}; observed {observed!r}"
        )

    bcftools_lines = [
        line.strip()
        for line in requirements.splitlines()
        if line.strip().startswith("bcftools:")
    ]
    bcftools_status = bcftools_lines[0] if len(bcftools_lines) == 1 else ""
    if "✓" not in bcftools_status:
        observed = bcftools_status or "missing bcftools status"
        raise RuntimeError(
            f"bcftools requirement check did not pass: observed {observed!r}"
        )


def preflight(config: PilotConfig, *, require_clean: bool = True) -> dict[str, object]:
    """Validate the exact tools and immutable inputs used by the pilot."""
    required_paths = [
        config.input_vcf,
        Path(f"{config.input_vcf}.tbi"),
        config.cache_dir / "annotation.yaml",
        config.cache_dir / "vcfcache_annotated.bcf",
        config.cache_dir / "vcfcache_annotated.bcf.csi",
        config.params_file,
        VCFCACHE_CMD,
        TIME_CMD,
    ]
    bundled_provenance: dict[str, object] | None = None
    if config.cache_dir == DEFAULT_CACHE:
        required_paths.append(DEFAULT_CACHE_PROVENANCE)
    missing = [str(path) for path in required_paths if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing pilot prerequisites: {missing}")
    if config.cache_dir == DEFAULT_CACHE:
        bundled_provenance = validate_default_cache_provenance(DEFAULT_CACHE_PROVENANCE)
    for command in ("bcftools", "apptainer"):
        if shutil.which(command) is None:
            raise RuntimeError(f"Required command not found: {command}")
    if require_clean and not tracked_tree_is_clean():
        raise RuntimeError(
            "Tracked worktree changes must be committed before benchmarking"
        )

    requirements = run_checked(
        [
            VCFCACHE_CMD,
            "annotate",
            "-a",
            config.cache_dir,
            "-y",
            config.params_file,
            "--requirements",
        ]
    ).stdout
    annotation = yaml.safe_load((config.cache_dir / "annotation.yaml").read_text())
    if not isinstance(annotation, dict):
        raise RuntimeError("Cache annotation recipe must be a YAML mapping")
    validate_requirements_output(
        requirements,
        expected_tool_version=str(annotation.get("required_tool_version", "")),
    )

    input_records = int(
        run_checked(["bcftools", "index", "--nrecords", config.input_vcf]).stdout
    )
    memory_kib = 0
    for line in Path("/proc/meminfo").read_text().splitlines():
        if line.startswith("MemTotal:"):
            memory_kib = int(line.split()[1])
            break
    manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "commit": config.commit,
        "tracked_tree_clean": tracked_tree_is_clean(),
        "hostname": socket.gethostname(),
        "platform": platform.platform(),
        "cpu_count": os.cpu_count(),
        "memory_bytes": memory_kib * 1024,
        "cpu_model": next(
            (
                line.split(":", maxsplit=1)[1].strip()
                for line in Path("/proc/cpuinfo").read_text().splitlines()
                if line.startswith("model name")
            ),
            "unknown",
        ),
        "slurm": {
            key: value for key, value in os.environ.items() if key.startswith("SLURM_")
        },
        "input_vcf": str(config.input_vcf),
        "input_records": input_records,
        "input_bytes": config.input_vcf.stat().st_size,
        "input_sha256": sha256sum(config.input_vcf),
        "cache_dir": str(config.cache_dir),
        "bundled_cache_provenance": bundled_provenance,
        "cache_bcf_bytes": (config.cache_dir / "vcfcache_annotated.bcf").stat().st_size,
        "annotation_yaml_sha256": sha256sum(config.cache_dir / "annotation.yaml"),
        "params_yaml_sha256": sha256sum(config.params_file),
        "vcfcache_version": _version([VCFCACHE_CMD, "--version"]),
        "bcftools_version": _version(["bcftools", "--version"]),
        "apptainer_version": _version(["apptainer", "--version"]),
        "requirements_output": requirements,
    }
    config.pilot_root.mkdir(parents=True, exist_ok=True)
    write_json_atomic(config.pilot_root / "environment.json", manifest)
    print(
        f"Pilot preflight OK: {config.sample}, {input_records:,} records, "
        f"commit {config.commit}"
    )
    return manifest


def annotation_command(config: PilotConfig, mode: str, run_dir: Path) -> list[str]:
    """Build the paired command; only --uncached differs between modes."""
    if mode not in MODES:
        raise ValueError(f"Unknown mode: {mode}")
    command = [
        str(VCFCACHE_CMD),
        "annotate",
        "-a",
        str(config.cache_dir),
        "-i",
        str(config.input_vcf),
        "-o",
        str(run_dir / "output.bcf"),
        "--stats-dir",
        str(run_dir / "stats"),
        "--statistics",
        PUBLICATION_STATISTICS_MODE,
        "-y",
        str(config.params_file),
        "--force",
    ]
    if mode == "uncached":
        command.append("--uncached")
    return command


def parse_elapsed(value: str) -> float:
    """Parse GNU time's h:mm:ss or m:ss elapsed format."""
    fields = value.split(":")
    if len(fields) == 3:
        hours, minutes, seconds = fields
        return int(hours) * 3600 + int(minutes) * 60 + float(seconds)
    if len(fields) == 2:
        minutes, seconds = fields
        return int(minutes) * 60 + float(seconds)
    return float(value)


def parse_gnu_time(path: Path) -> dict[str, float | int]:
    """Extract stable resource fields from GNU time --verbose output."""
    values: dict[str, str] = {}
    for line in path.read_text().splitlines():
        if ": " not in line:
            continue
        key, value = line.strip().rsplit(": ", maxsplit=1)
        values[key] = value.strip()
    elapsed_key = next(
        key for key in values if key.startswith("Elapsed (wall clock) time")
    )
    return {
        "wall_seconds_gnu": parse_elapsed(values[elapsed_key]),
        "user_seconds": float(values["User time (seconds)"]),
        "system_seconds": float(values["System time (seconds)"]),
        "cpu_percent": int(values["Percent of CPU this job got"].rstrip("%")),
        "max_rss_kib": int(values["Maximum resident set size (kbytes)"]),
        "filesystem_inputs": int(values["File system inputs"]),
        "filesystem_outputs": int(values["File system outputs"]),
        "voluntary_context_switches": int(values["Voluntary context switches"]),
        "involuntary_context_switches": int(values["Involuntary context switches"]),
    }


def find_stats_dir(run_dir: Path) -> Path:
    """Locate the single stats directory produced by VCFcache."""
    matches = list((run_dir / "stats").glob("*_vcstats"))
    if len(matches) != 1:
        raise RuntimeError(f"Expected one VCFcache stats directory, found {matches}")
    return matches[0]


def calculate_cache_hit_rate(
    mode: str, variant_counts: dict[str, object], output_records: int
) -> float | None:
    """Calculate hit rate using the normalized total when input count is absent."""
    if mode != "cached":
        return None
    total = (
        variant_counts.get("input_variants")
        or variant_counts.get("total_output")
        or output_records
    )
    tool_annotated = variant_counts.get("tool_annotated")
    if not isinstance(total, int) or not isinstance(tool_annotated, int) or total <= 0:
        return None
    return 1 - (tool_annotated / total)


def write_json_atomic(path: Path, value: object) -> None:
    """Write JSON without exposing a partially-written result."""
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_suffix(path.suffix + ".partial")
    partial.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    partial.replace(path)


def run_one(config: PilotConfig, mode: str) -> dict[str, object]:
    """Execute and instrument one immutable pilot run."""
    cgroup_before = cgroup_v2_snapshot()
    run_dir = config.run_dir(mode)
    if run_dir.exists():
        raise FileExistsError(f"Pilot run already exists: {run_dir}")
    run_dir.mkdir(parents=True)
    command = annotation_command(config, mode, run_dir)
    time_file = run_dir / "resource_usage.txt"
    log_file = run_dir / "command.log"
    write_json_atomic(
        run_dir / "command.json",
        {
            "mode": mode,
            "command": command,
            "commit": config.commit,
            "started_at": datetime.now(timezone.utc).isoformat(),
        },
    )
    wrapped = [
        str(TIME_CMD),
        "--verbose",
        "--output",
        str(time_file),
        *command,
    ]
    started_at = datetime.now(timezone.utc)
    started = time.monotonic()
    environment = os.environ.copy()
    environment.update({"LC_ALL": "C", "LANG": "C", "TMPDIR": str(run_dir / "tmp")})
    (run_dir / "tmp").mkdir()
    with log_file.open("w") as log:
        completed = subprocess.run(
            wrapped,
            stdout=log,
            stderr=subprocess.STDOUT,
            env=environment,
            check=False,
        )
    wall_seconds = time.monotonic() - started
    if completed.returncode:
        write_json_atomic(
            run_dir / "failure.json",
            {
                "returncode": completed.returncode,
                "wall_seconds": wall_seconds,
                "log": str(log_file),
            },
        )
        raise RuntimeError(f"{mode} pilot failed; see {log_file}")

    output_bcf = run_dir / "output.bcf"
    output_index = Path(f"{output_bcf}.csi")
    if not output_bcf.exists() or not output_index.exists():
        raise RuntimeError(f"{mode} output or index missing")
    output_records = int(
        run_checked(["bcftools", "index", "--nrecords", output_bcf]).stdout
    )
    header = run_checked(["bcftools", "view", "--header-only", output_bcf]).stdout
    if "##INFO=<ID=CSQ" not in header:
        raise RuntimeError(f"{mode} output has no CSQ header")

    stats_dir = find_stats_dir(run_dir)
    compare_stats_path = stats_dir / "compare_stats.yaml"
    compare_stats = yaml.safe_load(compare_stats_path.read_text()) or {}
    variant_counts = compare_stats.get("variant_counts", {}) or {}
    hit_rate = calculate_cache_hit_rate(mode, variant_counts, output_records)

    cgroup_after = cgroup_v2_snapshot()
    metrics: dict[str, object] = {
        "mode": mode,
        "replicate": config.replicate,
        "sample": config.sample,
        "commit": config.commit,
        "started_at": started_at.isoformat(),
        "completed_at": datetime.now(timezone.utc).isoformat(),
        "wall_seconds": wall_seconds,
        **parse_gnu_time(time_file),
        "input_vcf": str(config.input_vcf),
        "output_bcf": str(output_bcf),
        "output_bytes": output_bcf.stat().st_size,
        "output_sha256": sha256sum(output_bcf),
        "output_index_sha256": sha256sum(output_index),
        "output_records": output_records,
        "stats_dir": str(stats_dir),
        "variant_counts": variant_counts,
        "cache_hit_rate": hit_rate,
        "command": command,
        "statistics_mode": PUBLICATION_STATISTICS_MODE,
        "cgroup_v2": {"before": cgroup_before, "after": cgroup_after},
    }
    if os.environ.get("VCFCACHE_REQUIRE_CGROUP_METRICS") == "1" and not cgroup_after:
        raise RuntimeError("Required Slurm cgroup-v2 metrics are unavailable")
    write_json_atomic(run_dir / "metrics.json", metrics)
    print(f"{mode.capitalize()} pilot completed in {wall_seconds:.1f}s")
    return metrics


def _query_process(path: Path) -> subprocess.Popen[str]:
    query_format = (
        r"%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/AC\t%INFO/AN"
        r"[\t%GT]\t%INFO/CSQ\n"
    )
    contigs = sorted(
        line.split("\t", maxsplit=1)[0]
        for line in run_checked(["bcftools", "index", "--stats", path])
        .stdout.strip()
        .splitlines()
        if line
    )
    command = ["bcftools", "query", "--allow-undef-tags"]
    if contigs:
        # VEP --fork can preserve coordinates while emitting contigs in a
        # different header order between cached and uncached runs. Query the
        # indexed contigs in one canonical order before streaming comparison.
        command.extend(["--regions", ",".join(contigs)])
    command.extend(["--format", query_format, str(path)])
    return subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )


def _csq_header(path: Path) -> str:
    header = run_checked(["bcftools", "view", "--header-only", path]).stdout
    return next(
        (line for line in header.splitlines() if line.startswith("##INFO=<ID=CSQ")),
        "",
    )


def _csq_fields(header: str) -> tuple[str, ...]:
    """Extract the ordered VEP CSQ field names from its INFO header."""
    match = re.search(r'Format: ([^"]+)', header)
    return tuple(match.group(1).split("|")) if match else ()


def _locus_groups(
    stream: Iterable[str],
) -> Iterable[tuple[tuple[str, str], list[list[str]]]]:
    """Yield records grouped by CHROM/POS, preserving order within each locus."""
    locus: tuple[str, str] | None = None
    records: list[list[str]] = []
    for record_number, line in enumerate(stream, start=1):
        fields = line.rstrip("\n").split("\t", maxsplit=8)
        if len(fields) != 9:
            raise RuntimeError(f"Unexpected bcftools query output at {record_number}")
        next_locus = (fields[0], fields[1])
        if locus is not None and next_locus != locus:
            yield locus, records
            records = []
        locus = next_locus
        records.append(fields)
    if locus is not None:
        yield locus, records


def _canonical_csq(
    value: str,
    ignored_indices: tuple[int, ...] = (),
    unordered_indices: tuple[int, ...] = (),
) -> str:
    items: list[str] = []
    for item in value.split(","):
        values = item.split("|")
        for index in unordered_indices:
            if index < len(values) and values[index]:
                values[index] = "&".join(sorted(values[index].split("&")))
        for index in ignored_indices:
            if index < len(values):
                values[index] = ""
        items.append("|".join(values))
    return ",".join(sorted(items))


def _canonical_record(
    fields: list[str],
    ignored_indices: tuple[int, ...] = (),
    unordered_indices: tuple[int, ...] = (),
) -> tuple[tuple[str, ...], str]:
    return tuple(fields[:8]), _canonical_csq(
        fields[8], ignored_indices, unordered_indices
    )


def semantic_compare(
    cached_bcf: Path, uncached_bcf: Path, *, mismatch_limit: int = 20
) -> dict[str, object]:
    """Compare semantic records, with narrowly documented VEP exceptions."""
    cached_header = _csq_header(cached_bcf)
    uncached_header = _csq_header(uncached_bcf)
    csq_fields = _csq_fields(cached_header) if cached_header == uncached_header else ()
    ignored_csq_fields = tuple(
        field for field in KNOWN_VEP_IGNORED_CSQ_FIELDS if field in csq_fields
    )
    ignored_indices = tuple(csq_fields.index(field) for field in ignored_csq_fields)
    unordered_csq_fields = tuple(
        field for field in KNOWN_VEP_UNORDERED_CSQ_FIELDS if field in csq_fields
    )
    unordered_indices = tuple(csq_fields.index(field) for field in unordered_csq_fields)
    cached = _query_process(cached_bcf)
    uncached = _query_process(uncached_bcf)
    assert cached.stdout is not None
    assert uncached.stdout is not None
    sentinel = object()
    records = 0
    key_mismatches = 0
    annotation_mismatches = 0
    ignored_annotation_mismatches = 0
    annotation_order_only = 0
    record_order_only = 0
    examples: list[dict[str, object]] = []
    cached_digest = hashlib.sha256()
    uncached_digest = hashlib.sha256()

    cached_groups = _locus_groups(cached.stdout)
    uncached_groups = _locus_groups(uncached.stdout)
    for group_number, pair in enumerate(
        zip_longest(cached_groups, uncached_groups, fillvalue=sentinel), start=1
    ):
        cached_group, uncached_group = pair
        if cached_group is sentinel or uncached_group is sentinel:
            key_mismatches += 1
            if len(examples) < mismatch_limit:
                examples.append({"group": group_number, "kind": "locus_count_mismatch"})
            continue
        assert isinstance(cached_group, tuple)
        assert isinstance(uncached_group, tuple)
        cached_locus, cached_fields = cached_group
        uncached_locus, uncached_fields = uncached_group
        records += max(len(cached_fields), len(uncached_fields))
        if cached_locus != uncached_locus:
            key_mismatches += 1
            if len(examples) < mismatch_limit:
                examples.append(
                    {
                        "group": group_number,
                        "kind": "locus",
                        "cached": cached_locus,
                        "uncached": uncached_locus,
                    }
                )
        cached_original_keys = [tuple(fields[:8]) for fields in cached_fields]
        uncached_original_keys = [tuple(fields[:8]) for fields in uncached_fields]
        cached_fields.sort(key=lambda fields: tuple(fields[:8]))
        uncached_fields.sort(key=lambda fields: tuple(fields[:8]))
        if cached_original_keys != uncached_original_keys and sorted(
            cached_original_keys
        ) == sorted(uncached_original_keys):
            record_order_only += 1

        for occurrence, record_pair in enumerate(
            zip_longest(cached_fields, uncached_fields, fillvalue=sentinel), start=1
        ):
            cached_record, uncached_record = record_pair
            if cached_record is sentinel or uncached_record is sentinel:
                key_mismatches += 1
                if len(examples) < mismatch_limit:
                    examples.append(
                        {
                            "group": group_number,
                            "locus": cached_locus,
                            "occurrence": occurrence,
                            "kind": "record_count_mismatch",
                        }
                    )
                continue
            assert isinstance(cached_record, list)
            assert isinstance(uncached_record, list)
            cached_key, cached_raw_canonical = _canonical_record(
                cached_record, unordered_indices=unordered_indices
            )
            uncached_key, uncached_raw_canonical = _canonical_record(
                uncached_record, unordered_indices=unordered_indices
            )
            _, cached_canonical = _canonical_record(
                cached_record, ignored_indices, unordered_indices
            )
            _, uncached_canonical = _canonical_record(
                uncached_record, ignored_indices, unordered_indices
            )
            if cached_key != uncached_key:
                key_mismatches += 1
                if len(examples) < mismatch_limit:
                    examples.append(
                        {
                            "group": group_number,
                            "kind": "variant_or_input_info",
                            "cached": cached_key,
                            "uncached": uncached_key,
                        }
                    )
            if cached_record[8] != uncached_record[8]:
                if cached_raw_canonical == uncached_raw_canonical:
                    annotation_order_only += 1
                elif cached_canonical == uncached_canonical:
                    ignored_annotation_mismatches += 1
                    if len(examples) < mismatch_limit:
                        examples.append(
                            {
                                "group": group_number,
                                "kind": "CSQ_ignored_known_vep_issue",
                                "key": cached_key[:4],
                                "ignored_fields": list(ignored_csq_fields),
                                "reference": KNOWN_VEP_ISSUE_URL,
                            }
                        )
                else:
                    annotation_mismatches += 1
                    if len(examples) < mismatch_limit:
                        examples.append(
                            {
                                "group": group_number,
                                "kind": "CSQ",
                                "key": cached_key[:4],
                                "cached": cached_record[8][:500],
                                "uncached": uncached_record[8][:500],
                            }
                        )
            cached_digest.update(
                ("\t".join(cached_key) + "\t" + cached_canonical + "\n").encode()
            )
            uncached_digest.update(
                ("\t".join(uncached_key) + "\t" + uncached_canonical + "\n").encode()
            )

    cached_stderr = cached.communicate()[1]
    uncached_stderr = uncached.communicate()[1]
    if cached.returncode or uncached.returncode:
        raise RuntimeError(
            "bcftools query failed: "
            f"cached={cached.returncode} {cached_stderr}; "
            f"uncached={uncached.returncode} {uncached_stderr}"
        )
    return {
        "comparator": "vep_semantic_with_documented_exceptions_v1",
        "semantic_pass": (
            key_mismatches == 0
            and annotation_mismatches == 0
            and cached_header == uncached_header
            and cached_digest.hexdigest() == uncached_digest.hexdigest()
        ),
        "records_compared": records,
        "key_mismatches": key_mismatches,
        "annotation_mismatches": annotation_mismatches,
        "raw_annotation_mismatches": (
            annotation_mismatches + ignored_annotation_mismatches
        ),
        "ignored_annotation_mismatches": ignored_annotation_mismatches,
        "ignored_csq_fields": list(ignored_csq_fields),
        "unordered_csq_fields": list(unordered_csq_fields),
        "ignored_difference_reference": (
            KNOWN_VEP_ISSUE_URL if ignored_annotation_mismatches else None
        ),
        "annotation_order_only": annotation_order_only,
        "record_order_only_loci": record_order_only,
        "csq_headers_equal": cached_header == uncached_header,
        "cached_semantic_sha256": cached_digest.hexdigest(),
        "uncached_semantic_sha256": uncached_digest.hexdigest(),
        "examples": examples,
    }


def _strict_relevant_header(path: Path) -> tuple[str, ...]:
    """Return semantic header definitions in an order-independent form."""
    header = run_checked(["bcftools", "view", "--header-only", path]).stdout
    prefixes = ("##contig=<", "##FILTER=<", "##INFO=<", "##FORMAT=<")
    definitions = sorted(
        line for line in header.splitlines() if line.startswith(prefixes)
    )
    columns = [line for line in header.splitlines() if line.startswith("#CHROM\t")]
    if len(columns) != 1:
        raise RuntimeError(f"Expected one #CHROM header in {path}")
    return (*definitions, columns[0])


def _canonical_strict_info(value: str) -> str:
    """Canonicalize INFO tag order and CSQ entry order, but no values."""
    if value in {"", "."}:
        return value
    fields: list[tuple[str, str]] = []
    for item in value.split(";"):
        key, separator, content = item.partition("=")
        if key == "CSQ" and separator:
            content = ",".join(sorted(content.split(",")))
        fields.append((key, item if not separator else f"{key}={content}"))
    return ";".join(item for _key, item in sorted(fields))


def _full_record_process(path: Path) -> subprocess.Popen[str]:
    contigs = sorted(
        line.split("\t", maxsplit=1)[0]
        for line in run_checked(["bcftools", "index", "--stats", path])
        .stdout.strip()
        .splitlines()
        if line
    )
    command = ["bcftools", "view"]
    if contigs:
        # Both BCFs can be validly indexed while assigning different numeric
        # header IDs to contigs. Traverse their indexes by the same contig-name
        # order so header layout cannot masquerade as a record difference.
        command.extend(["--regions", ",".join(contigs)])
    command.extend(["--no-header", str(path)])
    return subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )


def _strict_locus_groups(
    stream: Iterable[str],
) -> Iterable[tuple[tuple[str, str], list[str]]]:
    """Yield complete canonical records grouped by their coordinate."""
    locus: tuple[str, str] | None = None
    records: list[str] = []
    for record_number, line in enumerate(stream, start=1):
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 8:
            raise RuntimeError(f"Malformed VCF record {record_number}")
        next_locus = (fields[0], fields[1])
        if locus is not None and next_locus != locus:
            yield locus, records
            records = []
        locus = next_locus
        fields[7] = _canonical_strict_info(fields[7])
        records.append("\t".join(fields) + "\n")
    if locus is not None:
        yield locus, records


def strict_semantic_compare(
    cached_bcf: Path, uncached_bcf: Path, *, mismatch_limit: int = 20
) -> dict[str, object]:
    """Compare complete fastVEP records without VEP-specific exceptions."""
    cached_header = _strict_relevant_header(cached_bcf)
    uncached_header = _strict_relevant_header(uncached_bcf)
    cached = _full_record_process(cached_bcf)
    uncached = _full_record_process(uncached_bcf)
    assert cached.stdout is not None
    assert uncached.stdout is not None
    sentinel = object()
    records = 0
    mismatches = 0
    locus_mismatches = 0
    record_order_only_loci = 0
    examples: list[dict[str, object]] = []
    cached_digest = hashlib.sha256()
    uncached_digest = hashlib.sha256()
    cached_groups = _strict_locus_groups(cached.stdout)
    uncached_groups = _strict_locus_groups(uncached.stdout)
    for group_number, pair in enumerate(
        zip_longest(cached_groups, uncached_groups, fillvalue=sentinel), start=1
    ):
        cached_group, uncached_group = pair
        if cached_group is sentinel or uncached_group is sentinel:
            mismatches += 1
            locus_mismatches += 1
            if len(examples) < mismatch_limit:
                examples.append({"group": group_number, "kind": "locus_count_mismatch"})
            continue
        assert isinstance(cached_group, tuple)
        assert isinstance(uncached_group, tuple)
        cached_locus, cached_records = cached_group
        uncached_locus, uncached_records = uncached_group
        records += max(len(cached_records), len(uncached_records))
        if cached_locus != uncached_locus:
            locus_mismatches += 1
            if len(examples) < mismatch_limit:
                examples.append(
                    {
                        "group": group_number,
                        "kind": "locus",
                        "cached": cached_locus,
                        "uncached": uncached_locus,
                    }
                )
        if cached_records != uncached_records and sorted(cached_records) == sorted(
            uncached_records
        ):
            record_order_only_loci += 1
        cached_records.sort()
        uncached_records.sort()
        for occurrence, record_pair in enumerate(
            zip_longest(cached_records, uncached_records, fillvalue=sentinel), start=1
        ):
            cached_record, uncached_record = record_pair
            if cached_record is not sentinel:
                assert isinstance(cached_record, str)
                cached_digest.update(cached_record.encode())
            if uncached_record is not sentinel:
                assert isinstance(uncached_record, str)
                uncached_digest.update(uncached_record.encode())
            if cached_record != uncached_record:
                mismatches += 1
                if len(examples) < mismatch_limit:
                    examples.append(
                        {
                            "group": group_number,
                            "locus": cached_locus,
                            "occurrence": occurrence,
                            "kind": "complete_record",
                            "cached": (
                                cached_record[:1000]
                                if isinstance(cached_record, str)
                                else None
                            ),
                            "uncached": (
                                uncached_record[:1000]
                                if isinstance(uncached_record, str)
                                else None
                            ),
                        }
                    )
    cached_stderr = cached.communicate()[1]
    uncached_stderr = uncached.communicate()[1]
    if cached.returncode or uncached.returncode:
        raise RuntimeError(
            "bcftools view failed: "
            f"cached={cached.returncode} {cached_stderr}; "
            f"uncached={uncached.returncode} {uncached_stderr}"
        )
    cached_sha256 = cached_digest.hexdigest()
    uncached_sha256 = uncached_digest.hexdigest()
    headers_equal = cached_header == uncached_header
    return {
        "comparator": "fastvep_complete_record_and_header_v3",
        "semantic_pass": (
            mismatches == 0
            and locus_mismatches == 0
            and headers_equal
            and cached_sha256 == uncached_sha256
        ),
        "records_compared": records,
        "record_mismatches": mismatches,
        "locus_mismatches": locus_mismatches,
        "record_order_only_loci": record_order_only_loci,
        "raw_annotation_mismatches": mismatches,
        "ignored_annotation_mismatches": 0,
        "relevant_headers_equal": headers_equal,
        "cached_header_sha256": hashlib.sha256(
            "\n".join(cached_header).encode()
        ).hexdigest(),
        "uncached_header_sha256": hashlib.sha256(
            "\n".join(uncached_header).encode()
        ).hexdigest(),
        "cached_semantic_sha256": cached_sha256,
        "uncached_semantic_sha256": uncached_sha256,
        "canonicalization": [
            "INFO tag order",
            "CSQ entry order",
            "complete-record order within identical CHROM/POS loci",
            "indexed contig traversal order by contig name",
        ],
        "ignored_fields": [],
        "examples": examples,
    }


def compare_outputs(config: PilotConfig) -> dict[str, object]:
    """Run the full semantic comparison and save its report."""
    cached = config.run_dir("cached") / "output.bcf"
    uncached = config.run_dir("uncached") / "output.bcf"
    report = semantic_compare(cached, uncached)
    write_json_atomic(config.comparison_path, report)
    status = "PASS" if report["semantic_pass"] else "FAIL"
    print(f"Semantic comparison {status}: {report['records_compared']:,} records")
    if not report["semantic_pass"]:
        raise RuntimeError(
            f"Cached and uncached outputs differ; see " f"{config.comparison_path}"
        )
    return report


def summarize(config: PilotConfig) -> dict[str, object]:
    """Create the paired pilot summary consumed by later plotting code."""
    metrics = {
        mode: json.loads((config.run_dir(mode) / "metrics.json").read_text())
        for mode in MODES
    }
    uncached_seconds = float(metrics["uncached"]["wall_seconds"])
    cached_seconds = float(metrics["cached"]["wall_seconds"])
    summary = {
        "sample": config.sample,
        "commit": config.commit,
        "replicate": config.replicate,
        "input_records": (
            metrics["uncached"]["variant_counts"].get("input_variants")
            or metrics["uncached"]["variant_counts"].get("total_output")
            or metrics["uncached"]["output_records"]
        ),
        "uncached_wall_seconds": uncached_seconds,
        "cached_wall_seconds": cached_seconds,
        "speedup": uncached_seconds / cached_seconds,
        "wall_seconds_saved": uncached_seconds - cached_seconds,
        "fraction_wall_time_saved": 1 - (cached_seconds / uncached_seconds),
        "cache_hit_rate": metrics["cached"].get("cache_hit_rate"),
        "uncached_user_seconds": metrics["uncached"]["user_seconds"],
        "cached_user_seconds": metrics["cached"]["user_seconds"],
        "uncached_max_rss_kib": metrics["uncached"]["max_rss_kib"],
        "cached_max_rss_kib": metrics["cached"]["max_rss_kib"],
        "uncached_cgroup_memory_peak_bytes": (
            metrics["uncached"].get("cgroup_v2", {}).get("after") or {}
        ).get("memory_peak_bytes"),
        "cached_cgroup_memory_peak_bytes": (
            metrics["cached"].get("cgroup_v2", {}).get("after") or {}
        ).get("memory_peak_bytes"),
    }
    comparison_path = config.comparison_path
    legacy_comparison_path = config.pilot_root / "semantic_comparison.json"
    if (
        not comparison_path.exists()
        and config.replicate == 1
        and legacy_comparison_path.exists()
    ):
        comparison_path = legacy_comparison_path
    if comparison_path.exists():
        summary["semantic_comparison"] = json.loads(comparison_path.read_text())
    write_json_atomic(config.summary_path, summary)
    print(
        f"Pilot speedup: {summary['speedup']:.2f}x; "
        f"wall time saved: {summary['wall_seconds_saved']:.1f}s"
    )
    return summary


def build_parser() -> argparse.ArgumentParser:
    """Build the pilot CLI."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "command", choices=("preflight", "run", "compare", "summarize", "all")
    )
    parser.add_argument("--mode", choices=MODES)
    parser.add_argument("--replicate", type=int, default=1)
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--params", type=Path, default=DEFAULT_PARAMS)
    return parser


def main() -> None:
    """Command-line entry point."""
    args = build_parser().parse_args()
    config = PilotConfig(
        data_root=args.data_root.expanduser().resolve(),
        input_vcf=args.input.expanduser().resolve(),
        cache_dir=args.cache.expanduser().resolve(),
        params_file=args.params.expanduser().resolve(),
        replicate=args.replicate,
    )
    if config.replicate < 1:
        raise ValueError("--replicate must be positive")
    if args.command == "preflight":
        preflight(config)
    elif args.command == "run":
        if not args.mode:
            raise ValueError("--mode is required for run")
        preflight(config)
        run_one(config, args.mode)
    elif args.command == "compare":
        compare_outputs(config)
    elif args.command == "summarize":
        summarize(config)
    else:
        preflight(config)
        run_one(config, "uncached")
        run_one(config, "cached")
        compare_outputs(config)
        summarize(config)


if __name__ == "__main__":
    main()
