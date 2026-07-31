#!/usr/bin/env python3
"""Run and validate a paired cached/uncached VCFcache publication pilot."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
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
    "/mnt/data/vcfcache_data/caches/gnomad_v4.1_GRCh38_joint_af001/"
    "cache/vep115.2_everything"
)
DEFAULT_PARAMS = DEFAULT_CACHE / "params.snapshot.yaml"
VCFCACHE_CMD = REPO_ROOT / ".venv/bin/vcfcache"
TIME_CMD = Path("/usr/bin/time")
MODES = ("uncached", "cached")


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


def preflight(config: PilotConfig, *, require_clean: bool = True) -> dict[str, object]:
    """Validate the exact tools and immutable inputs used by the pilot."""
    required_paths = (
        config.input_vcf,
        Path(f"{config.input_vcf}.tbi"),
        config.cache_dir / "annotation.yaml",
        config.cache_dir / "vcfcache_annotated.bcf",
        config.cache_dir / "vcfcache_annotated.bcf.csi",
        config.params_file,
        VCFCACHE_CMD,
        TIME_CMD,
    )
    missing = [str(path) for path in required_paths if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing pilot prerequisites: {missing}")
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
    if "annotation tool: 115.2  ✓" not in requirements:
        raise RuntimeError("VEP 115.2 requirement check did not pass")
    if "bcftools:" not in requirements or "✓" not in requirements:
        raise RuntimeError("bcftools requirement check did not pass")

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
        "input_vcf": str(config.input_vcf),
        "input_records": input_records,
        "input_bytes": config.input_vcf.stat().st_size,
        "input_sha256": sha256sum(config.input_vcf),
        "cache_dir": str(config.cache_dir),
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
        "output_records": output_records,
        "stats_dir": str(stats_dir),
        "variant_counts": variant_counts,
        "cache_hit_rate": hit_rate,
        "command": command,
    }
    write_json_atomic(run_dir / "metrics.json", metrics)
    print(f"{mode.capitalize()} pilot completed in {wall_seconds:.1f}s")
    return metrics


def _query_process(path: Path) -> subprocess.Popen[str]:
    query_format = (
        r"%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/AC\t%INFO/AN"
        r"[\t%GT]\t%INFO/CSQ\n"
    )
    return subprocess.Popen(
        ["bcftools", "query", "--format", query_format, str(path)],
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


def _canonical_record(fields: list[str]) -> tuple[tuple[str, ...], str]:
    return tuple(fields[:8]), ",".join(sorted(fields[8].split(",")))


def semantic_compare(
    cached_bcf: Path, uncached_bcf: Path, *, mismatch_limit: int = 20
) -> dict[str, object]:
    """Compare semantic records, ignoring within-locus and CSQ item order."""
    cached = _query_process(cached_bcf)
    uncached = _query_process(uncached_bcf)
    assert cached.stdout is not None
    assert uncached.stdout is not None
    sentinel = object()
    records = 0
    key_mismatches = 0
    annotation_mismatches = 0
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
            cached_key, cached_canonical = _canonical_record(cached_record)
            uncached_key, uncached_canonical = _canonical_record(uncached_record)
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
                if cached_canonical == uncached_canonical:
                    annotation_order_only += 1
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
    cached_header = _csq_header(cached_bcf)
    uncached_header = _csq_header(uncached_bcf)
    return {
        "semantic_pass": (
            key_mismatches == 0
            and annotation_mismatches == 0
            and cached_header == uncached_header
            and cached_digest.hexdigest() == uncached_digest.hexdigest()
        ),
        "records_compared": records,
        "key_mismatches": key_mismatches,
        "annotation_mismatches": annotation_mismatches,
        "annotation_order_only": annotation_order_only,
        "record_order_only_loci": record_order_only,
        "csq_headers_equal": cached_header == uncached_header,
        "cached_semantic_sha256": cached_digest.hexdigest(),
        "uncached_semantic_sha256": uncached_digest.hexdigest(),
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
