#!/usr/bin/env python3
"""Micro-profile the VCFcache 100%-hit fastVEP path on one WGS."""

# This is an operational benchmark runner, not an installed public API.
# ruff: noqa: D103,E501

from __future__ import annotations

import argparse
import json
import statistics
import subprocess
import time
from pathlib import Path
from typing import Sequence


def parse_time(path: Path) -> dict[str, float | int]:
    values: dict[str, str] = {}
    for line in path.read_text().splitlines():
        if ": " in line:
            key, value = line.strip().rsplit(": ", maxsplit=1)
            values[key] = value
    return {
        "user_seconds": float(values["User time (seconds)"]),
        "system_seconds": float(values["System time (seconds)"]),
        "max_rss_kib": int(values["Maximum resident set size (kbytes)"]),
        "filesystem_inputs": int(values["File system inputs"]),
        "filesystem_outputs": int(values["File system outputs"]),
    }


def timed(
    root: Path,
    name: str,
    command: Sequence[str | Path],
    *,
    output: Path | None = None,
) -> dict[str, object]:
    metrics_path = root / f"{name}.json"
    if metrics_path.exists():
        return json.loads(metrics_path.read_text())
    time_path = root / f"{name}.time.txt"
    log_path = root / f"{name}.log"
    started = time.monotonic()
    with log_path.open("w") as log:
        completed = subprocess.run(
            ["/usr/bin/time", "--verbose", "--output", time_path, *map(str, command)],
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
            text=True,
        )
    result: dict[str, object] = {
        "command": [str(value) for value in command],
        "returncode": completed.returncode,
        "wall_seconds": time.monotonic() - started,
        **parse_time(time_path),
    }
    if output and output.exists():
        result["output_bytes"] = output.stat().st_size
    metrics_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    if completed.returncode:
        raise RuntimeError(f"{name} failed; see {log_path}")
    return result


def startup_profile(
    command: Sequence[str | Path], repeats: int = 25
) -> dict[str, object]:
    observations = []
    for _ in range(repeats):
        started = time.perf_counter()
        subprocess.run(
            [str(value) for value in command],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True,
        )
        observations.append(time.perf_counter() - started)
    return {
        "command": [str(value) for value in command],
        "repeats": repeats,
        "median_seconds": statistics.median(observations),
        "min_seconds": min(observations),
        "max_seconds": max(observations),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pilot-root",
        type=Path,
        default=Path("/mnt/data/vcfcache_benchmarks/fastvep_pilot"),
    )
    parser.add_argument("--threads", type=int, default=16)
    args = parser.parse_args()

    pilot = args.pilot_root.resolve()
    root = pilot / "profile_100_hit"
    root.mkdir(parents=True, exist_ok=True)
    environment = json.loads((pilot / "environment.json").read_text())
    input_vcf = Path(environment["input_vcf"])
    cache = pilot / "cache/full/vcfcache_annotated.bcf"
    fastvep = pilot / "tools/fastVEP/target/release/fastvep"
    filtered = root / "input.filtered.bcf"
    lookup = root / "lookup.bcf"
    final = root / "final.bcf"
    missing = root / "missing.bcf"
    threads = str(args.threads)

    stages: dict[str, dict[str, object]] = {}
    stages["input_parse_filter"] = timed(
        root,
        "input_parse_filter",
        ["bcftools", "view", "-e", 'ALT="*"', "-Ou", "-o", "/dev/null", input_vcf],
    )
    stages["input_filter_compress"] = timed(
        root,
        "input_filter_compress",
        [
            "bcftools",
            "view",
            "-e",
            'ALT="*"',
            "-Ob",
            "-o",
            filtered,
            "--threads",
            threads,
            input_vcf,
        ],
        output=filtered,
    )
    stages["input_index"] = timed(
        root, "input_index", ["bcftools", "index", "--force", filtered]
    )
    stages["decode_control"] = timed(
        root,
        "decode_control",
        ["bcftools", "view", "-Ou", "-o", "/dev/null", filtered],
    )
    stages["cache_lookup_no_compression"] = timed(
        root,
        "cache_lookup_no_compression",
        [
            "bcftools",
            "annotate",
            "-a",
            cache,
            "-c",
            "INFO",
            "-Ou",
            "-o",
            "/dev/null",
            filtered,
        ],
    )
    stages["cache_lookup_compress"] = timed(
        root,
        "cache_lookup_compress",
        [
            "bcftools",
            "annotate",
            "-a",
            cache,
            "-c",
            "INFO",
            "-Ob",
            "-o",
            lookup,
            "--threads",
            threads,
            filtered,
        ],
        output=lookup,
    )
    stages["cache_lookup_index"] = timed(
        root, "cache_lookup_index", ["bcftools", "index", "--force", lookup]
    )
    stages["miss_filter_no_compression"] = timed(
        root,
        "miss_filter_no_compression",
        ["bcftools", "filter", "-i", 'INFO/CSQ==""', "-Ou", "-o", "/dev/null", lookup],
    )
    stages["miss_filter_compress"] = timed(
        root,
        "miss_filter_compress",
        [
            "bcftools",
            "filter",
            "-i",
            'INFO/CSQ==""',
            "-Ob",
            "-o",
            missing,
            "--threads",
            threads,
            lookup,
        ],
        output=missing,
    )
    stages["miss_index"] = timed(
        root, "miss_index", ["bcftools", "index", "--force", missing]
    )
    stages["output_filter_no_compression"] = timed(
        root,
        "output_filter_no_compression",
        ["bcftools", "view", "-i", 'INFO/CSQ!=""', "-Ou", "-o", "/dev/null", lookup],
    )
    stages["output_filter_compress"] = timed(
        root,
        "output_filter_compress",
        [
            "bcftools",
            "view",
            "-i",
            'INFO/CSQ!=""',
            "-Ob",
            "-o",
            final,
            "--threads",
            threads,
            lookup,
        ],
        output=final,
    )
    stages["output_index"] = timed(
        root, "output_index", ["bcftools", "index", "--force", final]
    )

    startup = {
        "shell_true": startup_profile(["/bin/true"]),
        "bcftools": startup_profile(["bcftools", "--version"]),
        "fastvep": startup_profile([fastvep, "--version"]),
    }
    estimates = {
        "cache_lookup_increment_over_decode_seconds": max(
            0.0,
            float(stages["cache_lookup_no_compression"]["wall_seconds"])
            - float(stages["decode_control"]["wall_seconds"]),
        ),
        "lookup_compression_and_write_increment_seconds": max(
            0.0,
            float(stages["cache_lookup_compress"]["wall_seconds"])
            - float(stages["cache_lookup_no_compression"]["wall_seconds"]),
        ),
        "final_compression_and_write_increment_seconds": max(
            0.0,
            float(stages["output_filter_compress"]["wall_seconds"])
            - float(stages["output_filter_no_compression"]["wall_seconds"]),
        ),
    }
    summary = {
        "input": str(input_vcf),
        "cache": str(cache),
        "threads": args.threads,
        "method": (
            "Single exploratory run per full-WGS stage. Null-output controls separate "
            "decode/filter/lookup from compressed-file writing; deltas are diagnostic and "
            "not strictly additive. Startup uses 25 repetitions."
        ),
        "stages": stages,
        "startup": startup,
        "derived_estimates": estimates,
    }
    (root / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    with (root / "summary.tsv").open("w") as handle:
        handle.write(
            "stage\twall_seconds\tuser_seconds\tsystem_seconds\tmax_rss_kib\toutput_bytes\n"
        )
        for name, metrics in stages.items():
            handle.write(
                f"{name}\t{metrics['wall_seconds']}\t{metrics['user_seconds']}\t"
                f"{metrics['system_seconds']}\t{metrics['max_rss_kib']}\t"
                f"{metrics.get('output_bytes', '')}\n"
            )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
