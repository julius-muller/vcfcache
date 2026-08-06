#!/usr/bin/env python3
"""Benchmark direct fastVEP versus 90%-hit VCFcache across CPU counts."""

# Operational benchmark runner, not an installed public API.
# ruff: noqa: D101,D103,E501

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shlex
import subprocess
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from itertools import zip_longest
from pathlib import Path

import yaml  # type: ignore[import-untyped]


@dataclass(frozen=True)
class Config:
    pilot: Path
    repo: Path
    source: Path
    core_counts: tuple[int, ...]

    @property
    def root(self) -> Path:
        return self.pilot / "core_scaling"

    @property
    def input_vcf(self) -> Path:
        environment = json.loads((self.pilot / "environment.json").read_text())
        return Path(environment["input_vcf"])

    @property
    def python(self) -> Path:
        return self.repo / ".venv/bin/python"

    @property
    def placeholder(self) -> Path:
        return self.pilot / "cache/placeholder"

    @property
    def cache_90(self) -> Path:
        return self.pilot / "cache/hit-090"


def environment(config: Config, cores: int) -> dict[str, str]:
    result = os.environ.copy()
    result.update(
        {
            "PYTHONPATH": str(config.source),
            "RAYON_NUM_THREADS": str(cores),
            "LC_ALL": "C",
            "LANG": "C",
            "TMPDIR": str(config.root / "tmp"),
        }
    )
    return result


def prepare(config: Config) -> None:
    for path in (config.root, config.root / "tmp", config.root / "config"):
        path.mkdir(parents=True, exist_ok=True)
    required = (
        config.source / "vcfcache/__init__.py",
        config.placeholder / "vcfcache_annotated.bcf.csi",
        config.cache_90 / "vcfcache_annotated.bcf.csi",
        config.pilot / "config/params.yaml",
    )
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(
            "Missing scaling prerequisite(s): " + ", ".join(missing)
        )

    base_params = yaml.safe_load((config.pilot / "config/params.yaml").read_text())
    for cores in config.core_counts:
        params = dict(base_params)
        params["threads"] = cores
        params["rayon_threads"] = cores
        params["temp_dir"] = str(config.root / "tmp")
        (config.root / f"config/params.cores-{cores}.yaml").write_text(
            yaml.safe_dump(params, sort_keys=False)
        )
    manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "host_logical_cpus": os.cpu_count(),
        "core_counts": config.core_counts,
        "input": str(config.input_vcf),
        "cache": str(config.cache_90),
        "statistics": "light",
        "repeats": 1,
        "execution": "serial",
    }
    (config.root / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )


def command(config: Config, cores: int, mode: str, run_dir: Path) -> list[str | Path]:
    uncached = mode == "direct"
    cache = config.placeholder if uncached else config.cache_90
    result: list[str | Path] = [
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
        config.root / f"config/params.cores-{cores}.yaml",
        "--skip-split-multiallelic",
        "--force",
    ]
    if uncached:
        result.append("--uncached")
    return result


def parse_resource_usage(path: Path) -> dict[str, float | int]:
    values: dict[str, str] = {}
    for line in path.read_text().splitlines():
        if ": " in line:
            key, value = line.strip().rsplit(": ", maxsplit=1)
            values[key] = value
    return {
        "user_seconds": float(values["User time (seconds)"]),
        "system_seconds": float(values["System time (seconds)"]),
        "cpu_percent": int(values["Percent of CPU this job got"].rstrip("%")),
        "max_rss_kib": int(values["Maximum resident set size (kbytes)"]),
    }


def workflow_seconds(log: Path) -> float | None:
    match = re.findall(
        r"Workflow completed successfully in ([0-9.]+)s", log.read_text()
    )
    return float(match[-1]) if match else None


def run_cell(config: Config, cores: int, mode: str) -> dict[str, object]:
    run_dir = config.root / f"runs/cores-{cores}/{mode}"
    metrics_path = run_dir / "metrics.json"
    if metrics_path.exists():
        previous = json.loads(metrics_path.read_text())
        if previous.get("returncode") == 0:
            return previous
        failed = config.root / "runs_failed"
        failed.mkdir(parents=True, exist_ok=True)
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        run_dir.rename(failed / f"cores-{cores}_{mode}_{stamp}")
    if run_dir.exists():
        raise FileExistsError(f"Incomplete run without metrics: {run_dir}")
    run_dir.mkdir(parents=True)
    args = command(config, cores, mode, run_dir)
    resource_path = run_dir / "resource_usage.txt"
    log_path = run_dir / "command.log"
    started = time.monotonic()
    with log_path.open("w") as log:
        completed = subprocess.run(
            ["/usr/bin/time", "--verbose", "--output", resource_path, *map(str, args)],
            cwd=config.root,
            env=environment(config, cores),
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
            text=True,
        )
    result: dict[str, object] = {
        "cores": cores,
        "mode": mode,
        "command": shlex.join(map(str, args)),
        "returncode": completed.returncode,
        "wall_seconds": time.monotonic() - started,
        "workflow_seconds": workflow_seconds(log_path),
        "output": str(run_dir / "output.bcf"),
        **parse_resource_usage(resource_path),
    }
    metrics_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    if completed.returncode:
        raise RuntimeError(f"Scaling cell failed: {cores} cores {mode}; see {log_path}")
    return result


def canonical_info(value: str) -> str:
    if value in ("", "."):
        return value
    fields = []
    for item in value.split(";"):
        if item.startswith("CSQ="):
            item = "CSQ=" + ",".join(sorted(item[4:].split(",")))
        fields.append(item)
    return ";".join(sorted(fields))


def compare_outputs(direct: Path, cached: Path, destination: Path) -> dict[str, object]:
    if destination.exists():
        return json.loads(destination.read_text())
    processes = [
        subprocess.Popen(
            ["bcftools", "view", "-H", str(path)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        for path in (direct, cached)
    ]
    assert processes[0].stdout is not None and processes[1].stdout is not None
    digests = [hashlib.sha256(), hashlib.sha256()]
    mismatches = 0
    records = 0
    for pair in zip_longest(processes[0].stdout, processes[1].stdout):
        records += 1
        normalized: list[str | None] = []
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
    stderr = [process.communicate()[1] for process in processes]
    returncodes = [process.returncode for process in processes]
    result = {
        "records": records,
        "mismatches": mismatches,
        "sha256": [digest.hexdigest() for digest in digests],
        "returncodes": returncodes,
        "stderr": stderr,
    }
    result["pass"] = (
        mismatches == 0
        and result["sha256"][0] == result["sha256"][1]  # type: ignore[index]
        and returncodes == [0, 0]
    )
    destination.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def collect(config: Config) -> dict[str, object]:
    rows = []
    for cores in config.core_counts:
        direct = json.loads(
            (config.root / f"runs/cores-{cores}/direct/metrics.json").read_text()
        )
        cached = json.loads(
            (config.root / f"runs/cores-{cores}/cached/metrics.json").read_text()
        )
        equality_path = config.root / f"runs/cores-{cores}/equality.json"
        equality = compare_outputs(
            Path(direct["output"]), Path(cached["output"]), equality_path
        )
        if not equality["pass"]:
            raise RuntimeError(f"Output equality failed at {cores} cores")
        for values in (direct, cached):
            wall = float(values["wall_seconds"])
            rows.append(
                {
                    "cores": cores,
                    "mode": values["mode"],
                    "wall_seconds": wall,
                    "workflow_seconds": values["workflow_seconds"],
                    "max_rss_kib": values["max_rss_kib"],
                    "cpu_percent": values["cpu_percent"],
                    "speedup_vs_direct_same_cores": (
                        float(direct["wall_seconds"]) / wall
                    ),
                    "semantic_pass": True,
                }
            )
    summary = {"exploratory": True, "rows": rows}
    (config.root / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    columns = tuple(rows[0])
    with (config.root / "summary.tsv").open("w") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write("\t".join(str(row[column]) for column in columns) + "\n")
    print(json.dumps(summary, indent=2))
    return summary


def benchmark(config: Config) -> None:
    prepare(config)
    for cores in config.core_counts:
        for mode in ("direct", "cached"):
            print(f"Running {mode} at {cores} core(s)", flush=True)
            run_cell(config, cores, mode)
        print(f"Validating outputs at {cores} core(s)", flush=True)
        direct = config.root / f"runs/cores-{cores}/direct/output.bcf"
        cached = config.root / f"runs/cores-{cores}/cached/output.bcf"
        equality = config.root / f"runs/cores-{cores}/equality.json"
        if not compare_outputs(direct, cached, equality)["pass"]:
            raise RuntimeError(f"Output equality failed at {cores} cores")
    collect(config)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("run", "collect", "status"))
    parser.add_argument(
        "--pilot-root",
        type=Path,
        default=Path("/mnt/data/vcfcache_benchmarks/fastvep_pilot"),
    )
    parser.add_argument(
        "--repo", type=Path, default=Path("/mnt/data/home/appuser-projects/vcfcache")
    )
    parser.add_argument("--source", type=Path)
    parser.add_argument("--cores", type=int, nargs="+", default=None)
    args = parser.parse_args()

    available = os.cpu_count() or 1
    counts = tuple(args.cores or (1, 10, available))
    if any(value < 1 or value > available for value in counts):
        parser.error(f"--cores values must be between 1 and {available}")
    source = args.source or args.pilot_root / "heavy/source"
    config = Config(
        args.pilot_root.resolve(), args.repo.resolve(), source.resolve(), counts
    )
    if args.command == "run":
        benchmark(config)
    elif args.command == "collect":
        collect(config)
    else:
        print(
            json.dumps(
                {
                    str(cores): {
                        mode: (
                            config.root / f"runs/cores-{cores}/{mode}/metrics.json"
                        ).exists()
                        for mode in ("direct", "cached")
                    }
                    for cores in counts
                },
                indent=2,
            )
        )


if __name__ == "__main__":
    main()
