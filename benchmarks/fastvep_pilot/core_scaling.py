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
import shutil
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
    asset_subdir: str = ""
    input_override: Path | None = None
    repeats: int = 3
    repeat_control: int = 10
    run_name: str = "core_scaling_affinity"

    @property
    def root(self) -> Path:
        return self.pilot / self.run_name

    @property
    def assets(self) -> Path:
        return self.pilot / self.asset_subdir if self.asset_subdir else self.pilot

    @property
    def allowed_cpus(self) -> tuple[int, ...]:
        return tuple(sorted(os.sched_getaffinity(0)))

    def repetitions(self, cores: int) -> range:
        return (
            range(1, self.repeats + 1) if cores == self.repeat_control else range(1, 2)
        )

    @property
    def input_vcf(self) -> Path:
        if self.input_override is not None:
            return self.input_override
        environment = json.loads((self.pilot / "environment.json").read_text())
        return Path(environment["input_vcf"])

    @property
    def python(self) -> Path:
        return self.repo / ".venv/bin/python"

    @property
    def placeholder(self) -> Path:
        return self.assets / "cache/placeholder"

    @property
    def cache_90(self) -> Path:
        controlled = self.root / "cache/hit-090"
        if controlled.joinpath("vcfcache_annotated.bcf.csi").exists():
            return controlled
        return self.assets / "cache/hit-090"


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


def prepare(config: Config, *, require_cache: bool = True) -> None:
    for path in (config.root, config.root / "tmp", config.root / "config"):
        path.mkdir(parents=True, exist_ok=True)
    required = (
        config.source / "vcfcache/__init__.py",
        config.placeholder / "vcfcache_annotated.bcf.csi",
        config.assets / "config/params.yaml",
    )
    if require_cache:
        required += (config.cache_90 / "vcfcache_annotated.bcf.csi",)
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(
            "Missing scaling prerequisite(s): " + ", ".join(missing)
        )

    base_params = yaml.safe_load((config.assets / "config/params.yaml").read_text())
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
        "process_allowed_cpus": config.allowed_cpus,
        "core_counts": config.core_counts,
        "input": str(config.input_vcf),
        "cache": str(config.cache_90),
        "asset_root": str(config.assets),
        "statistics": "light",
        "repeat_control_thread_limit": config.repeat_control,
        "repeat_control_repeats": config.repeats,
        "cpu_quota": "enforced with taskset process affinity",
        "execution": "serial",
    }
    (config.root / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )


def write_controlled_cache_metadata(
    cache: Path, input_vcf: Path, records: int
) -> None:
    """Write metadata required for a controlled cache to validate."""
    (cache / "cache.json").write_text(
        json.dumps(
            {
                "construction": "POS modulo 100 < 90",
                "input": str(input_vcf),
                "records": records,
                "target_hit_rate": 0.9,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )


def copy_required_cache_structure(source_root: Path, destination_root: Path) -> None:
    """Copy immutable root structure required by frozen VCFcache snapshots."""
    for name in ("blueprint", "workflow"):
        source = source_root / name
        if not source.is_dir():
            raise FileNotFoundError(f"Missing cache structure source: {source}")
        shutil.copytree(source, destination_root / name, dirs_exist_ok=True)


def build_controlled_cache(config: Config) -> None:
    """Build a deterministic 90%-hit cache from an independent input genome."""
    prepare(config, require_cache=False)
    cores = max(config.core_counts)
    direct = run_cell(config, cores, "direct", repeat=1)
    direct_output = Path(str(direct["output"]))
    cache = config.root / "cache/hit-090"
    output = cache / "vcfcache_annotated.bcf"
    index = Path(f"{output}.csi")
    cache.mkdir(parents=True, exist_ok=True)
    copy_required_cache_structure(config.assets, config.root)
    if not index.exists():
        shutil.copy2(
            config.assets / "cache/hit-090/annotation.yaml", cache / "annotation.yaml"
        )
        shutil.copy2(
            config.assets / "config/params.yaml", cache / "params.snapshot.yaml"
        )
        keep = "INFO/CSQ,INFO/FV_CLINVAR,INFO/FV_GNOMAD,INFO/FV_1KG,INFO/FV_TOPMED"
        shell = (
            f"bcftools view -G -Ou {shlex.quote(str(direct_output))} | "
            f"bcftools annotate -x {shlex.quote('^' + keep)} -Ou | "
            f"bcftools view -i {shlex.quote('POS%100<90')} -Ob "
            f"-o {shlex.quote(str(output))}"
        )
        subprocess.run(["bash", "-euo", "pipefail", "-c", shell], check=True)
        subprocess.run(["bcftools", "index", "--force", str(output)], check=True)
    records = int(
        subprocess.check_output(
            ["bcftools", "index", "-n", str(output)], text=True
        ).strip()
    )
    write_controlled_cache_metadata(cache, config.input_vcf, records)
    provenance = {
        "input": str(config.input_vcf),
        "source_direct_output": str(direct_output),
        "selection_rule": "POS%100<90",
        "purpose": "controlled approximately 90 percent hit cache",
        "annotation_recipe": str(cache / "annotation.yaml"),
    }
    (cache / "controlled_cache_provenance.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n"
    )
    prepare(config)


def command(config: Config, cores: int, mode: str, run_dir: Path) -> list[str | Path]:
    uncached = mode == "direct"
    cache = config.placeholder if uncached else config.cache_90
    cpu_list = ",".join(map(str, config.allowed_cpus[:cores]))
    result: list[str | Path] = [
        "taskset",
        "--cpu-list",
        cpu_list,
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


def run_cell(
    config: Config, cores: int, mode: str, repeat: int = 1
) -> dict[str, object]:
    run_dir = config.root / f"runs/cores-{cores}/repeat-{repeat}/{mode}"
    metrics_path = run_dir / "metrics.json"
    if metrics_path.exists():
        previous = json.loads(metrics_path.read_text())
        if previous.get("returncode") == 0:
            return previous
        failed = config.root / "runs_failed"
        failed.mkdir(parents=True, exist_ok=True)
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        run_dir.rename(failed / f"cores-{cores}_repeat-{repeat}_{mode}_{stamp}")
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
        "configured_thread_limit": cores,
        "enforced_cpu_affinity": list(config.allowed_cpus[:cores]),
        "repeat": repeat,
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
        for repeat in config.repetitions(cores):
            cell = config.root / f"runs/cores-{cores}/repeat-{repeat}"
            direct = json.loads((cell / "direct/metrics.json").read_text())
            cached = json.loads((cell / "cached/metrics.json").read_text())
            equality_path = cell / "equality.json"
            equality = compare_outputs(
                Path(direct["output"]), Path(cached["output"]), equality_path
            )
            if not equality["pass"]:
                raise RuntimeError(f"Output equality failed at {cores} threads")
            for values in (direct, cached):
                wall = float(values["wall_seconds"])
                rows.append(
                    {
                        "configured_thread_limit": cores,
                        "repeat": repeat,
                        "enforced_cpu_affinity": ",".join(
                            map(str, config.allowed_cpus[:cores])
                        ),
                        "mode": values["mode"],
                        "wall_seconds": wall,
                        "workflow_seconds": values["workflow_seconds"],
                        "max_rss_kib": values["max_rss_kib"],
                        "cpu_percent": values["cpu_percent"],
                        "speedup_vs_direct_same_threads": (
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
        for repeat in config.repetitions(cores):
            for mode in ("direct", "cached"):
                print(
                    f"Running {mode} at {cores} configured threads, repeat {repeat}",
                    flush=True,
                )
                run_cell(config, cores, mode, repeat)
            cell = config.root / f"runs/cores-{cores}/repeat-{repeat}"
            if not compare_outputs(
                cell / "direct/output.bcf",
                cell / "cached/output.bcf",
                cell / "equality.json",
            )["pass"]:
                raise RuntimeError(f"Output equality failed at {cores} threads")
    collect(config)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "command", choices=("prepare-cache", "run", "collect", "status")
    )
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
        "--input-vcf",
        type=Path,
        help="Override the pilot input, for example with an independent WGS",
    )
    parser.add_argument("--cores", type=int, nargs="+", default=None)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--repeat-control", type=int, default=10)
    parser.add_argument("--run-name", default="core_scaling_affinity")
    parser.add_argument(
        "--asset-subdir",
        default="",
        help=(
            "Optional directory below --pilot-root containing the cache and base "
            "params, for example 'heavy' for the ten-database stress recipe"
        ),
    )
    args = parser.parse_args()

    available = len(os.sched_getaffinity(0))
    counts = tuple(args.cores or (1, 10, available))
    if any(value < 1 or value > available for value in counts):
        parser.error(f"--cores values must be between 1 and {available}")
    source = args.source or args.pilot_root / "heavy/source"
    if args.repeats < 1 or args.repeat_control not in counts:
        parser.error(
            "--repeats must be positive and --repeat-control must be in --cores"
        )
    config = Config(
        args.pilot_root.resolve(),
        args.repo.resolve(),
        source.resolve(),
        counts,
        args.asset_subdir,
        args.input_vcf.resolve() if args.input_vcf else None,
        args.repeats,
        args.repeat_control,
        args.run_name,
    )
    if args.command == "prepare-cache":
        build_controlled_cache(config)
    elif args.command == "run":
        benchmark(config)
    elif args.command == "collect":
        collect(config)
    else:
        print(
            json.dumps(
                {
                    str(cores): {
                        mode: (
                            config.root
                            / f"runs/cores-{cores}/repeat-1/{mode}/metrics.json"
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
