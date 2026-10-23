#!/usr/bin/env python3
"""Full-CLI adaptive-routing benchmark; preparation/audits are outside timing.

Use a fresh output directory. Engineering controls are cache-derived, never
biological hit-rate evidence. A Slurm allocation can give each timed CLI its
own step with --slurm-steps, keeping condition memory peaks independent.
"""

import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

import pysam
import yaml

from benchmarks.run_pilot import strict_semantic_compare


def controls(cache, destination):
    """Draw exact-size, exact-span all-hit engineering controls from a cache."""
    destination.mkdir(parents=True, exist_ok=False)
    inputs = []
    with pysam.VariantFile(str(cache / "vcfcache_annotated.bcf")) as reader:
        names = [
            name
            for name in reader.header.contigs
            if name.removeprefix("chr").isdigit()
            and 1 <= int(name.removeprefix("chr")) <= 22
        ]
        specs = [(f"size-{n}", n, 1 if n == 1 else 22) for n in (1, 250, 2000, 75000)]
        specs += [(f"span-{k}", 250, k) for k in (1, 2, 4, 8)]
        for name, n, k in specs:
            path = destination / f"{name}.bcf"
            with pysam.VariantFile(str(path), "wb", header=reader.header) as output:
                for i, contig in enumerate(names[:k]):
                    need = n // k + (i < n % k)
                    for j, record in enumerate(reader.fetch(contig, 20_000_000)):
                        if j == need:
                            break
                        record.info.clear()
                        output.write(record)
                    else:
                        raise RuntimeError(f"not enough cache records on {contig}")
            subprocess.run(["bcftools", "index", str(path)], check=True)
            assert (
                int(subprocess.check_output(["bcftools", "index", "-n", str(path)]))
                == n
            )
            inputs.append(path)
    return inputs


def digest(path):
    """Hash a source or executable without loading it all into memory."""
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def main():
    """Execute isolated timed invocations and untimed complete-record audits."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--params", type=Path)
    parser.add_argument("--input", type=Path, nargs="+")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cpus", type=int, default=8)
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--warmup", action="store_true")
    parser.add_argument("--generate-controls", action="store_true")
    parser.add_argument("--control-names", nargs="+")
    parser.add_argument(
        "--arms",
        nargs="+",
        choices=["legacy", "original-fused", "auto", "regions", "contigs"],
        default=["legacy", "auto"],
    )
    parser.add_argument("--original-root", type=Path)
    parser.add_argument("--slurm-steps", action="store_true")
    args = parser.parse_args()
    if args.cpus < 1 or args.repeats < 1:
        parser.error("CPUs and repetitions must be positive")
    if "original-fused" in args.arms and not args.original_root:
        parser.error("original-fused needs --original-root")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    inputs = (
        controls(args.cache, output / "inputs")
        if args.generate_controls
        else args.input
    )
    if not inputs:
        parser.error("provide inputs or --generate-controls")
    if args.control_names:
        inputs = [path for path in inputs if path.stem in args.control_names]
        if not inputs:
            parser.error("no controls matched --control-names")
    params = yaml.safe_load(
        (args.params or args.cache / "params.snapshot.yaml").read_text()
    )
    params["threads"] = args.cpus
    for name in ("vep_forks", "fastvep_threads"):
        if name in params:
            params[name] = min(int(params[name]), args.cpus)
    params_path = output / "params.yaml"
    params_path.write_text(yaml.safe_dump(params))
    root = Path(__file__).resolve().parents[1]
    manifest = dict(
        platform=platform.platform(),
        python=sys.version,
        cpus=args.cpus,
        repeats=args.repeats,
        warmup=args.warmup,
        inputs=[str(path) for path in inputs],
        arms=args.arms,
        cache=str(args.cache),
        source_sha256={
            str(path.relative_to(root)): digest(path)
            for path in (
                root / "vcfcache/database/fused.py",
                root / "vcfcache/database/adaptive.py",
                root / "native/src/main.rs",
                root / "native/src/adaptive.rs",
            )
        },
        native_sha256=digest(root / "native/target/release/vcfcache-engine"),
        bcftools=subprocess.check_output(["bcftools", "--version"], text=True),
        original_root=str(args.original_root),
        slurm_steps=args.slurm_steps,
        note="cache-derived controls are engineering inputs, not cohort coverage results",
    )
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    measurements = []
    for source in inputs:
        reference = None
        for arm in args.arms:
            for rep in range(0 if args.warmup else 1, args.repeats + 1):
                directory = output / f"{source.stem}-{arm}-{rep}"
                directory.mkdir()
                result = directory / "output.bcf"
                command = [
                    sys.executable,
                    "-c",
                    "from vcfcache.cli import main; main()",
                    "annotate",
                    "-a",
                    str(args.cache),
                    "-i",
                    str(source),
                    "-o",
                    str(result),
                    "-y",
                    str(params_path),
                    "--stats-dir",
                    str(directory),
                    "--statistics",
                    "light",
                    "--engine",
                    "legacy" if arm == "legacy" else "fused",
                ]
                if arm not in {"legacy", "original-fused"}:
                    command += ["--engine-lookup", arm]
                environment = dict(os.environ)
                deployment = args.original_root if arm == "original-fused" else root
                environment["PYTHONPATH"] = str(deployment)
                environment["VCFCACHE_ENGINE"] = str(
                    deployment / "native/target/release/vcfcache-engine"
                )
                command = [sys.executable, str(root / "benchmarks/adaptive_cell.py"),
                           "--metrics", str(directory / "resources.json"), "--", *command]
                if args.slurm_steps:
                    command = [
                        "srun",
                        "--exclusive",
                        "--exact",
                        "--cpu-bind=cores",
                        "-n1",
                        f"-c{args.cpus}",
                        *command,
                    ]
                else:
                    available = sorted(os.sched_getaffinity(0))
                    if len(available) < args.cpus:
                        parser.error("requested CPUs exceed current affinity")
                    command = [
                        "taskset",
                        "-c",
                        ",".join(map(str, available[: args.cpus])),
                        *command,
                    ]
                start = time.perf_counter()
                with (directory / "run.log").open("w") as log:
                    subprocess.run(
                        command,
                        env=environment,
                        stdout=log,
                        stderr=subprocess.STDOUT,
                        check=True,
                    )
                wall = time.perf_counter() - start
                if rep == 0:
                    continue
                if reference is None:
                    reference = result
                comparison = strict_semantic_compare(result, reference)
                (directory / "comparison.json").write_text(
                    json.dumps(comparison, indent=2) + "\n"
                )
                row = dict(
                    input=str(source),
                    arm=arm,
                    rep=rep,
                    wall_seconds=wall,
                    semantic_pass=comparison["semantic_pass"],
                )
                measurements.append(row)
                (output / "measurements.json").write_text(
                    json.dumps(measurements, indent=2) + "\n"
                )
                print(json.dumps(row), flush=True)
                if not row["semantic_pass"]:
                    raise RuntimeError(
                        f"complete-record comparison failed: {directory}"
                    )


if __name__ == "__main__":
    main()
