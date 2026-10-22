#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-only
# Copyright (c) 2026 Julius Müller
"""Local synthetic engineering benchmark; never a biological/publication cohort.

Generate data outside timing, run one CPU-affinity-matched cell per condition,
then compare complete records outside timing. Existing outputs are not replaced.
"""

import argparse
import csv
import hashlib
import json
import os
import platform
import resource
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path

import pysam
import yaml

from benchmarks.run_pilot import strict_semantic_compare
from vcfcache.database.fused import native_binary
from vcfcache.database.workflow_manager import WorkflowManager


def checked(*command):
    """Run an inexpensive diagnostic or data-preparation command."""
    return subprocess.run(
        [str(value) for value in command], capture_output=True, text=True, check=True
    ).stdout.strip()


def generate(directory, records, hit_rate, payload):
    """Build sorted WGS-shaped records, a partial cache and a mock recipe."""
    directory.mkdir(parents=True)
    header = pysam.VariantHeader()
    chromosomes = [str(value) for value in (2, 1, *range(3, 23))]
    for chromosome in chromosomes:
        header.contigs.add(chromosome, length=300000000)
    header.info.add("CSQ", ".", "String", "Synthetic consequences")
    header.info.add("AF", "A", "Float", "Allele frequency")
    header.info.add("AC", "A", "Integer", "Allele count")
    header.info.add("AN", 1, "Integer", "Allele number")
    for number in range(10):
        header.info.add(f"DB{number}", 1, "String", "Synthetic supplementary lookup")
    header.formats.add("GT", 1, "String", "Genotype")
    header.formats.add("DP", 1, "Integer", "Depth")
    header.add_sample("engineering-sample")
    with (
        pysam.VariantFile(str(directory / "sample.bcf"), "wb", header=header) as source,
        pysam.VariantFile(str(directory / "cache.bcf"), "wb", header=header) as cache,
        pysam.VariantFile(str(directory / "truth.bcf"), "wb", header=header) as truth,
    ):
        for index in range(records):
            chromosome = chromosomes[min(21, index * 22 // records)]
            position = index * 10
            alleles = ("AT", "A") if index % 50 == 0 else ("A", "C")
            record = source.new_record(
                contig=chromosome,
                start=position,
                alleles=alleles,
                qual=50,
                filter="PASS",
            )
            record.info["AF"], record.info["AC"], record.info["AN"] = (0.5,), (1,), 2
            record.samples[0]["GT"] = (0, 1)
            record.samples[0]["DP"] = 30 + index % 20
            source.write(record)
            token = hashlib.blake2b(str(index).encode(), digest_size=32).hexdigest()
            record.info["CSQ"] = (
                "C|missense_variant|ENST_"
                + str(index % 20000)
                + "|"
                + (token * ((payload + 63) // 64))[:payload]
            )
            for number in range(10):
                record.info[f"DB{number}"] = token[number : number + 8]
            truth.write(record)
            if index % 1000 < round(hit_rate * 1000):
                cache.write(record)
    for filename in ("sample.bcf", "cache.bcf", "truth.bcf"):
        checked("bcftools", "index", directory / filename)
    recipe = {
        "annotation_cmd": f"bcftools annotate -a {directory / 'truth.bcf'} ${{INPUT_BCF}} -c INFO -Ob -W --threads ${{params.threads}} -o ${{OUTPUT_BCF}}",
        "must_contain_info_tag": "CSQ",
        "required_tool_version": "synthetic",
        "genome_build": "GRCh38",
    }
    (directory / "annotation.yaml").write_text(yaml.safe_dump(recipe))


def cell(args):
    """Run only the measured workflow; the parent performs the equality audit."""
    os.sched_setaffinity(0, set(args.affinity))
    output = args.cell_output
    output.mkdir(parents=True)
    params = output / "params.yaml"
    params.write_text(
        yaml.safe_dump(
            {
                "bcftools_cmd": "bcftools",
                "annotation_tool_cmd": "bcftools",
                "genome_build": "GRCh38",
                "threads": len(args.affinity),
                "temp_dir": str(output),
            }
        )
    )
    manager = WorkflowManager(
        input_file=args.data / "sample.bcf",
        output_dir=output,
        name="synthetic",
        params_file=params,
        anno_config_file=args.data / "annotation.yaml",
        output_file=output / "output.bcf",
    )
    start_usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = time.perf_counter()
    manager.run("annotate", db_bcf=args.data / "cache.bcf", engine=args.engine)
    wall = time.perf_counter() - started
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    result = {
        "engine": args.engine,
        "configured_cpus": len(args.affinity),
        "affinity": args.affinity,
        "wall_seconds": wall,
        "cpu_seconds": usage.ru_utime
        + usage.ru_stime
        - start_usage.ru_utime
        - start_usage.ru_stime,
        "max_child_process_rss_kib": usage.ru_maxrss,
        "memory_note": "largest child-process RSS, NOT simultaneous aggregate peak",
        "artifact_bytes": sum(
            path.stat().st_size for path in output.rglob("*") if path.is_file()
        ),
        "stats": manager.last_run_stats,
    }
    (output / "measurement.json").write_text(json.dumps(result, indent=2) + "\n")


def main():
    """Generate fresh conditions and execute paired measurements serially."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--work-dir", type=Path)
    parser.add_argument("--records", type=int, default=1000000)
    parser.add_argument("--hit-rates", type=float, nargs="+", default=[0.8, 0.9])
    parser.add_argument("--payloads", type=int, nargs="+", default=[32, 1024])
    parser.add_argument("--cpus", type=int, nargs="+", default=[1, 6, 12])
    parser.add_argument("--cell", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--engine", choices=["legacy", "fused"])
    parser.add_argument("--affinity", type=int, nargs="+")
    parser.add_argument("--data", type=Path)
    parser.add_argument("--cell-output", type=Path)
    args = parser.parse_args()
    if args.cell:
        cell(args)
        return
    if args.output is None:
        parser.error("--output is required")
    if args.records < 22 or any(not 0 < rate < 1 for rate in args.hit_rates):
        parser.error("use at least 22 records and hit rates strictly between 0 and 1")
    if any(value < 1 for value in args.cpus + args.payloads):
        parser.error("CPUs and payload sizes must be positive")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    work = (
        args.work_dir.resolve()
        if args.work_dir
        else Path(tempfile.mkdtemp(prefix="vcfcache-fused-bench-"))
    )
    work.mkdir(parents=True, exist_ok=True)
    available = sorted(os.sched_getaffinity(0))
    cpus = sorted({min(value, len(available)) for value in args.cpus})
    binary = Path(native_binary())
    manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "synthetic": True,
        "records": args.records,
        "payloads": args.payloads,
        "hit_rates": args.hit_rates,
        "cpu_limits": cpus,
        "available_affinity": available,
        "work_directory": str(work),
        "platform": platform.platform(),
        "python": sys.version,
        "bcftools": checked("bcftools", "--version"),
        "native_version": checked(binary, "--version"),
        "native_sha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
        "commit": checked("git", "rev-parse", "HEAD"),
        "tracked_diff_sha256": hashlib.sha256(
            checked("git", "diff", "HEAD").encode()
        ).hexdigest(),
        "source_sha256": {
            str(path): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in [
                Path(__file__),
                Path("vcfcache/database/fused.py"),
                Path("native/src/main.rs"),
                Path("native/Cargo.lock"),
            ]
        },
        "measurement": "one run per condition; no technical repeats; generation and complete-record comparisons untimed",
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    rows = []
    for payload in args.payloads:
        for hit_rate in args.hit_rates:
            name = f"payload-{payload}-hit-{hit_rate:g}"
            data = work / name / "data"
            generate(data, args.records, hit_rate, payload)
            for count in cpus:
                results = {}
                for engine in ("legacy", "fused"):
                    label = f"{name}-cpu-{count}-{engine}"
                    destination = work / label
                    with (
                        (output / f"{label}.stdout.txt").open("wb") as stdout,
                        (output / f"{label}.stderr.txt").open("wb") as stderr,
                    ):
                        subprocess.run(
                            [
                                sys.executable,
                                "-m",
                                "benchmarks.fused_engine",
                                "--cell",
                                "--engine",
                                engine,
                                "--data",
                                str(data),
                                "--cell-output",
                                str(destination),
                                "--affinity",
                                *[str(value) for value in available[:count]],
                            ],
                            stdout=stdout,
                            stderr=stderr,
                            check=True,
                        )
                    result = json.loads((destination / "measurement.json").read_text())
                    result.update(payload_bytes=payload, cache_hit_rate=hit_rate)
                    results[engine] = result
                    (output / f"{label}.json").write_text(
                        json.dumps(result, indent=2) + "\n"
                    )
                audit = strict_semantic_compare(
                    work / f"{name}-cpu-{count}-fused/output.bcf",
                    work / f"{name}-cpu-{count}-legacy/output.bcf",
                )
                (output / f"{name}-cpu-{count}-equality.json").write_text(
                    json.dumps(audit, indent=2) + "\n"
                )
                if not audit["semantic_pass"]:
                    raise RuntimeError(
                        f"complete-record comparison failed: {name} cpu={count}"
                    )
                row = {
                    "records": args.records,
                    "payload_bytes": payload,
                    "cache_hit_rate": hit_rate,
                    "cpus": count,
                    "legacy_seconds": results["legacy"]["wall_seconds"],
                    "fused_seconds": results["fused"]["wall_seconds"],
                    "semantic_pass": True,
                }
                row["speedup"] = row["legacy_seconds"] / row["fused_seconds"]
                row["reduction_percent"] = 100 * (1 - 1 / row["speedup"])
                rows.append(row)
                with (output / "summary.tsv").open("w") as stream:
                    writer = csv.DictWriter(
                        stream,
                        fieldnames=list(row),
                        delimiter="\t",
                        lineterminator="\n",
                    )
                    writer.writeheader()
                    writer.writerows(rows)
                print(json.dumps(row), flush=True)
    print(f"Synthetic results: {output}\nRaw generated data retained at: {work}")


if __name__ == "__main__":
    main()
