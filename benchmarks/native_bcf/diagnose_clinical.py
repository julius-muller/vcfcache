#!/usr/bin/env python3
"""Recheck two discrepant clinical loci in isolation; never ignore fields."""

import argparse
import json
import subprocess
import sys
from pathlib import Path

import pysam

from benchmarks.adaptive_routing import claim_directory
from benchmarks.run_pilot import strict_semantic_compare


def subset(source, destination, loci):
    """Copy complete records at selected loci in the source header order."""
    with (
        pysam.VariantFile(str(source)) as inp,
        pysam.VariantFile(str(destination), "wb", header=inp.header) as out,
    ):
        for contig in inp.header.contigs:
            for name, pos in sorted(loci):
                if contig == name:
                    for record in inp.fetch(name, pos - 1, pos):
                        if record.pos == pos:
                            out.write(record)
    subprocess.run(["bcftools", "index", str(destination)], check=True)


def main():
    """Distinguish full-input versus small-input results with the same recipe."""
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("cache", "input", "cached", "direct", "comparison", "output", "work"):
        parser.add_argument("--" + key, type=Path, required=True)
    args = parser.parse_args()
    claim_directory(args.output)
    claim_directory(args.work)
    comparison = json.loads(args.comparison.read_text())
    loci = [
        (row["locus"][0], int(row["locus"][1])) for row in comparison["examples"][:2]
    ]
    if not loci:
        raise RuntimeError("no discrepant loci to investigate")
    for source, name in (
        (args.input, "input"),
        (args.cached, "cached-reference"),
        (args.direct, "direct-reference"),
    ):
        subset(source, args.work / f"{name}.bcf", loci)
    log = args.output / "isolated-direct.out"
    output = args.work / "isolated-direct.bcf"
    with log.open("w") as stream:
        subprocess.run(
            [
                sys.executable,
                "-m",
                "vcfcache",
                "annotate",
                "-i",
                str(args.work / "input.bcf"),
                "-a",
                str(args.cache),
                "-o",
                str(output),
                "-y",
                str(args.cache / "params.snapshot.yaml"),
                "--uncached",
                "--statistics",
                "light",
                "--stats-dir",
                str(args.work / "stats"),
            ],
            stdout=stream,
            stderr=subprocess.STDOUT,
            check=True,
        )
    result = dict(
        loci=loci,
        isolated_vs_full_direct=strict_semantic_compare(
            output, args.work / "direct-reference.bcf"
        ),
        isolated_vs_cached=strict_semantic_compare(
            output, args.work / "cached-reference.bcf"
        ),
        note="Diagnostic only: same recipe, different input context; no field exceptions or cache edits.",
    )
    (args.output / "comparison.json").write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
