#!/usr/bin/env python3
"""Prepare matched miss/minimal inputs outside timing, then measure direct CLI."""

import argparse
import json
import shlex
import subprocess
import sys
from pathlib import Path

import pysam

from benchmarks.adaptive_routing import claim_directory


def main():
    """Keep genomic span in the minimal input; do not assume a constant startup."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cpus", type=int, default=16)
    parser.add_argument("--work-dir", type=Path, required=True)
    args = parser.parse_args()
    claim_directory(args.output)
    claim_directory(args.work_dir)
    misses = args.work_dir / "misses.bcf"
    normalized = args.work_dir / "normalized.bcf"
    subprocess.run(
        [
            "bash",
            "-euo",
            "pipefail",
            "-c",
            f"bcftools norm -m- {shlex.quote(str(args.input))} -Ou | "
            f"bcftools view -e 'ALT=\"*\"' -Ob -W -o {shlex.quote(str(normalized))}",
        ],
        check=True,
    )
    command = (
        f"bcftools annotate -a {shlex.quote(str(args.cache/'vcfcache_annotated.bcf'))} {shlex.quote(str(normalized))} -c INFO -Ou | "
        f"bcftools view -i 'INFO/CSQ==\"\"' -Ob -W -o {shlex.quote(str(misses))}"
    )
    subprocess.run(["bash", "-euo", "pipefail", "-c", command], check=True)
    minimal = args.work_dir / "minimal-same-contigs.bcf"
    seen = set()
    with (
        pysam.VariantFile(str(misses)) as source,
        pysam.VariantFile(str(minimal), "wb", header=source.header) as out,
    ):
        for record in source:
            if record.contig not in seen:
                seen.add(record.contig)
                out.write(record)
    subprocess.run(["bcftools", "index", str(minimal)], check=True)
    (args.output / "preparation.json").write_text(
        json.dumps(
            dict(
                full_input=str(args.input),
                source_cache=str(args.cache),
                preparation_command=command,
                minimal_contigs=sorted(seen),
                note="Input preparation is untimed; measurements are direct full-CLI costs, not isolated startup. Compare with the preceding full-input direct cell.",
            ),
            indent=2,
        )
        + "\n"
    )
    subprocess.run(
        [
            sys.executable,
            "-m",
            "benchmarks.adaptive_routing",
            "--cache",
            str(args.cache),
            "--input",
            str(misses),
            str(minimal),
            "--output",
            str(args.output / "direct"),
            "--work-dir",
            str(args.work_dir / "runs"),
            "--cpus",
            str(args.cpus),
            "--repeats",
            "1",
            "--arms",
            "direct",
            "--slurm-steps",
        ],
        check=True,
    )


if __name__ == "__main__":
    main()
