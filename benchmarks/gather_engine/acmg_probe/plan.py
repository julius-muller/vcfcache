"""Submit the chr1 probe: plain and cached fastVEP 0.5.0 for each recipe tier.

Two GRCh38 genomes, restricted to chr1, each in one exclusive unit so that all
tiers of a sample are timed on the same node. Uses the cohort runner's cells.
"""

import argparse
import json
import subprocess
from pathlib import Path

from benchmarks.gather_engine import cohort
from benchmarks.run_pilot import sha256sum

COHORT = "/results/gather-engine/cohort-v116-fv050-20261111/manifest.json"
PROBE = "/results/fastvep05/acmg"
SAMPLES = ["KPGP-00128", "LP6005443-DNA_B07"]
TIERS = ["A", "B0", "B", "C"]


def main():
    """Build one unit per sample and submit them with the cohort runner."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--deployment", required=True)
    parser.add_argument("--commit", required=True)
    args = parser.parse_args()
    data = json.loads(Path(COHORT).read_text())
    bed = Path(PROBE, "chr1.bed")
    bed.write_text("chr1\t0\t248956422\n")
    units = []
    for sample in SAMPLES:
        template = next(
            c
            for u in data["units"]
            for c in u["cells"]
            if (c["tool"], c["level"], c["assembly"], c["assay"], c["sample"])
            == ("fastvep", "clinical", "GRCh38", "wgs", sample)
        )
        records = int(
            subprocess.check_output(
                f"{cohort.BCFTOOLS} view -H -R {bed} {template['input']['path']} | wc -l",
                shell=True,
                text=True,
            )
        )
        cells = []
        for tier in TIERS:
            name = f"fastvep050_{tier}_chr1"
            cache = f"{PROBE}/caches/{name}/cache/{name}"
            cells.append(
                dict(
                    template,
                    level=f"tier{tier}",
                    bed=str(bed),
                    input=dict(template["input"], records=records),
                    cache=cache,
                    cache_sha256=sha256sum(Path(cache, "vcfcache_annotated.bcf.csi")),
                    recipe_level=6,
                    repeats=2,
                    arms=["plain", "cached"],
                    expected=60,
                )
            )
        units.append(dict(id=f"probe-{sample}", tier="A", cells=cells))
    cohort.submit(args.root, args.deployment, args.commit, units, probe=PROBE)


if __name__ == "__main__":
    main()
