#!/usr/bin/env python3
"""Reduce an input-context-sensitive clinical discrepancy without ignoring fields."""

import argparse
import hashlib
import json
import shlex
import subprocess
from pathlib import Path

import pysam
import yaml


def key(record):
    """Identify the complete variant allele, not merely its position."""
    return record.contig, record.pos, tuple(record.alleles)


def main():
    """Delta-debug a measured context into a local, independently runnable case."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    base = Path("/mnt/data/tmp/integration-clinical-1329")
    frozen = Path("/results/branch-integration/clinical-1329")
    manifest = json.loads((frozen / "manifest.json").read_text())
    target_raw = manifest["keys"][0]
    target = target_raw[0], target_raw[1], tuple(target_raw[2])
    with pysam.VariantFile(str(base / "blueprint-context-100000-sorted.bcf")) as source:
        header = source.header.copy()
        records = [r.copy() for r in source]
    reference = {
        (r["key"][0], r["key"][1], tuple(r["key"][2])): tuple(r["csq"])
        for r in json.loads((frozen / "blueprint_original.csq.json").read_text())
    }[target]
    target_record = next(r for r in records if key(r) == target)
    candidates = [r for r in records if key(r) != target]
    recipe = manifest["recipe"]["annotation_cmd"]
    for name, value in manifest["params"].items():
        recipe = recipe.replace("${params." + name + "}", shlex.quote(str(value)))
    rows = []

    def differs(context):
        index = len(rows)
        inp = args.output / "input.bcf"
        out = args.output / "output.bcf"
        # Preserve the exact input header, alleles and coordinate order.
        contig_order = {name: i for i, name in enumerate(header.contigs)}
        ordered = sorted([target_record, *context], key=lambda r: (contig_order[r.contig], r.pos))
        with pysam.VariantFile(str(inp), "wb", header=header) as destination:
            for record in ordered:
                destination.write(record)
        command = recipe.replace("${INPUT_BCF}", shlex.quote(str(inp))).replace(
            "${OUTPUT_BCF}", shlex.quote(str(out))
        )
        # Each iteration replaces only this probe's own scratch output/index.
        for previous in (out, Path(str(out) + ".csi")):
            previous.unlink(missing_ok=True)
        with (args.output / f"probe-{index}.log").open("w") as log:
            subprocess.run(["bash", "-euo", "pipefail", "-c", command], stdout=log, stderr=log, check=True)
        with pysam.VariantFile(str(out)) as result:
            observed = next(tuple(sorted(r.info["CSQ"])) for r in result if key(r) == target)
        changed = observed != reference
        rows.append(dict(probe=index, context_records=len(context), changed=changed))
        (args.output / "iterations.json").write_text(json.dumps(rows, indent=2) + "\n")
        print(json.dumps(rows[-1]), flush=True)
        return changed, observed, command

    if not differs(candidates)[0] or differs([])[0]:
        raise RuntimeError("context-only precondition did not reproduce")
    # Classic ddmin, using complements: keep a deletion whenever it preserves
    # the discrepancy. It handles interactions, not just a monotone threshold.
    granularity = 2
    while candidates:
        block = max(1, (len(candidates) + granularity - 1) // granularity)
        reduced = False
        for start in range(0, len(candidates), block):
            complement = candidates[:start] + candidates[start + block:]
            if differs(complement)[0]:
                candidates = complement
                granularity = max(2, granularity - 1)
                reduced = True
                break
        if not reduced:
            if granularity >= len(candidates):
                break
            granularity = min(len(candidates), granularity * 2)
        if len(rows) > 120:
            raise RuntimeError("probe budget exceeded; retain last measured context")
    changed, observed, command = differs(candidates)
    assert changed
    # Final BCF retains the exact reproducing input. Include readable VCF too.
    subprocess.run(["bcftools", "view", str(args.output / "input.bcf"), "-Ov", "-o", str(args.output / "minimal.vcf")], check=True)
    (args.output / "reproduce.sh").write_text("#!/bin/bash\nset -euo pipefail\n" + command + "\n")
    (args.output / "result.json").write_text(json.dumps(dict(
        target=target, other_keys=[key(r) for r in candidates],
        isolated_csq=reference, contextual_csq=observed,
        tool_sha256=manifest["tool_sha256"], recipe=manifest["recipe"],
        params=manifest["params"], input_sha256=hashlib.sha256((args.output / "input.bcf").read_bytes()).hexdigest(),
        note="Diagnostic CSQ-only reduction, not a relaxed benchmark comparator; minimal with respect to single-record deletion."
    ), indent=2) + "\n")
    # Keep all commands local; no issue or external report is submitted.
    (args.output / "params.yaml").write_text(yaml.safe_dump(manifest["params"]))


if __name__ == "__main__":
    main()
