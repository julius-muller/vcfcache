#!/usr/bin/env python3
"""Reproduce clinical discrepancies by varying record metadata, not engines."""

import argparse
import hashlib
import json
import shlex
import subprocess
import time
from pathlib import Path

import pysam
import yaml


def key(record):
    """Use the full allele key; never match positions alone."""
    return record.contig, record.pos, tuple(record.alleles)


def csqs(path, wanted=None):
    """Canonicalize CSQ-entry ordering only, retaining every value/occurrence."""
    result = {}
    with pysam.VariantFile(str(path)) as source:
        for record in source:
            k = key(record)
            if wanted is not None and k not in wanted:
                continue
            if k in result:
                raise RuntimeError(f"duplicate key requires explicit handling: {k}")
            result[k] = tuple(sorted(record.info.get("CSQ", ())))
    return result


def extract(source_path, destination, wanted):
    """Fetch exact matching records without modifying the original cache/input."""
    found = set()
    with (
        pysam.VariantFile(str(source_path)) as source,
        pysam.VariantFile(str(destination), "wb", header=source.header) as out,
    ):
        for contig in source.header.contigs:
            positions = sorted({pos for name, pos, _ in wanted if name == contig})
            for pos in positions:
                for record in source.fetch(contig, pos - 1, pos):
                    if key(record) in wanted:
                        if key(record) in found:
                            raise RuntimeError("ambiguous extraction")
                        found.add(key(record))
                        out.write(record)
    if found != wanted:
        raise RuntimeError(f"missing exact source keys: {wanted-found}")
    subprocess.run(["bcftools", "index", str(destination)], check=True)


def main():
    """Probe all WES discrepancies and preserve raw commands/results."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--context-spans", type=int, nargs="+", default=[0, 1000])
    parser.add_argument("--whole-contig", action="store_true")
    args = parser.parse_args()
    for directory in (args.output, args.work):
        directory.mkdir(parents=True, exist_ok=True)
        if any(directory.iterdir()):
            raise FileExistsError(directory)
        (directory / ".owner").touch(exist_ok=False)
    root = Path("/mnt/data/tmp/branch-eval-veprs-clinical-1323")
    cached = csqs(root / "runs/wes-native-1/output.bcf")
    direct = csqs(root / "runs/wes-direct-1/output.bcf")
    if set(cached) != set(direct):
        raise RuntimeError("input-key discrepancy is not an annotation-only case")
    wanted = {k for k in direct if cached[k] != direct[k]}
    sample = args.work / "sample.bcf"
    blueprint = args.work / "blueprint.bcf"
    extract(Path("/results/branch-eval/inputs/wes.bcf"), sample, wanted)
    extract(
        Path(
            "/mnt/data/vcfcache_benchmarks/bundled_zenodo_blueprints/gnomad_v4.1_GRCh38_joint_af001/blueprint/vcfcache.bcf"
        ),
        blueprint,
        wanted,
    )
    cache_subset = args.work / "cache-payload.bcf"
    extract(root / "cache/veprs_clinical/vcfcache_annotated.bcf", cache_subset, wanted)
    payload = csqs(cache_subset)
    manifest_payload = dict(
        records=len(payload),
        matches_cached=sum(payload[k] == cached[k] for k in wanted),
        matches_direct=sum(payload[k] == direct[k] for k in wanted),
    )
    params = yaml.safe_load(
        (root / "cache/veprs_clinical/params.snapshot.yaml").read_text()
    )
    recipe = yaml.safe_load((root / "cache/veprs_clinical/annotation.yaml").read_text())
    manifest = dict(
        discrepancies=len(wanted),
        keys=sorted(wanted),
        recipe=recipe,
        params=params,
        tool_sha256=hashlib.sha256(
            Path(params["annotation_tool_cmd"]).read_bytes()
        ).hexdigest(),
        note="CSQ-only diagnosis; this does not replace or relax full-record benchmark audits.",
    )
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    manifest["cache_payload_check"] = manifest_payload
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    # Context probes retain every blueprint allele, not just the target key.
    reference = Path(
        "/mnt/data/vcfcache_benchmarks/bundled_zenodo_blueprints/gnomad_v4.1_GRCh38_joint_af001/blueprint/vcfcache.bcf"
    )
    context_inputs = []
    for span in args.context_spans:
        destination = args.work / f"blueprint-context-{span}.bcf"
        positions = {(name, pos) for name, pos, _ in wanted}
        with (
            pysam.VariantFile(str(reference)) as source,
            pysam.VariantFile(str(destination), "wb", header=source.header) as out,
        ):
            seen = set()
            for contig in source.header.contigs:
                for name, pos in sorted(positions):
                    if name != contig:
                        continue
                    for record in source.fetch(
                        name, max(0, pos - 1 - span), pos + span
                    ):
                        k = key(record)
                        if k in seen:
                            continue
                        seen.add(k)
                        out.write(record)
        # Nearby intervals can overlap; sorting preserves the recipe contract.
        sorted_file = args.work / f"blueprint-context-{span}-sorted.bcf"
        subprocess.run(
            ["bcftools", "sort", str(destination), "-Ob", "-o", str(sorted_file)],
            check=True,
        )
        context_inputs.append((f"blueprint_context_{span}", sorted_file, [], 16))
    if args.whole_contig:
        destination = args.work / "blueprint-full-contigs.bcf"
        contigs = sorted({name for name, _, _ in wanted})
        subprocess.run([
            "bcftools", "view", str(reference), "-r", ",".join(contigs),
            "-Ob", "-o", str(destination)
        ], check=True)
        context_inputs.append(("blueprint_full_contigs", destination, [], 16))
    rows = []
    variants = [
        ("sample_original", sample, [], 16),
        ("sample_repeat", sample, [], 16),
        ("sample_single_thread", sample, [], 1),
        ("blueprint_original", blueprint, [], 16),
        ("sample_sites_only", sample, ["view", "-G"], 16),
        ("sample_no_info", sample, ["annotate", "-x", "INFO"], 16),
        ("sample_no_id", sample, ["annotate", "-x", "ID"], 16),
        ("sample_no_samples_or_info", sample, ["annotate", "-x", "INFO"], 16),
        *context_inputs,
    ]
    for name, source, transform, threads in variants:
        if transform:
            transformed = args.work / f"{name}.bcf"
            subprocess.run(
                ["bcftools", *transform, str(source), "-Ob", "-o", str(transformed)],
                check=True,
            )
            source = transformed
        if name == "sample_no_samples_or_info":
            transformed = args.work / f"{name}-sites.bcf"
            subprocess.run(
                ["bcftools", "view", "-G", str(source), "-Ob", "-o", str(transformed)],
                check=True,
            )
            source = transformed
        output = args.work / f"{name}-annotated.bcf"
        command = recipe["annotation_cmd"]
        for field, value in dict(params, threads=threads).items():
            command = command.replace(
                "${params." + field + "}", shlex.quote(str(value))
            )
        command = command.replace("${INPUT_BCF}", shlex.quote(str(source))).replace(
            "${OUTPUT_BCF}", shlex.quote(str(output))
        )
        (args.output / f"{name}.command.sh").write_text(
            "set -euo pipefail\n" + command + "\n"
        )
        start = time.perf_counter()
        with (args.output / f"{name}.out").open("w") as log:
            subprocess.run(
                ["bash", "-euo", "pipefail", "-c", command],
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
            )
        observed = csqs(output, wanted)
        rows.append(
            dict(
                condition=name,
                threads=threads,
                seconds=time.perf_counter() - start,
                records=int(subprocess.check_output(["bcftools", "index", "-n", str(output)])),
                matches_direct=sum(observed.get(k) == direct[k] for k in wanted),
                matches_cached=sum(observed.get(k) == cached[k] for k in wanted),
            )
        )
        (args.output / f"{name}.csq.json").write_text(
            json.dumps([dict(key=k, csq=v) for k, v in observed.items()], indent=2)
            + "\n"
        )
        (args.output / "results.json").write_text(json.dumps(rows, indent=2) + "\n")
        print(json.dumps(rows[-1]), flush=True)


if __name__ == "__main__":
    main()
