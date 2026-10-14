#!/usr/bin/env python3
"""Untimed SGDP coverage/QC screen using the publication normalization helpers.

Keeps all new data below --output. Existing samples and caches are read-only.
The common-region analysis is a GQ/DP proxy, not a validated callable mask.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import tempfile
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path

from benchmarks.prepare_external_wgs import Selected, _rename_map, prepare_one

NEW_SAMPLE = "LP6005677-DNA_D03"
SAMPLES = {
    "LP6005592-DNA_C05": "Khomani_San",
    "SS6004480": "Dinka",
    "LP6005677-DNA_G01": "Luo",
    "LP6005442-DNA_B04": "Georgian",
    NEW_SAMPLE: "Khomani_San",
}
AUTOSOMES = {str(n) for n in range(1, 23)}


def write_table(path, rows):
    rows = list(rows)
    if not rows:
        raise ValueError(f"Empty results: {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def digest(path):
    checksum = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            checksum.update(chunk)
    return checksum.hexdigest()


def stream(commands):
    """Stream a checked pipeline without risking a stderr pipe deadlock."""
    processes = []
    previous = None
    with tempfile.TemporaryFile(mode="w+") as errors:
        for command in commands:
            process = subprocess.Popen(
                command, stdin=previous, stdout=subprocess.PIPE,
                stderr=errors, text=True,
            )
            if previous is not None:
                previous.close()
            previous = process.stdout
            processes.append(process)
        yield from previous
        codes = [process.wait() for process in processes]
        if any(codes):
            errors.seek(0)
            raise RuntimeError(f"Pipeline failed {codes}: {errors.read()}")


def intersect(a, b):
    """Intersect sorted, disjoint zero-based half-open intervals."""
    result = []
    i = j = 0
    while i < len(a) and j < len(b):
        lo, hi = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if lo < hi:
            result.append((lo, hi))
        if a[i][1] <= b[j][1]:
            i += 1
        else:
            j += 1
    return result


def append_interval(intervals, lo, hi):
    if intervals and lo <= intervals[-1][1]:
        intervals[-1] = (intervals[-1][0], max(hi, intervals[-1][1]))
    else:
        intervals.append((lo, hi))


def quality_and_regions(source, sample, out):
    """Bounded chr20 source-quality/control-region diagnostic."""
    intervals = {}
    bed = out / f"{sample}.chr20.genotyped_regions.bed"
    stats = out / f"{sample}.source_quality.tsv"
    if stats.exists() and bed.exists():
        with stats.open() as handle:
            row = next(csv.DictReader(handle, delimiter="\t"))
        with bed.open() as handle:
            for line in handle:
                chrom, lo, hi = line.split()
                intervals.setdefault(chrom, []).append((int(lo), int(hi)))
        return row, intervals
    format_string = r"%CHROM\t%POS\t%END\t%REF\t%ALT\t%QUAL\t%FILTER[\t%GT\t%DP\t%GQ\t%MIN_DP]\n"
    contigs = {line.split("\t")[0] for line in subprocess.check_output(
        ["bcftools", "index", "-s", str(source)], text=True
    ).splitlines()}
    region = "chr20" if "chr20" in contigs else "20"
    if region not in contigs:
        raise RuntimeError(f"No chromosome 20 in {source}")
    awk = Path(__file__).with_name("source_quality.awk")
    lines = list(stream([
        ["bcftools", "query", "-r", region, "-f", format_string, str(source)],
        ["awk", "-v", f"bed={bed}", "-f", str(awk)],
    ]))
    row = {"sample": sample, "population": SAMPLES[sample], "source_scope": "chr20"}
    names = ["carried_source_records", "gq_missing", "gq_below20", "dp_missing", "dp_below10", "qual_missing", "qual_below30"]
    row.update(dict(zip(names, map(int, lines[-1].split()))))
    with bed.open() as handle:
        for line in handle:
            chrom, lo, hi = line.split()
            intervals.setdefault(chrom, []).append((int(lo), int(hi)))
    row["high_confidence_genotyped_bases"] = sum(hi-lo for runs in intervals.values() for lo, hi in runs)
    write_table(out / f"{sample}.source_quality.tsv", [row])
    return row, intervals


def count_variants(path, cache=None, regions=None):
    commands = []
    if cache:
        commands.append(["bcftools", "isec", "-c", "none", "-n=2", "-w1", "-Ou", str(path), str(cache)])
    fmt = r"%REF\t%ALT\t%QUAL\n"
    query = ["bcftools", "query", "-f", fmt]
    if regions:
        query += ["-R", str(regions), "--regions-overlap", "0"]
    if cache:
        # Region restrictions must be applied to isec's indexed inputs.
        if regions:
            commands[0][2:2] = ["-R", str(regions), "--regions-overlap", "0"]
        query = ["bcftools", "query", "-f", fmt]
    else:
        query.append(str(path))
    commands.append(query)
    commands.append(["awk", "-f", str(Path(__file__).with_name("variant_counts.awk"))])
    counts = Counter()
    for line in stream(commands):
        kind, quality, count = line.rstrip().split("\t")
        counts[(kind, quality)] = int(count)
    return counts


def screen(path, sample, cache, out, scope="whole_autosomes", regions=None):
    target = out / f"{sample}.{scope}.tsv"
    if target.exists():
        with target.open() as handle:
            return list(csv.DictReader(handle, delimiter="\t"))
    totals = count_variants(path, regions=regions)
    hits = count_variants(path, cache=cache, regions=regions)
    rows = []
    for (kind, quality), count in sorted(totals.items()):
        if count == 0:
            continue
        hit_count = hits[(kind, quality)]
        rows.append({
            "sample": sample, "population": SAMPLES[sample],
            "scope": scope, "variant_type": kind, "quality_subset": quality,
            "input_records": count, "cache_hits": hit_count,
            "cache_misses": count-hit_count, "hit_rate": hit_count/count,
            "cache_strategy": "gnomad_any_stratum_af_0.01",
            "timed_annotation": False,
        })
    write_table(target, rows)
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--cache", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    out = args.output / "results"
    out.mkdir(exist_ok=True)
    stem = NEW_SAMPLE + ".hs38DH"
    name = stem + ".BQSR.autosome.g.vcf.gz"
    url = f"https://ddbj.nig.ac.jp/public/public-human-genomes/GRCh38/SGDP/GVCF/{stem}/{name}"
    source = args.output / "sources/sgdp" / NEW_SAMPLE / name
    source.parent.mkdir(parents=True, exist_ok=True)
    for address, target in [(url, source), (url + ".tbi", Path(str(source)+".tbi"))]:
        if not target.exists():
            print(f"Downloading {address}", flush=True)
            subprocess.run(["curl", "-fL", "--retry", "3", "--connect-timeout", "30", "--output", str(target)+".partial", address], check=True)
            Path(str(target)+".partial").replace(target)
    row = Selected(
        cohort="sgdp", sample=NEW_SAMPLE, role="screen_only", provider="DDBJ_NIG",
        population="Khomani_San", region="Africa", sex="", assembly="GRCh38",
        source_kind="gVCF", url=url, index_url=url+".tbi", source_name=name,
        source_bytes=source.stat().st_size, upstream_md5="", documented_overlap="none",
        landing_url=url.rsplit("/", 1)[0], selection_seed="lower_hit_followup_20261014",
        selection_key="Khomani_San_followup",
    )
    print("Preparing second Khomani San using existing publication helper", flush=True)
    # gVCF reference blocks cannot produce a carried allele after splitting.
    # Remove only hom-reference/missing records ahead of the unchanged helper
    # to avoid reference-normalizing millions of non-variant blocks. Keep the
    # complete source for genotype quality and common-region diagnostics.
    carried_source = source.with_name(source.stem + ".carried.vcf.gz")
    if not carried_source.exists():
        subprocess.run([
            "bcftools", "view", "-i", 'GT="alt"', "--threads", "8",
            "-Oz", "-o", str(carried_source)+".partial", str(source),
        ], check=True)
        Path(str(carried_source)+".partial").replace(carried_source)
        subprocess.run(["bcftools", "index", "-t", str(carried_source)], check=True)
    prepared = prepare_one(
        args.output, replace(row, source_name=carried_source.name),
        args.reference, _rename_map(args.output),
    )
    paths = {
        s: prepared if s == NEW_SAMPLE else args.root / "samples/GRCh38/sgdp" / f"{s}.GRCh38.small_variants.vcf.gz"
        for s in SAMPLES
    }
    provenance = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "purpose": "untimed coverage and QC diagnostic; no new annotation performance claims",
        "source_url": url, "cache": str(args.cache), "reference": str(args.reference),
        "preselection": 'GT="alt" before the existing normalization helper; complete source retained',
        "bcftools_version": subprocess.check_output(["bcftools", "--version"], text=True).splitlines()[0],
        "common_region_definition": "chr20 PASS/dot gVCF records; complete GT; GQ>=20; MIN_DP>=10 when present, otherwise DP>=10; intersection of all five samples; variant start in interval; compare with matched chr20 baseline",
        "samples": SAMPLES,
        "source_sha256": digest(source),
        "reference_sha256": digest(args.reference),
        "cache_sha256": digest(args.cache),
        "normalized_input_sha256": {sample: digest(path) for sample, path in paths.items()},
    }
    (out / "SCREEN_PROVENANCE.json").write_text(json.dumps(provenance, indent=2)+"\n")
    print("Screening whole-autosome allele coverage", flush=True)
    with ThreadPoolExecutor(max_workers=3) as pool:
        results = list(pool.map(lambda item: screen(item[1], item[0], args.cache, out), paths.items()))
    write_table(out / "coverage_screen.tsv", [r for rows in results for r in rows])
    print("Constructing chr20 common high-confidence genotyped regions", flush=True)
    common = None
    quality_rows = []
    def source_qc(sample):
        raw = source if sample == NEW_SAMPLE else args.root / "sources/sgdp" / sample / f"{sample}.hs38DH.BQSR.autosome.g.vcf.gz"
        result = quality_and_regions(raw, sample, out)
        print(f"Source QC complete: {sample}", flush=True)
        return result

    with ThreadPoolExecutor(max_workers=3) as pool:
        quality_results = list(pool.map(source_qc, SAMPLES))
    for qrow, intervals in quality_results:
        quality_rows.append(qrow)
        common = intervals if common is None else {chrom: intersect(runs, intervals.get(chrom, [])) for chrom, runs in common.items()}
    write_table(out / "source_quality.tsv", quality_rows)
    bed = out / "common_genotyped_regions.bed"
    with bed.open("w") as handle:
        for chrom in sorted(common, key=lambda c: int(c.removeprefix("chr"))):
            for lo, hi in common[chrom]:
                handle.write(f"{chrom}\t{lo}\t{hi}\n")
    common_bases = sum(hi-lo for runs in common.values() for lo, hi in runs)
    if common_bases <= 0:
        raise RuntimeError("Empty common genotyped-region intersection")
    provenance["common_genotyped_bases"] = common_bases
    (out / "SCREEN_PROVENANCE.json").write_text(json.dumps(provenance, indent=2)+"\n")
    chromosome_bed = out / "chr20_control.bed"
    chromosome_bed.write_text("chr20\t0\t64444167\n")
    with ThreadPoolExecutor(max_workers=3) as pool:
        results = list(pool.map(lambda item: screen(item[1], item[0], args.cache, out, "chr20_all_prepared_alleles", chromosome_bed), paths.items()))
    write_table(out / "coverage_chr20.tsv", [r for rows in results for r in rows])
    with ThreadPoolExecutor(max_workers=3) as pool:
        results = list(pool.map(lambda item: screen(item[1], item[0], args.cache, out, "common_genotyped_regions", bed), paths.items()))
    write_table(out / "coverage_common_regions.tsv", [r for rows in results for r in rows])
    (out / "COMPLETE").write_text(datetime.now(timezone.utc).isoformat()+"\n")
    print("COMPLETE", flush=True)


if __name__ == "__main__":
    main()
