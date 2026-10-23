#!/usr/bin/env python3
"""Build immutable standard-BCF vep-rs pilot caches from the bundled blueprint."""

import argparse
import hashlib
import json
import shlex
import shutil
import subprocess
import time
from pathlib import Path

import yaml

from benchmarks.adaptive_routing import claim_directory, digest


def run(command):
    """Run one checked setup pipeline outside benchmark timing."""
    print("$", command, flush=True)
    subprocess.run(["bash", "-euo", "pipefail", "-c", command], check=True)


def main():
    """Freeze resources and build each annotator-specific cache once."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--transcripts", type=Path, required=True)
    args = parser.parse_args()
    claim_directory(args.output)
    if not (args.transcripts / "COMPLETE").exists():
        raise RuntimeError("genome-wide transcript conversion is incomplete")

    def quote(path):
        return shlex.quote(str(path))

    source = Path(
        "/mnt/data/vcfcache_benchmarks/bundled_zenodo_blueprints/gnomad_v4.1_GRCh38_joint_af001/blueprint/vcfcache.bcf"
    )
    tool = Path("/results/veprs/dist/vep")
    gnomad = Path("/results/veprs/data/gnomad_joint_gw_GRCh38.vcf.gz")
    clinvar = Path("/results/clinical_sources/clinvar_GRCh38.vcf.gz")
    tsv = args.output / "clinvar.tsv.gz"
    # Use a native local htslib bgzip/tabix; the existing container supplies
    # these tools on workers where they are not installed on PATH.
    container = "apptainer exec -B /mnt/data -B /results /results/veprs/dist/ensembl-vep_116.2.sif"
    raw = args.output / "clinvar.tsv"
    raw.write_text("#chr\tpos(1-based)\tref\talt\tCV_CLNSIG\tCV_REVSTAT\tCV_CLNDN\n")
    with raw.open("ab") as stream:
        subprocess.run(
            [
                "bcftools",
                "query",
                "-f",
                "%CHROM\t%POS\t%REF\t%ALT\t%INFO/CLNSIG\t%INFO/CLNREVSTAT\t%INFO/CLNDN\n",
                str(clinvar),
            ],
            stdout=stream,
            check=True,
        )
    run(f"{container} bgzip {quote(raw)} && {container} tabix -s1 -b2 -e2 {quote(tsv)}")
    provenance = dict(
        kind="locally_built_veprs_from_bundled_blueprint",
        original_kind="bundled_zenodo_blueprint",
        blueprint=str(source),
        blueprint_sha256=digest(source),
        tool_sha256=digest(tool),
        tool_version=subprocess.check_output(
            [str(tool), "--version"], text=True
        ).strip(),
        transcript_root=str(args.transcripts),
        transcript_conversion_source_hashes=(
            args.transcripts / "SOURCE_SHA256SUMS"
        ).read_text(),
        resources={
            str(p): digest(p)
            for p in (
                gnomad,
                Path(str(gnomad) + ".tbi"),
                clinvar,
                Path(str(clinvar) + ".tbi"),
                tsv,
                Path(str(tsv) + ".tbi"),
            )
        },
        bcftools=subprocess.check_output(["bcftools", "--version"], text=True),
    )
    transcript_hash = hashlib.sha256()
    transcript_files = sorted(args.transcripts.rglob("*.json"))
    for path in transcript_files:
        transcript_hash.update(
            str(path.relative_to(args.transcripts)).encode()
            + b"\0"
            + digest(path).encode()
            + b"\n"
        )
    provenance["transcript_tree_sha256"] = transcript_hash.hexdigest()
    provenance["transcript_files"] = len(transcript_files)
    params = dict(
        genome_build="GRCh38",
        bcftools_cmd="bcftools",
        annotation_tool_cmd=str(tool),
        tool_version_command=f"{tool} --version",
        threads=16,
        temp_dir="/mnt/data/tmp",
        transcript_cache=str(args.transcripts),
        optional_checks={},
    )
    for name, extra in (
        ("veprs_core", ""),
        (
            "veprs_clinical",
            f"--plugin {quote('gnomADc,'+str(gnomad))} --plugin {quote('dbNSFP,'+str(tsv)+',CV_CLNSIG,CV_REVSTAT,CV_CLNDN')}",
        ),
    ):
        base = args.output / name
        cache = base / "cache" / name
        cache.mkdir(parents=True)
        (base / "blueprint").mkdir()
        (base / "workflow").mkdir()
        shutil.copyfile(
            "/results/richness/caches/fastvep_clinical_grch38_af001/workflow/init.yaml",
            base / "workflow/init.yaml",
        )
        command = (
            "${params.bcftools_cmd} view ${INPUT_BCF} -Ov | "
            "${params.annotation_tool_cmd} --assembly GRCh38 --json_cache ${params.transcript_cache} "
            "--fork ${params.threads} --vcf --no_stats --input_file STDIN --output_file STDOUT "
            + extra
            + " | ${params.bcftools_cmd} view -Ob -W --threads ${params.threads} -o ${OUTPUT_BCF}"
        )
        recipe = dict(
            annotation_cmd=command,
            must_contain_info_tag="CSQ",
            genome_build="GRCh38",
            required_tool_version="116.2",
        )
        (cache / "annotation.yaml").write_text(yaml.safe_dump(recipe, sort_keys=False))
        (cache / "params.snapshot.yaml").write_text(
            yaml.safe_dump(params, sort_keys=False)
        )
        target = cache / "vcfcache_annotated.bcf"
        build = command.replace("${INPUT_BCF}", quote(source)).replace(
            "${OUTPUT_BCF}", quote(target)
        )
        for key, value in params.items():
            build = build.replace("${params." + key + "}", str(value))
        start = time.perf_counter()
        run(build)
        duration = time.perf_counter() - start
        data = dict(
            provenance,
            build_seconds=duration,
            recipe_sha256=digest(cache / "annotation.yaml"),
            cache_sha256=digest(target),
            cache_bytes=target.stat().st_size,
            input_records=int(
                subprocess.check_output(["bcftools", "index", "-n", str(source)])
            ),
            output_records=int(
                subprocess.check_output(["bcftools", "index", "-n", str(target)])
            ),
        )
        (cache / "build_provenance.json").write_text(json.dumps(data, indent=2) + "\n")
        print(
            json.dumps(
                dict(name=name, build_seconds=duration, records=data["output_records"])
            ),
            flush=True,
        )
    (args.output / "COMPLETE").touch()


if __name__ == "__main__":
    main()
