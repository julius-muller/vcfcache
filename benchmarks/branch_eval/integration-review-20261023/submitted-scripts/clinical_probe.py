#!/usr/bin/env python3
"""Reproduce clinical discrepancies by varying record metadata, not engines."""

import argparse
import hashlib
import json
import shlex
import subprocess
import time
from collections import Counter
from pathlib import Path

import pysam
import yaml


def key(record):
    """Use the full allele key; never match positions alone."""
    return record.contig, record.pos, tuple(record.alleles)


def csqs(path):
    """Canonicalize CSQ-entry ordering only, retaining every value/occurrence."""
    result = {}
    with pysam.VariantFile(str(path)) as source:
        for record in source:
            k = key(record)
            if k in result:
                raise RuntimeError(f"duplicate key requires explicit handling: {k}")
            result[k] = tuple(sorted(record.info.get("CSQ", ())))
    return result


def extract(source_path, destination, wanted):
    """Fetch exact matching records without modifying the original cache/input."""
    found = set()
    with pysam.VariantFile(str(source_path)) as source, pysam.VariantFile(
        str(destination), "wb", header=source.header
    ) as out:
        for contig in source.header.contigs:
            positions = sorted({pos for name,pos,_ in wanted if name == contig})
            for pos in positions:
                for record in source.fetch(contig,pos-1,pos):
                    if key(record) in wanted:
                        if key(record) in found:
                            raise RuntimeError("ambiguous extraction")
                        found.add(key(record))
                        out.write(record)
    if found != wanted:
        raise RuntimeError(f"missing exact source keys: {wanted-found}")
    subprocess.run(["bcftools","index",str(destination)],check=True)


def main():
    """Probe all WES discrepancies and preserve raw commands/results."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--work",type=Path,required=True)
    args=parser.parse_args()
    for directory in (args.output,args.work):
        directory.mkdir(parents=True,exist_ok=True)
        if any(directory.iterdir()):
            raise FileExistsError(directory)
        (directory/".owner").touch(exist_ok=False)
    root=Path("/mnt/data/tmp/branch-eval-veprs-clinical-1323")
    cached=csqs(root/"runs/wes-native-1/output.bcf")
    direct=csqs(root/"runs/wes-direct-1/output.bcf")
    if set(cached)!=set(direct):
        raise RuntimeError("input-key discrepancy is not an annotation-only case")
    wanted={k for k in direct if cached[k]!=direct[k]}
    sample=args.work/"sample.bcf"
    blueprint=args.work/"blueprint.bcf"
    extract(Path("/results/branch-eval/inputs/wes.bcf"),sample,wanted)
    extract(Path("/mnt/data/vcfcache_benchmarks/bundled_zenodo_blueprints/gnomad_v4.1_GRCh38_joint_af001/blueprint/vcfcache.bcf"),blueprint,wanted)
    params=yaml.safe_load((root/"cache/veprs_clinical/params.snapshot.yaml").read_text())
    recipe=yaml.safe_load((root/"cache/veprs_clinical/annotation.yaml").read_text())
    manifest=dict(discrepancies=len(wanted),keys=sorted(wanted),recipe=recipe,
                  params=params,tool_sha256=hashlib.sha256(Path(params["annotation_tool_cmd"]).read_bytes()).hexdigest(),
                  note="CSQ-only diagnosis; this does not replace or relax full-record benchmark audits.")
    (args.output/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    rows=[]
    variants=[("sample_original",sample,[],16),("sample_repeat",sample,[],16),
              ("sample_single_thread",sample,[],1),("blueprint_original",blueprint,[],16),
              ("sample_sites_only",sample,["view","-G"],16),
              ("sample_no_info",sample,["annotate","-x","INFO"],16),
              ("sample_no_id",sample,["annotate","-x","ID"],16),
              ("sample_no_samples_or_info",sample,["annotate","-x","INFO"],16)]
    for name,source,transform,threads in variants:
        if transform:
            transformed=args.work/f"{name}.bcf"
            subprocess.run(["bcftools",*transform,str(source),"-Ob","-o",str(transformed)],check=True)
            source=transformed
        if name=="sample_no_samples_or_info":
            transformed=args.work/f"{name}-sites.bcf"
            subprocess.run(["bcftools","view","-G",str(source),"-Ob","-o",str(transformed)],check=True)
            source=transformed
        output=args.work/f"{name}-annotated.bcf"
        command=recipe["annotation_cmd"]
        for field,value in dict(params,threads=threads).items():
            command=command.replace("${params."+field+"}",shlex.quote(str(value)))
        command=command.replace("${INPUT_BCF}",shlex.quote(str(source))).replace("${OUTPUT_BCF}",shlex.quote(str(output)))
        (args.output/f"{name}.command.sh").write_text("set -euo pipefail\n"+command+"\n")
        start=time.perf_counter()
        with (args.output/f"{name}.out").open("w") as log:
            subprocess.run(["bash","-euo","pipefail","-c",command],stdout=log,stderr=subprocess.STDOUT,check=True)
        observed=csqs(output)
        rows.append(dict(condition=name,threads=threads,seconds=time.perf_counter()-start,
                         records=len(observed),matches_direct=sum(observed.get(k)==direct[k] for k in wanted),
                         matches_cached=sum(observed.get(k)==cached[k] for k in wanted)))
        (args.output/f"{name}.csq.json").write_text(json.dumps([dict(key=k,csq=v) for k,v in observed.items()],indent=2)+"\n")
        (args.output/"results.json").write_text(json.dumps(rows,indent=2)+"\n")
        print(json.dumps(rows[-1]),flush=True)


if __name__ == "__main__":
    main()
