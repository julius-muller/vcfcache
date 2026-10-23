#!/usr/bin/env python3
"""Verify the reduced clinical case across thread counts and plugin selection."""

import argparse
import json
import shlex
import subprocess
from pathlib import Path

import pysam


def main():
    """Keep the external binary unchanged and retain every diagnostic output."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if any(args.output.iterdir()):
        raise FileExistsError(args.output)
    (args.output / ".owner").touch(exist_ok=False)
    source = Path("/results/branch-integration/minimal-clinical-1331")
    manifest = json.loads((source / "result.json").read_text())
    target = manifest["target"][1]
    rows = []
    with pysam.VariantFile(str(source / "input.bcf")) as inp:
        header = inp.header.copy()
        records = [r.copy() for r in inp]
    for mode in ("both", "dbnsfp-only"):
        for threads in (1, 16):
            for subset in (
                "all",
                "target-only",
                "first-and-target",
                "second-and-target",
            ):
                name = f"{mode}-{threads}-{subset}"
                inp = args.output / f"{name}.bcf"
                selected = {
                    "all": records,
                    "target-only": records[-1:],
                    "first-and-target": [records[0], records[-1]],
                    "second-and-target": records[1:],
                }[subset]
                with pysam.VariantFile(str(inp), "wb", header=header) as out:
                    for record in selected:
                        out.write(record)
                output = args.output / f"{name}-output.bcf"
                command = manifest["recipe"]["annotation_cmd"]
                if mode == "dbnsfp-only":
                    command = command.replace(
                        "--plugin gnomADc,/results/veprs/data/gnomad_joint_gw_GRCh38.vcf.gz ",
                        "",
                    )
                for key, value in dict(manifest["params"], threads=threads).items():
                    command = command.replace(
                        "${params." + key + "}", shlex.quote(str(value))
                    )
                command = command.replace(
                    "${INPUT_BCF}", shlex.quote(str(inp))
                ).replace("${OUTPUT_BCF}", shlex.quote(str(output)))
                (args.output / f"{name}.sh").write_text(
                    "set -euo pipefail\n" + command + "\n"
                )
                with (args.output / f"{name}.log").open("w") as log:
                    subprocess.run(
                        ["bash", "-euo", "pipefail", "-c", command],
                        stdout=log,
                        stderr=log,
                        check=True,
                    )
                with pysam.VariantFile(str(output)) as out:
                    description = out.header.info["CSQ"].description
                    fields = description.split("Format: ", 1)[1].split("|")
                    indices = [
                        fields.index("dbNSFP_" + field)
                        for field in ("CV_CLNSIG", "CV_REVSTAT", "CV_CLNDN")
                    ]
                    record = next(r for r in out if r.pos == target)
                    values = [
                        tuple(csq.split("|")[i] for i in indices)
                        for csq in record.info["CSQ"]
                    ]
                rows.append(
                    dict(
                        mode=mode,
                        threads=threads,
                        subset=subset,
                        clinical_values=values,
                    )
                )
                (args.output / "results.json").write_text(
                    json.dumps(rows, indent=2) + "\n"
                )
                print(json.dumps(rows[-1]), flush=True)


if __name__ == "__main__":
    main()
