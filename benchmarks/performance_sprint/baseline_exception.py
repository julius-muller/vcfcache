"""Independent direct controls only when an old digest cannot be reproduced."""

import argparse
import json
import os
import shlex
import subprocess
import tempfile
from pathlib import Path

from benchmarks.engine_profile.run import PYTHON, read
from benchmarks.performance_sprint.manuscript import annotation_digest
from benchmarks.run_pilot import strict_semantic_compare, write_json_atomic


def direct_command(task, source, output, auxiliary):
    """Use the original C2 bare-wrapper contract, not an extra cached annotation."""
    params = task["parameters"]
    if task["tool"] == "fastvep":
        return [
            "/results/richness/scripts/fastvep_run.sh",
            str(source),
            str(output),
            str(auxiliary),
            params["annotation_tool_cmd"],
            params["transcript_cache"],
            params["fasta"],
            f"/results/richness/sa/clinical_{task['assembly']}",
        ]
    script = (
        params["snpeff_clinical_run"]
        if task["level"] == "clinical"
        else params["snpeff_parallel"]
    )
    return [
        script,
        str(source),
        str(output),
        str(auxiliary),
        params["java_cmd"],
        params["snpeff_jar"],
        params["snpeff_config"],
        params["snpeff_data"],
        params["snpeff_genome"],
    ]


def run_direct(task, source, output, directory, deployment, environment):
    """Measure one new bare direct exception in its own fresh Slurm-step cgroup."""
    directory.mkdir(parents=True, exist_ok=False)
    auxiliary = output.parent / "direct-aux"
    auxiliary.mkdir(exist_ok=True)
    command = direct_command(task, source, output, auxiliary)
    launch = [
        "srun",
        "--exclusive",
        "--exact",
        "--cpu-bind=cores",
        "-n1",
        f"-c{task['cpus']}",
        PYTHON,
        str(Path(deployment) / "benchmarks/engine_profile/measure.py"),
        "--directory",
        str(directory),
        "--",
        *command,
    ]
    env = dict(
        environment,
        RAYON_NUM_THREADS=str(task["parameters"].get("rayon_threads", task["cpus"])),
        SNPEFF_BCFTOOLS_THREADS=str(task["cpus"]),
    )
    write_json_atomic(
        directory / "command.json",
        dict(
            command=launch,
            reason="Historical annotation/key digest was not reproduced; independent direct control required",
            contract="original C2 bare wrapper; no cache/CLI in the timed direct path",
        ),
    )
    with (directory / "run.log").open("w") as log:
        subprocess.run(launch, env=env, stdout=log, stderr=log, check=True)
    return read(directory / "resources.json")


def main():
    """Adjudicate retained failed small inputs without repeating cached timings."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("campaign", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cell", action="append")
    parser.add_argument("--java-command")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    manifest = read(args.campaign / "manifest.json")
    for identifier in args.cell or ["d3a4e19d85523300", "3693ce45be7d4941"]:
        task = next(t for t in manifest["tasks"] if t["id"] == identifier)
        if args.java_command:
            task["parameters"]["java_cmd"] = args.java_command
        directory = args.output / identifier
        directory.mkdir()
        with tempfile.TemporaryDirectory(
            prefix="direct-exception-", dir="/mnt/data/tmp"
        ) as temporary:
            scratch = Path(temporary)
            source = scratch / "input.bcf"
            pipeline = [
                [
                    "bcftools",
                    "view",
                    task["input"]["path"],
                    "-R",
                    task["preprocess"]["bed"]["path"],
                    "-Ou",
                ],
                ["bcftools", "norm", "-m-", "-Ou"],
                ["bcftools", "view", "-e", 'ALT="*"', "-Ou"],
                [
                    "bcftools",
                    "norm",
                    "--rm-dup",
                    "exact",
                    "-Ob",
                    "-W",
                    "-o",
                    str(source),
                ],
            ]
            subprocess.run(
                [
                    "bash",
                    "-o",
                    "pipefail",
                    "-c",
                    " | ".join(shlex.join(p) for p in pipeline),
                ],
                check=True,
            )
            output = directory / "direct.bcf"
            metrics = run_direct(
                task,
                source,
                output,
                directory / "measurement",
                manifest["freeze"]["build"]["root"],
                dict(os.environ, LC_ALL="C", **task["environment"]),
            )
            comparison = strict_semantic_compare(
                args.campaign / "attempts/a001" / identifier / "output.bcf", output
            )
            comparison["new_direct_digest"] = annotation_digest(
                output, task["baseline"]
            )
            comparison["historical_direct_digest"] = task["baseline"]["direct_digest"]
            comparison["direct_seconds"] = metrics["wall_seconds"]
            write_json_atomic(directory / "comparison.json", comparison)
            print(
                json.dumps(
                    dict(
                        cell=identifier,
                        semantic_pass=comparison["semantic_pass"],
                        direct_seconds=metrics["wall_seconds"],
                        examples=comparison["examples"][:3],
                    )
                ),
                flush=True,
            )


if __name__ == "__main__":
    main()
