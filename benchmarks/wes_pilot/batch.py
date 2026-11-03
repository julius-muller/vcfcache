"""Experimental local resident-cache session and equally batched direct control.

This is an API/session experiment, NOT single-sample CLI latency. Loading is
reported explicitly. Aggregate memory is session-level, never per-sample peak.
"""

import argparse
import json
import os
import shlex
import subprocess
import tempfile
import time
import traceback
from pathlib import Path

import yaml

from benchmarks.engine_profile.run import PYTHON
from benchmarks.performance_sprint.cache_stage import node_cache
from benchmarks.performance_sprint.manuscript import verified_file
from benchmarks.run_pilot import sha256sum, strict_semantic_compare, write_json_atomic
from vcfcache.database.workflow_manager import BcftoolsCommand, WorkflowManager


def manager(item, directory, output):
    """Construct the same recipe/parameters for either batch arm."""
    directory.mkdir(parents=True, exist_ok=True)
    return WorkflowManager(
        input_file=Path(item["prepared"]),
        output_dir=directory,
        name="batch",
        anno_config_file=Path(item["cache"]) / "annotation.yaml",
        params_file=Path(item["params"]),
        output_file=output,
        verbosity=0,
    )


def phase(prepared, kind, binary):
    """Run one whole session inside one fresh Slurm step."""
    items = json.loads(prepared.read_text())
    root = prepared.parent
    rows = []
    process = None
    started = time.perf_counter()
    load = None
    if kind == "resident-batch":
        jobs = []
        for item in items:
            folder = root / kind / item["id"]
            folder.mkdir(parents=True)
            os.mkfifo(folder / "input.pipe")
            jobs.append(
                dict(
                    input=str(folder / "input.pipe"),
                    cache=str(Path(item["cache"]) / "vcfcache_annotated.bcf"),
                    misses=str(folder / "misses.bcf"),
                    annotated=str(folder / "annotated.bcf"),
                    output=str(folder / "output.bcf"),
                    stats=str(folder / "native_stats.json"),
                    tag=item["tag"],
                    workers=4,
                    cpus=8,
                    memory_bytes=2 * 1024**3,
                    preserve=False,
                    lookup="regions",
                    window_bp=16384,
                )
            )
        write_json_atomic(root / "jobs.json", jobs)
        stderr = (root / "resident-batch.stderr").open("w")
        process = subprocess.Popen(
            [binary, "resident-batch", "--jobs", str(root / "jobs.json")],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=stderr,
            text=True,
        )

        def receive():
            line = process.stdout.readline()
            if not line:
                raise RuntimeError(
                    f"resident session ended ({process.wait()}): "
                    + (root / "resident-batch.stderr").read_text()[-2000:]
                )
            value = json.loads(line)
            with (root / "events.jsonl").open("a") as stream:
                stream.write(json.dumps(dict(unix=time.time(), **value)) + "\n")
            return value

        event = receive()
        assert event["event"] == "cache-ready"
        load = dict(**event, controller_elapsed_seconds=time.perf_counter() - started)
        write_json_atomic(root / "load.json", load)
    try:
        for index, item in enumerate(items):
            folder = root / kind / item["id"]
            if process:
                event = receive()
                assert event["event"] == "sample-ready" and event["index"] == index
            begin = time.perf_counter()
            output = folder / "output.bcf"
            obj = manager(item, folder / "manager", output)
            if kind == "direct-batch":
                obj.run(
                    "annotate-nocache",
                    db_bcf=Path(item["cache"]) / "vcfcache_annotated.bcf",
                )
            else:
                job = jobs[index]
                # Normalization is counted for both batch arms.
                command = (
                    "bcftools norm -m- "
                    + shlex.quote(item["prepared"])
                    + " -Ou > "
                    + shlex.quote(job["input"])
                )
                normalizer = subprocess.Popen(
                    ["bash", "-euo", "pipefail", "-c", command],
                    stderr=(folder / "normalize.stderr").open("w"),
                )
                try:
                    event = receive()
                    assert event["event"] == "misses-ready"
                    assert normalizer.wait() == 0
                    assert event["input_variants"] == item["records"]
                    if event["missing_variants"]:
                        auxiliary = folder / "auxiliary"
                        auxiliary.mkdir()
                        recipe = obj._substitute_variables(
                            obj.nfa_config_content["annotation_cmd"],
                            extra_vars={
                                "INPUT_BCF": job["misses"],
                                "OUTPUT_BCF": job["annotated"],
                                "AUXILIARY_DIR": str(auxiliary),
                            },
                        )
                        BcftoolsCommand(recipe, obj.logger, folder).run()
                        subprocess.run(
                            ["bcftools", "index", "-f", job["annotated"]], check=True
                        )
                        count = int(
                            subprocess.check_output(
                                ["bcftools", "index", "-n", job["annotated"]]
                            )
                        )
                        assert (
                            count == event["missing_variants"]
                        ), "annotation truncation"
                    process.stdin.write('{"action":"resume"}\n')
                    process.stdin.flush()
                    event = receive()
                    assert (
                        event["event"] == "complete"
                    ), "unsupported pairing/session output"
                finally:
                    if normalizer.poll() is None:
                        normalizer.terminate()
                        normalizer.wait()
            seconds = time.perf_counter() - begin
            assert (
                int(subprocess.check_output(["bcftools", "index", "-n", str(output)]))
                == item["records"]
            )
            rows.append(
                dict(
                    id=item["id"],
                    sample=item["sample"],
                    mode=kind,
                    wall_seconds=seconds,
                    records=item["records"],
                )
            )
            write_json_atomic(
                root / (kind + ".json"),
                dict(
                    rows=rows,
                    load=load,
                    scope="batch API; per-sample time includes manager, normalization, annotation and output; memory peak belongs to full session",
                ),
            )
        if process:
            assert receive()["event"] == "batch-complete" and process.wait() == 0
    finally:
        if process and process.poll() is None:
            process.terminate()
            process.wait()
    return rows


def driver(root, tool, deployment, binary):
    """Prepare outside timing; run both batch steps; compare all fields afterward."""
    manifest = json.loads((root / "manifest.json").read_text())
    targets = sorted(
        [t for t in manifest["all_targets"] if t["tool"] == tool],
        key=lambda t: t["sample"],
    )
    assert len(targets) == 12
    evidence = root / ("batch-" + tool)
    evidence.mkdir(exist_ok=False)
    state = dict(
        state="running",
        job_id=os.environ.get("SLURM_JOB_ID"),
        tool=tool,
        started_unix=time.time(),
    )
    write_json_atomic(evidence / "state.json", state)
    with tempfile.TemporaryDirectory(
        prefix="wes-batch-", dir="/mnt/data/tmp"
    ) as temporary:
        scratch = Path(temporary)
        try:
            items = []
            for task in targets:
                source = verified_file(task["input"], manifest["asset_attestations"])
                bed = verified_file(
                    task["preprocess"]["bed"], manifest["asset_attestations"]
                )
                cache = node_cache(task["cache"])
                for name, digest in task["cache"]["files"].items():
                    if task["cache"].get("storage") != "node_local":
                        verified_file(
                            dict(path=str(cache / name), sha256=digest),
                            manifest["asset_attestations"],
                        )
                for asset in (
                    task["annotation_assets"]
                    + manifest["asset_groups"][task["annotation_asset_group"]]
                ):
                    verified_file(asset, manifest["asset_attestations"])
                folder = scratch / task["id"]
                folder.mkdir()
                prepared = folder / "input.bcf"
                pipeline = f"bcftools view -R {shlex.quote(str(bed))} {shlex.quote(str(source))} -Ou | bcftools norm -m- -Ou | bcftools view -e 'ALT=\"*\"' -Ou | bcftools norm --rm-dup exact -Ob -o {prepared} -W --threads 8"
                subprocess.run(["bash", "-euo", "pipefail", "-c", pipeline], check=True)
                assert (
                    int(
                        subprocess.check_output(
                            ["bcftools", "index", "-n", str(prepared)]
                        )
                    )
                    == task["input"]["records"]
                )
                params = dict(
                    task["parameters"], threads=8, temp_dir=str(folder / "tmp")
                )
                (folder / "tmp").mkdir()
                (folder / "params.yaml").write_text(yaml.safe_dump(params))
                tag = yaml.safe_load((cache / "annotation.yaml").read_text())[
                    "must_contain_info_tag"
                ]
                items.append(
                    dict(
                        id=task["id"],
                        sample=task["sample"],
                        cache=str(cache),
                        prepared=str(prepared),
                        params=str(folder / "params.yaml"),
                        records=task["input"]["records"],
                        tag=tag,
                    )
                )
            assert len({i["cache"] for i in items}) == 1
            write_json_atomic(scratch / "prepared.json", items)
            for kind in ["direct-batch", "resident-batch"]:
                directory = evidence / kind
                directory.mkdir()
                command = [
                    "srun",
                    "--exclusive",
                    "--exact",
                    "--cpu-bind=cores",
                    "-n1",
                    "-c8",
                    PYTHON,
                    str(Path(deployment) / "benchmarks/engine_profile/measure.py"),
                    "--directory",
                    str(directory),
                    "--",
                    PYTHON,
                    "-m",
                    "benchmarks.wes_pilot.batch",
                    str(scratch / "prepared.json"),
                    "--phase",
                    kind,
                    "--binary",
                    binary,
                ]
                env = dict(
                    os.environ,
                    PYTHONPATH=deployment,
                    VCFCACHE_ROOT=deployment,
                    VCFCACHE_ENGINE=binary,
                )
                for key in list(env):
                    if key.startswith("VCFCACHE_EXPERIMENT_"):
                        env.pop(key)
                write_json_atomic(
                    directory / "command.json",
                    dict(
                        command=command,
                        tasks=targets,
                        binary_sha256=sha256sum(Path(binary)),
                    ),
                )
                with (directory / "run.log").open("w") as log:
                    subprocess.run(command, env=env, stdout=log, stderr=log, check=True)
                metrics = json.loads((directory / "resources.json").read_text())
                assert (
                    metrics["exit_code"] == 0
                    and len(metrics["affinity"]) == 8
                    and metrics["cgroup_after"]
                )
            for item, task in zip(items, targets, strict=True):
                direct = scratch / "direct-batch" / item["id"] / "output.bcf"
                cached = scratch / "resident-batch" / item["id"] / "output.bcf"
                compare = strict_semantic_compare(cached, direct)
                write_json_atomic(evidence / (item["id"] + "-comparison.json"), compare)
                assert compare["semantic_pass"]
                original = (
                    Path(task["baseline"]["direct_reuse"]["attempt"]) / "direct.bcf"
                )
                compare = strict_semantic_compare(direct, original)
                write_json_atomic(
                    evidence / (item["id"] + "-historical-direct.json"), compare
                )
                assert compare["semantic_pass"]
            state["state"] = "complete"
        except Exception as error:
            state.update(
                state="failed", error=str(error), traceback=traceback.format_exc()
            )
        finally:
            for name in [
                "prepared.json",
                "jobs.json",
                "load.json",
                "direct-batch.json",
                "resident-batch.json",
                "events.jsonl",
                "resident-batch.stderr",
            ]:
                if (scratch / name).exists():
                    (evidence / name).write_bytes((scratch / name).read_bytes())
            state["finished_unix"] = time.time()
            write_json_atomic(evidence / "state.json", state)
    return state


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--tool", choices=["veprs", "fastvep"])
    parser.add_argument("--phase", choices=["direct-batch", "resident-batch"])
    parser.add_argument("--deployment")
    parser.add_argument("--binary", required=True)
    args = parser.parse_args()
    if args.phase:
        phase(args.root, args.phase, args.binary)
    else:
        result = driver(args.root, args.tool, args.deployment, args.binary)
        print(json.dumps(result, indent=2))
        raise SystemExit(0 if result["state"] == "complete" else 1)
