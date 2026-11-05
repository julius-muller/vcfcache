"""Pinned Intel-only pilot; durable state and untimed strict semantic audits."""

import argparse
import csv
import json
import os
import shlex
import shutil
import socket
import subprocess
import tarfile
import tempfile
import time
import traceback
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import yaml

from benchmarks.run_pilot import sha256sum, strict_semantic_compare, write_json_atomic

PYTHON = "/results/deployment/vcfcache-native-pilot-20261023c/.venv/bin/python"
PARENT = Path(
    "/mnt/data/slurm-results/python-adaptive-manuscript-20261101-ee2162473/campaign/manifest.json"
)
MODES = ["old-direct", "old-cached", "new-direct", "new-cached", "new-regions"]
INTEL = ["sl-w1", "sl-w2", "sl-w3", "sl-w4", "sl-w6"]


def mapped(path):
    """Controller and worker spellings of the shared mount."""
    return str(path).replace("/mnt/data/slurm-results/", "/results/", 1)


def controller(path):
    """Resolve a worker-visible path on the controller."""
    return Path(str(path).replace("/results/", "/mnt/data/slurm-results/", 1))


def definition(root):
    """Require the frozen definition, never an edited live job recipe."""
    assert sha256sum(root / "manifest.json") == (root / "READY").read_text().strip()
    return json.loads((root / "manifest.json").read_text())


def verify(path, digest):
    """Verify immutable inputs/cache payloads without editing them."""
    path = Path(path)
    assert path.exists() and sha256sum(path) == digest, f"Asset changed: {path}"
    return path


def prepare(args):
    """Submit five end-to-end blocks and one reader control, excluding sl-w5."""
    parent = json.loads(PARENT.read_text())
    selected = [
        ("veprs", "core", "wes"),
        ("veprs", "core", "wgs"),
        ("fastvep", "clinical", "wes"),
        ("fastvep", "clinical", "wgs"),
        ("veprs", "clinical", "wes"),
    ]
    tasks, cases = [], []
    for index, ((tool, level, assay), node) in enumerate(
        zip(selected, INTEL, strict=True)
    ):
        matches = [
            t
            for t in parent["tasks"]
            if t["sample"] == "KPGP-00139"
            and (t["tool"], t["level"], t["assay"]) == (tool, level, assay)
            and t["assembly"] == "GRCh38"
        ]
        assert len(matches) == 1
        tasks.append(matches[0])
        cases.append(
            dict(
                id=f"{tool}-{level}-{assay}",
                task_id=matches[0]["id"],
                node=node,
                offset=index,
            )
        )
    cases.append(dict(id="reader", task_id=tasks[0]["id"], node="sl-w1", offset=0))
    args.root.mkdir(parents=True, exist_ok=False)
    (args.root / "logs").mkdir()
    data = dict(
        tasks=tasks,
        cases=cases,
        cpus=8,
        repetitions=1,
        modes=MODES,
        deployment=args.deployment,
        baseline=args.baseline,
        tools=args.tools,
        commit=args.commit,
        baseline_commit="136bafa90",
        parent_sha256=sha256sum(PARENT),
        excluded_nodes=["sl-w5"],
        asset_groups=parent["freeze"]["annotation_asset_groups"],
        tools_sha256=sha256sum(controller(args.tools) / "bin/bcftools"),
        old_tools_sha256=sha256sum(Path("/usr/local/bin/bcftools")),
    )
    write_json_atomic(args.root / "manifest.json", data)
    (args.root / "READY").write_text(sha256sum(args.root / "manifest.json") + "\n")
    submitted, jobs = [], []
    for case in cases:
        launch = [
            PYTHON,
            "-m",
            "benchmarks.input_preserving_pilot.run",
            "run",
            mapped(args.root),
            "--case",
            case["id"],
        ]
        command = [
            "sbatch",
            "--parsable",
            "--partition=benchmark",
            "--exclusive",
            "--nodes=1",
            "--ntasks=1",
            "--cpus-per-task=8",
            "--mem=28000M",
            "--time=02:00:00",
            "--chdir=/tmp",
            "--nodelist=" + case["node"],
            "--job-name=input-preserving-" + case["id"],
            "--output=" + mapped(args.root / "logs") + "/" + case["id"] + "-%j.out",
        ]
        if case["id"] == "reader":
            command.append("--dependency=afterany:" + jobs[0])
        command += [
            "--wrap",
            "PYTHONPATH=" + shlex.quote(args.deployment) + " " + shlex.join(launch),
        ]
        job = subprocess.check_output(command, text=True).strip().split(";")[0]
        jobs.append(job)
        submitted.append(dict(case=case, job_id=job, command=command))
        write_json_atomic(args.root / "submissions.json", submitted)
    review = [
        "sbatch",
        "--parsable",
        "--partition=benchmark",
        "--ntasks=1",
        "--cpus-per-task=1",
        "--mem=3000M",
        "--time=00:10:00",
        "--chdir=/tmp",
        "--exclude=sl-w5",
        "--dependency=afterany:" + ":".join(jobs),
        "--job-name=input-preserving-review",
        "--output=" + mapped(args.root / "logs/review-%j.out"),
        "--wrap",
        "PYTHONPATH="
        + shlex.quote(args.deployment)
        + " "
        + shlex.join(
            [
                PYTHON,
                "-m",
                "benchmarks.input_preserving_pilot.run",
                "collect",
                mapped(args.root),
            ]
        ),
    ]
    job = subprocess.check_output(review, text=True).strip().split(";")[0]
    write_json_atomic(
        args.root / "review_submission.json", dict(job_id=job, command=review)
    )
    print(json.dumps(dict(jobs=jobs, review=job)))


def wrapper(directory, actual, latest):
    """Thread every supported stage, including commands inside frozen recipes."""
    directory.mkdir()
    destination = directory / "bcftools"
    script = f"""#!/bin/bash
set -euo pipefail
actual={shlex.quote(str(actual))}
subcommand=${{1-}}
case "$subcommand" in
  annotate|view|filter|concat|index{'|sort' if latest else ''})
    shift
    for arg in "$@"; do
      case "$arg" in --threads|--threads=*|-@*) exec "$actual" "$subcommand" "$@";; esac
    done
    exec "$actual" "$subcommand" "$@" --threads 8;;
esac
exec "$actual" "$@"
"""
    destination.write_text(script)
    destination.chmod(0o755)
    return destination


def prepared_input(task, scratch):
    """Assay extraction/BCF encoding only; never split or remove any variants."""
    source = verify(task["input"]["path"], task["input"]["sha256"])
    output = scratch / (task["assay"] + ".bcf")
    command = [
        "/usr/local/bin/bcftools",
        "view",
        str(source),
        "-Ob1",
        "-W",
        "--threads",
        "8",
        "-o",
        str(output),
    ]
    if task.get("preprocess", {}).get("bed"):
        bed = task["preprocess"]["bed"]
        command += ["-R", str(verify(bed["path"], bed["sha256"]))]
    subprocess.run(command, check=True)
    return output


def timed(data, evidence, mode, command, environment):
    """Each measured cell has a fresh eight-CPU cgroup/peak-memory counter."""
    directory = evidence / mode
    directory.mkdir()
    launch = [
        "srun",
        "--exclusive",
        "--exact",
        "--cpu-bind=cores",
        "-n1",
        "-c8",
        PYTHON,
        str(Path(data["deployment"]) / "benchmarks/input_preserving_pilot/measure.py"),
        "--directory",
        str(directory),
        "--",
        *map(str, command),
    ]
    write_json_atomic(
        directory / "command.json",
        dict(
            command=launch,
            environment={
                k: environment[k]
                for k in ("PATH", "PYTHONPATH", "VCFCACHE_PILOT_REGION_MODE")
            },
            commit=(
                data["baseline_commit"] if mode.startswith("old-") else data["commit"]
            ),
        ),
    )
    with (directory / "run.log").open("w") as log:
        subprocess.run(launch, env=environment, stdout=log, stderr=log, check=True)
    resources = json.loads((directory / "resources.json").read_text())
    assert resources["exit_code"] == 0 and len(resources["affinity"]) == 8
    return resources


def run(root, identifier):
    """One same-node block; no full-cohort or best-retry automation."""
    data = definition(root)
    case = next(c for c in data["cases"] if c["id"] == identifier)
    evidence = root / "cases" / identifier
    evidence.mkdir(parents=True, exist_ok=False)
    state = dict(
        state="running",
        node=socket.gethostname(),
        job_id=os.environ.get("SLURM_JOB_ID"),
        rows=[],
        started=time.time(),
    )
    write_json_atomic(evidence / "state.json", state)
    scratch = Path(tempfile.mkdtemp(prefix="input-preserving-", dir="/mnt/data/tmp"))
    state["scratch"] = str(scratch)  # Retain on failure/time limit for audit recovery.
    try:
        assert state["node"] == case["node"] and state["node"] in INTEL
        hardware = dict(
            lscpu=json.loads(subprocess.check_output(["lscpu", "-J"], text=True)),
            storage=json.loads(
                subprocess.check_output(
                    ["findmnt", "-J", "-T", str(scratch)], text=True
                )
            ),
            free_bytes=shutil.disk_usage(scratch).free,
        )
        assert any(
            x["field"].rstrip(":") == "Model name" and "Broadwell" in x["data"]
            for x in hardware["lscpu"]["lscpu"]
        )
        assert hardware["free_bytes"] > 30 * 1024**3
        write_json_atomic(evidence / "hardware.json", hardware)
        verify("/usr/local/bin/bcftools", data["old_tools_sha256"])
        verify(Path(data["tools"]) / "bin/bcftools", data["tools_sha256"])
        task = next(t for t in data["tasks"] if t["id"] == case["task_id"])
        cache = Path(task["cache"]["directory"])
        for name, digest in task["cache"]["files"].items():
            verify(cache / name, digest)
        state["phase"] = "verify_annotation_assets"
        write_json_atomic(evidence / "state.json", state)
        for asset in (
            task["annotation_assets"]
            + data["asset_groups"][task["annotation_asset_group"]]
        ):
            verify(asset["path"], asset["sha256"])
        write_json_atomic(evidence / "task.json", task)
        source = prepared_input(task, scratch)
        records = int(subprocess.check_output(["bcftools", "index", "-n", str(source)]))
        assert (
            records == task["input"]["records"]
        ), "Extraction differs from frozen assay; do not silently normalize it"
        state["input_sha256"] = sha256sum(source)
        state["records"] = records
        old = wrapper(scratch / "old-bin", Path("/usr/local/bin/bcftools"), False)
        new = wrapper(scratch / "new-bin", Path(data["tools"]) / "bin/bcftools", True)
        if identifier == "reader":
            reader(data, evidence, state, scratch, source, cache, old, new, "wes")
            whole = next(
                t
                for t in data["tasks"]
                if (t["tool"], t["level"], t["assay"]) == ("veprs", "core", "wgs")
            )
            whole_input = prepared_input(whole, scratch)
            whole_hash = sha256sum(whole_input)
            reader(data, evidence, state, scratch, whole_input, cache, old, new, "wgs")
            assert sha256sum(whole_input) == whole_hash
            state["reader_wgs_sha256"] = whole_hash
        else:
            outputs = {}
            shift = case["offset"] % len(MODES)
            order = MODES[shift:] + MODES[:shift]
            for mode in order:
                state["phase"] = mode
                write_json_atomic(evidence / "state.json", state)
                original = mode.startswith("old-")
                binary = old if original else new
                deployment = data["baseline"] if original else data["deployment"]
                output = scratch / (mode + ".bcf")
                outputs[mode] = output
                stats = scratch / (mode + "-stats")
                parameters = dict(
                    task["parameters"],
                    bcftools_cmd=str(binary),
                    threads=8,
                    temp_dir=str(scratch / "tmp"),
                )
                (scratch / "tmp").mkdir(exist_ok=True)
                params = scratch / (mode + ".yaml")
                params.write_text(yaml.safe_dump(parameters))
                environment = dict(
                    os.environ,
                    PYTHONPATH=deployment,
                    VCFCACHE_ROOT=deployment,
                    PATH=str(binary.parent)
                    + ":"
                    + str(Path(data["tools"]) / "bin")
                    + ":"
                    + os.environ["PATH"],
                    VCFCACHE_PILOT_REGION_MODE=(
                        "windows" if mode == "new-regions" else "off"
                    ),
                    LC_ALL="C",
                )
                command = [
                    PYTHON,
                    "-c",
                    "from vcfcache.cli import main; main()",
                    "annotate",
                    "-a",
                    cache,
                    "-i",
                    source,
                    "-o",
                    output,
                    "-y",
                    params,
                    "--statistics",
                    "light",
                    "--stats-dir",
                    stats,
                    "--preserve-unannotated",
                    "-v",
                ]
                if mode.endswith("direct"):
                    command.append("--uncached")
                resources = timed(data, evidence, mode, command, environment)
                found = list(stats.rglob("compare_stats.yaml"))
                assert len(found) == 1
                counters = yaml.safe_load(found[0].read_text())
                assert counters["variant_counts"]["total_output"] == records
                if mode in {"new-cached", "new-regions"}:
                    assert (
                        counters["cache_matching"] == "biallelic-exact"
                    ), "Do not count a duplicate-key direct bypass as cached timing"
                write_json_atomic(evidence / mode / "statistics.json", counters)
                state["rows"].append(
                    dict(
                        mode=mode,
                        seconds=resources["wall_seconds"],
                        records=records,
                        missing=counters["variant_counts"].get("missing_variants"),
                    )
                )
                write_json_atomic(evidence / "state.json", state)
            # No comparison reads between timed conditions. Known ClinVar
            # discrepancies are reported, not ignored or promoted as validated.
            state["phase"] = "untimed_semantic_audit"
            write_json_atomic(evidence / "state.json", state)
            with ProcessPoolExecutor(max_workers=2) as pool:
                pending = {
                    mode: pool.submit(
                        strict_semantic_compare,
                        outputs[mode],
                        outputs["old-direct" if mode == "old-cached" else "new-direct"],
                    )
                    for mode in ("old-cached", "new-cached", "new-regions")
                }
                for mode, future in pending.items():
                    comparison = future.result()
                    write_json_atomic(
                        evidence / mode / "against_direct.json", comparison
                    )
            write_json_atomic(
                evidence / "direct_versions.json",
                strict_semantic_compare(outputs["new-direct"], outputs["old-direct"]),
            )
        assert sha256sum(source) == state["input_sha256"]
        state.update(state="complete", finished=time.time())
    except Exception as error:
        state.update(
            state="failed",
            error=str(error),
            traceback=traceback.format_exc(),
            finished=time.time(),
        )
        raise
    finally:
        write_json_atomic(evidence / "state.json", state)
    # Outputs remain available locally until explicitly cleaned, even on success.


def reader(data, evidence, state, scratch, source, cache, old, new, assay):
    """Same-CPU reader experiment, no annotator and no output compression."""
    from vcfcache.database.preserving import inspect_input

    bed = scratch / f"{assay}-windows.bed"
    inspect_input(str(new), source, bed)
    occupied = scratch / f"{assay}-contigs.bed"
    text = subprocess.check_output(
        [str(new), "index", "--stats", str(source)], text=True
    )
    occupied.write_text(
        "".join(
            f"{v[0]}\t0\t{v[1]}\n"
            for line in text.splitlines()
            if (v := line.split("\t"))
        )
    )
    for version, binary in (("1.21", old), ("1.25", new)):
        for regions in (None, occupied, bed):
            mode = (
                f"reader-{assay}-"
                + version
                + "-"
                + (regions.stem if regions else "unrestricted")
            )
            command = [
                str(binary),
                "annotate",
                "--pair-logic",
                "exact",
                "-a",
                str(cache / "vcfcache_annotated.bcf"),
                "-c",
                "INFO",
                str(source),
                "-Ou",
                "-o",
                os.devnull,
            ]
            if regions:
                command += ["-R", str(regions), "--regions-overlap", "0"]
            env = dict(
                os.environ,
                PYTHONPATH=data["deployment"],
                PATH=str(binary.parent) + ":" + os.environ["PATH"],
                VCFCACHE_PILOT_REGION_MODE="off",
            )
            metrics = timed(data, evidence, mode, command, env)
            state["rows"].append(
                dict(
                    mode=mode,
                    seconds=metrics["wall_seconds"],
                    records=int(
                        subprocess.check_output(
                            [str(binary), "index", "-n", str(source)]
                        )
                    ),
                )
            )
            write_json_atomic(evidence / "state.json", state)


def collect(root):
    """Review all terminal and partial evidence; never enqueue another campaign."""
    data = definition(root)
    cases, rows = [], []
    for case in data["cases"]:
        path = root / "cases" / case["id"] / "state.json"
        state = (
            json.loads(path.read_text())
            if path.exists()
            else dict(state="not_started", rows=[])
        )
        cases.append(dict(definition=case, **state))
        for row in state["rows"]:
            comparison = (
                root / "cases" / case["id"] / row["mode"] / "against_direct.json"
            )
            semantic = json.loads(comparison.read_text()) if comparison.exists() else {}
            rows.append(
                dict(
                    case=case["id"],
                    node=case["node"],
                    state=state["state"],
                    **row,
                    semantic_pass=semantic.get("semantic_pass"),
                    mismatches=semantic.get("record_mismatches"),
                )
            )
    write_json_atomic(
        root / "REVIEW.json",
        dict(
            cases=cases,
            reviewed=time.time(),
            complete=sum(c["state"] == "complete" for c in cases),
        ),
    )
    if rows:
        with (root / "conditions.tsv").open("w") as stream:
            writer = csv.DictWriter(
                stream,
                fieldnames=list(dict.fromkeys(k for r in rows for k in r)),
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerows(rows)
    report = [
        "# Input-preserving main / HTSlib 1.25 pilot",
        "",
        "Single-run engineering observations; not a new cohort estimate.",
        "All pairs are same-worker, eight-CPU measurements. sl-w5 is excluded.",
        "Code and toolchain changes are combined in old/new end-to-end comparisons.",
        "A completed block is not necessarily annotation-equivalent; consult its comparisons.",
        "",
        "| Case | State | Original direct/cached (s) | New direct/cached/regions (s) | Strict new cached / regions |",
        "|---|---|---|---|---|",
    ]
    for case in cases:
        name = case["definition"]["id"]
        if name == "reader":
            continue
        values = {r["mode"]: r for r in rows if r["case"] == name}

        def elapsed(mode, values=values):
            return f"{values[mode]['seconds']:.2f}" if mode in values else "NA"

        def passed(mode, values=values):
            return str(values.get(mode, {}).get("semantic_pass", "pending"))

        report.append(
            f"| {name} | {case['state']} | {elapsed('old-direct')} / {elapsed('old-cached')} | "
            f"{elapsed('new-direct')} / {elapsed('new-cached')} / {elapsed('new-regions')} | "
            f"{passed('new-cached')} / {passed('new-regions')} |"
        )
    report += [
        "",
        "Known vep-rs clinical annotation discrepancies are not ignored.",
        "Input preprocessing is absent from both new arms; validation overhead is timed.",
        "Reader-only cells have no annotator or compressed-output cost and are not end-to-end speedups.",
        "See conditions.tsv, cases/*/state.json, comparisons, command/resource records and the frozen manifest.",
        "No follow-up campaign is submitted automatically.",
    ]
    (root / "RESULTS.md").write_text("\n".join(report) + "\n")
    if all(c["state"] in {"complete", "failed"} for c in cases):
        archive = root / "terminal-evidence.tar.gz"
        if not archive.exists():
            with tarfile.open(archive, "x:gz") as output:
                for name in (
                    "manifest.json",
                    "READY",
                    "submissions.json",
                    "review_submission.json",
                    "cases",
                    "logs",
                ):
                    if (root / name).exists():
                        output.add(root / name, arcname=name)
        write_json_atomic(
            root / "SNAPSHOT.json",
            dict(
                archive_sha256=sha256sum(archive),
                manifest_sha256=sha256sum(root / "manifest.json"),
            ),
        )
        if shutil.which("Rscript") and not (root / "figures-v1").exists():
            rendering = subprocess.run(
                [
                    "Rscript",
                    "--vanilla",
                    str(
                        Path(data["deployment"])
                        / "benchmarks/input_preserving_pilot/render.R"
                    ),
                    str(root),
                    "figures-v1",
                ],
                capture_output=True,
                text=True,
                env=dict(
                    os.environ, R_LIBS_USER=str(Path(data["tools"]) / "R-library")
                ),
            )
            (root / "render.log").write_text(rendering.stdout + rendering.stderr)
            if rendering.returncode:
                print(
                    "R rendering needs review; timings and comparisons remain available"
                )
    print(json.dumps({c["definition"]["id"]: c["state"] for c in cases}))


def main():
    """Prepare, execute or review only the requested bounded experiment."""
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    prepare_parser = sub.add_parser("prepare")
    prepare_parser.add_argument("root", type=Path)
    for name in ("deployment", "baseline", "tools", "commit"):
        prepare_parser.add_argument("--" + name, required=True)
    run_parser = sub.add_parser("run")
    run_parser.add_argument("root", type=Path)
    run_parser.add_argument("--case", required=True)
    collect_parser = sub.add_parser("collect")
    collect_parser.add_argument("root", type=Path)
    args = parser.parse_args()
    if args.action == "prepare":
        prepare(args)
    elif args.action == "run":
        run(args.root, args.case)
    else:
        collect(args.root)


if __name__ == "__main__":
    main()
