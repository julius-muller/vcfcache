"""Bounded Python-adaptive HTSlib 1.25 factorial engineering pilot."""

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
MODES = [
    "baseline-121",
    "baseline-125",
    "prepared-auto125",
    "prepared-contigs2",
    "prepared-windows2",
    "prepared-windows4",
    "direct-125",
]
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
    """Submit five end-to-end blocks and their review, excluding sl-w5."""
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
                offset=[0, 3, 6, 2, 5][index],
            )
        )
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
        baseline_commit="8ec84894f27ede92676c3d4a5d193f765086b77d",
        parent_sha256=sha256sum(PARENT),
        excluded_nodes=["sl-w5"],
        asset_groups={
            k: v
            for k, v in parent["freeze"]["annotation_asset_groups"].items()
            if k in {"veprs", "fastvep"}
        },
        tools_sha256=sha256sum(controller(args.tools) / "bin/bcftools"),
        old_tools_sha256=sha256sum(Path("/usr/local/bin/bcftools")),
        tool_versions={
            label: subprocess.check_output(
                [str(binary), "--version"], text=True
            ).strip()
            for label, binary in (
                ("baseline", Path("/usr/local/bin/bcftools")),
                ("latest", controller(args.tools) / "bin/bcftools"),
            )
        },
    )
    data["deployment_hashes"] = {}
    for label, directory in (
        ("prototype", args.deployment),
        ("reference", args.baseline),
    ):
        folder = controller(directory)
        data["deployment_hashes"][label] = {
            str(p.relative_to(folder)): sha256sum(p)
            for p in sorted(folder.rglob("*.py"))
        }
    write_json_atomic(args.root / "manifest.json", data)
    (args.root / "READY").write_text(sha256sum(args.root / "manifest.json") + "\n")
    submitted, jobs = [], []
    for case in cases:
        launch = [
            PYTHON,
            "-m",
            "benchmarks.adaptive_hts125_pilot.run",
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
            "--job-name=adaptive-hts125-" + case["id"],
            "--output=" + mapped(args.root / "logs") + "/" + case["id"] + "-%j.out",
        ]
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
        "--job-name=adaptive-hts125-review",
        "--output=" + mapped(args.root / "logs/review-%j.out"),
        "--wrap",
        "PYTHONPATH="
        + shlex.quote(args.deployment)
        + " "
        + shlex.join(
            [
                PYTHON,
                "-m",
                "benchmarks.adaptive_hts125_pilot.run",
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


def preflight_prepared(bcftools, source, tag):
    """Verify the experimental trusted-input assertion outside all timings."""
    header = subprocess.check_output([bcftools, "view", "-h", str(source)], text=True)
    has_tag = f"##INFO=<ID={tag}," in header
    fmt = "%CHROM\t%POS\t%REF\t%ALT" + (f"\t%INFO/{tag}" if has_tag else "") + "\n"
    locus, seen, records = None, set(), 0
    with tempfile.TemporaryFile() as errors:
        process = subprocess.Popen(
            [bcftools, "query", "-f", fmt, str(source)],
            stdout=subprocess.PIPE,
            stderr=errors,
            text=True,
        )
        assert process.stdout is not None
        try:
            with process.stdout:
                for line in process.stdout:
                    fields = line.rstrip("\n").split("\t")
                    contig, pos, ref, alt = fields[:4]
                    if "," in alt or alt in {"*", ".", ""}:
                        raise ValueError(
                            "Pilot requires already-biallelic, star-free input; do not silently transform it"
                        )
                    if has_tag and fields[4] not in {"", "."}:
                        raise ValueError(
                            "Existing input annotation could masquerade as a cache hit"
                        )
                    if locus != (contig, pos):
                        locus, seen = (contig, pos), set()
                    if (ref, alt) in seen:
                        raise ValueError(
                            "Duplicate variant keys are outside this prepared-input pilot"
                        )
                    seen.add((ref, alt))
                    records += 1
            if process.wait():
                errors.seek(0)
                raise RuntimeError(errors.read().decode(errors="replace"))
        finally:
            if process.poll() is None:
                process.terminate()
                process.wait()
    return dict(
        records=records,
        biallelic=True,
        star_free=True,
        unique_keys=True,
        unannotated=True,
    )


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
        str(Path(data["deployment"]) / "benchmarks/adaptive_hts125_pilot/measure.py"),
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
                for k in (
                    "PATH",
                    "PYTHONPATH",
                    "VCFCACHE_PILOT_PREPARED",
                    "VCFCACHE_PILOT_LOOKUP",
                    "VCFCACHE_PILOT_THREAD_BUDGET",
                    "VCFCACHE_PILOT_MERGE_WINDOWS",
                )
            },
            commit=(
                data["baseline_commit"]
                if mode.startswith("baseline-")
                else data["commit"]
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
    scratch = Path(tempfile.mkdtemp(prefix="adaptive-hts125-", dir="/mnt/data/tmp"))
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
        for label, folder in (
            ("prototype", data["deployment"]),
            ("reference", data["baseline"]),
        ):
            for relative, digest in data["deployment_hashes"][label].items():
                verify(Path(folder) / relative, digest)
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
        config = yaml.safe_load((cache / "annotation.yaml").read_text())
        preflight = preflight_prepared(
            str(Path(data["tools"]) / "bin/bcftools"),
            source,
            config["must_contain_info_tag"],
        )
        assert preflight["records"] == records
        write_json_atomic(evidence / "prepared_preflight.json", preflight)
        source.chmod(0o444)
        old = wrapper(scratch / "old-bin", Path("/usr/local/bin/bcftools"), False)
        new = wrapper(scratch / "new-bin", Path(data["tools"]) / "bin/bcftools", True)
        outputs = {}
        shift = case["offset"] % len(MODES)
        order = MODES[shift:] + MODES[:shift]
        for mode in order:
            state["phase"] = mode
            write_json_atomic(evidence / "state.json", state)
            original = mode.startswith("baseline-")
            binary = old if mode == "baseline-121" else new
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
                VCFCACHE_PILOT_REGION_MODE="off",
                VCFCACHE_PILOT_PREPARED="0" if original else "1",
                VCFCACHE_PILOT_LOOKUP=(
                    "windows"
                    if "windows" in mode
                    else ("contigs" if "contigs" in mode else "auto")
                ),
                VCFCACHE_PILOT_THREAD_BUDGET="1" if mode[-1:] in {"2", "4"} else "0",
                VCFCACHE_PILOT_MERGE_WINDOWS="1" if "windows" in mode else "0",
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
            command += ["--engine", "legacy" if mode == "direct-125" else "adaptive"]
            if mode == "direct-125":
                command.append("--uncached")
            elif mode[-1:] in {"2", "4"}:
                command += ["--engine-workers", mode[-1]]
            resources = timed(data, evidence, mode, command, environment)
            found = list(stats.rglob("compare_stats.yaml"))
            assert len(found) == 1
            counters = yaml.safe_load(found[0].read_text())
            assert counters["variant_counts"]["total_output"] == records
            if not original and mode != "direct-125":
                assert counters["pilot_prepared"] and counters["engine"] == "adaptive"
                if mode[-1:] in {"2", "4"}:
                    assert (
                        counters["shards"]
                        * (1 + counters["pilot_codec_threads_per_shard"])
                        <= 8
                    )
            write_json_atomic(evidence / mode / "statistics.json", counters)
            state["rows"].append(
                dict(
                    mode=mode,
                    seconds=resources["wall_seconds"],
                    records=records,
                    missing=counters["variant_counts"].get("missing_variants"),
                    shards=counters.get("shards"),
                    cache_access=counters.get("cache_access"),
                    codec_threads=counters.get("pilot_codec_threads_per_shard"),
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
                    outputs["direct-125"],
                )
                for mode in MODES
                if mode != "direct-125"
            }
            for mode, future in pending.items():
                comparison = future.result()
                write_json_atomic(evidence / mode / "against_direct.json", comparison)
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
    for row in rows:
        peers = {r["mode"]: r for r in rows if r["case"] == row["case"]}
        row["speedup_vs_prepared_direct"] = (
            peers["direct-125"]["seconds"] / row["seconds"]
            if "direct-125" in peers
            else None
        )
        row["gain_over_baseline121"] = (
            peers["baseline-121"]["seconds"] / row["seconds"]
            if "baseline-121" in peers
            else None
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
        "# Python-adaptive / HTSlib 1.25 pilot",
        "",
        "Same-worker, eight-CPU, single-run engineering cells; sl-w5 excluded.",
        "Baseline 1.21/1.25 isolates the external toolchain. Prepared variants are preflight verified.",
        "Full record/header audits are untimed. Known clinical discrepancies remain unvalidated.",
        "",
        "| Case | State | Direct | Baseline 1.21 | Baseline 1.25 | Prepared auto | Contigs 2 | Windows 2 | Windows 4 |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for case in cases:
        name = case["definition"]["id"]
        values = {r["mode"]: r for r in rows if r["case"] == name}
        times = [
            f"{values[m]['seconds']:.2f}" if m in values else "NA"
            for m in ["direct-125"] + MODES[:-1]
        ]
        report.append(f"| {name} | {case['state']} | " + " | ".join(times) + " |")
    report += [
        "",
        "Times are seconds, not fold speedups. Consult semantic_pass/mismatches in conditions.tsv before quoting gains.",
        "Prepared-input assertions and full external comparisons are outside timing. Region construction and normal cache work are timed.",
        "These trusted-input experimental switches are NOT a generic production replacement for main.",
        "No full-cohort or follow-up jobs are submitted automatically.",
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
                        / "benchmarks/adaptive_hts125_pilot/render.R"
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
