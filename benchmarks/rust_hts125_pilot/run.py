"""Freeze, execute and collect a bounded Rust/HTSlib comparison, never a cohort rerun."""

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
from benchmarks.rust_hts125_pilot.common import (
    controller,
    definition,
    mapped,
    preflight_prepared,
    prepared_input,
    verify,
    wrapper,
)

PYTHON = "/results/deployment/vcfcache-native-pilot-20261023c/.venv/bin/python"
PARENT = Path(
    "/mnt/data/slurm-results/python-adaptive-manuscript-20261101-ee2162473/campaign/manifest.json"
)
INTEL = ["sl-w1", "sl-w2", "sl-w3", "sl-w4", "sl-w6"]
# (engine, external bcftools, linked helper, input preserved). Defaults unchanged.
CONDITIONS = {
    "adaptive-control121": ("adaptive", "old", "old", False),
    "adaptive-external125": ("adaptive", "latest", "old", False),
    "adaptive-linked125": ("adaptive", "latest", "latest", False),
    "adaptive-prepared125": ("adaptive", "latest", "latest", True),
    "resident-control121": ("native", "old", "old", False),
    "resident-external125": ("native", "latest", "old", False),
    "resident-linked125": ("native", "latest", "latest", False),
    "resident-prepared125": ("native", "latest", "latest", True),
    "direct-prepared125": ("legacy", "latest", "latest", True),
}


def modes_for(assay):
    """Resident input is WES-only here; never mislabel WGS fallback as native."""
    return [m for m in CONDITIONS if assay == "wes" or not m.startswith("resident-")]


def prepare(args):
    """Submit exactly six immutable blocks and an after-any review."""
    parent = json.loads(PARENT.read_text())
    selection = [
        ("veprs", "core", "KPGP-00139", "wes", "sl-w1"),
        ("fastvep", "clinical", "KPGP-00139", "wes", "sl-w2"),
        ("veprs", "core", "LP6005592-DNA_C05", "wes", "sl-w3"),
        ("fastvep", "clinical", "LP6005592-DNA_C05", "wes", "sl-w4"),
        ("veprs", "core", "KPGP-00139", "wgs", "sl-w6"),
        ("fastvep", "clinical", "KPGP-00139", "wgs", "sl-w1"),
    ]
    cases, tasks = [], []
    for i, (tool, level, sample, assay, node) in enumerate(selection):
        found = [
            t
            for t in parent["tasks"]
            if (t["tool"], t["level"], t["sample"], t["assay"], t["assembly"])
            == (tool, level, sample, assay, "GRCh38")
        ]
        assert len(found) == 1
        task = found[0]
        tasks.append(task)
        modes = modes_for(assay)
        offset = (i * 2) % len(modes)
        cases.append(
            dict(
                id=f"{tool}-{assay}-{sample}",
                task_id=task["id"],
                node=node,
                modes=modes[offset:] + modes[:offset],
            )
        )
    args.root.mkdir(parents=True, exist_ok=False)
    (args.root / "logs").mkdir()
    helpers = {}
    for label in ("old", "latest"):
        path = Path(args.helpers) / ("vcfcache-engine-" + label)
        actual = controller(path)
        info = json.loads(
            subprocess.check_output([str(actual), "build-info"], text=True)
        )
        assert info["htslib_version"] == ("1.19.1" if label == "old" else "1.25")
        assert info["libdeflate"]
        helpers[label] = dict(path=str(path), sha256=sha256sum(actual), info=info)
    assert (
        helpers["old"]["info"]["htslib_features"]
        == helpers["latest"]["info"]["htslib_features"]
    )
    build_file = controller(Path(args.helpers) / "BUILD.json")
    build = json.loads(build_file.read_text())
    for relative, digest in build["source_hashes"].items():
        verify(controller(args.deployment) / "native" / relative, digest)
    data = dict(
        tasks=tasks,
        cases=cases,
        cpus=8,
        repetitions=1,
        deployment=args.deployment,
        tools=args.tools,
        helpers=helpers,
        commit=args.commit,
        source_manifest_sha256=sha256sum(PARENT),
        excluded_nodes=["sl-w5"],
        conditions=CONDITIONS,
        native_build=build,
        native_build_sha256=sha256sum(build_file),
        asset_groups={
            k: v
            for k, v in parent["freeze"]["annotation_asset_groups"].items()
            if k in {"veprs", "fastvep"}
        },
        bcftools={
            label: dict(
                path=str(path),
                sha256=sha256sum(controller(path)),
                version=subprocess.check_output(
                    [str(controller(path)), "--version"], text=True
                ),
            )
            for label, path in (
                ("old", Path("/usr/local/bin/bcftools")),
                ("latest", Path(args.tools) / "bin/bcftools"),
            )
        },
        source_hashes={
            str(p.relative_to(controller(args.deployment))): sha256sum(p)
            for p in sorted(controller(args.deployment).rglob("*.py"))
        },
    )
    write_json_atomic(args.root / "manifest.json", data)
    (args.root / "READY").write_text(sha256sum(args.root / "manifest.json") + "\n")
    jobs, submitted = [], []
    for i, case in enumerate(cases):
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
            "--job-name=rust125-" + case["id"],
            "--output=" + mapped(args.root / "logs") + "/" + case["id"] + "-%j.out",
        ]
        if i == 5:
            command += ["--dependency=afterany:" + jobs[0]]
        command += [
            "--wrap",
            "PYTHONPATH="
            + shlex.quote(args.deployment)
            + " "
            + shlex.join(
                [
                    PYTHON,
                    "-m",
                    "benchmarks.rust_hts125_pilot.run",
                    "run",
                    mapped(args.root),
                    "--case",
                    case["id"],
                ]
            ),
        ]
        job = subprocess.check_output(command, text=True).strip().split(";")[0]
        jobs.append(job)
        submitted.append(dict(case=case, job_id=job, command=command))
        write_json_atomic(args.root / "submissions.json", submitted)
    command = [
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
        "--job-name=rust125-review",
        "--output=" + mapped(args.root / "logs/review-%j.out"),
        "--wrap",
        "PYTHONPATH="
        + shlex.quote(args.deployment)
        + " "
        + shlex.join(
            [
                PYTHON,
                "-m",
                "benchmarks.rust_hts125_pilot.run",
                "collect",
                mapped(args.root),
                "--terminal",
            ]
        ),
    ]
    review = subprocess.check_output(command, text=True).strip().split(";")[0]
    write_json_atomic(
        args.root / "review_submission.json", dict(job_id=review, command=command)
    )
    print(json.dumps(dict(jobs=jobs, review=review)))


def timed(data, evidence, mode, command, environment):
    """Time one complete CLI in a fresh eight-CPU Slurm step."""
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
        str(Path(data["deployment"]) / "benchmarks/rust_hts125_pilot/measure.py"),
        "--directory",
        str(directory),
        "--",
        *map(str, command),
    ]
    write_json_atomic(
        directory / "command.json",
        dict(
            command=launch,
            commit=data["commit"],
            environment={
                k: v
                for k, v in environment.items()
                if k
                in {
                    "PATH",
                    "PYTHONPATH",
                    "VCFCACHE_ROOT",
                    "VCFCACHE_ENGINE",
                    "VCFCACHE_PILOT_PREPARED",
                }
                or k.startswith("VCFCACHE_EXPERIMENT_")
            },
        ),
    )
    with (directory / "run.log").open("w") as log:
        subprocess.run(launch, env=environment, stdout=log, stderr=log, check=True)
    resources = json.loads((directory / "resources.json").read_text())
    assert resources["exit_code"] == 0 and len(resources["affinity"]) == 8
    return resources


def run(root, identifier):
    """Execute a same-worker block, retaining scratch and evidence on failure."""
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
    scratch = Path(tempfile.mkdtemp(prefix="rust-hts125-", dir="/mnt/data/tmp"))
    state["scratch"] = str(scratch)
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
        for relative, digest in data["source_hashes"].items():
            verify(Path(data["deployment"]) / relative, digest)
        for relative, digest in data["native_build"]["source_hashes"].items():
            verify(Path(data["deployment"]) / "native" / relative, digest)
        for tool in list(data["helpers"].values()) + list(data["bcftools"].values()):
            verify(tool["path"], tool["sha256"])
        task = next(t for t in data["tasks"] if t["id"] == case["task_id"])
        cache = Path(task["cache"]["directory"])
        for name, digest in task["cache"]["files"].items():
            verify(cache / name, digest)
        state["phase"] = "untimed_asset_verification"
        write_json_atomic(evidence / "state.json", state)
        for asset in (
            task["annotation_assets"]
            + data["asset_groups"][task["annotation_asset_group"]]
        ):
            verify(asset["path"], asset["sha256"])
        state["phase"] = "untimed_input_extraction"
        write_json_atomic(evidence / "state.json", state)
        source = prepared_input(task, scratch)
        config = yaml.safe_load((cache / "annotation.yaml").read_text())
        proof = preflight_prepared(
            data["bcftools"]["latest"]["path"], source, config["must_contain_info_tag"]
        )
        assert proof["records"] == task["input"]["records"]
        write_json_atomic(evidence / "prepared_preflight.json", proof)
        write_json_atomic(evidence / "task.json", task)
        source.chmod(0o444)
        state.update(input_sha256=sha256sum(source), records=proof["records"])
        binaries = {
            label: wrapper(
                scratch / (label + "-bin"), Path(tool["path"]), label == "latest"
            )
            for label, tool in data["bcftools"].items()
        }
        outputs = {}
        for mode in case["modes"]:
            state["phase"] = mode
            write_json_atomic(evidence / "state.json", state)
            engine, external, linked, prepared = CONDITIONS[mode]
            binary = binaries[external]
            outputs[mode] = scratch / (mode + ".bcf")
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
            environment = {
                k: v
                for k, v in os.environ.items()
                if not k.startswith("VCFCACHE_EXPERIMENT_")
                and k != "VCFCACHE_PROFILE_DIR"
            }
            environment.update(
                PYTHONPATH=data["deployment"],
                VCFCACHE_ROOT=data["deployment"],
                VCFCACHE_ENGINE=data["helpers"][linked]["path"],
                VCFCACHE_PILOT_PREPARED="1" if prepared else "0",
                PATH=str(binary.parent)
                + ":"
                + str(Path(data["tools"]) / "bin")
                + ":"
                + os.environ["PATH"],
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
                outputs[mode],
                "-y",
                params,
                "--engine",
                engine,
                "--statistics",
                "light",
                "--stats-dir",
                stats,
                "--preserve-unannotated",
                "-v",
            ]
            if mode.startswith("direct-"):
                command.append("--uncached")
            elif engine == "native":
                command += ["--engine-workers", "4", "--engine-memory-mb", "2048"]
            resources = timed(data, evidence, mode, command, environment)
            found = list(stats.rglob("compare_stats.yaml"))
            assert len(found) == 1
            counters = yaml.safe_load(found[0].read_text())
            actual = int(
                subprocess.check_output(
                    [str(binary), "index", "-n", str(outputs[mode])]
                )
            )
            assert (
                actual == proof["records"] == counters["variant_counts"]["total_output"]
            )
            if engine != "legacy":
                assert counters["effective_engine"] == (
                    "native" if engine == "native" else "adaptive"
                ), "Unexpected backend/fallback"
                assert bool(counters["pilot_prepared"]) == prepared
                assert not counters.get("fallback_reason")
            write_json_atomic(evidence / mode / "statistics.json", counters)
            state["rows"].append(
                dict(
                    mode=mode,
                    seconds=resources["wall_seconds"],
                    records=actual,
                    missing=counters["variant_counts"].get("missing_variants"),
                    effective_engine=counters.get("effective_engine", "legacy"),
                    workers=counters.get("workers"),
                    lookup_codec_threads=counters.get("lookup_codec_threads"),
                )
            )
            write_json_atomic(evidence / "state.json", state)
        state["phase"] = "untimed_semantic_audit"
        write_json_atomic(evidence / "state.json", state)
        with ProcessPoolExecutor(max_workers=2) as pool:
            pending = {
                mode: pool.submit(
                    strict_semantic_compare, output, outputs["direct-prepared125"]
                )
                for mode, output in outputs.items()
                if mode != "direct-prepared125"
            }
            for mode, future in pending.items():
                comparison = future.result()
                write_json_atomic(evidence / mode / "against_direct.json", comparison)
                assert comparison["semantic_pass"], f"Strict comparison failed: {mode}"
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


def collect(root, terminal=False):
    """Collect partial/terminal evidence without submitting more annotation."""
    data = definition(root)
    cases, rows = [], []
    for case in data["cases"]:
        file = root / "cases" / case["id"] / "state.json"
        state = (
            json.loads(file.read_text())
            if file.exists()
            else dict(state="not_started", rows=[])
        )
        if terminal and state["state"] not in {"complete", "failed"}:
            state = dict(state, state_before_terminal_review=state["state"])
            state["state"] = "interrupted"
            state["terminal_review_reason"] = (
                "After-any dependency confirms all jobs ended; runner did not record completion"
            )
        cases.append(dict(definition=case, **state))
        for entry in state["rows"]:
            file = root / "cases" / case["id"] / entry["mode"] / "against_direct.json"
            comparison = json.loads(file.read_text()) if file.exists() else {}
            rows.append(
                dict(
                    case=case["id"],
                    node=case["node"],
                    state=state["state"],
                    **entry,
                    semantic_pass=comparison.get("semantic_pass"),
                    mismatches=comparison.get("record_mismatches"),
                )
            )
    for row in rows:
        peers = {r["mode"]: r for r in rows if r["case"] == row["case"]}
        row["speedup_vs_prepared_direct"] = (
            peers["direct-prepared125"]["seconds"] / row["seconds"]
            if "direct-prepared125" in peers
            else None
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
        "# Rust / HTSlib upgrade pilot",
        "",
        "Intel workers only, fixed eight-CPU affinity, single-run engineering cells.",
        "External upgrade, linked-library upgrade and input preservation are separate conditions.",
        "Strict complete-record/header audits have no ignored fields. Ratios use prepared direct; control input preprocessing differs.",
        "",
        "| Case | State | " + " | ".join(CONDITIONS) + " |",
        "|---|---|" + "---:|" * len(CONDITIONS),
    ]
    for case in cases:
        identifier = case["definition"]["id"]
        values = {r["mode"]: r for r in rows if r["case"] == identifier}
        report.append(
            "| "
            + identifier
            + " | "
            + case["state"]
            + " | "
            + " | ".join(
                f"{values[m]['seconds']:.2f}" if m in values else "NA"
                for m in CONDITIONS
            )
            + " |"
        )
    report += [
        "",
        "Times are seconds, not fold speedups. Require semantic_pass and actual backend before quoting performance.",
        "Resident WGS is deliberately not attempted: an adaptive fallback would not measure resident performance.",
        "Prepared input is an experimental externally verified domain, not a generic product assertion.",
        "No full-cohort jobs or product-default changes are triggered automatically.",
    ]
    (root / "RESULTS.md").write_text("\n".join(report) + "\n")
    if all(c["state"] in {"complete", "failed", "interrupted"} for c in cases):
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
                    "conditions.tsv",
                    "REVIEW.json",
                    "RESULTS.md",
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
    print(json.dumps({c["definition"]["id"]: c["state"] for c in cases}))


def main():
    """Prepare, run or review only the explicitly requested bounded pilot."""
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    p = sub.add_parser("prepare")
    p.add_argument("root", type=Path)
    for name in ("deployment", "tools", "helpers", "commit"):
        p.add_argument("--" + name, required=True)
    p = sub.add_parser("run")
    p.add_argument("root", type=Path)
    p.add_argument("--case", required=True)
    p = sub.add_parser("collect")
    p.add_argument("root", type=Path)
    p.add_argument("--terminal", action="store_true")
    args = parser.parse_args()
    if args.action == "prepare":
        prepare(args)
    elif args.action == "run":
        run(args.root, args.case)
    else:
        collect(args.root, terminal=args.terminal)


if __name__ == "__main__":
    main()
