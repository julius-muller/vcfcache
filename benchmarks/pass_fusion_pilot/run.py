"""Freeze and run a bounded same-node A/B, with fresh direct controls."""

import argparse
import json
import os
import shlex
import shutil
import socket
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

MODES = ["direct", "control-auto", "control4", "fused4", "fused2"]
PARENT = Path(
    "/mnt/data/slurm-results/python-adaptive-manuscript-20261101-ee2162473/campaign/manifest.json"
)


def mapped(path):
    """Map the controller's shared mount to its worker-visible spelling."""
    return str(path).replace("/mnt/data/slurm-results/", "/results/", 1)


def manifest(root):
    """Never accept modified experiment definitions after submission."""
    document = root / "manifest.json"
    assert sha256sum(document) == (root / "READY").read_text().strip()
    return json.loads(document.read_text())


def prepare(args):
    """Pin typical runs to all nodes; queue low-hit validation on both CPU classes."""
    original = json.loads(PARENT.read_text())
    samples = ["KPGP-00139", "LP6005592-DNA_C05"]
    tasks = [
        t
        for t in original["tasks"]
        if t["tool"] == "veprs"
        and t["level"] == "core"
        and t["assay"] == "wgs"
        and t["assembly"] == "GRCh38"
        and t["sample"] in samples
    ]
    assert len(tasks) == 2
    cases = []
    for index, node in enumerate([f"sl-w{i}" for i in range(1, 7)]):
        cases.append(
            dict(
                id=f"typical-{node}",
                sample=samples[0],
                node=node,
                requires=None,
                offset=index,
            )
        )
    for node in ["sl-w2", "sl-w5"]:
        cases.append(
            dict(
                id=f"lower-hit-{node}",
                sample=samples[1],
                node=node,
                requires=f"typical-{node}",
                offset=len(cases),
            )
        )
    args.root.mkdir(parents=True, exist_ok=False)
    (args.root / "logs").mkdir()
    data = dict(
        tasks=tasks,
        cases=cases,
        repetitions=2,
        cpus=8,
        modes=MODES,
        deployment=args.deployment,
        commit=args.commit,
        parent_sha256=sha256sum(PARENT),
        asset_attestations=original["freeze"]["asset_attestations"],
        asset_groups=original["freeze"]["annotation_asset_groups"],
        design="fresh same-node direct each replicate; reverse condition order in replicate 2; no technical best-picking",
    )
    write_json_atomic(args.root / "manifest.json", data)
    (args.root / "READY").write_text(sha256sum(args.root / "manifest.json") + "\n")
    jobs = {}
    submissions = []
    for case in cases:
        command = [
            PYTHON,
            "-m",
            "benchmarks.pass_fusion_pilot.run",
            "run",
            mapped(args.root),
            "--case",
            case["id"],
        ]
        submit = [
            "sbatch",
            "--parsable",
            "--partition=benchmark",
            "--exclusive",
            "--nodes=1",
            "--ntasks=1",
            "--cpus-per-task=8",
            "--mem=28000M",
            "--time=01:00:00",
            "--nodelist=" + case["node"],
            "--job-name=pass-fusion-" + case["id"],
            "--output=" + mapped(args.root / "logs") + "/" + case["id"] + "-%j.out",
        ]
        if case["requires"]:
            submit.append("--dependency=afterany:" + jobs[case["requires"]])
        submit += [
            "--wrap",
            "cd /tmp; PYTHONPATH="
            + shlex.quote(args.deployment)
            + " "
            + shlex.join(command),
        ]
        job = subprocess.check_output(submit, text=True).strip().split(";")[0]
        jobs[case["id"]] = job
        submissions.append(
            dict(case=case["id"], node=case["node"], job_id=job, command=submit)
        )
        write_json_atomic(args.root / "submissions.json", submissions)
    # A short durable review is automatic even after a benchmark failure.
    wrap = (
        "cd /tmp; PYTHONPATH="
        + shlex.quote(args.deployment)
        + " "
        + shlex.join(
            [
                PYTHON,
                "-m",
                "benchmarks.pass_fusion_pilot.run",
                "collect",
                mapped(args.root),
            ]
        )
    )
    review = [
        "sbatch",
        "--parsable",
        "--partition=benchmark",
        "--nodes=1",
        "--ntasks=1",
        "--cpus-per-task=1",
        "--mem=2000M",
        "--time=00:05:00",
        "--dependency=afterany:" + ":".join(jobs.values()),
        "--output=" + mapped(args.root / "logs/review-%j.out"),
        "--wrap",
        wrap,
    ]
    review_job = subprocess.check_output(review, text=True).strip().split(";")[0]
    write_json_atomic(
        args.root / "review_submission.json", dict(job_id=review_job, command=review)
    )
    print(json.dumps(dict(jobs=jobs, review=review_job), indent=2))


def run(root, identifier):
    """One exclusive node block; compare only fresh controls on that same node."""
    data = manifest(root)
    case = next(c for c in data["cases"] if c["id"] == identifier)
    evidence = root / "cases" / identifier
    evidence.mkdir(parents=True, exist_ok=False)
    state = dict(
        state="running",
        case=case,
        job_id=os.environ.get("SLURM_JOB_ID"),
        node=socket.gethostname(),
        rows=[],
        started_unix=time.time(),
    )
    write_json_atomic(evidence / "state.json", state)
    try:
        assert state["node"] == case["node"]
        if case["requires"]:
            previous = json.loads(
                (root / "cases" / case["requires"] / "state.json").read_text()
            )
            if previous["state"] != "complete":
                state.update(
                    state="skipped", reason="same-node typical block did not pass"
                )
                return
        hardware = dict(
            lscpu=json.loads(subprocess.check_output(["lscpu", "-J"], text=True)),
            storage=json.loads(
                subprocess.check_output(
                    [
                        "findmnt",
                        "-J",
                        "-T",
                        "/mnt/data/tmp",
                        "-o",
                        "SOURCE,FSTYPE,TARGET",
                    ],
                    text=True,
                )
            ),
            scratch_free_bytes=shutil.disk_usage("/mnt/data/tmp").free,
            kernel=subprocess.check_output(["uname", "-a"], text=True).strip(),
        )
        write_json_atomic(evidence / "hardware.json", hardware)
        assert (
            hardware["scratch_free_bytes"] > 25 * 1024**3
        ), "not enough safe local scratch"
        task = next(t for t in data["tasks"] if t["sample"] == case["sample"])
        attest = data["asset_attestations"]
        source = verified_file(task["input"], attest)
        cache = node_cache(task["cache"])
        for name, digest in task["cache"]["files"].items():
            if task["cache"].get("storage") != "node_local":
                verified_file(dict(path=str(cache / name), sha256=digest), attest)
        for asset in (
            task["annotation_assets"]
            + data["asset_groups"][task["annotation_asset_group"]]
        ):
            verified_file(asset, attest)
        with tempfile.TemporaryDirectory(
            prefix="pass-fusion-", dir="/mnt/data/tmp"
        ) as temporary:
            scratch = Path(temporary)
            prepared = scratch / "input.bcf"
            # Exactly the frozen WGS extraction contract; both timed arms still
            # perform their ordinary, unchanged normalization afterwards.
            commands = [
                ["bcftools", "view", str(source), "-Ou"],
                ["bcftools", "norm", "-m-", "-Ou"],
                ["bcftools", "view", "-e", 'ALT="*"', "-Ou"],
                [
                    "bcftools",
                    "norm",
                    "--rm-dup",
                    "exact",
                    "-Ob",
                    "-o",
                    str(prepared),
                    "-W",
                    "--threads",
                    "8",
                ],
            ]
            if task.get("preprocess", {}).get("bed"):
                commands[0] += [
                    "-R",
                    str(verified_file(task["preprocess"]["bed"], attest)),
                ]
            subprocess.run(
                [
                    "bash",
                    "-euo",
                    "pipefail",
                    "-c",
                    " | ".join(shlex.join(c) for c in commands),
                ],
                check=True,
            )
            assert (
                int(subprocess.check_output(["bcftools", "index", "-n", str(prepared)]))
                == task["input"]["records"]
            )
            state["prepared_sha256"] = sha256sum(prepared)
            parameters = dict(
                task["parameters"], threads=8, temp_dir=str(scratch / "tmp")
            )
            (scratch / "tmp").mkdir()
            params = scratch / "params.yaml"
            params.write_text(yaml.safe_dump(parameters))
            for repeat in range(data["repetitions"]):
                shift = case["offset"] % len(MODES)
                order = MODES[shift:] + MODES[:shift]
                if repeat:
                    order = list(reversed(order))
                with tempfile.TemporaryDirectory(
                    prefix=f"repeat-{repeat+1}-", dir=scratch
                ) as rep_temp:
                    outputs = {}
                    for mode in order:
                        directory = evidence / f"r{repeat+1:02d}" / mode
                        directory.mkdir(parents=True)
                        output = Path(rep_temp) / (mode + ".bcf")
                        outputs[mode] = output
                        stats_dir = Path(rep_temp) / (mode + "-stats")
                        command = [
                            PYTHON,
                            "-c",
                            "from vcfcache.cli import main; main()",
                            "annotate",
                            "-a",
                            str(cache),
                            "-i",
                            str(prepared),
                            "-o",
                            str(output),
                            "-y",
                            str(params),
                            "--engine",
                            "adaptive",
                            "--statistics",
                            "light",
                            "--stats-dir",
                            str(stats_dir),
                            "-v",
                        ]
                        if mode == "direct":
                            command += ["--uncached"]
                        elif mode != "control-auto":
                            command += ["--engine-workers", mode[-1]]
                        environment = {
                            k: v
                            for k, v in os.environ.items()
                            if not k.startswith("VCFCACHE_EXPERIMENT_")
                            and k != "VCFCACHE_PROFILE_DIR"
                        }
                        environment.update(
                            PYTHONPATH=data["deployment"],
                            VCFCACHE_ROOT=data["deployment"],
                            LC_ALL="C",
                            VCFCACHE_EXPERIMENT_FUSED_PUBLISH=(
                                "1" if mode.startswith("fused") else "0"
                            ),
                        )
                        launch = [
                            "srun",
                            "--exclusive",
                            "--exact",
                            "--cpu-bind=cores",
                            "-n1",
                            "-c8",
                            PYTHON,
                            str(
                                Path(data["deployment"])
                                / "benchmarks/engine_profile/measure.py"
                            ),
                            "--directory",
                            str(directory),
                            "--",
                            *command,
                        ]
                        write_json_atomic(
                            directory / "command.json",
                            dict(
                                command=launch,
                                task=task,
                                parameters=parameters,
                                mode=mode,
                                experimental_environment={
                                    "VCFCACHE_EXPERIMENT_FUSED_PUBLISH": environment[
                                        "VCFCACHE_EXPERIMENT_FUSED_PUBLISH"
                                    ]
                                },
                                source_commit=data["commit"],
                            ),
                        )
                        with (directory / "run.log").open("w") as stream:
                            subprocess.run(
                                launch,
                                env=environment,
                                stdout=stream,
                                stderr=stream,
                                check=True,
                            )
                        metrics = json.loads((directory / "resources.json").read_text())
                        assert (
                            metrics["exit_code"] == 0 and len(metrics["affinity"]) == 8
                        )
                        paths = list(stats_dir.rglob("compare_stats.yaml"))
                        assert len(paths) == 1
                        stats = yaml.safe_load(paths[0].read_text())
                        write_json_atomic(directory / "statistics.json", stats)
                        counts = stats["variant_counts"]
                        assert (
                            counts["input_variants"]
                            == counts["total_output"]
                            == task["input"]["records"]
                        )
                        if mode != "direct":
                            assert stats["engine"] == "adaptive" and stats[
                                "fused_publish"
                            ] == mode.startswith("fused")
                            if mode.startswith("fused"):
                                assert (
                                    stats["shards"]
                                    * (2 + stats["publish_threads_per_shard"])
                                    <= 8
                                )
                        row = dict(
                            repeat=repeat + 1,
                            mode=mode,
                            seconds=metrics["wall_seconds"],
                            records=counts["total_output"],
                            missing=counts.get("missing_variants"),
                            shards=stats.get("shards"),
                        )
                        state["rows"].append(row)
                        write_json_atomic(evidence / "state.json", state)
                    # No comparison I/O between timed conditions in a repeat.
                    for mode in order:
                        if mode == "direct":
                            continue
                        comparison = strict_semantic_compare(
                            outputs[mode], outputs["direct"]
                        )
                        write_json_atomic(
                            evidence
                            / f"r{repeat+1:02d}"
                            / mode
                            / "against_direct.json",
                            comparison,
                        )
                        assert (
                            comparison["semantic_pass"]
                            and comparison["ignored_fields"] == []
                        ), "strict annotation discrepancy"
            state["state"] = "complete"
    except Exception as error:
        state.update(state="failed", error=str(error), traceback=traceback.format_exc())
        raise
    finally:
        state["finished_unix"] = time.time()
        write_json_atomic(evidence / "state.json", state)


def collect(root):
    """Persist terminal or explicitly incomplete states; Slurm accounting is optional."""
    data = manifest(root)
    rows = []
    for case in data["cases"]:
        path = root / "cases" / case["id"] / "state.json"
        state = (
            json.loads(path.read_text()) if path.exists() else dict(state="not_started")
        )
        rows.append(dict(case=case, **state))
    result = dict(
        reviewed_unix=time.time(),
        cases=rows,
        complete=sum(r["state"] == "complete" for r in rows),
        failed=sum(r["state"] == "failed" for r in rows),
        not_terminal=sum(
            r["state"] not in {"complete", "failed", "skipped"} for r in rows
        ),
    )
    write_json_atomic(root / "REVIEW.json", result)
    print(json.dumps({k: v for k, v in result.items() if k != "cases"}, indent=2))


def main():
    """Submit only explicit experiments; a review never triggers more runs."""
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    prepare_parser = sub.add_parser("prepare")
    prepare_parser.add_argument("root", type=Path)
    prepare_parser.add_argument("--deployment", required=True)
    prepare_parser.add_argument("--commit", required=True)
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
