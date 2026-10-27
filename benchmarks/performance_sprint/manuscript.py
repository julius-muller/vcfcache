"""Gated cached-only manuscript replay with durable, independently timed cells.

This runner deliberately cannot infer approval from a passing asset audit.
Tasks must bind exact inputs, cache files, parameters and historical evidence;
the candidate freeze and reuse approval are separate, explicit inputs.
"""

import argparse
import fcntl
import hashlib
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

from benchmarks.engine_profile.run import PYTHON, persist_run, read, verify
from benchmarks.run_pilot import (
    semantic_compare,
    sha256sum,
    strict_semantic_compare,
    write_json_atomic,
)


def validate_task(task):
    """Reject incomplete reuse bindings and implicit direct annotation requests."""
    if not task["id"].replace("-", "").isalnum():
        raise ValueError("unsafe task id")
    if task["cpus"] != task["parameters"]["threads"] or not 1 <= task["cpus"] <= 16:
        raise ValueError("CPU settings do not match frozen parameters")
    reference = task["baseline"]
    if reference["cpus"] != task["cpus"]:
        raise ValueError("CPU allocation differs from historical direct baseline")
    paired = reference["kind"] == "new_paired_reference"
    if paired and (task["family"] != "veprs_extension" or task["tool"] != "veprs"):
        raise ValueError(
            "fresh direct measurements are restricted to the new vep-rs extension"
        )
    if not paired and reference["reuse_decision"] != "approved":
        raise ValueError("baseline reuse has not been approved")
    if (not paired and reference["direct_seconds"] <= 0) or not reference.get(
        "approval_evidence"
    ):
        raise ValueError("baseline timing/evidence missing")
    if reference["kind"] in {"complete_records", "composed_vep_reference"}:
        for name in ["direct", "cached"]:
            if set(reference[name]) != {"path", "sha256"}:
                raise ValueError("full output reference requires exact SHA256 binding")
    elif reference["kind"] == "historical_annotation_digest":
        for name in ["direct_digest", "cached_digest", "query_format", "blocksort"]:
            if not reference.get(name):
                raise ValueError(f"historical digest binding missing: {name}")
    elif not paired:
        raise ValueError("unsupported baseline evidence")
    if task["engine"] not in {"legacy", "adaptive", "native"}:
        raise ValueError("engine must be explicitly frozen")
    required = {
        "annotation.yaml",
        "vcfcache_annotated.bcf",
        "vcfcache_annotated.bcf.csi",
    }
    if not required.issubset(task["cache"]["files"]):
        raise ValueError("cache recipe/payload/index fingerprints required")
    if not task["input"].get("sha256") or not task["annotation_assets"]:
        raise ValueError("input and annotation resource fingerprints required")
    if set(task.get("environment", {})) - {
        "SNPEFF_HEAP",
        "SNPSIFT_JAR",
        "CLINVAR_VCF",
        "GNOMAD_VCF",
    }:
        raise ValueError("unsupported annotator environment override")
    for key in ["vep_forks", "rayon_threads"]:
        if task["parameters"].get(key, 1) > task["cpus"]:
            raise ValueError("annotator threads exceed matched CPU allocation")


def minutes(task):
    """Historical cached runtime plus explicit audit allowance, capped at 24 h."""
    return task.get(
        "wall_minutes",
        min(1440, max(30, int(task["historical_cached_seconds"] * 2 / 60 + 90))),
    )


def exclusion_args(freeze):
    """Use every eligible worker unless a reviewed exclusion records its reason."""
    exclusions = freeze.get("node_exclusions", {})
    for node, reason in exclusions.items():
        if (
            not node.replace("-", "").isalnum()
            or not isinstance(reason, str)
            or not reason.strip()
        ):
            raise ValueError(
                "each node exclusion requires a literal hostname and an explicit reason"
            )
    return ["--exclude=" + ",".join(sorted(exclusions))] if exclusions else []


def check_file(record):
    """Hash outside timing; a missing historical artifact is never silently rebuilt."""
    path = Path(record["path"])
    if sha256sum(path) != record["sha256"]:
        raise RuntimeError(f"frozen file changed: {path}")
    return path


def verified_file(record, attestations):
    """Reuse a frozen SHA audit only while size and nanosecond mtime agree."""
    path = Path(record["path"])
    pinned = attestations.get(str(path))
    if pinned and pinned["sha256"] == record["sha256"]:
        stat = path.stat()
        if stat.st_size == pinned["bytes"] and stat.st_mtime_ns == pinned["mtime_ns"]:
            return path
    return check_file(record)


def annotation_digest(path, reference):
    """Reproduce C2's exact blocksort/per-contig MD5 protocol in streaming memory.

    It intentionally covers only the original query columns, not all INFO/FORMAT.
    No extra CSQ canonicalization or undocumented exception is introduced.
    """
    blocksort = check_file(reference["blocksort"])
    query = subprocess.Popen(
        ["bcftools", "query", "-f", reference["query_format"], str(path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert query.stdout is not None
    sort = subprocess.Popen(
        ["awk", "-f", str(blocksort)],
        stdin=query.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    query.stdout.close()
    assert sort.stdout is not None
    digests = {}
    for line in sort.stdout:
        contig = line.split(b"\t", 1)[0]
        digests.setdefault(contig, hashlib.md5()).update(line)
    _, sorting_error = sort.communicate()
    _, query_error = query.communicate()
    if query.returncode or sort.returncode:
        raise RuntimeError(f"digest query failed: {query_error!r} {sorting_error!r}")
    result = hashlib.md5()
    for contig, digest in sorted(digests.items()):
        result.update(contig + b" " + digest.hexdigest().encode() + b"\n")
    return result.hexdigest()


def prepare(root, document, freeze):
    """Freeze only a reviewed candidate and a complete approved historical matrix."""
    if not freeze.get("approved") or not freeze.get("tests_passed"):
        raise ValueError("candidate freeze is not approved")
    exclusion_args(freeze)
    if freeze["tasks_sha256"] != sha256sum(document):
        raise ValueError("freeze approval does not bind this task document")
    tasks = read(document)
    if len({t["id"] for t in tasks}) != len(tasks):
        raise ValueError("duplicate tasks")
    if not tasks:
        raise ValueError("empty manuscript matrix")
    for task in tasks:
        validate_task(task)
    expected_ids = set(freeze["required_task_ids"])
    if expected_ids != {t["id"] for t in tasks}:
        raise ValueError(
            "required experiment families/samples were dropped or added without review"
        )
    root.mkdir(parents=True, exist_ok=False)
    for directory in ["attempts", "logs", "reviews"]:
        (root / directory).mkdir()
    manifest = dict(
        root=str(root).replace("/mnt/data/slurm-results/", "/results/", 1),
        tasks=tasks,
        freeze=freeze,
        created_unix=time.time(),
        mode="cached_only_no_implicit_direct_reruns",
        deadline="2026-11-02T18:00:00+01:00",
    )
    write_json_atomic(root / "manifest.json", manifest)
    (root / "READY").write_text(sha256sum(root / "manifest.json") + "\n")


def latest(root, identifier):
    """Read the last explicit attempt; never overwrite a successful cell."""
    paths = sorted((root / "attempts").glob(f"*/{identifier}/state.json"))
    return read(paths[-1]) if paths else {}


def collect(root):
    """Persist a complete snapshot even when a job failed before writing state."""
    manifest = verify(root)
    history = (
        read(root / "submissions.json") if (root / "submissions.json").exists() else []
    )
    jobs = {j["cell"]: j["job_id"] for batch in history for j in batch["jobs"]}
    active = None
    if jobs and shutil.which("squeue"):
        probe = subprocess.run(
            ["squeue", "-h", "-j", ",".join(jobs.values()), "-o", "%i"],
            text=True,
            capture_output=True,
        )
        if probe.returncode == 0:
            active = set(probe.stdout.split())
    rows = []
    for task in manifest["tasks"]:
        row = dict(id=task["id"], family=task["family"], **latest(root, task["id"]))
        row.setdefault("state", "not_started")
        if task["id"] in jobs:
            row["submitted_job_id"] = jobs[task["id"]]
            row["live_scheduler"] = (
                "unavailable"
                if active is None
                else ("active" if jobs[task["id"]] in active else "not_active")
            )
            if (
                row["state"] in {"running", "not_started"}
                and active is not None
                and jobs[task["id"]] not in active
            ):
                row.update(
                    state="incomplete",
                    reason="job is not active and has no completed checkpoint; inspect logs before explicit retry",
                )
        rows.append(row)
    snapshot = dict(
        created_unix=time.time(),
        rows=rows,
        complete=sum(r.get("state") == "complete" for r in rows),
    )
    destination = root / "reviews" / f"review-{time.time_ns()}.json"
    write_json_atomic(destination, snapshot)
    return snapshot


def submit(root, selected=None):
    """One job/step per cached condition, explicit retries, afterany collector."""
    manifest = verify(root)
    selected = set(selected or [])
    if selected - {t["id"] for t in manifest["tasks"]}:
        raise ValueError("unknown retry task")
    with (root / ".submit.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        history = (
            read(root / "submissions.json")
            if (root / "submissions.json").exists()
            else []
        )
        if history and not selected:
            raise ValueError("already submitted: resume requires explicit task ids")
        old = [j["job_id"] for b in history for j in b["jobs"]]
        old.extend(b["collector"] for b in history if b.get("collector"))
        if (
            old
            and subprocess.check_output(
                ["squeue", "-h", "-j", ",".join(old), "-o", "%i"], text=True
            ).strip()
        ):
            raise ValueError("previous submission is still active")
        tasks = [t for t in manifest["tasks"] if not selected or t["id"] in selected]
        if any(
            latest(root, t["id"]).get("state")
            in {"complete", "complete_with_direct_difference"}
            for t in tasks
        ):
            raise ValueError("refusing to rerun a successful condition")
        attempt = len(history) + 1
        batch = dict(attempt=attempt, jobs=[], created_unix=time.time())
        history.append(batch)
        write_json_atomic(root / "submissions.json", history)
        build = manifest["freeze"]["build"]
        exclude = exclusion_args(manifest["freeze"])
        worker = manifest["root"]
        prefix = f"cd /tmp; PYTHONPATH={shlex.quote(build['root'])} {PYTHON} -m benchmarks.performance_sprint.manuscript"
        try:
            # Longest historical cached conditions first to reduce the tail.
            for task in sorted(tasks, key=lambda t: -t["historical_cached_seconds"]):
                command = [
                    "sbatch",
                    "--parsable",
                    "--partition=benchmark",
                    "--nodes=1",
                    "--ntasks=1",
                    f"--cpus-per-task={task['cpus']}",
                    "--mem=28000M",
                    "--exclusive",
                    *exclude,
                    f"--time={minutes(task)}",
                    f"--output={worker}/logs/{task['id']}-a{attempt}-%j.out",
                    "--wrap",
                    f"{prefix} run {shlex.quote(worker)} --cell {task['id']} --attempt {attempt}",
                ]
                job = subprocess.check_output(command, text=True).strip().split(";")[0]
                if not job.isdigit():
                    raise RuntimeError(f"unexpected sbatch response: {job}")
                batch["jobs"].append(dict(cell=task["id"], job_id=job, command=command))
                write_json_atomic(root / "submissions.json", history)
        finally:
            if batch["jobs"]:
                collector = [
                    "sbatch",
                    "--parsable",
                    "--partition=benchmark",
                    "--nodes=1",
                    "--ntasks=1",
                    "--cpus-per-task=1",
                    "--mem=2000M",
                    *exclude,
                    "--time=00:30:00",
                    "--dependency=afterany:"
                    + ":".join(j["job_id"] for j in batch["jobs"]),
                    f"--output={worker}/logs/collect-a{attempt}-%j.out",
                    "--wrap",
                    f"{prefix} collect {shlex.quote(worker)}",
                ]
                batch["collector"] = (
                    subprocess.check_output(collector, text=True).strip().split(";")[0]
                )
                write_json_atomic(root / "submissions.json", history)
    return batch


def dispatch(root):
    """Release remaining primary cells only after the six first cells pass."""
    manifest = verify(root)
    initial = manifest["freeze"]["first_cells"]
    bad = {
        identifier: latest(root, identifier)
        for identifier in initial
        if latest(root, identifier).get("state")
        not in {"complete", "complete_with_direct_difference"}
    }
    if bad:
        write_json_atomic(
            root / "DISPATCH_BLOCKED.json",
            dict(
                reason="first primary cells failed validation; inspect failures before explicit retry",
                cells=bad,
                created_unix=time.time(),
            ),
        )
        raise RuntimeError("first-cell validation blocked automatic dispatch")
    remaining = [t["id"] for t in manifest["tasks"] if t["id"] not in initial]
    result = submit(root, remaining)
    write_json_atomic(root / "AUTO_DISPATCH.json", result)
    return result


def start(root):
    """Submit first primary cells plus a durable gated continuation, not repeats."""
    manifest = verify(root)
    batch = submit(root, manifest["freeze"]["first_cells"])
    deployment = manifest["freeze"]["build"]["root"]
    command = [
        "sbatch",
        "--parsable",
        "--partition=benchmark",
        "--cpus-per-task=1",
        "--mem=2000M",
        "--time=00:30:00",
        f"--dependency=afterany:{batch['collector']}",
        f"--output={manifest['root']}/logs/dispatch-%j.out",
        "--wrap",
        f"cd /tmp; PYTHONPATH={shlex.quote(deployment)} {PYTHON} -m benchmarks.performance_sprint.manuscript dispatch {shlex.quote(manifest['root'])}",
    ]
    job = subprocess.check_output(command, text=True).strip().split(";")[0]
    receipt = dict(first_batch=batch, dispatch_job=job, command=command)
    write_json_atomic(root / "START.json", receipt)
    return receipt


def without_csq_compare(first, second, scratch, cpus):
    """Check complete non-CSQ records before classifying known clinical CSQ differences."""
    outputs = []
    for number, source in enumerate([first, second]):
        output = scratch / f"non-csq-{number}.bcf"
        subprocess.run(
            [
                "bcftools",
                "annotate",
                "-x",
                "INFO/CSQ",
                str(source),
                "-Ob",
                "-W",
                "--threads",
                str(cpus),
                "-o",
                str(output),
            ],
            check=True,
        )
        outputs.append(output)
    return strict_semantic_compare(*outputs)


def run(root, identifier, attempt):
    """Prepare untimed, run one isolated timed step, audit, then checkpoint."""
    manifest = verify(root)
    task = next(t for t in manifest["tasks"] if t["id"] == identifier)
    attestations = manifest["freeze"].get("asset_attestations", {})

    def checked_file(record):
        return verified_file(record, attestations)

    validate_task(task)
    evidence = root / "attempts" / f"a{attempt:03d}" / identifier
    evidence.mkdir(parents=True, exist_ok=False)
    state = dict(
        state="running",
        started_unix=time.time(),
        node=socket.gethostname(),
        job_id=os.environ.get("SLURM_JOB_ID"),
    )
    write_json_atomic(evidence / "state.json", state)
    with tempfile.TemporaryDirectory(
        prefix=f"manuscript-{identifier}-", dir="/mnt/data/tmp"
    ) as directory:
        scratch = Path(directory)
        output = scratch / "output.bcf"
        try:
            for artifact in task["baseline"]["approval_evidence"]:
                checked_file(artifact)
            for path, digest in manifest["freeze"]["build"]["files"].items():
                check_file(
                    dict(
                        path=str(Path(manifest["freeze"]["build"]["root"]) / path),
                        sha256=digest,
                    )
                )
            source = checked_file(task["input"])
            for name, digest in task["cache"]["files"].items():
                checked_file(
                    dict(
                        path=str(Path(task["cache"]["directory"]) / name), sha256=digest
                    )
                )
            for asset in task["annotation_assets"]:
                checked_file(asset)
            group = task.get("annotation_asset_group")
            if group:
                for asset in manifest["freeze"]["annotation_asset_groups"][group]:
                    checked_file(asset)
            if task.get("preprocess"):
                preparation = task["preprocess"]
                view = ["bcftools", "view", str(source), "-Ou"]
                if preparation.get("bed"):
                    view += ["-R", str(checked_file(preparation["bed"]))]
                source = scratch / "input.bcf"
                # All paths are shlex-quoted; the fixed pipeline exactly matches C2.
                pipeline = [
                    view,
                    ["bcftools", "norm", "-m-", "-Ou"],
                    ["bcftools", "view", "-e", 'ALT="*"', "-Ou"],
                    [
                        "bcftools",
                        "norm",
                        "--rm-dup",
                        "exact",
                        "-Ob",
                        "-o",
                        str(source),
                        "-W",
                        "--threads",
                        str(task["cpus"]),
                    ],
                ]
                subprocess.run(
                    [
                        "bash",
                        "-o",
                        "pipefail",
                        "-c",
                        " | ".join(shlex.join(c) for c in pipeline),
                    ],
                    check=True,
                )
            count = int(
                subprocess.check_output(
                    ["bcftools", "index", "-n", str(source)], text=True
                )
            )
            if count != task["input"]["records"]:
                raise RuntimeError("prepared record count differs from historical cell")
            state["prepared_input_sha256"] = sha256sum(source)
            parameters = dict(task["parameters"], temp_dir=str(scratch / "tmp"))
            (scratch / "tmp").mkdir()
            (scratch / "params.yaml").write_text(yaml.safe_dump(parameters))
            build = manifest["freeze"]["build"]
            command = [
                PYTHON,
                "-c",
                "from vcfcache.cli import main; main()",
                "annotate",
                "-a",
                task["cache"]["directory"],
                "-i",
                str(source),
                "-o",
                str(output),
                "-y",
                str(scratch / "params.yaml"),
                "--engine",
                task["engine"],
                "--statistics",
                "light",
                "--stats-dir",
                str(scratch / "stats"),
            ]
            if task["engine"] == "native":
                command += ["--engine-memory-mb", str(task.get("memory_mb", 512))]
            environment = dict(
                os.environ,
                PYTHONPATH=build["root"],
                VCFCACHE_ROOT=build["root"],
                VCFCACHE_ENGINE=build["root"]
                + "/native/target/release/vcfcache-engine",
                LC_ALL="C",
                **task.get("environment", {}),
            )
            for key in list(environment):
                if (
                    key.startswith("VCFCACHE_EXPERIMENT_")
                    or key == "VCFCACHE_PROFILE_DIR"
                ):
                    environment.pop(key)
            launch = [
                "srun",
                "--exclusive",
                "--exact",
                "--cpu-bind=cores",
                "-n1",
                f"-c{task['cpus']}",
                PYTHON,
                str(Path(build["root"]) / "benchmarks/engine_profile/measure.py"),
                "--directory",
                str(evidence),
                "--",
                *command,
            ]
            write_json_atomic(
                evidence / "command.json",
                dict(command=launch, task=task, parameters=parameters),
            )
            direct_output = None
            if task["baseline"]["kind"] == "new_paired_reference":
                direct_output = scratch / "direct.bcf"
                direct_directory = evidence / "direct"
                direct_directory.mkdir()
                direct_command = list(command)
                direct_command[direct_command.index("--engine") + 1] = "legacy"
                direct_command[direct_command.index("-o") + 1] = str(direct_output)
                direct_command[direct_command.index("--stats-dir") + 1] = str(
                    scratch / "direct-stats"
                )
                direct_command.append("--uncached")
                direct_launch = list(launch[: launch.index("--") + 1]) + direct_command
                direct_launch[direct_launch.index("--directory") + 1] = str(
                    direct_directory
                )
                write_json_atomic(
                    direct_directory / "command.json",
                    dict(command=direct_launch, parameters=parameters),
                )
                with (direct_directory / "run.log").open("w") as log:
                    subprocess.run(
                        direct_launch,
                        env=environment,
                        stdout=log,
                        stderr=log,
                        check=True,
                    )
                state["direct_metrics"] = read(direct_directory / "resources.json")
                if not state["direct_metrics"]["cgroup_after"]:
                    raise RuntimeError(
                        "new direct step lacks independent cgroup metrics"
                    )
                state["direct_seconds"] = state["direct_metrics"]["wall_seconds"]
            with (evidence / "run.log").open("w") as log:
                subprocess.run(
                    launch, env=environment, stdout=log, stderr=log, check=True
                )
            state["metrics"] = read(evidence / "resources.json")
            if not state["metrics"]["cgroup_after"]:
                raise RuntimeError("condition-specific cgroup metrics unavailable")
            statistics = list((scratch / "stats").glob("*_vcstats/compare_stats.yaml"))
            if len(statistics) != 1:
                raise RuntimeError("missing or ambiguous annotation statistics")
            counts = yaml.safe_load(statistics[0].read_text())
            write_json_atomic(evidence / "statistics.json", counts)
            state["input_records"] = count
            state["cache_hit_rate"] = (
                1 - counts["variant_counts"]["missing_variants"] / count
                if count
                else None
            )
            state["effective_engine"] = counts.get("effective_engine", task["engine"])
            if (
                int(
                    subprocess.check_output(
                        ["bcftools", "index", "-n", str(output)], text=True
                    )
                )
                != count
            ):
                raise RuntimeError("output count mismatch")
            reference = task["baseline"]
            state["direct_seconds"] = state.get(
                "direct_seconds", reference.get("direct_seconds")
            )
            if reference["kind"] in {"complete_records", "composed_vep_reference"}:
                cached = checked_file(reference["cached"])
                direct = checked_file(reference["direct"])
                exact = strict_semantic_compare(output, cached)
                comparator = (
                    semantic_compare
                    if task["tool"].lower() == "vep"
                    else strict_semantic_compare
                )
                against_direct = comparator(output, direct)
                comparison = dict(
                    against_cached=exact,
                    against_direct=against_direct,
                    validation_scope="complete cached-reference records; annotator-specific direct comparator",
                )
                write_json_atomic(evidence / "comparison.json", comparison)
                if reference["kind"] == "composed_vep_reference":
                    # This is the same input/annotation recipe plus a plugin
                    # that emits no fields. Do not call it a retained per-delay
                    # full-record validation when those BCFs no longer exist.
                    composed = semantic_compare(output, cached)
                    comparison["composed_cached_reference"] = composed
                    comparison["validation_scope"] = (
                        "composed VEP semantic reference; per-delay historical timing retained, not full-record per-delay validation"
                    )
                    write_json_atomic(evidence / "comparison.json", comparison)
                    if not composed["semantic_pass"]:
                        raise RuntimeError(
                            "synthetic output differs from composed annotation reference"
                        )
                elif not exact["semantic_pass"]:
                    if reference.get("allow_documented_vep_exception"):
                        documented = semantic_compare(output, cached)
                        other_fields = without_csq_compare(
                            output, cached, scratch, task["cpus"]
                        )
                        comparison.update(
                            documented_vep_cached=documented,
                            non_csq_cached=other_fields,
                        )
                        write_json_atomic(evidence / "comparison.json", comparison)
                        if not (
                            documented["semantic_pass"]
                            and other_fields["semantic_pass"]
                            and exact["relevant_headers_equal"]
                        ):
                            raise RuntimeError(
                                "new output differs beyond the documented VEP exception"
                            )
                        comparison["validation_scope"] = (
                            "complete non-CSQ records and documented VEP CSQ equivalence; raw differences retained"
                        )
                    else:
                        raise RuntimeError(
                            "new output differs from frozen cached reference"
                        )
                if not against_direct["semantic_pass"] and not reference.get(
                    "inherited_direct_difference"
                ):
                    raise RuntimeError("unexpected direct-output difference")
            elif reference["kind"] == "new_paired_reference":
                comparison = strict_semantic_compare(output, direct_output)
                comparison["validation_scope"] = (
                    "fresh paired complete-record/header comparison; no ignored fields"
                )
                if not comparison["semantic_pass"]:
                    if task["level"] != "clinical":
                        write_json_atomic(evidence / "comparison.json", comparison)
                        raise RuntimeError("new core vep-rs paired output differs")
                    non_csq = without_csq_compare(
                        output, direct_output, scratch, task["cpus"]
                    )
                    comparison["non_csq_comparison"] = non_csq
                    if not non_csq["semantic_pass"]:
                        write_json_atomic(evidence / "comparison.json", comparison)
                        raise RuntimeError(
                            "clinical vep-rs differs outside its predeclared CSQ limitation"
                        )
                    comparison["validation_scope"] = (
                        "pinned clinical vep-rs timing only; CSQ discrepancies retained and NOT validated speedup evidence"
                    )
                    state["timing_only"] = True
            else:
                digest = annotation_digest(output, reference)
                comparison = dict(
                    digest=digest,
                    cached_match=digest == reference["cached_digest"],
                    direct_match=digest == reference["direct_digest"],
                    validation_scope="historical annotation/key digest only; NOT complete-record direct validation",
                )
                write_json_atomic(evidence / "comparison.json", comparison)
                if not comparison["cached_match"]:
                    if (
                        str(reference.get("historical_acceptance", "")).startswith(
                            "FAIL_CONTENT"
                        )
                        and comparison["direct_match"]
                    ):
                        state["timing_only"] = True
                        comparison["validation_scope"] = (
                            "predeclared SnpEff context exclusion now matches historical direct digest; raw cached-reference difference retained"
                        )
                    else:
                        raise RuntimeError(
                            "new output differs from accepted historical cached digest"
                        )
                elif not comparison["direct_match"]:
                    state["timing_only"] = True
            write_json_atomic(evidence / "comparison.json", comparison)
            state.update(
                state="complete_with_direct_difference"
                if state.get("timing_only")
                else "complete",
                output_sha256=sha256sum(output),
                output_index_sha256=sha256sum(Path(str(output) + ".csi")),
                speedup=state["direct_seconds"] / state["metrics"]["wall_seconds"],
                validation_scope=comparison["validation_scope"],
            )
        except Exception as error:
            state.update(
                state="failed", error=str(error), traceback=traceback.format_exc()
            )
            if output.exists():
                shutil.copy2(output, evidence / "output.bcf")
                if Path(str(output) + ".csi").exists():
                    shutil.copy2(str(output) + ".csi", evidence / "output.bcf.csi")
        finally:
            persist_run(scratch, evidence)
            if (
                task["baseline"]["kind"] == "new_paired_reference"
                and (scratch / "direct.bcf").exists()
            ):
                for name in ["direct.bcf", "direct.bcf.csi"]:
                    if (scratch / name).exists():
                        shutil.copy2(scratch / name, evidence / name)
            state["finished_unix"] = time.time()
            write_json_atomic(evidence / "state.json", state)
    return state


def main():
    """Expose explicit preparation, submission and days-later review commands."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "command",
        choices=[
            "prepare",
            "submit",
            "resume",
            "run",
            "collect",
            "status",
            "start",
            "dispatch",
        ],
    )
    parser.add_argument("campaign", type=Path)
    parser.add_argument("--tasks", type=Path)
    parser.add_argument("--freeze", type=Path)
    parser.add_argument("--cell", action="append")
    parser.add_argument("--attempt", type=int, default=1)
    args = parser.parse_args()
    if args.command == "prepare":
        if not args.tasks or not args.freeze:
            parser.error("prepare requires --tasks and --freeze")
        prepare(args.campaign, args.tasks, read(args.freeze))
        result = dict(state="prepared_not_submitted")
    elif args.command in {"submit", "resume"}:
        if args.command == "resume" and not args.cell:
            parser.error("resume requires explicit --cell values")
        result = submit(args.campaign, args.cell)
    elif args.command == "run":
        if not args.cell or len(args.cell) != 1:
            parser.error("run requires exactly one --cell")
        result = run(args.campaign, args.cell[0], args.attempt)
    elif args.command in {"start", "dispatch"}:
        result = (
            start(args.campaign) if args.command == "start" else dispatch(args.campaign)
        )
    else:
        result = collect(args.campaign)
    print(json.dumps(result, indent=2))
    if args.command == "run" and result["state"] not in {
        "complete",
        "complete_with_direct_difference",
    }:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
