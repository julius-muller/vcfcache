"""Scoped header-guard recovery; preserve the parent freeze and direct controls."""

import argparse
import copy
import fcntl
import hashlib
import json
import os
import shlex
import shutil
import subprocess
import time
from pathlib import Path

from benchmarks.engine_profile.run import PYTHON, read, verify
from benchmarks.run_pilot import sha256sum, write_json_atomic


def filesystem(path):
    """Controller and worker views refer to the same shared filesystem."""
    value = str(path)
    if value.startswith("/results/") and not Path("/results").exists():
        value = value.replace("/results/", "/mnt/data/slurm-results/", 1)
    return Path(value)


def worker_path(path):
    """Express a shared path using the worker mount prefix."""
    return str(path).replace("/mnt/data/slurm-results/", "/results/", 1)


def direct_binding(task):
    """Exact direct condition, excluding routing/build and historical durations."""
    keys = [
        "id",
        "family",
        "sample",
        "assembly",
        "tool",
        "level",
        "assay",
        "cpus",
        "input",
        "preprocess",
        "parameters",
        "environment",
        "annotation_assets",
        "annotation_asset_group",
    ]
    return hashlib.sha256(
        json.dumps({k: task.get(k) for k in keys}, sort_keys=True).encode()
    ).hexdigest()


def reuse_direct(task, directory):
    """Freeze a finished original control, without rerunning its annotator."""
    reuse = task["baseline"]["direct_reuse"]
    source = filesystem(reuse["attempt"])
    command_file = source / "command.json"
    if sha256sum(command_file) != reuse["command_sha256"]:
        raise RuntimeError("reused direct condition command changed")
    original = read(command_file)["task"]
    if direct_binding(original) != direct_binding(task):
        raise RuntimeError("reused direct input/recipe/CPU/resource binding differs")
    state = read(source / "state.json")
    if state["state"] == "running" or str(state["job_id"]) != reuse["job_id"]:
        raise RuntimeError("original control has not finished under its expected job")
    metrics_file = source / "direct/resources.json"
    metrics = read(metrics_file)
    if metrics["exit_code"] != 0 or not metrics["cgroup_after"]:
        raise RuntimeError(
            "original direct control did not succeed with aggregate metrics"
        )
    if len(metrics["affinity"]) != task["cpus"] or metrics != state.get(
        "direct_metrics"
    ):
        raise RuntimeError(
            "original direct metrics do not match its checkpoint/allocation"
        )
    output = source / "direct.bcf"
    files = {
        name: sha256sum(source / name)
        for name in ["direct.bcf", "direct.bcf.csi", "direct/resources.json"]
    }
    for name, expected in reuse.get("files", {}).items():
        if files.get(name) != expected:
            raise RuntimeError(f"reused direct artifact changed: {name}")
    directory.mkdir()
    for name in ["resources.json", "command.json"]:
        shutil.copy2(source / "direct" / name, directory / name)
    write_json_atomic(
        directory / "REUSED.json",
        dict(
            source=worker_path(source),
            job_id=reuse["job_id"],
            files=files,
            binding=direct_binding(task),
            frozen_unix=time.time(),
            reason="Only receiving-header compatibility changed; exact original direct configuration retained",
        ),
    )
    return output, metrics


def live_jobs():
    """Read scheduler states without changing any job."""
    return {
        line.split()[0]: line.split()[1]
        for line in subprocess.check_output(
            ["squeue", "-h", "-o", "%i %T"], text=True
        ).splitlines()
    }


def register_replacements(root, manifest, batch):
    """Register new job IDs in the parent; cancel only their held predecessors."""
    repair = manifest["freeze"]["repair_of"]
    parent = filesystem(repair["root"])
    if sha256sum(parent / "manifest.json") != repair["manifest_sha256"]:
        raise RuntimeError("parent candidate freeze changed")
    record = dict(
        batch,
        repair_root=worker_path(root),
        repair_manifest_sha256=sha256sum(root / "manifest.json"),
    )
    with (parent / ".submit.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        history = read(parent / "submissions.json")
        existing = next(
            (
                b
                for b in history
                if b.get("repair_root") == record["repair_root"]
                and b["attempt"] == batch["attempt"]
            ),
            None,
        )
        if existing is None:
            history.append(record)
        else:
            existing.update(record)
        write_json_atomic(parent / "submissions.json", history)
    live = live_jobs()
    old = [repair["original_jobs"][j["cell"]] for j in batch["jobs"]]
    cancellable = [job for job in old if live.get(job) == "PENDING"]
    if cancellable:
        subprocess.run(["scancel", *cancellable], check=True)
    write_json_atomic(
        root / f"replacements-a{batch['attempt']}.json",
        dict(
            jobs=batch["jobs"],
            cancelled_held_predecessors=cancellable,
            running_predecessors_untouched=[j for j in old if live.get(j) == "RUNNING"],
        ),
    )


def finish_dispatch(root, manifest, batch):
    """Collect parent and repairs together, even if one attempt fails."""
    repair = manifest["freeze"]["repair_of"]
    build = manifest["freeze"]["build"]
    command = [
        "sbatch",
        "--parsable",
        "--partition=benchmark",
        "--cpus-per-task=1",
        "--mem=2000M",
        "--time=00:30:00",
        f"--dependency=afterany:{repair['original_collector']}:{batch['collector']}",
        f"--output={manifest['root']}/logs/final-parent-review-%j.out",
        "--wrap",
        f"cd /tmp; PYTHONPATH={shlex.quote(build['root'])} {PYTHON} -m benchmarks.performance_sprint.manuscript collect {shlex.quote(repair['root'])}",
    ]
    job = subprocess.check_output(command, text=True).strip().split(";")[0]
    write_json_atomic(
        root / "FINAL_PARENT_COLLECTOR.json", dict(job_id=job, command=command)
    )


def prepare_repair(parent, output, build):
    """Replace exactly the affected 36 slots; never repeat a successful slot."""
    from benchmarks.performance_sprint.manuscript import latest, prepare, start

    parent = filesystem(parent)
    output = filesystem(output)
    original = verify(parent)
    selected = [
        copy.deepcopy(t)
        for t in original["tasks"]
        if t["tool"] == "snpeff"
        and t["level"] == "clinical"
        and t["assembly"] == "GRCh37"
    ]
    if len(selected) != 36:
        raise RuntimeError("unexpected repair scope; refuse implicit cohort changes")
    preserved = [
        task["id"]
        for task in selected
        if latest(parent, task["id"]).get("state")
        in {"complete", "complete_with_direct_difference"}
    ]
    selected = [task for task in selected if task["id"] not in preserved]
    if not selected:
        raise RuntimeError("no unsuccessful conditions remain to repair")
    history = read(parent / "submissions.json")
    jobs = {j["cell"]: j["job_id"] for b in history for j in b["jobs"]}
    live = live_jobs()
    parent_proof = dict(
        path=worker_path(parent / "manifest.json"),
        sha256=sha256sum(parent / "manifest.json"),
    )
    for task in selected:
        state = latest(parent, task["id"])
        task["baseline"]["approval_evidence"] = task["baseline"][
            "approval_evidence"
        ] + [parent_proof]
        attempts = sorted((parent / "attempts").glob(f"*/{task['id']}/state.json"))
        if attempts:
            source = attempts[-1].parent
            if (source / "direct/command.json").exists():
                # For an in-flight original, the successor waits for its job.
                # Its final content hashes are frozen before cached timing.
                reuse = dict(
                    attempt=worker_path(source),
                    job_id=str(state["job_id"]),
                    command_sha256=sha256sum(source / "command.json"),
                )
                if state["state"] != "running":
                    reuse["files"] = {
                        name: sha256sum(source / name)
                        for name in [
                            "direct.bcf",
                            "direct.bcf.csi",
                            "direct/resources.json",
                        ]
                    }
                task["baseline"]["direct_reuse"] = reuse
                if live.get(reuse["job_id"]) in {"RUNNING", "COMPLETING"}:
                    task["after_jobs"] = [reuse["job_id"]]
    output.mkdir(parents=True, exist_ok=False)
    document = output / "tasks.json"
    write_json_atomic(document, selected)
    freeze = copy.deepcopy(original["freeze"])
    gate = min(
        (
            t
            for t in selected
            if t["assay"] == "panel"
            and t["baseline"].get("historical_acceptance") == "PASS"
        ),
        key=lambda t: t["historical_cached_seconds"],
    )["id"]
    freeze.update(
        build=build,
        tasks_sha256=sha256sum(document),
        required_task_ids=[t["id"] for t in selected],
        first_cells=[gate],
        aliases={},
        attempt_base=2,
        created_unix=time.time(),
        approved_decision="User approved targeted INFO-cardinality repair and reuse of completed direct controls",
        repair_of=dict(
            root=worker_path(parent),
            manifest_sha256=parent_proof["sha256"],
            original_jobs={t["id"]: jobs[t["id"]] for t in selected},
            original_collector=next(
                b["collector"] for b in reversed(history) if not b.get("repair_root")
            ),
            preserved_successful_ids=preserved,
        ),
    )
    write_json_atomic(output / "freeze.json", freeze)
    root = output / "campaign"
    prepare(root, document, freeze)
    # Only our newly created empty directory is removed. The original attempt
    # tree is preserved; the repair writes new, higher-numbered attempts there.
    (root / "attempts").rmdir()
    (root / "attempts").symlink_to(os.path.relpath(parent / "attempts", root))
    receipt = start(root)
    write_json_atomic(output / "SUBMITTED.json", receipt)
    print(json.dumps(dict(root=worker_path(root), gate=gate, start=receipt), indent=2))


def main():
    """Prepare and submit only the approved header-guard replacement cohort."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("parent", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--build", type=Path, required=True)
    args = parser.parse_args()
    prepare_repair(args.parent, args.output, read(filesystem(args.build)))


if __name__ == "__main__":
    main()
