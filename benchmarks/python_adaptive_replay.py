"""Freeze a cached-only Python-adaptive replay; never remeasure direct arms.

Run on the Slurm controller from an immutable deployment. Existing campaign
manifests, results and figures remain untouched. The existing six-cell gate
automatically releases the remainder; an afterany collector persists status.
"""

import argparse
import copy
import json
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from benchmarks.engine_profile.run import read, verify
from benchmarks.performance_sprint import manuscript
from benchmarks.performance_sprint.repair_headers import (
    direct_binding,
    filesystem,
    worker_path,
)
from benchmarks.run_pilot import sha256sum, write_json_atomic


def freeze_replay(parent, repairs, destination, deployment, commit, tests):
    """Bind every retained direct control and freeze a new independent matrix."""
    original = verify(parent)
    tasks = {task["id"]: copy.deepcopy(task) for task in original["tasks"]}
    source_roots = [parent]
    for repair in repairs:
        replacement = verify(repair)
        for task in replacement["tasks"]:
            if task["id"] not in tasks:
                raise ValueError("repair adds an unknown condition")
            tasks[task["id"]] = copy.deepcopy(task)
        source_roots.append(repair)
    source_proofs = [
        dict(
            path=worker_path(root / "manifest.json"),
            sha256=sha256sum(root / "manifest.json"),
        )
        for root in source_roots
    ]

    def bind(task):
        # Existing repairs can share attempt trees. Select the latest explicit
        # attempt, not the fastest measurement or the most favorable outcome.
        attempts = sorted(
            {
                p.resolve()
                for root in source_roots
                for p in (root / "attempts").glob(f"*/{task['id']}/state.json")
            },
            key=lambda p: (p.parent.parent.name, str(p)),
        )
        if not attempts:
            raise RuntimeError(f"missing prior checkpoint: {task['id']}")
        source = attempts[-1].parent
        state = read(source / "state.json")
        if state["state"] == "running":
            raise RuntimeError("prior condition is still running")
        task["engine"] = "adaptive"
        task.pop("after_jobs", None)
        task["submission_priority"] = (
            1 if task["family"] == "synthetic_pipeline_spectrum" else 0
        )
        task["previous_result"] = dict(
            path=worker_path(source / "state.json"),
            sha256=sha256sum(source / "state.json"),
            state=state["state"],
        )
        task["historical_cached_seconds"] = state.get("metrics", {}).get(
            "wall_seconds", task["historical_cached_seconds"]
        )
        reference = task["baseline"]
        reference["approval_evidence"] += source_proofs
        if reference["kind"] == "new_paired_reference":
            # A repaired cached attempt can itself have reused an earlier direct
            # arm; bind the original retained output, not an absent local copy.
            if reference.get("direct_reuse"):
                source = filesystem(reference["direct_reuse"]["attempt"])
            command = read(source / "command.json")
            if direct_binding(command["task"]) != direct_binding(task):
                raise RuntimeError(f"direct configuration changed: {task['id']}")
            prior = read(source / "state.json")
            metrics = read(source / "direct/resources.json")
            if (
                metrics["exit_code"] != 0
                or metrics != prior.get("direct_metrics")
                or len(metrics["affinity"]) != task["cpus"]
            ):
                raise RuntimeError(f"invalid retained direct measurement: {task['id']}")
            reference["direct_reuse"] = dict(
                attempt=worker_path(source),
                job_id=str(prior["job_id"]),
                command_sha256=sha256sum(source / "command.json"),
                files={
                    name: sha256sum(source / name)
                    for name in [
                        "direct.bcf",
                        "direct.bcf.csi",
                        "direct/resources.json",
                    ]
                },
            )
            reference["direct_seconds"] = metrics["wall_seconds"]
        manuscript.validate_task(task)
        return task

    # Only four read-only hashing workers on the controller; timed annotation
    # starts afterward on exclusive workers with the original CPU allocation.
    with ThreadPoolExecutor(max_workers=4) as pool:
        bound = list(pool.map(bind, tasks.values()))
    destination.mkdir(parents=True, exist_ok=False)
    document = destination / "tasks.json"
    write_json_atomic(document, bound)
    freeze = copy.deepcopy(original["freeze"])
    for key in ["repair_of", "attempt_base"]:
        freeze.pop(key, None)
    files = {}
    for directory in ["vcfcache", "resources", "benchmarks"]:
        for path in sorted((deployment / directory).rglob("*")):
            if (
                path.is_file()
                and "__pycache__" not in path.parts
                and path.suffix in {".py", ".yaml", ".sh", ".awk"}
            ):
                files[str(path.relative_to(deployment))] = sha256sum(path)
    freeze.update(
        approved=True,
        tests_passed=tests,
        approved_decision="User-approved Python/bcftools default and cached-only replay; urgency prioritizes primary figures",
        build=dict(
            commit=commit,
            root=worker_path(deployment),
            files=files,
            implementation="python-bcftools-adaptive; no native helper",
        ),
        engine_policy="Python/bcftools adaptive only; never relabel historical Rust timings",
        tasks_sha256=sha256sum(document),
        required_task_ids=[t["id"] for t in bound],
        created_unix=time.time(),
        node_exclusions={},
    )
    write_json_atomic(destination / "freeze.json", freeze)
    manuscript.prepare(destination / "campaign", document, freeze)
    estimate = dict(
        cells=len(bound),
        direct_reruns=0,
        reused_paired_controls=sum(
            t["baseline"]["kind"] == "new_paired_reference" for t in bound
        ),
        cached_node_hours=sum(t["historical_cached_seconds"] for t in bound) / 3600,
        note="Historical routing implementation; Python runtime and untimed audits may differ. Six workers; synthetic delay cells submitted last.",
    )
    write_json_atomic(destination / "PLAN.json", estimate)
    return estimate


def main():
    """Prepare the replay and optionally submit its automatically gated launch."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parent", type=Path, required=True)
    parser.add_argument("--repair", type=Path, action="append", default=[])
    parser.add_argument("--destination", type=Path, required=True)
    parser.add_argument("--deployment", type=Path, required=True)
    parser.add_argument("--commit", required=True)
    parser.add_argument("--tests", required=True)
    parser.add_argument("--start", action="store_true")
    args = parser.parse_args()
    result = freeze_replay(
        args.parent,
        args.repair,
        args.destination,
        args.deployment,
        args.commit,
        args.tests,
    )
    if args.start:
        result["submission"] = manuscript.start(args.destination / "campaign")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
