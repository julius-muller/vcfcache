"""Freeze and submit a bounded screen from the exact completed WES conditions."""

import argparse
import json
import shlex
import subprocess
from datetime import datetime, timedelta, timezone
from pathlib import Path

from benchmarks.engine_profile.run import PYTHON
from benchmarks.run_pilot import sha256sum, write_json_atomic


def main():
    """Retain all 24 targets while scheduling only six deterministic screen cases."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--deployment", required=True)
    parser.add_argument("--commit", required=True)
    parser.add_argument("--binary", required=True)
    parser.add_argument(
        "--modes", nargs="+", default=["direct", "adaptive", "resident"]
    )
    parser.add_argument("--reference-screen")
    parser.add_argument(
        "--selection", choices=["screen", "remaining"], default="screen"
    )
    parser.add_argument("--balance-order", action="store_true")
    parser.add_argument("--exclude-ids", nargs="*", default=[])
    parser.add_argument("--no-submit", action="store_true")
    parser.add_argument("--deadline-utc")
    args = parser.parse_args()
    parent = Path(
        "/mnt/data/slurm-results/python-adaptive-manuscript-20261101-ee2162473/campaign/manifest.json"
    )
    original = json.loads(parent.read_text())
    tasks = [
        t
        for t in original["tasks"]
        if t["assay"] == "wes"
        and t["assembly"] == "GRCh38"
        and ((t["tool"], t["level"]) in {("veprs", "core"), ("fastvep", "clinical")})
    ]
    assert len(tasks) == 24
    selected = []
    for tool in ["veprs", "fastvep"]:
        order = sorted(
            (t for t in tasks if t["tool"] == tool),
            key=lambda t: (
                t["baseline"]["direct_seconds"] / t["historical_cached_seconds"],
                t["sample"],
            ),
        )
        selected.append(next(t["sample"] for t in order if t["sample"] not in selected))

    # The completed source state records the observed hit rate.
    def hit(task):
        paths = sorted(
            (parent.parent / "attempts").glob("*/" + task["id"] + "/state.json")
        )
        return json.loads(paths[-1].read_text())["cache_hit_rate"]

    order = sorted(
        (t for t in tasks if t["tool"] == "veprs"), key=lambda t: (hit(t), t["sample"])
    )
    selected.append(next(t["sample"] for t in order if t["sample"] not in selected))
    chosen = [
        t for t in tasks if (t["sample"] in selected) == (args.selection == "screen")
    ]
    assert len(chosen) == (6 if args.selection == "screen" else 18)
    # A confirmed configuration is never rerun on already completed targets.
    known = {t["id"] for t in chosen}
    assert set(args.exclude_ids) <= known
    chosen = [t for t in chosen if t["id"] not in args.exclude_ids]
    orders = {}
    for index, task in enumerate(chosen):
        shift = index % len(args.modes) if args.balance_order else 0
        orders[task["id"]] = args.modes[shift:] + args.modes[:shift]
    args.root.mkdir(parents=True, exist_ok=False)
    worker = str(args.root).replace("/mnt/data/slurm-results/", "/results/", 1)
    (args.root / "logs").mkdir()
    binary_local = Path(args.binary.replace("/results/", "/mnt/data/slurm-results/", 1))
    freeze = dict(
        tasks=chosen,
        all_targets=tasks,
        selected_samples=selected,
        source_manifest_sha256=sha256sum(parent),
        commit=args.commit,
        deployment=args.deployment,
        binary=args.binary,
        binary_sha256=sha256sum(binary_local),
        asset_attestations=original["freeze"]["asset_attestations"],
        asset_groups=original["freeze"]["annotation_asset_groups"],
        deadline_utc=args.deadline_utc
        or (datetime.now(timezone.utc) + timedelta(minutes=30)).isoformat(),
        modes=args.modes,
        excluded_ids=args.exclude_ids,
        reference_screen=args.reference_screen,
        modes_by_id=orders,
    )
    write_json_atomic(args.root / "manifest.json", freeze)
    (args.root / "READY").write_text(sha256sum(args.root / "manifest.json") + "\n")
    if args.no_submit:
        print(json.dumps(dict(state="prepared_not_submitted", root=worker)))
        return
    jobs = []
    for task in chosen:
        command = [
            PYTHON,
            "-m",
            "benchmarks.wes_pilot.run",
            worker,
            "--cell",
            task["id"],
            "--deployment",
            args.deployment,
            "--binary",
            args.binary,
            "--modes",
            *orders[task["id"]],
        ]
        wrap = (
            "cd /tmp; PYTHONPATH="
            + shlex.quote(args.deployment)
            + " "
            + shlex.join(command)
        )
        submit = [
            "sbatch",
            "--parsable",
            "--partition=benchmark",
            "--exclusive",
            "--nodes=1",
            "--ntasks=1",
            "--cpus-per-task=8",
            "--mem=28000M",
            "--time=00:20:00",
            "--output=" + worker + "/logs/" + task["id"] + "-%j.out",
            "--wrap",
            wrap,
        ]
        job = subprocess.check_output(submit, text=True).strip().split(";")[0]
        jobs.append(
            dict(
                cell=task["id"],
                sample=task["sample"],
                tool=task["tool"],
                job_id=job,
                command=submit,
            )
        )
        write_json_atomic(args.root / "submissions.json", jobs)
    print(json.dumps(dict(selected=selected, jobs=jobs), indent=2))


if __name__ == "__main__":
    main()
