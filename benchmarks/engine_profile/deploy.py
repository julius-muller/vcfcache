"""Deploy frozen source/helper pairs and submit only the preparation job.

The real matrix is submitted explicitly after READY and preflight inspection.
This script does not build binaries or modify existing deployments/campaigns.
"""

import argparse
import datetime
import json
import shlex
import subprocess
import tarfile
import tempfile
from pathlib import Path

from benchmarks.run_pilot import REPO_ROOT, sha256sum, write_json_atomic

HOST = "ubuntu@10.133.241.128"
PYTHON = "/results/deployment/vcfcache-native-pilot-20261023c/.venv/bin/python"


def remote(command):
    """Execute a scoped controller command without an interactive session."""
    return subprocess.check_output(
        ["ssh", "-o", "BatchMode=yes", HOST, command], text=True
    )


def main():
    """Freeze and deploy two snapshots, then durably record preparation acceptance."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-helper", type=Path, required=True)
    parser.add_argument("--profile-helper", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    if args.receipt.exists():
        parser.error("receipt already exists")
    if subprocess.check_output(
        ["git", "status", "--porcelain", "--untracked-files=no"],
        cwd=REPO_ROOT,
        text=True,
    ).strip():
        parser.error("commit tracked changes before freezing deployments")
    head = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, text=True
    ).strip()
    baseline = subprocess.check_output(
        ["git", "rev-parse", "296855d72"], cwd=REPO_ROOT, text=True
    ).strip()
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    label = f"{stamp}-{head[:9]}"
    campaign = f"/results/engine-profile/{label}"
    request = dict(builds={}, created_utc=stamp, campaign=campaign)
    with tempfile.TemporaryDirectory(prefix="vcfcache-profile-deploy-") as temp:
        directory = Path(temp)
        for name, commit, helper in (
            ("baseline", baseline, args.baseline_helper),
            ("profile", head, args.profile_helper),
        ):
            target = f"/results/deployment/engine-profile-{label}-{name}"
            mapped = target.replace("/results/", "/mnt/data/slurm-results/", 1)
            local = directory / name
            local.mkdir()
            archive = directory / f"{name}.tar"
            paths = ["vcfcache", "resources", "pyproject.toml"]
            if name == "profile":
                paths += [
                    "benchmarks/__init__.py",
                    "benchmarks/run_pilot.py",
                    "benchmarks/engine_profile",
                ]
            with archive.open("wb") as stream:
                subprocess.run(
                    ["git", "archive", commit, *paths],
                    cwd=REPO_ROOT,
                    stdout=stream,
                    check=True,
                )
            with tarfile.open(archive) as bundle:
                bundle.extractall(local, filter="data")
            binary = local / "native/target/release/vcfcache-engine"
            binary.parent.mkdir(parents=True)
            import shutil

            shutil.copy2(helper, binary)
            request["builds"][name] = dict(
                root=target,
                commit=commit,
                helper_version=subprocess.check_output(
                    [str(helper.resolve()), "--version"], text=True
                ).strip(),
                files={
                    str(p.relative_to(local)): sha256sum(p)
                    for p in sorted(local.rglob("*"))
                    if p.is_file()
                },
            )
            remote(f"cd /tmp && sudo -u appuser mkdir {shlex.quote(mapped)}")
            subprocess.run(
                [
                    "rsync",
                    "-a",
                    "--rsync-path=sudo -u appuser rsync",
                    str(local) + "/",
                    f"{HOST}:{mapped}/",
                ],
                check=True,
            )
        request_path = directory / "request.json"
        write_json_atomic(request_path, request)
        remote_request = (
            f"/mnt/data/slurm-results/deployment/engine-profile-{label}-request.json"
        )
        subprocess.run(
            [
                "rsync",
                "-a",
                "--rsync-path=sudo -u appuser rsync",
                str(request_path),
                f"{HOST}:{remote_request}",
            ],
            check=True,
        )
        worker_request = remote_request.replace(
            "/mnt/data/slurm-results/", "/results/", 1
        )
        profile = request["builds"]["profile"]["root"]
        wrap = f"cd /tmp; export PYTHONPATH={shlex.quote(profile)}; {PYTHON} -m benchmarks.engine_profile.run prepare {shlex.quote(campaign)} --request {shlex.quote(worker_request)}"
        command = [
            "sudo",
            "-u",
            "appuser",
            "sbatch",
            "--parsable",
            "--partition=benchmark",
            "--nodelist=sl-w2",
            "--nodes=1",
            "--ntasks=1",
            "--cpus-per-task=16",
            "--mem=28000M",
            "--exclusive",
            "--time=02:00:00",
            "--job-name=engine-profile-prepare",
            f"--output=/results/deployment/engine-profile-{label}-prepare-%j.out",
            "--wrap",
            wrap,
        ]
        # Save provenance before asking Slurm to avoid losing deployment paths.
        receipt = dict(request=request, prepare_command=command, controller=HOST)
        args.receipt.parent.mkdir(parents=True, exist_ok=True)
        write_json_atomic(args.receipt, receipt)
        job = remote("cd /tmp && " + shlex.join(command)).strip().split(";")[0]
        if not job.isdigit():
            raise RuntimeError(f"unexpected Slurm response: {job}")
        receipt["prepare_job_id"] = job
        write_json_atomic(args.receipt, receipt)
        print(
            json.dumps(
                dict(campaign=campaign, prepare_job_id=job, receipt=str(args.receipt)),
                indent=2,
            )
        )


if __name__ == "__main__":
    main()
