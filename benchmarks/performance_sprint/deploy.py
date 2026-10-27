"""Deploy a tested sprint snapshot and submit cached-only screening jobs."""

import argparse
import datetime
import json
import shlex
import shutil
import subprocess
import tarfile
import tempfile
from pathlib import Path

from benchmarks.engine_profile.deploy import HOST, remote
from benchmarks.engine_profile.run import PYTHON
from benchmarks.run_pilot import REPO_ROOT, sha256sum, write_json_atomic


def main():
    """Record deployment identity before scheduling the immutable pilot."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--helper", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument(
        "--screen",
        choices=("initial", "trusted-hits", "memory-preflight", "native-route"),
        default="initial",
    )
    args = parser.parse_args()
    if args.receipt.exists():
        parser.error("receipt already exists")
    if subprocess.check_output(
        ["git", "diff", "--name-only", "HEAD"], cwd=REPO_ROOT, text=True
    ).strip():
        parser.error("commit tracked changes first")
    commit = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, text=True
    ).strip()
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    label = f"{stamp}-{commit[:9]}"
    target = f"/results/deployment/performance-sprint-{label}"
    mapped = target.replace("/results/", "/mnt/data/slurm-results/", 1)
    campaign = f"/results/performance-sprint/{label}"
    with tempfile.TemporaryDirectory(prefix="vcfcache-sprint-deploy-") as temp:
        directory = Path(temp)
        paths = subprocess.check_output(
            [
                "git",
                "ls-files",
                "vcfcache",
                "resources",
                "pyproject.toml",
                "benchmarks/__init__.py",
                "benchmarks/run_pilot.py",
                "benchmarks/engine_profile/*.py",
                "benchmarks/engine_profile/job.sh",
                "benchmarks/engine_profile/plot.R",
                "benchmarks/performance_sprint/*.py",
            ],
            cwd=REPO_ROOT,
            text=True,
        ).splitlines()
        archive = directory / "source.tar"
        with archive.open("wb") as out:
            subprocess.run(
                ["git", "archive", commit, *paths],
                cwd=REPO_ROOT,
                stdout=out,
                check=True,
            )
        source = directory / "source"
        source.mkdir()
        with tarfile.open(archive) as bundle:
            bundle.extractall(source, filter="data")
        binary = source / "native/target/release/vcfcache-engine"
        binary.parent.mkdir(parents=True)
        shutil.copy2(args.helper, binary)
        request = dict(
            screen=args.screen,
            profile=dict(
                root=target,
                commit=commit,
                helper_version=subprocess.check_output(
                    [str(args.helper.resolve()), "--version"], text=True
                ).strip(),
                files={
                    str(p.relative_to(source)): sha256sum(p)
                    for p in sorted(source.rglob("*"))
                    if p.is_file()
                },
            ),
        )
        remote(f"cd /tmp && sudo -u appuser mkdir {shlex.quote(mapped)}")
        subprocess.run(
            [
                "rsync",
                "-a",
                "--rsync-path=sudo -u appuser rsync",
                str(source) + "/",
                f"{HOST}:{mapped}/",
            ],
            check=True,
        )
        request_file = directory / "request.json"
        write_json_atomic(request_file, request)
        subprocess.run(
            [
                "rsync",
                "-a",
                "--rsync-path=sudo -u appuser rsync",
                str(request_file),
                f"{HOST}:{mapped}-request.json",
            ],
            check=True,
        )
        parent = "/mnt/data/slurm-results/engine-profile/20261023T190736Z-97545d739"
        mapped_campaign = campaign.replace("/results/", "/mnt/data/slurm-results/", 1)
        python = PYTHON.replace("/results/", "/mnt/data/slurm-results/", 1)
        prefix = ["sudo", "-u", "appuser", "env", f"PYTHONPATH={mapped}", python]
        prepare = [
            *prefix,
            "-m",
            "benchmarks.performance_sprint.pilot",
            "prepare",
            mapped_campaign,
            "--parent",
            parent,
            "--request",
            mapped + "-request.json",
        ]
        submit = [
            *prefix,
            "-m",
            "benchmarks.engine_profile.run",
            "submit",
            mapped_campaign,
        ]
        receipt = dict(
            campaign=campaign,
            deployment=target,
            request=request,
            prepare=prepare,
            submit=submit,
        )
        args.receipt.parent.mkdir(parents=True, exist_ok=True)
        write_json_atomic(args.receipt, receipt)
        remote("cd /tmp && " + shlex.join(prepare))
        result = remote("cd /tmp && " + shlex.join(submit))
        receipt["submission"] = json.loads(result)
        write_json_atomic(args.receipt, receipt)
        print(json.dumps(receipt["submission"], indent=2))


if __name__ == "__main__":
    main()
