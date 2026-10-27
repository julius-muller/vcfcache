"""Freeze a verified candidate and submit autonomous manuscript preparation/launch."""

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
    """Persist the deployment and preparation receipt before any jobs are accepted."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--helper", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--reuse-launch-checks", type=Path)
    parser.add_argument("--reuse-cache-payloads", type=Path)
    args = parser.parse_args()
    if args.receipt.exists():
        raise FileExistsError(args.receipt)
    if subprocess.check_output(
        ["git", "diff", "--name-only", "HEAD"], cwd=REPO_ROOT, text=True
    ).strip():
        raise RuntimeError("commit tracked changes before freezing the candidate")
    commit = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, text=True
    ).strip()
    label = (
        datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        + "-"
        + commit[:9]
    )
    deployment = f"/results/deployment/manuscript-final-{label}"
    output = f"/results/manuscript-final/{label}"
    mapped = deployment.replace("/results/", "/mnt/data/slurm-results/", 1)
    paths = [
        "vcfcache",
        "resources",
        "pyproject.toml",
        "benchmarks/__init__.py",
        "benchmarks/run_pilot.py",
        "benchmarks/run_wgs_pipeline_spectrum_task.py",
        "benchmarks/engine_profile/measure.py",
        "benchmarks/engine_profile/run.py",
        "benchmarks/engine_profile/matrix.py",
        "benchmarks/engine_profile/__init__.py",
        "benchmarks/engine_profile/deploy.py",
        "benchmarks/engine_profile/review.py",
        "benchmarks/performance_sprint",
        "benchmarks/reproduction/final_clinical_campaign/cache_recipes",
        "benchmarks/reproduction/final_clinical_campaign/blocksort.frozen.awk",
        "benchmarks/figures/output/final/2026-08-25T224623Z-native-annotation-final-v3/source/campaign_metrics.tsv",
    ]
    # Avoid archiving repeated engineering traces and large result snapshots.
    paths.remove("benchmarks/performance_sprint")
    paths += subprocess.check_output(
        [
            "git",
            "ls-files",
            "benchmarks/performance_sprint/*.py",
            "benchmarks/performance_sprint/*.sh",
            "benchmarks/performance_sprint/reference-request-v2-20261027.json",
        ],
        cwd=REPO_ROOT,
        text=True,
    ).splitlines()
    with tempfile.TemporaryDirectory(prefix="vcfcache-final-deploy-") as temporary:
        root = Path(temporary)
        archive = root / "source.tar"
        with archive.open("wb") as stream:
            subprocess.run(
                ["git", "archive", commit, *paths],
                cwd=REPO_ROOT,
                stdout=stream,
                check=True,
            )
        source = root / "source"
        source.mkdir()
        with tarfile.open(archive) as bundle:
            bundle.extractall(source, filter="data")
        helper = source / "native/target/release/vcfcache-engine"
        helper.parent.mkdir(parents=True)
        shutil.copy2(args.helper, helper)
        build = dict(
            root=deployment,
            commit=commit,
            helper_version=subprocess.check_output(
                [str(args.helper.resolve()), "--version"], text=True
            ).strip(),
            files={
                str(p.relative_to(source)): sha256sum(p)
                for p in sorted(source.rglob("*"))
                if p.is_file()
            },
        )
        build_path = root / "build.json"
        write_json_atomic(build_path, build)
        remote("cd /tmp; sudo -u appuser mkdir -p " + shlex.quote(mapped))
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
        subprocess.run(
            [
                "rsync",
                "-a",
                "--rsync-path=sudo -u appuser rsync",
                str(build_path),
                f"{HOST}:{mapped}-build.json",
            ],
            check=True,
        )
        wrap = f"cd /tmp; export PYTHONPATH={shlex.quote(deployment)}; {PYTHON} -m benchmarks.performance_sprint.assemble --request {deployment}/benchmarks/performance_sprint/reference-request-v2-20261027.json --source {deployment} --output {output} --build {deployment}-build.json --submit"
        if args.reuse_launch_checks:
            wrap += " --reuse-launch-checks " + shlex.quote(
                str(args.reuse_launch_checks)
            )
        if args.reuse_cache_payloads:
            wrap += " --reuse-cache-payloads " + shlex.quote(
                str(args.reuse_cache_payloads)
            )
        command = [
            "sudo",
            "-u",
            "appuser",
            "sbatch",
            "--parsable",
            "--partition=benchmark",
            "--nodes=1",
            "--nodelist=sl-w1",
            "--ntasks=1",
            "--cpus-per-task=8",
            "--mem=28000M",
            "--time=04:00:00",
            "--chdir=/tmp",
            f"--output=/results/manuscript-final-prepare-{label}-%j.out",
            "--wrap",
            wrap,
        ]
        receipt = dict(
            build=build,
            root=output,
            preparation_command=command,
            decision="User selected A: freeze supported changes and execute full benchmark; no more optimization experiments",
            workflow="prepare all historical conditions and 72 new vep-rs pairs; six primary cells first, automatic gated remainder; durable afterany collectors",
        )
        args.receipt.parent.mkdir(parents=True, exist_ok=True)
        write_json_atomic(args.receipt, receipt)
        job = remote("cd /tmp; " + shlex.join(command)).strip()
        if not job.isdigit():
            raise RuntimeError(f"unexpected sbatch response: {job}")
        receipt["preparation_job_id"] = job
        write_json_atomic(args.receipt, receipt)
        print(
            json.dumps(
                dict(root=output, preparation_job_id=job, deployment=deployment),
                indent=2,
            )
        )


if __name__ == "__main__":
    main()
