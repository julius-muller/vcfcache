"""Immutable, separately timed WES pilot cases; preparation and comparison untimed."""

import argparse
import json
import os
import shlex
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


def run(root, identifier, deployment, modes, binary):
    """One case owns isolated steps and never overwrites an existing result."""
    manifest = json.loads((root / "manifest.json").read_text())
    assert sha256sum(root / "manifest.json") == (root / "READY").read_text().strip()
    task = next(t for t in manifest["tasks"] if t["id"] == identifier)
    evidence = root / "cases" / identifier
    evidence.mkdir(parents=True, exist_ok=False)
    state = dict(
        state="running",
        job_id=os.environ.get("SLURM_JOB_ID"),
        started_unix=time.time(),
        modes=modes,
        rows=[],
    )
    write_json_atomic(evidence / "state.json", state)
    with tempfile.TemporaryDirectory(
        prefix="wes-pilot-", dir="/mnt/data/tmp"
    ) as temporary:
        scratch = Path(temporary)
        try:
            attestations = manifest["asset_attestations"]
            source = verified_file(task["input"], attestations)
            cache = node_cache(task["cache"])
            for name, digest in task["cache"]["files"].items():
                if task["cache"].get("storage") != "node_local":
                    verified_file(
                        dict(path=str(cache / name), sha256=digest), attestations
                    )
            for asset in (
                task["annotation_assets"]
                + manifest["asset_groups"][task["annotation_asset_group"]]
            ):
                verified_file(asset, attestations)
            bed = verified_file(task["preprocess"]["bed"], attestations)
            prepared = scratch / "input.bcf"
            pipeline = f"bcftools view -R {shlex.quote(str(bed))} {shlex.quote(str(source))} -Ou | bcftools norm -m- -Ou | bcftools view -e 'ALT=\"*\"' -Ou | bcftools norm --rm-dup exact -Ob -o {prepared} -W --threads 8"
            subprocess.run(["bash", "-euo", "pipefail", "-c", pipeline], check=True)
            n = int(subprocess.check_output(["bcftools", "index", "-n", str(prepared)]))
            assert n == task["input"]["records"]
            params = dict(task["parameters"], threads=8, temp_dir=str(scratch / "tmp"))
            (scratch / "tmp").mkdir()
            (scratch / "params.yaml").write_text(yaml.safe_dump(params))
            reference = None
            if manifest.get("reference_screen"):
                origin = Path(manifest["reference_screen"]) / "cases" / identifier
                previous = json.loads((origin / "state.json").read_text())
                assert previous["state"] == "complete"
                metrics = json.loads((origin / "direct/resources.json").read_text())
                assert metrics["exit_code"] == 0 and len(metrics["affinity"]) == 8
                state["reused_direct"] = dict(
                    path=str(origin),
                    job_id=previous["job_id"],
                    wall_seconds=metrics["wall_seconds"],
                )
                binding = task["baseline"]["direct_reuse"]
                reference = Path(binding["attempt"]) / "direct.bcf"
                assert sha256sum(reference) == binding["files"]["direct.bcf"]
            for mode in modes:
                directory = evidence / mode
                directory.mkdir()
                output = scratch / f"{mode}.bcf"
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
                    str(scratch / "params.yaml"),
                    "--engine",
                    (
                        "legacy"
                        if mode == "direct"
                        else ("adaptive" if mode == "adaptive" else "native")
                    ),
                    "--statistics",
                    "light",
                    "--stats-dir",
                    str(scratch / f"stats-{mode}"),
                ]
                if mode == "direct":
                    command += ["--uncached"]
                elif mode != "adaptive":
                    command += ["--engine-memory-mb", "2048"]
                environment = dict(
                    os.environ,
                    PYTHONPATH=deployment,
                    VCFCACHE_ROOT=deployment,
                    VCFCACHE_ENGINE=binary,
                    LC_ALL="C",
                )
                for key in list(environment):
                    if (
                        key.startswith("VCFCACHE_EXPERIMENT_")
                        or key == "VCFCACHE_PROFILE_DIR"
                    ):
                        environment.pop(key)
                if mode == "coalesced":
                    environment["VCFCACHE_EXPERIMENT_COALESCE_CHUNKS"] = "1"
                launch = [
                    "srun",
                    "--exclusive",
                    "--exact",
                    "--cpu-bind=cores",
                    "-n1",
                    "-c8",
                    PYTHON,
                    str(Path(deployment) / "benchmarks/engine_profile/measure.py"),
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
                        parameters=params,
                        deployment=deployment,
                        binary=binary,
                        mode=mode,
                    ),
                )
                with (directory / "run.log").open("w") as log:
                    subprocess.run(
                        launch, env=environment, stdout=log, stderr=log, check=True
                    )
                metrics = json.loads((directory / "resources.json").read_text())
                assert (
                    metrics["exit_code"] == 0
                    and len(metrics["affinity"]) == 8
                    and metrics["cgroup_after"]
                )
                assert (
                    int(
                        subprocess.check_output(
                            ["bcftools", "index", "-n", str(output)]
                        )
                    )
                    == n
                )
                stats_paths = list(
                    (scratch / f"stats-{mode}").rglob("compare_stats.yaml")
                )
                assert len(stats_paths) == 1
                stats = yaml.safe_load(stats_paths[0].read_text())
                write_json_atomic(directory / "statistics.json", stats)
                if mode == "direct":
                    reference = output
                    binding = task["baseline"].get("direct_reuse")
                    if binding:
                        historical = Path(binding["attempt"]) / "direct.bcf"
                        comparison = strict_semantic_compare(output, historical)
                        assert comparison[
                            "semantic_pass"
                        ], "fresh direct differs from retained direct"
                        write_json_atomic(
                            directory / "against_historical_direct.json", comparison
                        )
                else:
                    assert reference is not None
                    comparison = strict_semantic_compare(output, reference)
                    write_json_atomic(directory / "against_direct.json", comparison)
                    assert comparison["semantic_pass"], "strict output mismatch"
                    if mode != "adaptive":
                        assert (
                            stats.get("effective_engine") == "native"
                        ), "resident candidate fell back"
                    adaptive_output = scratch / "adaptive.bcf"
                    if mode == "coalesced":
                        assert stats.get("coalesced_chunks") is True
                    if mode != "adaptive" and adaptive_output.exists():
                        comparison = strict_semantic_compare(output, adaptive_output)
                        write_json_atomic(
                            directory / "against_adaptive.json", comparison
                        )
                        assert comparison["semantic_pass"]
                state["rows"].append(
                    dict(
                        mode=mode,
                        wall_seconds=metrics["wall_seconds"],
                        output_sha256=sha256sum(output),
                        records=n,
                    )
                )
                write_json_atomic(evidence / "state.json", state)
            state["state"] = "complete"
        except Exception as error:
            state.update(
                state="failed", error=str(error), traceback=traceback.format_exc()
            )
        finally:
            for path in scratch.rglob("*"):
                if path.is_file() and path.suffix in {
                    ".json",
                    ".yaml",
                    ".txt",
                    ".sh",
                    ".stderr",
                    ".stdout",
                }:
                    # Preserve small execution evidence, not variant payloads.
                    destination = evidence / "cli" / path.relative_to(scratch)
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    if path.stat().st_size < 8 * 1024 * 1024:
                        destination.write_bytes(path.read_bytes())
            state["finished_unix"] = time.time()
            write_json_atomic(evidence / "state.json", state)
    return state


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--cell", required=True)
    parser.add_argument("--deployment", required=True)
    parser.add_argument("--binary", required=True)
    parser.add_argument(
        "--modes", nargs="+", default=["direct", "adaptive", "resident"]
    )
    args = parser.parse_args()
    result = run(args.root, args.cell, args.deployment, args.modes, args.binary)
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result["state"] == "complete" else 1)
