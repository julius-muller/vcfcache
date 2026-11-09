"""Freeze, submit and review six same-worker Python/Rust/direct blocks."""

import argparse
import csv
import json
import os
import shlex
import shutil
import socket
import statistics
import subprocess
import tarfile
import tempfile
import time
import traceback
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import yaml

from benchmarks.run_pilot import sha256sum, strict_semantic_compare, write_json_atomic

PYTHON = "/results/deployment/vcfcache-native-pilot-20261023c/.venv/bin/python"
PARENT = "/mnt/data/slurm-results/python-adaptive-manuscript-20261101-ee2162473/campaign/manifest.json"
MODES = ("direct", "python", "rust")
SELECTION = (
    ("veprs", "core", "KPGP-00139", "wgs", "sl-w1"),
    ("fastvep", "clinical", "KPGP-00139", "wgs", "sl-w2"),
    ("veprs", "core", "LP6005519-DNA_A06", "wgs", "sl-w3"),
    ("fastvep", "clinical", "LP6005519-DNA_A06", "wgs", "sl-w4"),
    ("veprs", "core", "KPGP-00139", "wes", "sl-w5"),
    ("fastvep", "clinical", "KPGP-00139", "wes", "sl-w6"),
)
CANDIDATE_SELECTION = (
    ("veprs", "core", "KPGP-00139", "wes", "sl-w1"),
    ("fastvep", "clinical", "KPGP-00139", "wes", "sl-w2"),
    ("veprs", "core", "LP6005592-DNA_C05", "wes", "sl-w3"),
    ("fastvep", "clinical", "LP6005592-DNA_C05", "wes", "sl-w4"),
    ("veprs", "core", "KPGP-00139", "wgs", "sl-w5"),
    ("fastvep", "clinical", "KPGP-00139", "wgs", "sl-w6"),
)
CANDIDATES = {
    "direct": dict(backend="python", engine="legacy", prepared=True),
    "python-windows2": dict(
        backend="python",
        engine="adaptive",
        prepared=True,
        workers=2,
        lookup="windows",
        thread_budget=True,
    ),
    "python-auto": dict(
        backend="python",
        engine="adaptive",
        prepared=True,
        lookup="auto",
        thread_budget=True,
    ),
    "python-contigs2": dict(
        backend="python",
        engine="adaptive",
        prepared=True,
        workers=2,
        lookup="contigs",
        thread_budget=True,
    ),
    "rust-adaptive": dict(backend="rust", engine="adaptive", prepared=True),
    "rust-native-w4-c0": dict(
        backend="rust", engine="native", prepared=True, workers=4, codec=0
    ),
    "rust-native-w4-c4": dict(
        backend="rust", engine="native", prepared=True, workers=4, codec=4
    ),
    "rust-native-w8-c0": dict(
        backend="rust", engine="native", prepared=True, workers=8, codec=0
    ),
    "rust-native-coalesced-w4-c4": dict(
        backend="rust",
        engine="native",
        prepared=True,
        workers=4,
        codec=4,
        coalesced=True,
    ),
}


def candidate_modes(assay):
    """Prespecify promising configurations, never silently label WGS fallback native."""
    if assay == "wes":
        return [
            "direct",
            "python-windows2",
            "rust-native-w4-c0",
            "rust-native-w4-c4",
            "rust-native-w8-c0",
            "rust-native-coalesced-w4-c4",
        ]
    if assay == "wgs":
        return ["direct", "python-auto", "python-contigs2", "rust-adaptive"]
    raise ValueError("Candidate pilot supports WES/WGS only")


def mode_spec(data, mode):
    """Read an explicit frozen condition; retain the original pilot's schema."""
    if "configurations" in data:
        return data["configurations"][mode]
    if mode not in MODES:
        raise ValueError("Unknown confirmation mode")
    return dict(
        backend="rust" if mode == "rust" else "python",
        engine="legacy" if mode == "direct" else "adaptive",
        prepared=False,
    )


def mapped(path):
    """Map the controller mount to its worker-visible spelling."""
    return str(path).replace("/mnt/data/slurm-results/", "/results/", 1)


def controller(path):
    """Map worker-visible paths for controller-only preparation."""
    return Path(str(path).replace("/results/", "/mnt/data/slurm-results/", 1))


def verify(path, digest):
    """Fail rather than silently use changed code, inputs or cache payloads."""
    path = Path(path)
    if not path.is_file() or sha256sum(path) != digest:
        raise RuntimeError(f"Frozen file missing or changed: {path}")
    return path


def definition(root):
    """Require the immutable definition before run, submit or review."""
    verify(root / "manifest.json", (root / "READY").read_text().strip())
    return json.loads((root / "manifest.json").read_text())


def order_for(repetition, offset, modes=MODES):
    """Rotate conditions; three-mode designs balance every position exactly."""
    shifts = (0, len(modes) // 2, len(modes) - 1) if len(modes) > 3 else (0, 1, 2)
    shift = (shifts[repetition - 1] + offset) % len(modes)
    return list(modes[shift:] + modes[:shift])


def command_for(mode, cache, source, output, params, stats, configuration=None):
    """Only the cached backend differs; direct uses the Python legacy CLI."""
    spec = configuration or mode_spec({}, mode)
    command = [
        PYTHON,
        "-c",
        "from vcfcache.cli import main; main()",
        "annotate",
        "-a",
        str(cache),
        "-i",
        str(source),
        "-o",
        str(output),
        "-y",
        str(params),
        "--engine",
        spec["engine"],
        "--statistics",
        "light",
        "--stats-dir",
        str(stats),
        "--preserve-unannotated",
        "-v",
    ]
    if mode == "direct":
        command.append("--uncached")
    if spec.get("workers"):
        command += ["--engine-workers", str(spec["workers"])]
    if spec["engine"] == "native":
        command += ["--engine-memory-mb", "2048"]
    return command


def freeze(args):
    """Bind known candidates and identical dependencies before any submission."""
    parent = json.loads(Path(PARENT).read_text())
    py = json.loads(
        controller("/results/adaptive-hts125-20261107/manifest.json").read_text()
    )
    rs = json.loads(
        controller("/results/rust-hts125-20261107/manifest.json").read_text()
    )
    backends = {
        "python": dict(
            path=py["deployment"],
            commit=py["commit"],
            files=py["deployment_hashes"]["prototype"],
        ),
        "rust": dict(
            path=rs["deployment"], commit=rs["commit"], files=rs["source_hashes"]
        ),
    }
    for backend in backends.values():
        for path, digest in backend["files"].items():
            verify(controller(backend["path"]) / path, digest)
    tools = dict(bcftools=rs["bcftools"]["latest"], helper=rs["helpers"]["latest"])
    for tool in tools.values():
        verify(controller(tool["path"]), tool["sha256"])
    cases, tasks = [], {}
    profile = getattr(args, "profile", "default")
    selection = CANDIDATE_SELECTION if profile == "candidates" else SELECTION
    for i, (tool, level, sample, assay, node) in enumerate(selection):
        matches = [
            t
            for t in parent["tasks"]
            if (t["tool"], t["level"], t["sample"], t["assay"], t["assembly"])
            == (tool, level, sample, assay, "GRCh38")
        ]
        if len(matches) != 1:
            raise RuntimeError(f"Nonunique source task: {selection[i]}")
        task = matches[0]
        tasks[task["id"]] = task
        cases.append(
            dict(
                id=f"{tool}-{assay}-{sample}",
                task_id=task["id"],
                node=node,
                offset=i,
                modes=(
                    candidate_modes(assay) if profile == "candidates" else list(MODES)
                ),
            )
        )
    runner_files = {
        str(p.relative_to(controller(args.deployment))): sha256sum(p)
        for p in sorted(controller(args.deployment).rglob("*.py"))
    }
    args.root.mkdir(parents=True, exist_ok=False)
    (args.root / "logs").mkdir()
    data = dict(
        profile=profile,
        created_unix=time.time(),
        deployment=args.deployment,
        commit=args.commit,
        runner_files=runner_files,
        backends=backends,
        tools=tools,
        parent_sha256=sha256sum(Path(PARENT)),
        tasks=list(tasks.values()),
        cases=cases,
        asset_groups=parent["freeze"]["annotation_asset_groups"],
        cpus=8,
        repetitions=3,
        modes=MODES,
        design="Same-worker fresh direct; three Latin-square repeats; unchanged default preprocessing in all arms; shared version 1.25; private node-local cache and scratch; no OS cache flushing",
        no_automatic_cohort_launch=True,
    )
    if profile == "candidates":
        data.update(
            configurations=CANDIDATES,
            modes=list(CANDIDATES),
            design="Prespecified best-candidate subset; three rotated repeats; shared prepared-input eligibility; exact biallelic matching; fixed eight CPUs and external/linked HTSlib 1.25; local cache/scratch; no native WGS or silent fallback; no OS cache flushing",
        )
    write_json_atomic(args.root / "manifest.json", data)
    (args.root / "READY").write_text(sha256sum(args.root / "manifest.json") + "\n")
    print(
        json.dumps(
            dict(
                root=str(args.root),
                cases=len(cases),
                timed_invocations=3 * sum(len(c["modes"]) for c in cases),
            )
        )
    )


def submit(root):
    """Submit each exclusive block once, then an automatic after-any review."""
    data = definition(root)
    if (root / "submission_intent.json").exists():
        raise RuntimeError(
            "Submission already attempted; inspect saved receipts rather than duplicating jobs"
        )
    write_json_atomic(
        root / "submission_intent.json",
        dict(started_unix=time.time(), cases=data["cases"]),
    )
    receipts = []
    for case in data["cases"]:
        launch = [
            PYTHON,
            "-m",
            "benchmarks.engine_confirmation.run",
            "run",
            mapped(root),
            "--case",
            case["id"],
        ]
        command = [
            "sbatch",
            "--parsable",
            "--partition=benchmark",
            "--exclusive",
            "--nodes=1",
            "--ntasks=1",
            "--cpus-per-task=8",
            "--mem=28000M",
            "--time=03:00:00",
            "--chdir=/tmp",
            "--nodelist=" + case["node"],
            "--job-name=confirm-" + case["id"],
            "--output=" + mapped(root / "logs") + "/" + case["id"] + "-%j.out",
            "--wrap",
            "PYTHONPATH=" + shlex.quote(data["deployment"]) + " " + shlex.join(launch),
        ]
        job = subprocess.check_output(command, text=True).strip().split(";")[0]
        if not job.isdigit():
            raise RuntimeError("Unrecognized sbatch receipt: " + job)
        receipts.append(
            dict(case=case["id"], node=case["node"], job_id=job, command=command)
        )
        write_json_atomic(root / "submissions.json", receipts)
    launch = [
        PYTHON,
        "-m",
        "benchmarks.engine_confirmation.run",
        "collect",
        mapped(root),
        "--terminal",
    ]
    command = [
        "sbatch",
        "--parsable",
        "--partition=benchmark",
        "--ntasks=1",
        "--cpus-per-task=1",
        "--mem=3000M",
        "--time=00:15:00",
        "--chdir=/tmp",
        "--nodelist=sl-w1",
        "--dependency=afterany:" + ":".join(r["job_id"] for r in receipts),
        "--job-name=confirm-review",
        "--output=" + mapped(root / "logs/review-%j.out"),
        "--wrap",
        "PYTHONPATH=" + shlex.quote(data["deployment"]) + " " + shlex.join(launch),
    ]
    review = subprocess.check_output(command, text=True).strip().split(";")[0]
    write_json_atomic(
        root / "review_submission.json", dict(job_id=review, command=command)
    )
    print(json.dumps(dict(jobs=receipts, review=review)))


def wrapper(directory, actual):
    """Use one pinned bcftools and thread supported recipe commands consistently."""
    directory.mkdir()
    target = directory / "bcftools"
    target.write_text("""#!/bin/bash
set -euo pipefail
actual=ACTUAL
case "${1-}" in
 annotate|view|filter|concat|index|sort)
  subcommand=$1; shift
  for arg in "$@"; do
   case "$arg" in --threads|--threads=*|-@*) exec "$actual" "$subcommand" "$@";; esac
  done
  exec "$actual" "$subcommand" "$@" --threads 8;;
esac
exec "$actual" "$@"
""".replace("ACTUAL", shlex.quote(str(actual))))
    target.chmod(0o755)
    return target


def local_cache(task, scratch):
    """Stage the same immutable cache bytes locally, avoiding concurrent NFS lookups."""
    original = Path(task["cache"]["directory"])
    root = scratch / "cache-root"
    destination = root / "cache" / original.name
    destination.mkdir(parents=True)
    (root / "workflow").mkdir()
    (root / "blueprint").mkdir()
    shutil.copy2(
        original.parents[1] / "workflow/init.yaml", root / "workflow/init.yaml"
    )
    for name, digest in task["cache"]["files"].items():
        shutil.copy2(original / name, destination / name)
        verify(destination / name, digest)
    return destination


def environment_for(data, mode, binary):
    """Select a frozen backend without inheriting experimental or profiling flags."""
    environment = {
        k: v
        for k, v in os.environ.items()
        if not k.startswith(("VCFCACHE_", "PYTHONPATH"))
    }
    spec = mode_spec(data, mode)
    backend = data["backends"][spec["backend"]]
    environment.update(
        PYTHONPATH=backend["path"],
        VCFCACHE_ROOT=backend["path"],
        PATH=str(binary.parent) + ":" + os.environ["PATH"],
        LC_ALL="C",
        PYTHONDONTWRITEBYTECODE="1",
        RAYON_NUM_THREADS="8",
        OMP_NUM_THREADS="8",
    )
    if spec["backend"] == "rust":
        environment["VCFCACHE_ENGINE"] = data["tools"]["helper"]["path"]
    if spec.get("prepared"):
        environment["VCFCACHE_PILOT_PREPARED"] = "1"
    if "lookup" in spec:
        environment["VCFCACHE_PILOT_LOOKUP"] = spec["lookup"]
    if spec.get("thread_budget"):
        environment["VCFCACHE_PILOT_THREAD_BUDGET"] = "1"
    if "codec" in spec:
        environment["VCFCACHE_EXPERIMENT_CODEC_THREADS"] = str(spec["codec"])
    if spec.get("coalesced"):
        environment["VCFCACHE_EXPERIMENT_COALESCE_CHUNKS"] = "1"
    return environment


def preflight_prepared(bcftools, source, tag):
    """Require the same eligible unmodified BCF domain for every prepared arm."""
    started = time.monotonic()
    header = subprocess.check_output([bcftools, "view", "-h", str(source)], text=True)
    has_tag = f"##INFO=<ID={tag}," in header
    fmt = "%CHROM\t%POS\t%REF\t%ALT" + (f"\t%INFO/{tag}" if has_tag else "") + "\n"
    locus, seen, records = None, set(), 0
    with tempfile.TemporaryFile() as errors:
        process = subprocess.Popen(
            [bcftools, "query", "-f", fmt, str(source)],
            stdout=subprocess.PIPE,
            stderr=errors,
            text=True,
        )
        try:
            with process.stdout:
                for line in process.stdout:
                    fields = line.rstrip("\n").split("\t")
                    contig, pos, ref, alt = fields[:4]
                    if "," in alt or alt in {"*", ".", ""}:
                        raise ValueError(
                            "Prepared pilot requires biallelic, star-free input; no silent conversion"
                        )
                    if has_tag and fields[4] not in {"", "."}:
                        raise ValueError(
                            "Prepared input already contains the required annotation tag"
                        )
                    if locus != (contig, pos):
                        locus, seen = (contig, pos), set()
                    if (ref, alt) in seen:
                        raise ValueError(
                            "Prepared input contains duplicate variant keys"
                        )
                    seen.add((ref, alt))
                    records += 1
            if process.wait():
                errors.seek(0)
                raise RuntimeError(errors.read().decode(errors="replace"))
        finally:
            if process.poll() is None:
                process.terminate()
                process.wait()
    return dict(
        records=records,
        biallelic=True,
        star_free=True,
        unique_keys=True,
        unannotated=True,
        wall_seconds=time.monotonic() - started,
    )


def audit_pair(pair):
    """Reuse the strict complete-record/header comparator, with no ignored fields."""
    return strict_semantic_compare(Path(pair[0]), Path(pair[1]))


def run(root, identifier):
    """Checkpoint all timings and audits; retain owned scratch even on failure."""
    data = definition(root)
    case = next(c for c in data["cases"] if c["id"] == identifier)
    evidence = root / "cases" / identifier
    evidence.mkdir(parents=True, exist_ok=False)
    state = dict(
        state="running",
        node=socket.gethostname(),
        job_id=os.environ.get("SLURM_JOB_ID"),
        phase="preflight",
        started_unix=time.time(),
        rows=[],
        audits={},
    )
    write_json_atomic(evidence / "state.json", state)
    scratch = Path(tempfile.mkdtemp(prefix="engine-confirmation-", dir="/mnt/data/tmp"))
    state["scratch"] = str(scratch)
    try:
        hardware = dict(
            lscpu=json.loads(subprocess.check_output(["lscpu", "-J"], text=True)),
            kernel=subprocess.check_output(["uname", "-a"], text=True).strip(),
            apptainer=subprocess.check_output(
                ["apptainer", "version"], text=True
            ).strip(),
            mounts=json.loads(
                subprocess.check_output(
                    ["findmnt", "-J", "-T", str(scratch), "-o", "SOURCE,FSTYPE,TARGET"],
                    text=True,
                )
            ),
            free_bytes=shutil.disk_usage(scratch).free,
            boot_id=Path("/proc/sys/kernel/random/boot_id").read_text().strip(),
            machine_id=Path("/etc/machine-id").read_text().strip(),
            load_average=os.getloadavg(),
        )
        write_json_atomic(evidence / "hardware.json", hardware)
        assert state["node"] == case["node"]
        assert any(
            "Broadwell" in x["data"]
            for x in hardware["lscpu"]["lscpu"]
            if x["field"].rstrip(":") == "Model name"
        )
        assert (
            "7.0.0-43-generic" in hardware["kernel"]
            and hardware["apptainer"] == "1.5.5"
        )
        assert hardware["free_bytes"] > 60 * 1024**3
        assert all(
            x["fstype"] not in {"nfs", "nfs4"}
            for x in hardware["mounts"]["filesystems"]
        )
        for name, digest in data["runner_files"].items():
            verify(Path(data["deployment"]) / name, digest)
        for backend in data["backends"].values():
            for name, digest in backend["files"].items():
                verify(Path(backend["path"]) / name, digest)
        for tool in data["tools"].values():
            verify(tool["path"], tool["sha256"])
        task = next(t for t in data["tasks"] if t["id"] == case["task_id"])
        write_json_atomic(evidence / "task.json", task)
        for asset in (
            task["annotation_assets"]
            + data["asset_groups"][task["annotation_asset_group"]]
        ):
            verify(asset["path"], asset["sha256"])
        state["phase"] = "stage_local_cache_and_input"
        write_json_atomic(evidence / "state.json", state)
        cache = local_cache(task, scratch)
        bcftools = data["tools"]["bcftools"]["path"]
        source = scratch / "input.bcf"
        original = verify(task["input"]["path"], task["input"]["sha256"])
        extract = [
            bcftools,
            "view",
            str(original),
            "-Ob1",
            "-W",
            "--threads",
            "8",
            "-o",
            str(source),
        ]
        if task.get("preprocess", {}).get("bed"):
            bed = task["preprocess"]["bed"]
            extract += ["-R", str(verify(bed["path"], bed["sha256"]))]
        subprocess.run(extract, check=True)
        assert (
            int(subprocess.check_output([bcftools, "index", "-n", str(source)]))
            == task["input"]["records"]
        )
        modes = case.get("modes", list(MODES))
        if any(mode_spec(data, mode).get("prepared") for mode in modes):
            state["phase"] = "shared_prepared_input_eligibility"
            write_json_atomic(evidence / "state.json", state)
            config = yaml.safe_load((cache / "annotation.yaml").read_text())
            proof = preflight_prepared(
                bcftools, source, config["must_contain_info_tag"]
            )
            assert proof["records"] == task["input"]["records"]
            write_json_atomic(evidence / "prepared_preflight.json", proof)
        state["input_sha256"] = sha256sum(source)
        source.chmod(0o444)
        binary = wrapper(scratch / "bin", bcftools)
        (scratch / "tmp").mkdir()
        parameters = dict(
            task["parameters"],
            bcftools_cmd=str(binary),
            threads=8,
            temp_dir=str(scratch / "tmp"),
        )
        params = scratch / "params.yaml"
        params.write_text(yaml.safe_dump(parameters))
        write_json_atomic(evidence / "parameters.json", parameters)
        outputs = {}
        for repetition in range(1, 4):
            for position, mode in enumerate(
                order_for(repetition, case["offset"], modes), 1
            ):
                spec = mode_spec(data, mode)
                cell = f"r{repetition:02d}-{mode}"
                directory = evidence / cell
                directory.mkdir()
                output, stats = scratch / (cell + ".bcf"), scratch / (cell + "-stats")
                command = command_for(mode, cache, source, output, params, stats, spec)
                environment = environment_for(data, mode, binary)
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
                        / "benchmarks/input_preserving_pilot/measure.py"
                    ),
                    "--directory",
                    str(directory),
                    "--",
                    "env",
                    "PYTHONPATH=" + environment["PYTHONPATH"],
                    *command,
                ]
                write_json_atomic(
                    directory / "command.json",
                    dict(
                        command=launch,
                        backend=mode,
                        configuration=spec,
                        backend_commit=data["backends"][spec["backend"]]["commit"],
                        environment={
                            k: v
                            for k, v in environment.items()
                            if k
                            in {
                                "PATH",
                                "PYTHONPATH",
                                "VCFCACHE_ROOT",
                                "VCFCACHE_ENGINE",
                                "RAYON_NUM_THREADS",
                                "OMP_NUM_THREADS",
                            }
                            or k.startswith(("VCFCACHE_PILOT_", "VCFCACHE_EXPERIMENT_"))
                        },
                    ),
                )
                state.update(phase=cell, last_progress_unix=time.time())
                write_json_atomic(evidence / "state.json", state)
                # The measuring script needs the frozen runner, while the CLI
                # child must import its selected backend, not runner/main code.
                measure_env = dict(
                    environment,
                    PYTHONPATH=data["deployment"] + ":" + environment["PYTHONPATH"],
                )
                with (directory / "run.log").open("w") as log:
                    subprocess.run(
                        launch, env=measure_env, stdout=log, stderr=log, check=True
                    )
                resources = json.loads((directory / "resources.json").read_text())
                assert resources["exit_code"] == 0 and len(resources["affinity"]) == 8
                count = int(
                    subprocess.check_output([bcftools, "index", "-n", str(output)])
                )
                if count != task["input"]["records"]:
                    raise RuntimeError(
                        f"Truncated output: {count} != {task['input']['records']}"
                    )
                files = list(stats.rglob("compare_stats.yaml"))
                assert len(files) == 1
                counters = yaml.safe_load(files[0].read_text())
                assert counters["variant_counts"]["total_output"] == count
                if spec["backend"] == "rust":
                    assert counters["effective_engine"] == spec[
                        "engine"
                    ] and not counters.get("fallback_reason")
                if mode != "direct" and spec.get("prepared"):
                    assert counters["pilot_prepared"]
                if spec["engine"] == "native":
                    assert counters["workers"] == spec["workers"]
                    assert counters["lookup_codec_threads"] == spec["codec"]
                    assert bool(counters["coalesced_chunks"]) == bool(
                        spec.get("coalesced")
                    )
                if spec["backend"] == "python" and "lookup" in spec:
                    expected_access = (
                        "windows" if spec["lookup"] == "windows" else "contigs"
                    )
                    assert counters["cache_access"] == expected_access
                    if spec.get("workers"):
                        assert counters["shards"] <= spec["workers"]
                write_json_atomic(directory / "statistics.json", counters)
                outputs[cell] = str(output)
                row = dict(
                    cell=cell,
                    repetition=repetition,
                    position=position,
                    mode=mode,
                    seconds=resources["wall_seconds"],
                    records=count,
                    missing=counters["variant_counts"].get("missing_variants"),
                    cpu_seconds=resources["cpu_seconds"],
                    peak_bytes=resources["cgroup_after"]["memory_peak_bytes"],
                    output=str(output),
                    cgroup=resources["cgroup_after"]["path"],
                    affinity=resources["affinity"],
                    effective_engine=counters.get("effective_engine", spec["engine"]),
                )
                state["rows"].append(row)
                write_json_atomic(evidence / "state.json", state)
        state.update(phase="untimed_strict_audits", last_progress_unix=time.time())
        write_json_atomic(evidence / "state.json", state)
        for row in state["rows"]:
            row["output_sha256"] = sha256sum(Path(row["output"]))
            row["output_index_sha256"] = sha256sum(Path(row["output"] + ".csi"))
        write_json_atomic(evidence / "state.json", state)
        # No audit reads between timed cells. Compare each cached output with
        # its own repeat's direct and both later directs with the first direct.
        os.environ["PATH"] = str(binary.parent) + ":" + os.environ["PATH"]
        pairs = {
            cell: (output, outputs[f"r{cell[1:3]}-direct"])
            for cell, output in outputs.items()
            if not cell.endswith("-direct")
        }
        pairs.update(
            {
                cell: (output, outputs["r01-direct"])
                for cell, output in outputs.items()
                if cell.endswith("-direct") and cell != "r01-direct"
            }
        )
        with ProcessPoolExecutor(max_workers=2) as pool:
            futures = {
                cell: pool.submit(audit_pair, pair) for cell, pair in pairs.items()
            }
            for cell, future in futures.items():
                comparison = future.result()
                write_json_atomic(evidence / cell / "against_direct.json", comparison)
                state["audits"][cell] = comparison["semantic_pass"]
                write_json_atomic(evidence / "state.json", state)
        assert sha256sum(source) == state["input_sha256"]
        passed = all(state["audits"].values())
        state.update(
            state="complete" if passed else "audit_failed", finished_unix=time.time()
        )
        if not passed:
            raise RuntimeError(
                "Complete-record/header equivalence failed; timings retained, not accepted"
            )
    except Exception as error:
        state.update(
            state="audit_failed" if state["state"] == "audit_failed" else "failed",
            error=str(error),
            traceback=traceback.format_exc(),
            finished_unix=time.time(),
        )
        raise
    finally:
        write_json_atomic(evidence / "state.json", state)


def collect(root, terminal=False):
    """Produce a days-later review even after failed, killed or incomplete jobs."""
    data = definition(root)
    states, rows, summary = [], [], []
    for case in data["cases"]:
        file = root / "cases" / case["id"] / "state.json"
        state = (
            json.loads(file.read_text())
            if file.exists()
            else dict(state="not_started", rows=[], audits={})
        )
        if terminal and state["state"] not in {"complete", "failed", "audit_failed"}:
            state = dict(
                state, state_before_terminal_review=state["state"], state="interrupted"
            )
        states.append(dict(case=case, state=state))
        direct = {
            r["repetition"]: r["seconds"]
            for r in state["rows"]
            if r["mode"] == "direct"
        }
        for row in state["rows"]:
            accepted = state["state"] == "complete"
            entry = dict(
                case=case["id"],
                node=case["node"],
                accepted=accepted,
                strict_pass=state["audits"].get(row["cell"]),
                **row,
            )
            entry["speedup_vs_direct"] = (
                direct.get(row["repetition"], 0) / row["seconds"] or None
            )
            rows.append(entry)
        for mode in case.get("modes", list(MODES)):
            selected = [
                r for r in rows if r["case"] == case["id"] and r["mode"] == mode
            ]
            if selected:
                gains = [
                    r["speedup_vs_direct"]
                    for r in selected
                    if r["speedup_vs_direct"] is not None
                ]
                summary.append(
                    dict(
                        case=case["id"],
                        node=case["node"],
                        mode=mode,
                        accepted=state["state"] == "complete",
                        repetitions=len(selected),
                        median_seconds=statistics.median(
                            r["seconds"] for r in selected
                        ),
                        min_seconds=min(r["seconds"] for r in selected),
                        max_seconds=max(r["seconds"] for r in selected),
                        median_speedup=statistics.median(gains) if gains else None,
                    )
                )
    for name, items in (("conditions.tsv", rows), ("summary.tsv", summary)):
        if items:
            with (root / name).open("w") as stream:
                writer = csv.DictWriter(
                    stream,
                    fieldnames=list(dict.fromkeys(k for item in items for k in item)),
                    delimiter="\t",
                    lineterminator="\n",
                    quoting=csv.QUOTE_NONNUMERIC,
                )
                writer.writeheader()
                writer.writerows(items)
    write_json_atomic(
        root / "REVIEW.json",
        dict(
            reviewed_unix=time.time(),
            terminal=terminal,
            cases=states,
            accepted_cases=sum(c["state"]["state"] == "complete" for c in states),
            accepted_timing_rows=sum(r["accepted"] for r in rows),
            no_further_campaign_submitted=True,
        ),
    )
    report = [
        "# Same-worker engine confirmation",
        "",
        "Fresh post-upgrade controls; fixed eight CPUs; three counterbalanced repeats.",
        "Private node-local cache; unchanged branch-default preprocessing; all bcftools/HTSlib at 1.25.",
        "Raw timing rows are not validated speedups unless the complete block passes strict audits.",
        "",
        "| Case | State | Direct s | Python s | Rust s | Python/direct speedup | Rust/direct speedup | Python/Rust cached ratio |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for item in states:
        case, state = item["case"], item["state"]
        selected = {r["mode"]: r for r in summary if r["case"] == case["id"]}
        cells = {r["cell"]: r for r in rows if r["case"] == case["id"]}
        gains = [
            cells[f"r{i:02d}-python"]["seconds"] / cells[f"r{i:02d}-rust"]["seconds"]
            for i in range(1, 4)
            if f"r{i:02d}-python" in cells and f"r{i:02d}-rust" in cells
        ]

        def fmt(x):
            return f"{x:.3f}" if x is not None else "NA"

        values = [fmt(selected.get(mode, {}).get("median_seconds")) for mode in MODES]
        values += [
            fmt(selected.get(mode, {}).get("median_speedup"))
            for mode in ("python", "rust")
        ]
        values += [fmt(statistics.median(gains) if gains else None)]
        report.append(
            "| "
            + case["id"]
            + " | "
            + state["state"]
            + " | "
            + " | ".join(values)
            + " |"
        )
    report += [
        "",
        "Technical repeats describe engineering variability, not independent biological replication.",
        (
            "Native and adaptive candidate rows are distinct; no silent backend fallback is accepted."
            if data.get("profile") == "candidates"
            else "Rust means adaptive routing, not resident/native WES; no silent backend fallback is accepted."
        ),
        "Cache transport, asset hashes and assay extraction are untimed; no OS cache flushing.",
        "No product defaults, manuscript figures or full-cohort campaigns are changed automatically.",
        "Local scratch/output paths and digests remain in case checkpoints for later audits.",
    ]
    if data.get("profile") == "candidates":
        report = candidate_report(data, states, rows, summary)
    (root / "RESULTS.md").write_text("\n".join(report) + "\n")
    if terminal:
        archive = root / "terminal-evidence.tar.gz"
        if not archive.exists():
            with tarfile.open(archive, "x:gz") as output:
                for name in (
                    "manifest.json",
                    "READY",
                    "submission_intent.json",
                    "submissions.json",
                    "review_submission.json",
                    "PREFLIGHT.json",
                    "cases",
                    "logs",
                    "conditions.tsv",
                    "summary.tsv",
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
    print(json.dumps({c["case"]["id"]: c["state"]["state"] for c in states}))


def candidate_report(data, states, rows, summary):
    """Display every candidate, not just a post-hoc fastest observation."""
    report = [
        "# Prepared-input best-candidate confirmation",
        "",
        "All conditions: fixed eight CPUs, external bcftools/HTSlib 1.25; Rust-linked HTSlib 1.25.",
        "Identical eligibility-verified input BCF; no normalization/conversion in timed arms; one shared untimed eligibility scan.",
        "Three rotated repeats; local cache and scratch; annotator databases retain their frozen locations.",
        "Only complete, strictly equivalent blocks are accepted. Other timings are diagnostic only.",
        "",
        "| Case | State | Configuration | Median s | Min s | Max s | Median speedup vs fresh direct | Python reference / candidate time |",
        "|---|---|---|---:|---:|---:|---:|---:|",
    ]

    def fmt(x):
        return f"{x:.3f}" if x is not None else "NA"

    for item in states:
        case, state = item["case"], item["state"]
        reference = "python-windows2" if "wes" in case["id"] else "python-auto"
        cells = {
            (r["repetition"], r["mode"]): r for r in rows if r["case"] == case["id"]
        }
        for entry in [r for r in summary if r["case"] == case["id"]]:
            paired = [
                cells[(i, reference)]["seconds"] / cells[(i, entry["mode"])]["seconds"]
                for i in range(1, 4)
                if (i, reference) in cells and (i, entry["mode"]) in cells
            ]
            values = [
                fmt(entry[k])
                for k in (
                    "median_seconds",
                    "min_seconds",
                    "max_seconds",
                    "median_speedup",
                )
            ]
            values.append(fmt(statistics.median(paired) if paired else None))
            report.append(
                "| "
                + " | ".join([case["id"], state["state"], entry["mode"], *values])
                + " |"
            )
    report += [
        "",
        "Candidate choices were prespecified from earlier pilots; this is not an exhaustive search or proof of a universal fastest implementation.",
        "Report all configurations; choosing and quoting a per-sample minimum would create selection bias.",
        "Prepared-input eligibility time is saved in cases/*/prepared_preflight.json; include appropriate validation/preparation costs before making end-to-end raw-VCF product claims.",
        "WES resident/native memory budget is 2048 MiB; requested workers/decoder settings and absence of fallback are verified.",
        "No native WGS is attempted; previous evidence shows a fallback, not a valid resident WGS measurement.",
        "Three repeats with more than three modes do not perfectly balance every order position. Read the min/max and every paired observation, not only medians.",
        "Prior default-routing results and historical direct baselines are not pooled with this pilot.",
        "Product defaults, manuscript figures and full-cohort scheduling remain unchanged.",
    ]
    return report


def main():
    """Expose explicit freeze/submit/run/review actions, never an automatic campaign."""
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    for action in ("freeze", "submit", "run", "collect"):
        p = sub.add_parser(action)
        p.add_argument("root", type=Path)
        if action == "freeze":
            p.add_argument("--deployment", required=True)
            p.add_argument("--commit", required=True)
            p.add_argument(
                "--profile", choices=["default", "candidates"], default="default"
            )
        if action == "run":
            p.add_argument("--case", required=True)
        if action == "collect":
            p.add_argument("--terminal", action="store_true")
    args = parser.parse_args()
    if args.action == "freeze":
        freeze(args)
    elif args.action == "submit":
        submit(args.root)
    elif args.action == "run":
        run(args.root, args.case)
    else:
        collect(args.root, args.terminal)


if __name__ == "__main__":
    main()
