#!/usr/bin/env python3
"""Prepare, submit and review an immutable profiling campaign without a live session."""

import argparse
import fcntl
import importlib.metadata
import json
import os
import shlex
import shutil
import socket
import subprocess
import sys
import tempfile
import time
import traceback
from pathlib import Path

import pysam
import yaml

from benchmarks.engine_profile.matrix import SAMPLES, design
from benchmarks.run_pilot import sha256sum, strict_semantic_compare, write_json_atomic

PYTHON = "/results/deployment/vcfcache-native-pilot-20261023c/.venv/bin/python"
WES_BED = Path(
    "/results/richness/regions/twist_human_core_exome_hg38.chr1-22-X.pad125.merged.bed"
)
ORIGINS = {
    "fastvep": Path(
        "/results/richness/caches/fastvep_clinical_grch38_af001/cache/fastvep_clinical_grch38_af001"
    ),
    "veprs-core": Path(
        "/results/native-bcf/veprs-caches-20261023b/veprs_core/cache/veprs_core"
    ),
    "vep": Path(
        "/mnt/data/vcfcache_benchmarks/bundled_zenodo_caches/gnomad_v4.1_GRCh38_joint_af001/cache/vep115.2_everything"
    ),
}


def read(path):
    """Read a complete atomic JSON artifact."""
    return json.loads(Path(path).read_text())


def checked(command, **options):
    """Run an explicit command; never silently ignore a setup failure."""
    return subprocess.run([str(v) for v in command], check=True, **options)


def claim(path):
    """Claim an empty directory safely, including NFS mkdir retry semantics."""
    path.mkdir(parents=True, exist_ok=True)
    if any(path.iterdir()):
        raise FileExistsError(path)
    with (path / ".owner").open("x") as stream:
        stream.write(f"{os.getpid()}\n")


def controller(path):
    """Map the cluster's shared mount only when running on its controller."""
    path = Path(path)
    if not Path("/results").exists() and str(path).startswith("/results/"):
        return Path("/mnt/data/slurm-results") / path.relative_to("/results")
    return path


def fingerprint(path):
    """Freeze regular files or sorted directory trees, outside timed cells."""
    path = Path(path)
    if path.is_file():
        return dict(path=str(path), bytes=path.stat().st_size, sha256=sha256sum(path))
    if not path.is_dir():
        raise FileNotFoundError(path)
    # Full tree fingerprints include transcript databases, not just recipes.
    files = {
        str(p.relative_to(path)): dict(bytes=p.stat().st_size, sha256=sha256sum(p))
        for p in sorted(path.rglob("*"))
        if p.is_file()
    }
    return dict(path=str(path), files=files)


def prepare(args):
    """Stage immutable fixtures on a worker; fail before scheduling if unsafe."""
    root = args.campaign.resolve()
    claim(root)
    request = read(args.request)
    (root / "logs").mkdir()
    (root / "attempts").mkdir()
    (root / "references").mkdir()
    groups = design()
    manifest = dict(
        schema_version=1,
        root=str(root),
        created_unix=time.time(),
        groups=groups,
        builds=request["builds"],
        recipes={},
        inputs={},
        max_cells=96,
        node_memory_mib=30000,
        memory_limit_mib=28000,
        safe_buffer_mib=8192,
        source_request=request,
        semantics="strict complete records; no ignored fields",
        storage_temperature="uncontrolled; first-use-after-copy is not cold",
    )
    # Reserve space for all retained references and conservatively failed cells.
    free = shutil.disk_usage(root).free
    if free < 300 * 1024**3:
        raise RuntimeError(
            f"shared storage requires 300 GiB free, found {free/1024**3:.1f}"
        )
    manifest["shared_free_bytes_at_prepare"] = free
    for build in manifest["builds"].values():
        for relative, expected in build["files"].items():
            if sha256sum(Path(build["root"]) / relative) != expected:
                raise RuntimeError(f"deployment drift: {build['root']}/{relative}")
    assets = root / "assets"
    assets.mkdir()
    resource_paths = set()
    for name, origin in ORIGINS.items():
        destination = assets / name / "cache" / origin.name
        destination.mkdir(parents=True)
        (destination.parents[1] / "workflow").mkdir()
        (destination.parents[1] / "blueprint").mkdir()
        shutil.copy2(
            origin.parents[1] / "workflow/init.yaml",
            destination.parents[1] / "workflow/init.yaml",
        )
        for file in origin.iterdir():
            if file.name in {
                "vcfcache_annotated.bcf",
                "vcfcache_annotated.bcf.csi",
                "annotation.yaml",
                "params.snapshot.yaml",
                "build_provenance.json",
            }:
                shutil.copy2(file, destination / file.name)
        recipe = yaml.safe_load((destination / "annotation.yaml").read_text())
        # The old vep-rs pilot included a params-only key in its recipe. Preserve
        # original bytes and document the metadata-only repair in the new bundle.
        if "tool_version_command" in recipe:
            shutil.copy2(
                destination / "annotation.yaml",
                destination / "annotation.original.yaml",
            )
            del recipe["tool_version_command"]
            (destination / "annotation.yaml").write_text(
                yaml.safe_dump(recipe, sort_keys=False)
            )
        params = yaml.safe_load((destination / "params.snapshot.yaml").read_text())
        recipe_resources = set()
        for key in (
            "fasta",
            "transcript_cache",
            "annotation_tool_cmd",
            "vep_cache",
            "vep_cache_dir",
            "vep_dir_cache",
        ):
            value = params.get(key)
            if value and Path(str(value)).exists():
                recipe_resources.add(str(value))
        command = recipe["annotation_cmd"].replace("\\\n", "")
        for key, value in params.items():
            command = command.replace("${params." + key + "}", str(value))
        tokens = shlex.split(command)
        for token in tokens:
            for piece in token.split(","):
                if piece.startswith("/") and Path(piece).exists():
                    if Path(piece) in {
                        Path("/"),
                        Path("/results"),
                        Path("/mnt/data"),
                        Path("/mnt/data/tmp"),
                    }:
                        raise RuntimeError(
                            f"refusing to fingerprint broad directory: {piece}"
                        )
                    recipe_resources.add(piece)
        for path in list(recipe_resources):
            for suffix in (".fai", ".gzi", ".tbi", ".csi"):
                if Path(path + suffix).is_file():
                    recipe_resources.add(path + suffix)
        if name == "vep":
            # Only GRCh38 is used; do not hash unrelated assemblies for this recipe.
            parent = str(params["vep_cache"])
            selected_cache = Path(parent) / "homo_sapiens/115_GRCh38"
            if not selected_cache.is_dir():
                raise FileNotFoundError(selected_cache)
            recipe_resources.discard(parent)
            recipe_resources.add(str(selected_cache))
        resource_paths.update(recipe_resources)
        manifest["recipes"][name] = dict(
            cache=str(destination),
            files={
                p.name: fingerprint(p) for p in destination.iterdir() if p.is_file()
            },
            params=params,
            source=str(origin),
            resources=sorted(recipe_resources),
            metadata_repair=(destination / "annotation.original.yaml").exists(),
        )
    # Reduce overlapping directory roots while still hashing every database.
    roots = [Path(p) for p in sorted(resource_paths)]
    roots = [p for p in roots if not any(q != p and q in p.parents for q in roots)]
    manifest["resources"] = [fingerprint(p) for p in roots]
    manifest["bed"] = fingerprint(WES_BED)
    for sample in SAMPLES:
        directory = assets / "inputs" / sample
        directory.mkdir(parents=True)
        source = Path(
            f"/mnt/data/vcfcache_benchmarks/external_wgs/samples/GRCh38/sgdp/{sample}.GRCh38.small_variants.vcf.gz"
        )
        if not source.is_file():
            raise FileNotFoundError(source)
        manifest["inputs"][sample] = dict(source=fingerprint(source))
        wgs = directory / "wgs.bcf"
        checked(["bcftools", "view", source, "-Ob", "-W", "--threads", "8", "-o", wgs])
        checked(
            [
                "bcftools",
                "view",
                "-R",
                WES_BED,
                wgs,
                "-Ob",
                "-W",
                "--threads",
                "8",
                "-o",
                directory / "wes.bcf",
            ]
        )
        for assay in ("wes", "wgs"):
            path = directory / f"{assay}.bcf"
            records, split = 0, True
            with pysam.VariantFile(str(path)) as stream:
                for record in stream:
                    records += 1
                    split = split and len(record.alleles) == 2
            manifest["inputs"][sample][assay] = dict(
                **fingerprint(path), records=records, already_split=split
            )
        checked(
            [
                "bcftools",
                "view",
                directory / "wes.bcf",
                "-Oz",
                "-W",
                "-o",
                directory / "wes.vcf.gz",
            ]
        )
        manifest["inputs"][sample]["wes_vcf"] = fingerprint(directory / "wes.vcf.gz")
    manifest["python"] = sys.version
    manifest["python_packages"] = sorted(
        (dist.metadata["Name"], dist.version)
        for dist in importlib.metadata.distributions()
    )
    manifest["bcftools"] = subprocess.check_output(["bcftools", "--version"], text=True)
    manifest["references_policy"] = (
        "Fresh references: new fixture snapshot has no previously verified direct/reference-output identity binding. References run once per sample/assay/recipe."
    )
    write_json_atomic(root / "manifest.json", manifest)
    (root / "READY").write_text(sha256sum(root / "manifest.json") + "\n")
    print(
        json.dumps(dict(root=str(root), cells=sum(len(g["cells"]) for g in groups))),
        flush=True,
    )


def verify(root):
    """Refuse a changed campaign manifest before submission or execution."""
    if (root / "READY").read_text().strip() != sha256sum(root / "manifest.json"):
        raise RuntimeError("campaign manifest was modified after preparation")
    return read(root / "manifest.json")


def submit(args):
    """Submit frozen groups and an afterany collector, recording each job immediately."""
    root = args.campaign.resolve()
    manifest = verify(root)
    if args.command == "resume" and not args.cells:
        raise ValueError(
            "resume requires explicit --cells; successful cells are never silently rerun"
        )
    selected = set(args.cells or [])
    all_ids = {cell["id"] for group in manifest["groups"] for cell in group["cells"]}
    if selected - all_ids:
        raise ValueError(f"unknown cells: {selected-all_ids}")
    with (root / ".submit.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        history = (
            read(root / "submissions.json")
            if (root / "submissions.json").exists()
            else []
        )
        if history and args.command == "submit":
            raise RuntimeError(
                "already submitted; inspect status and use explicit resume"
            )
        # No retry while a previous attempt could still be writing.
        if history:
            old_ids = [str(j["job_id"]) for batch in history for j in batch["jobs"]]
            old_ids += [
                str(batch["collector_job_id"])
                for batch in history
                if batch.get("collector_job_id")
            ]
            active = (
                subprocess.check_output(
                    ["squeue", "-h", "-j", ",".join(old_ids), "-o", "%i"], text=True
                ).strip()
                if old_ids
                else ""
            )
            if active:
                raise RuntimeError(f"previous jobs still active: {active}")
        attempt = len(history) + 1
        batch = dict(attempt=attempt, submitted_unix=time.time(), jobs=[])
        history.append(batch)
        write_json_atomic(root / "submissions.json", history)
        worker_root = manifest["root"]
        deployment = manifest.get("harness_root", manifest["builds"]["profile"]["root"])
        script = f"{deployment}/benchmarks/engine_profile/job.sh"
        try:
            for group in manifest["groups"]:
                cells = [
                    c["id"]
                    for c in group["cells"]
                    if not selected or c["id"] in selected
                ]
                if not cells:
                    continue
                for cell in cells:
                    previous = latest(root, cell)
                    if previous and previous.get("state") in {
                        "complete",
                        "complete_with_direct_difference",
                    }:
                        raise RuntimeError(f"refusing to rerun successful cell {cell}")
                command = [
                    "sbatch",
                    "--parsable",
                    f"--nodelist={group['node']}",
                    f"--output={worker_root}/logs/{group['id']}-a{attempt}-%j.out",
                    str(controller(script)),
                    worker_root,
                    group["id"],
                    str(attempt),
                    ",".join(cells),
                    deployment,
                ]
                job = subprocess.check_output(command, text=True).strip().split(";")[0]
                if not job.isdigit():
                    raise RuntimeError(f"unexpected sbatch response: {job}")
                batch["jobs"].append(
                    dict(group=group["id"], job_id=job, cells=cells, command=command)
                )
                write_json_atomic(root / "submissions.json", history)
        except Exception as error:
            batch["submission_error"] = str(error)
            write_json_atomic(root / "submissions.json", history)
            raise
        finally:
            if batch["jobs"]:
                dependency = "afterany:" + ":".join(j["job_id"] for j in batch["jobs"])
                command = [
                    "sbatch",
                    "--parsable",
                    "--partition=benchmark",
                    "--nodes=1",
                    "--ntasks=1",
                    "--cpus-per-task=1",
                    "--mem=4000M",
                    "--exclude=sl-w5",
                    "--time=01:00:00",
                    f"--dependency={dependency}",
                    f"--output={worker_root}/logs/collect-a{attempt}-%j.out",
                    "--wrap",
                    f"cd /tmp; PYTHONPATH={shlex.quote(deployment)} {PYTHON} -m benchmarks.engine_profile.run collect {shlex.quote(worker_root)}",
                ]
                batch["collector_job_id"] = (
                    subprocess.check_output(command, text=True).strip().split(";")[0]
                )
                write_json_atomic(root / "submissions.json", history)
        print(json.dumps(batch, indent=2))


def latest(root, cell):
    """Use the latest explicit attempt, never a silently overwritten result."""
    files = sorted((root / "attempts").glob(f"*/{cell}/state.json"))
    return read(files[-1]) if files else None


def persist_run(scratch, evidence):
    """Preserve logs and small metadata; work spools stay out of shared storage."""
    for path in scratch.rglob("*"):
        relative = path.relative_to(scratch)
        if path.is_file() and path.suffix in {
            ".yaml",
            ".json",
            ".log",
            ".txt",
            ".sh",
            ".stdout",
            ".stderr",
        }:
            dest = evidence / "cli" / relative
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, dest)


def resources_match(manifest, recipe):
    """Verify pinned resources before work; do not reuse stale cache provenance."""
    for resource in manifest["resources"]:
        needed = manifest["recipes"][recipe]["resources"]
        if not any(
            Path(resource["path"]) == Path(p)
            or Path(resource["path"]) in Path(p).parents
            for p in needed
        ):
            continue
        if fingerprint(resource["path"]) != resource:
            raise RuntimeError(f"annotation resource drift: {resource['path']}")


def cell_params(manifest, cell):
    """Keep engine CPU allocation separate from explicit VEP memory controls."""
    params = dict(manifest["recipes"][cell["recipe"]]["params"], threads=cell["cpus"])
    for name in ("vep_forks", "rayon_threads"):
        if name in params:
            params[name] = cell["cpus"]
    overrides = cell.get("annotation_overrides", {})
    if overrides and cell["recipe"] != "vep":
        raise ValueError("annotation overrides are restricted to VEP recovery")
    if set(overrides) - {"vep_forks", "vep_buffer"}:
        raise ValueError("unsupported annotation override")
    if any(type(value) is not int or value < 1 for value in overrides.values()):
        raise ValueError("annotation overrides must be positive integers")
    if overrides.get("vep_forks", 1) > cell["cpus"]:
        raise ValueError("VEP forks exceed CPU allocation")
    params.update(overrides)
    return params


def run_cell(manifest, cell, scratch, evidence, cache, refs):
    """Run one isolated Slurm step, audit outside timing, and durably finalize it."""
    state = dict(
        cell=cell,
        state="running",
        started_unix=time.time(),
        node=socket.gethostname(),
        job_id=os.environ.get("SLURM_JOB_ID"),
    )
    write_json_atomic(evidence / "state.json", state)
    try:
        if cell["role"] != "direct" and (cell["assay"], "direct") not in refs:
            state.update(state="skipped", reason="direct reference failed or missing")
            return
        if cell["role"] == "measure" and (cell["assay"], "cached") not in refs:
            state.update(state="skipped", reason="cached reference failed or missing")
            return
        if cell["memory_mb"] > manifest["safe_buffer_mib"]:
            state.update(state="skipped", reason="memory reservation policy")
            return
        source_info = manifest["inputs"][cell["sample"]][cell["assay"]]
        if cell["skip_split"] and not source_info["already_split"]:
            state.update(state="skipped", reason="input is not verified split")
            return
        if not cache.is_dir():
            state.update(
                state="skipped",
                reason="cache probe missing; resume first-use and repeated-access together",
            )
            return
        profile_root = evidence / "trace"
        profile_root.mkdir()
        source = scratch.parent / "inputs" / f"{cell['assay']}.{cell['form']}"
        output = scratch / "output.bcf"
        params = cell_params(manifest, cell)
        params["temp_dir"] = str(scratch / "tmp")
        Path(params["temp_dir"]).mkdir()
        params_file = scratch / "params.yaml"
        params_file.write_text(yaml.safe_dump(params))
        build = manifest["builds"][cell["build"]]
        entry = "from vcfcache.cli import main; main()"
        if cell["profile"]:
            entry = "from vcfcache.database.profiling import phase\nwith phase('cli.import_and_run'):\n from vcfcache.cli import main\n main()"
        command = [
            PYTHON,
            "-c",
            entry,
            "annotate",
            "-a",
            str(cache),
            "-i",
            str(source),
            "-o",
            str(output),
            "-y",
            str(params_file),
            "--stats-dir",
            str(scratch / "stats"),
            "--statistics",
            "light",
            "--engine",
            cell["engine"],
        ]
        if cell["role"] == "direct":
            command += ["--uncached"]
        if cell["engine"] in {"adaptive", "native"}:
            command += ["--engine-workers", str(cell["workers"])]
        if cell["engine"] == "native":
            command += ["--engine-memory-mb", str(cell["memory_mb"])]
        if cell["skip_split"]:
            command += ["--skip-split-multiallelic"]
        environment = dict(
            os.environ,
            PYTHONPATH=build["root"] + ":" + manifest["builds"]["profile"]["root"],
            VCFCACHE_ROOT=build["root"],
            VCFCACHE_ENGINE=build["root"] + "/native/target/release/vcfcache-engine",
        )
        environment.pop("VCFCACHE_PROFILE_DIR", None)
        experimental = cell.get("environment", {})
        allowed_experiments = {
            "VCFCACHE_EXPERIMENT_INDEXED_MISSES",
            "VCFCACHE_EXPERIMENT_WINDOW_BP",
        }
        if set(experimental) - allowed_experiments:
            raise ValueError("unrecognized experiment environment")
        environment.update(experimental)
        if cell["profile"]:
            environment["VCFCACHE_PROFILE_DIR"] = str(profile_root)
        if cell["perf"]:
            probe = (
                subprocess.run(
                    ["perf", "stat", "-e", "task-clock", "--", "true"],
                    capture_output=True,
                    text=True,
                )
                if shutil.which("perf")
                else None
            )
            if probe is None or probe.returncode:
                state.update(
                    state="skipped",
                    reason="perf unavailable or denied",
                    diagnostic=probe.stderr if probe else "not installed",
                )
                return
            command = [
                "perf",
                "record",
                "-F",
                "49",
                "--call-graph",
                "dwarf,4096",
                "-o",
                str(evidence / "perf.data"),
                "--",
                *command,
            ]
        # The sampler imports only benchmark code and runs in its own step cgroup.
        wrapper = [
            PYTHON,
            str(
                Path(manifest["builds"]["profile"]["root"])
                / "benchmarks/engine_profile/measure.py"
            ),
            "--directory",
            str(evidence),
            "--",
            *command,
        ]
        # Baseline deployment also receives benchmark-only harness sources; its
        # vcfcache package and helper remain the untouched frozen baseline.
        launch = [
            "srun",
            "--exclusive",
            "--exact",
            "--cpu-bind=cores",
            "-n1",
            f"-c{cell['cpus']}",
            *wrapper,
        ]
        write_json_atomic(
            evidence / "command.json",
            dict(
                command=launch,
                build=build,
                input=source_info,
                cache_files=manifest["recipes"][cell["recipe"]]["files"],
                profiling=cell["profile"],
                effective_parameters=params,
                experiment_environment=experimental,
                temperature=manifest["storage_temperature"],
            ),
        )
        with (evidence / "run.log").open("w") as log:
            completed = subprocess.run(
                launch, env=environment, cwd=scratch, stdout=log, stderr=log
            )
        persist_run(scratch, evidence)
        if completed.returncode:
            raise RuntimeError(
                f"Slurm step/CLI exit {completed.returncode}; see run.log"
            )
        if not output.is_file() or not Path(str(output) + ".csi").is_file():
            raise RuntimeError("missing output BCF/CSI")
        state["output_sha256"] = sha256sum(output)
        metrics = read(evidence / "resources.json")
        if not metrics["cgroup_after"]:
            raise RuntimeError("independent cgroup metrics unavailable")
        state["metrics"] = metrics
        stats_paths = list((scratch / "stats").glob("*_vcstats/compare_stats.yaml"))
        if len(stats_paths) != 1:
            raise RuntimeError("missing or ambiguous CLI statistics")
        stats = yaml.safe_load(stats_paths[0].read_text())
        write_json_atomic(evidence / "statistics.json", stats)
        state["effective_engine"] = stats.get(
            "effective_engine", stats.get("engine", "legacy")
        )
        if cell["profile"] and (
            not list(profile_root.glob("*.jsonl"))
            or "PROFILE_WRITE_FAILED" in (evidence / "run.log").read_text()
        ):
            raise RuntimeError("profiling trace incomplete or failed")
        if cell["perf"]:
            with (evidence / "perf-report.txt").open("w") as out:
                subprocess.run(
                    ["perf", "report", "--stdio", "-i", str(evidence / "perf.data")],
                    stdout=out,
                    stderr=subprocess.STDOUT,
                    check=True,
                )
        key = cell["assay"]
        if cell["role"] == "direct":
            comparison = strict_semantic_compare(output, output)
            state["audit_role"] = "reference self-check, not independent validation"
            refs[(key, "direct")] = output
        else:
            reference = refs.get((key, "direct"))
            if reference is None:
                raise RuntimeError(
                    "direct reference missing; never substitute another cached arm"
                )
            comparison = strict_semantic_compare(output, reference)
        write_json_atomic(evidence / "against_direct.json", comparison)
        state["direct_pass"] = comparison["semantic_pass"]
        if cell["role"] == "cached":
            refs[(key, "cached")] = output
        elif cell["role"] == "measure":
            baseline = refs.get((key, "cached"))
            if baseline is None:
                raise RuntimeError("uninstrumented cached reference missing")
            equivalence = strict_semantic_compare(output, baseline)
            write_json_atomic(evidence / "against_uninstrumented.json", equivalence)
            state["instrumentation_equivalent"] = equivalence["semantic_pass"]
            if not equivalence["semantic_pass"]:
                raise RuntimeError(
                    "output differs from frozen uninstrumented cached reference"
                )
        state["state"] = (
            "complete"
            if comparison["semantic_pass"]
            else "complete_with_direct_difference"
        )
    except Exception as error:
        state.update(state="failed", error=str(error), traceback=traceback.format_exc())
    finally:
        persist_run(scratch, evidence)
        output = scratch / "output.bcf"
        keep = cell["role"] != "measure" or state["state"] not in {
            "complete",
            "skipped",
        }
        if keep and output.exists():
            dest = evidence / "output.bcf"
            shutil.copy2(output, dest)
            if sha256sum(dest) != sha256sum(output):
                raise RuntimeError("durable output copy checksum mismatch")
            index = Path(str(output) + ".csi")
            if index.exists():
                shutil.copy2(index, Path(str(dest) + ".csi"))
            state["retained_output"] = str(dest)
        state["finished_unix"] = time.time()
        write_json_atomic(evidence / "state.json", state)
    return state


def group_run(args):
    """Run one pinned group without relying on the submitting shell."""
    root = args.campaign.resolve()
    manifest = verify(root)
    group = next(g for g in manifest["groups"] if g["id"] == args.group)
    if socket.gethostname().split(".")[0] != group["node"]:
        raise RuntimeError("group placed on a different node than frozen design")
    for relative, expected in manifest.get("harness_files", {}).items():
        if sha256sum(Path(manifest["harness_root"]) / relative) != expected:
            raise RuntimeError(f"recovery harness drift: {relative}")
    for build in manifest["builds"].values():
        for relative, expected in build["files"].items():
            if sha256sum(Path(build["root"]) / relative) != expected:
                raise RuntimeError(f"deployment drift: {build['root']}/{relative}")
    if (
        subprocess.check_output(["bcftools", "--version"], text=True)
        != manifest["bcftools"]
    ):
        raise RuntimeError("bcftools version drift")
    resources_match(manifest, group["recipe"])
    selected = set(args.cells.split(","))
    if shutil.disk_usage("/mnt/data").free < 60 * 1024**3:
        raise RuntimeError("worker requires 60 GiB scratch headroom")
    work = Path(
        tempfile.mkdtemp(prefix=f"engine-profile-{group['id']}-", dir="/mnt/data/tmp")
    )
    (work / ".owner").write_text(str(root))
    local_cache = (
        work
        / "bundle"
        / "cache"
        / Path(manifest["recipes"][group["recipe"]]["cache"]).name
    )
    shutil.copytree(
        Path(manifest["recipes"][group["recipe"]]["cache"]).parents[1], work / "bundle"
    )
    for name, info in manifest["recipes"][group["recipe"]]["files"].items():
        if sha256sum(local_cache / name) != info["sha256"]:
            raise RuntimeError("staged cache mismatch")
    input_dir = work / "inputs"
    input_dir.mkdir()
    for name in ("wes", "wgs", "wes_vcf"):
        info = manifest["inputs"][group["sample"]][name]
        suffix = "wes.vcf.gz" if name == "wes_vcf" else f"{name}.bcf"
        dest = input_dir / suffix
        shutil.copy2(info["path"], dest)
        shutil.copy2(str(info["path"]) + ".csi", str(dest) + ".csi")
        if sha256sum(dest) != info["sha256"]:
            raise RuntimeError("staged input mismatch")
    attempt = root / "attempts" / f"a{args.attempt:03d}"
    attempt.mkdir(parents=True, exist_ok=True)
    refs = {}
    reference_root = Path(manifest.get("reference_campaign", str(root)))
    reference_group = group
    if reference_root != root:
        reference_manifest = verify(reference_root)
        if (
            sha256sum(reference_root / "manifest.json")
            != manifest["reference_manifest_sha256"]
        ):
            raise RuntimeError("reference campaign drift")
        reference_group = next(
            g for g in reference_manifest["groups"] if g["id"] == group["id"]
        )
        if (reference_group["sample"], reference_group["recipe"]) != (
            group["sample"],
            group["recipe"],
        ):
            raise RuntimeError("reference group does not match pilot")
    for assay in sorted({cell["assay"] for cell in group["cells"]}):
        for role in ("direct", "cached"):
            refcell = next(
                c
                for c in reference_group["cells"]
                if c["assay"] == assay and c["role"] == role
            )
            previous = latest(reference_root, refcell["id"])
            if previous and previous.get("state") in {
                "complete",
                "complete_with_direct_difference",
            }:
                path = Path(previous["retained_output"])
                if sha256sum(path) != previous["output_sha256"]:
                    raise RuntimeError("reference output drift")
                refs[(assay, role)] = path
    for cell in group["cells"]:
        if cell["id"] not in selected:
            continue
        scratch = work / cell["id"]
        scratch.mkdir()
        evidence = attempt / cell["id"]
        claim(evidence)
        cache = local_cache
        if cell["purpose"] == "first-use-after-copy":
            storage = work / "storage-probe"
            shutil.copytree(work / "bundle", storage)
        if cell["purpose"] in {"first-use-after-copy", "repeated-access"}:
            cache = work / "storage-probe/cache" / local_cache.name
        run_cell(manifest, cell, scratch, evidence, cache, refs)
        # Only private, owned per-cell scratch is removed; retained references
        # are switched to their verified durable copies first.
        for key, value in list(refs.items()):
            if value == scratch / "output.bcf":
                refs[key] = evidence / "output.bcf"
        state = read(evidence / "state.json")
        if state["state"] in {"complete", "complete_with_direct_difference", "skipped"}:
            shutil.rmtree(scratch)
    write_json_atomic(
        root / f"group-{group['id']}-a{args.attempt}.json",
        dict(finished_unix=time.time(), scratch=str(work)),
    )


def main():
    """Expose submission and offline review without changing successful attempts."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "command", choices=("prepare", "submit", "status", "collect", "resume", "group")
    )
    parser.add_argument("campaign", type=Path)
    parser.add_argument("--request", type=Path)
    parser.add_argument("--cells", nargs="*")
    parser.add_argument("--group")
    parser.add_argument("--attempt", type=int, default=1)
    args = parser.parse_args()
    if args.command == "prepare":
        if not args.request:
            parser.error("prepare needs --request")
        prepare(args)
    elif args.command in {"submit", "resume"}:
        submit(args)
    elif args.command == "group":
        if not args.group or not args.cells or len(args.cells) != 1:
            parser.error("group needs --group and one comma-separated --cells argument")
        args.cells = args.cells[0]
        try:
            group_run(args)
        except BaseException:
            write_json_atomic(
                args.campaign / f"group-{args.group}-a{args.attempt}-failure.json",
                dict(
                    finished_unix=time.time(),
                    error=traceback.format_exc(),
                    job_id=os.environ.get("SLURM_JOB_ID"),
                ),
            )
            raise
    else:
        from benchmarks.engine_profile.review import collect

        collect(args.campaign.resolve(), write=args.command == "collect")


if __name__ == "__main__":
    main()
