"""Confirmation of the generic-input gather engine against plain annotation.

One Slurm job per node, one sample per node. Every cell is a fresh eight-CPU
step. Arms:

    plain   the recipe alone (reference deployment, no vcfcache checks)
    direct  this branch with --uncached (recipe plus output validation)
    pilot   trusted-input pilot engine (reference deployment), best pilot mode
    cached  this branch, no switches

An "unusual" WGS input (multiallelic records restored, duplicate records
injected on chr21) checks that irregular input still equals plain annotation.
"""

import argparse
import hashlib
import json
import os
import shlex
import shutil
import socket
import statistics
import subprocess
import tempfile
import time
import traceback
from pathlib import Path

import yaml

from benchmarks.run_pilot import strict_semantic_compare, write_json_atomic

PYTHON = "/results/deployment/vcfcache-native-pilot-20261023c/.venv/bin/python"
REFERENCE = "/results/deployment/veprs05-gather-20261110"
PARENT = "/results/python-adaptive-manuscript-20261101-ee2162473/campaign/manifest.json"
BCFTOOLS = "/results/input-preserving-20261105-tools/bin/bcftools"
CACHES = "/results/veprs05/caches"
VEP = "/results/veprs05/dist/vep-0.5.0-x86_64-unknown-linux-gnu/vep-0.5.0-x86_64-unknown-linux-gnu/vep"
PARAMETERS = dict(
    annotation_tool_cmd=VEP,
    genome_build="GRCh38",
    optional_checks={},
    tool_version_command=VEP + " --version",
    transcript_cache="/results/veprs/json_cache/GRCh38_native_20261023",
)
SAMPLES = {"sl-w1": "KPGP-00128", "sl-w2": "LP6005443-DNA_B07"}
# (recipe, assay, arms). Clinical WGS omits `direct`: validation cost does not
# depend on the recipe and a third 13-minute control adds nothing.
GROUPS = [
    ("core", "wgs", ["plain", "direct", "pilot", "cached"]),
    ("core", "unusual", ["plain", "cached"]),
    ("core", "wes", ["plain", "direct", "pilot", "cached"]),
    ("core", "panel", ["plain", "direct", "pilot", "cached"]),
    ("clinical", "wes", ["plain", "pilot", "cached"]),
    ("clinical", "panel", ["plain", "pilot", "cached"]),
    ("clinical", "wgs", ["plain", "pilot", "cached"]),
]
REPEATS = 2


def head(path):
    """Controller spelling of a worker-visible shared path."""
    return str(path).replace("/results/", "/mnt/data/slurm-results/", 1)


def wrapper(directory):
    """Thread every bcftools stage, including those inside frozen recipes."""
    directory.mkdir()
    script = directory / "bcftools"
    script.write_text(f"""#!/bin/bash
set -euo pipefail
actual={shlex.quote(BCFTOOLS)}
subcommand=${{1-}}
case "$subcommand" in
  annotate|view|filter|concat|index|sort)
    shift
    for arg in "$@"; do
      case "$arg" in --threads|--threads=*|-@*) exec "$actual" "$subcommand" "$@";; esac
    done
    exec "$actual" "$subcommand" "$@" --threads 8;;
esac
exec "$actual" "$@"
""")
    script.chmod(0o755)
    return script


def extract(task, folder):
    """Assay extraction and BCF encoding only; no record is split or removed."""
    output = folder / (task["assay"] + ".bcf")
    command = [BCFTOOLS, "view", task["input"]["path"], "-Ob1", "-W", "--threads", "8"]
    command += ["-o", str(output)]
    if task.get("preprocess", {}).get("bed"):
        command += ["-R", task["preprocess"]["bed"]["path"]]
    subprocess.run(command, check=True)
    return output


def unusual(source, folder):
    """Restore multiallelic records and repeat every 2000th chr21 record."""
    output = folder / "unusual.bcf"
    subprocess.run(
        f"{BCFTOOLS} norm -m+any --threads 8 {source} -Ov | "
        "awk '{print} !/^#/ && $1==\"chr21\" && ++n%2000==0 {print}' | "
        f"{BCFTOOLS} view -Ob1 -W --threads 8 -o {output}",
        shell=True,
        check=True,
        executable="/bin/bash",
    )
    return output


def record_digest(path):
    """Digest of all record lines; headers legitimately differ between runs."""
    process = subprocess.Popen(
        [BCFTOOLS, "view", "-H", "--threads", "4", str(path)], stdout=subprocess.PIPE
    )
    assert process.stdout is not None
    with process.stdout:
        digest = hashlib.file_digest(process.stdout, "md5").hexdigest()
    if process.wait():
        raise RuntimeError(f"Cannot read {path}")
    return digest


def cell(deployment, scratch, evidence, name, arm, assay, cache, source, binary):
    """Run one timed annotation as its own eight-CPU step."""
    directory = evidence / name
    directory.mkdir(parents=True)
    output = scratch / (name + ".bcf")
    stats = scratch / (name + "-stats")
    params = scratch / (name + ".yaml")
    (scratch / "tmp").mkdir(exist_ok=True)
    params.write_text(
        yaml.safe_dump(
            dict(
                PARAMETERS,
                bcftools_cmd=str(binary),
                threads=8,
                temp_dir=str(scratch / "tmp"),
            )
        )
    )
    code = REFERENCE if arm in {"plain", "pilot"} else deployment
    environment = dict(
        os.environ,
        PYTHONPATH=code,
        VCFCACHE_ROOT=code,
        PATH=str(binary.parent)
        + ":"
        + str(Path(BCFTOOLS).parent)
        + ":"
        + os.environ["PATH"],
        LC_ALL="C",
    )
    command = [PYTHON, "-c", "from vcfcache.cli import main; main()", "annotate"]
    command += ["-a", str(cache), "-i", str(source), "-o", str(output)]
    command += ["-y", str(params), "--statistics", "light", "--stats-dir", str(stats)]
    command += ["--preserve-unannotated", "-v"]
    if arm == "direct":
        command += ["--uncached"]
    elif arm in {"plain", "pilot"}:
        # The trusted-input switch also keeps the reference deployment's
        # uncached path from normalizing: `plain` is the recipe and nothing else.
        command += ["--engine", "legacy", "--uncached"] if arm == "plain" else []
        command += ["--engine", "adaptive"] if arm == "pilot" else []
        sparse = assay != "wgs"
        environment.update(
            VCFCACHE_PILOT_REGION_MODE="off",
            VCFCACHE_PILOT_PREPARED="1",
            VCFCACHE_PILOT_QUEUE="1",
            VCFCACHE_PILOT_LOOKUP="windows" if sparse else "contigs",
            VCFCACHE_PILOT_MERGE_WINDOWS="1" if sparse else "0",
            VCFCACHE_PILOT_THREAD_BUDGET="0",
        )
    launch = [
        "srun",
        "--exclusive",
        "--exact",
        "--cpu-bind=cores",
        "-n1",
        "-c8",
        PYTHON,
    ]
    launch += [str(Path(deployment) / "benchmarks/gather_engine/measure.py")]
    launch += ["--directory", str(directory), "--", *command]
    write_json_atomic(directory / "command.json", dict(command=launch))
    with (directory / "run.log").open("w") as log:
        subprocess.run(
            launch, env=environment, stdout=log, stderr=log, check=True, timeout=10800
        )
    resources = json.loads((directory / "resources.json").read_text())
    assert resources["exit_code"] == 0 and len(resources["affinity"]) == 8
    found = list(stats.rglob("compare_stats.yaml"))
    assert len(found) == 1
    counters = yaml.safe_load(found[0].read_text())
    shutil.rmtree(stats, ignore_errors=True)
    return output, dict(
        seconds=resources["wall_seconds"],
        cpu_seconds=resources["cpu_seconds"],
        peak_gib=resources["cgroup_after"]["memory_peak_bytes"] / 2**30,
        output_bytes=output.stat().st_size,
        output_records=counters["variant_counts"]["total_output"],
        missing=counters["variant_counts"].get("missing_variants"),
        shards=counters.get("shards"),
        duplicate_sites=counters.get("duplicate_sites"),
        cache_access=counters.get("cache_access"),
    )


def node(root, name, deployment):
    """All cells of one node."""
    parent = json.loads(Path(PARENT).read_text())
    sample = SAMPLES[name]
    tasks = {
        t["assay"]: t
        for t in parent["tasks"]
        if (t["tool"], t["level"], t["assembly"], t["sample"])
        == ("veprs", "core", "GRCh38", sample)
    }
    evidence = root / "nodes" / name
    evidence.mkdir(parents=True, exist_ok=False)
    state = dict(state="running", node=socket.gethostname(), sample=sample)
    state.update(rows=[], audits=[], failures=[], started=time.time())
    path = evidence / "state.json"
    scratch = Path(tempfile.mkdtemp(prefix="gather-", dir="/mnt/data/tmp"))
    binary = wrapper(scratch / "bin")
    try:
        assert state["node"] == name
        sources = {}
        for assay in ("wgs", "wes", "panel"):
            folder = scratch / assay
            folder.mkdir()
            sources[assay] = extract(tasks[assay], folder)
        sources["unusual"] = unusual(sources["wgs"], scratch / "wgs")
        state["records"] = {
            assay: int(subprocess.check_output([BCFTOOLS, "index", "-n", str(source)]))
            for assay, source in sources.items()
        }
        write_json_atomic(path, state)
        for recipe, assay, arms in GROUPS:
            group = f"{recipe}-{assay}"
            cache = Path(CACHES) / f"veprs05_{recipe}" / "cache" / f"veprs05_{recipe}"
            try:
                outputs = {}
                for repeat in range(1, REPEATS + 1):
                    for arm in arms if repeat % 2 else arms[::-1]:
                        output, row = cell(
                            deployment,
                            scratch,
                            evidence,
                            f"{group}-{arm}-r{repeat}",
                            arm,
                            assay,
                            cache,
                            sources[assay],
                            binary,
                        )
                        assert row["output_records"] == state["records"][assay]
                        outputs[(arm, repeat)] = output
                        state["rows"].append(
                            dict(
                                recipe=recipe,
                                assay=assay,
                                arm=arm,
                                repeat=repeat,
                                **row,
                            )
                        )
                        write_json_atomic(path, state)
                comparison = strict_semantic_compare(
                    outputs[("cached", 1)], outputs[("plain", 1)]
                )
                write_json_atomic(evidence / f"{group}-against_plain.json", comparison)
                state["audits"].append(
                    dict(
                        recipe=recipe,
                        assay=assay,
                        semantic_pass=comparison.get("semantic_pass"),
                        mismatches=comparison.get("record_mismatches"),
                        records_equal_plain=record_digest(outputs[("cached", 1)])
                        == record_digest(outputs[("plain", 1)]),
                        repeat_identical=record_digest(outputs[("cached", 1)])
                        == record_digest(outputs[("cached", 2)]),
                        index=Path(str(outputs[("cached", 1)]) + ".csi").exists(),
                    )
                )
            except Exception as error:
                state["failures"].append(
                    dict(
                        group=group,
                        error=str(error)[:2000],
                        traceback=traceback.format_exc()[-4000:],
                    )
                )
            finally:
                write_json_atomic(path, state)
                for stale in scratch.glob(f"{group}-*"):
                    if stale.is_dir():
                        shutil.rmtree(stale, ignore_errors=True)
                    else:
                        stale.unlink(missing_ok=True)
        state.update(state="complete", finished=time.time())
        shutil.rmtree(scratch, ignore_errors=True)
    except Exception as error:
        state.update(state="failed", error=str(error), traceback=traceback.format_exc())
        raise
    finally:
        write_json_atomic(path, state)


def ratio(mean, arm):
    """Speedup of plain annotation over one arm, when both were timed."""
    return f"{mean['plain'] / mean[arm]:.2f}x" if arm in mean else "-"


def collect(root):
    """Summarise whatever exists; safe to rerun while jobs are still going."""
    report = [
        "# Generic-input gather engine: confirmation",
        "",
        "Seconds are the mean of two repeats; one sample per node; eight CPUs.",
        "Speedup is plain / arm. `plain` is the recipe alone, `direct` is",
        "--uncached on this branch, `pilot` is the trusted-input pilot engine.",
        "",
        "| Node | Sample | Recipe | Input | Plain s | Direct | Pilot | Cached s | Cached | Misses | Shards (duplicate sites) | Access | Equal to plain | Repeat identical |",
        "|---|---|---|---|---:|---:|---:|---:|---:|---:|---|---|---|---|",
    ]
    failures = []
    for path in sorted((root / "nodes").glob("*/state.json")):
        state = json.loads(path.read_text())
        failures += [(state["node"], f) for f in state["failures"]]
        audits = {(a["recipe"], a["assay"]): a for a in state["audits"]}
        for recipe, assay, arms in GROUPS:
            rows = [
                r for r in state["rows"] if (r["recipe"], r["assay"]) == (recipe, assay)
            ]
            mean = {
                arm: statistics.mean(r["seconds"] for r in rows if r["arm"] == arm)
                for arm in arms
                if any(r["arm"] == arm for r in rows)
            }
            if "plain" not in mean or "cached" not in mean:
                continue
            last = [r for r in rows if r["arm"] == "cached"][-1]
            audit = audits.get((recipe, assay), {})
            report.append(
                f"| {state['node']} | {state['sample']} | {recipe} | {assay} | "
                f"{mean['plain']:.1f} | {ratio(mean, 'direct')} | {ratio(mean, 'pilot')} | "
                f"{mean['cached']:.1f} | {ratio(mean, 'cached')} | {last['missing']} | "
                f"{last['shards']} ({last['duplicate_sites']}) | {last['cache_access']} | "
                f"{audit.get('records_equal_plain', 'pending')} / strict "
                f"{audit.get('semantic_pass', 'pending')} | "
                f"{audit.get('repeat_identical', 'pending')} |"
            )
        report.append(f"<!-- {state['node']}: {state['state']} -->")
    report += ["", f"Failures: {len(failures)}."]
    report += [
        f"- {n} {f['group']}: {f['error'].splitlines()[-1][:300]}" for n, f in failures
    ]
    (root / "RESULTS.md").write_text("\n".join(report) + "\n")
    print("\n".join(report))


def submit(root, deployment):
    """One exclusive job per node; run on the controller."""
    Path(head(root), "logs").mkdir(parents=True)
    jobs = []
    for name in SAMPLES:
        launch = [PYTHON, "-m", "benchmarks.gather_engine.run", "node", str(root)]
        launch += ["--node", name, "--deployment", deployment]
        command = ["sbatch", "--parsable", "--partition=benchmark", "--exclusive"]
        command += ["--nodes=1", "--ntasks=1", "--cpus-per-task=8", "--mem=28000M"]
        command += ["--time=03:00:00", "--chdir=/tmp", "--nodelist=" + name]
        command += ["--job-name=gather-" + name, f"--output={root}/logs/{name}-%j.out"]
        command += [
            "--wrap",
            f"PYTHONPATH={shlex.quote(deployment)} " + shlex.join(launch),
        ]
        jobs.append(subprocess.check_output(command, text=True).strip().split(";")[0])
    print(json.dumps(dict(jobs=jobs)))


def main():
    """Submit, run one node or summarise."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=["submit", "node", "collect"])
    parser.add_argument("root", type=Path)
    parser.add_argument("--node")
    parser.add_argument("--deployment")
    args = parser.parse_args()
    if args.action == "submit":
        submit(args.root, args.deployment)
    elif args.action == "node":
        node(args.root, args.node, args.deployment)
    else:
        collect(args.root)


if __name__ == "__main__":
    main()
