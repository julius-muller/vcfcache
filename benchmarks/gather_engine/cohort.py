"""Cohort benchmark of the gather engine against the plain annotation recipe.

Every unit is one exclusive Slurm job that runs on whichever node is free and
times the plain recipe and the cached run of the same sample on that node, one
eight-CPU step at a time. System bcftools, no wrappers, no engine switches.

    plain      the recipe alone (benchmarks/gather_engine/plain.py)
    cached     vcfcache annotate at the recipe's own output compression level
    cached_l1  vcfcache annotate at its default level 1 (fast groups only)

Tier A: all twelve samples of every group whose plain run is short.
Tier B: four samples of the three WGS groups whose plain run takes hours.
"""

import argparse
import csv
import hashlib
import json
import os
import re
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

from benchmarks.run_pilot import sha256sum, strict_semantic_compare, write_json_atomic

PYTHON = "/results/deployment/vcfcache-native-pilot-20261023c/.venv/bin/python"
PARENT = "/results/python-adaptive-manuscript-20261101-ee2162473/campaign/manifest.json"
BCFTOOLS = "/usr/local/bin/bcftools"
VEPRS = "/results/veprs05/dist/vep-0.5.0-x86_64-unknown-linux-gnu/vep-0.5.0-x86_64-unknown-linux-gnu/vep"
VEPRS_PARAMETERS = dict(
    annotation_tool_cmd=VEPRS,
    bcftools_cmd="bcftools",
    genome_build="GRCh38",
    optional_checks={},
    threads=8,
    tool_version_command=VEPRS + " --version",
    transcript_cache="/results/veprs/json_cache/GRCh38_native_20261023",
)
GROUPS = [
    ("veprs", "core"),
    ("veprs", "clinical"),
    ("fastvep", "clinical"),
    ("snpeff", "core"),
    ("snpeff", "clinical"),
    ("vep", "everything"),
]
# Plain runs longer than this go to tier B with fewer samples.
LONG = 3000
TIER_B_SAMPLES = 4
ASSAY_ORDER = ["panel", "wes", "wgs"]


def head(path):
    """Controller spelling of a worker-visible shared path."""
    return str(path).replace("/results/", "/mnt/data/slurm-results/", 1)


def plan(root, deployment, commit, smoke):
    """Freeze the design and submit every unit plus a final summary job."""
    parent = json.loads(Path(head(PARENT)).read_text())
    cells = []
    for task in parent["tasks"]:
        if (task["tool"], task["level"]) not in GROUPS:
            continue
        cell = dict(
            tool=task["tool"],
            level=task["level"],
            assembly=task["assembly"],
            assay=task["assay"],
            sample=task["sample"],
            input=task["input"],
            bed=((task.get("preprocess") or {}).get("bed") or {}).get("path"),
            cache=task["cache"]["directory"],
            parameters=task["parameters"],
            environment=task["environment"],
            expected=task["baseline"]["direct_seconds"],
        )
        if task["tool"] == "veprs":
            name = f"veprs05_{task['level']}"
            cell.update(
                cache=f"/results/veprs05/caches/{name}/cache/{name}",
                parameters=VEPRS_PARAMETERS,
                environment={},
                expected=(
                    dict(core=60, clinical=800)[task["level"]]
                    if task["assay"] == "wgs"
                    else 60
                ),
            )
        recipe = yaml.safe_load(
            Path(head(cell["cache"]), "annotation.yaml").read_text()
        )
        level = re.search(r"-Ob(\d)", recipe["annotation_cmd"])
        cell["recipe_level"] = int(level.group(1)) if level else 6
        cell["cache_sha256"] = sha256sum(
            Path(head(cell["cache"]), "vcfcache_annotated.bcf.csi")
        )
        cell["repeats"] = 2 if cell["expected"] < 120 else 1
        cell["arms"] = ["plain", "cached"]
        if cell["recipe_level"] != 1 and cell["expected"] < 1000:
            cell["arms"].append("cached_l1")
        cells.append(cell)
    units = []
    short = [c for c in cells if c["expected"] < LONG]
    for assembly, sample in sorted({(c["assembly"], c["sample"]) for c in short}):
        members = [
            c for c in short if (c["assembly"], c["sample"]) == (assembly, sample)
        ]
        members.sort(key=lambda c: (c["expected"], ASSAY_ORDER.index(c["assay"])))
        units.append(dict(id=f"A-{assembly}-{sample}", tier="A", cells=members))
    long = [c for c in cells if c["expected"] >= LONG]
    kinds = sorted({(c["tool"], c["level"], c["assembly"]) for c in long})
    chosen = {}
    for kind in kinds:
        samples = sorted(
            c["sample"] for c in long if (c["tool"], c["level"], c["assembly"]) == kind
        )
        # Alternate cohorts: sorted names start with every KPGP sample.
        half = TIER_B_SAMPLES // 2
        chosen[kind] = samples[:half] + samples[-half:]
    for rank in range(TIER_B_SAMPLES):
        for kind in kinds:
            cell = next(
                c
                for c in long
                if (c["tool"], c["level"], c["assembly"], c["sample"])
                == (*kind, chosen[kind][rank])
            )
            units.append(
                dict(
                    id="B-" + "-".join(kind) + "-" + cell["sample"],
                    tier="B",
                    cells=[cell],
                )
            )
    if smoke:
        # Panel cells of one sample per assembly: every tool and recipe once.
        picked = {}
        for candidate in units:
            assembly = candidate["id"].split("-")[1]
            if candidate["tier"] == "A" and assembly not in picked:
                picked[assembly] = [
                    c for c in candidate["cells"] if c["assay"] == "panel"
                ]
        units = [
            dict(id="smoke", tier="A", cells=[c for v in picked.values() for c in v])
        ]
    Path(head(root), "logs").mkdir(parents=True)
    write_json_atomic(
        Path(head(root), "manifest.json"),
        dict(
            deployment=deployment,
            commit=commit,
            bcftools=BCFTOOLS,
            parent_sha256=sha256sum(Path(head(PARENT))),
            design=__doc__,
            units=units,
        ),
    )
    jobs = {}
    for order, unit in enumerate(units):
        launch = [PYTHON, "-m", "benchmarks.gather_engine.cohort", "unit", str(root)]
        launch += ["--unit", unit["id"]]
        command = ["sbatch", "--parsable", "--partition=benchmark", "--exclusive"]
        command += ["--nodes=1", "--ntasks=1", "--cpus-per-task=8", "--mem=28000M"]
        command += ["--time=" + ("08:00:00" if unit["tier"] == "B" else "05:00:00")]
        command += ["--chdir=/tmp", f"--nice={order}", "--job-name=" + unit["id"][:40]]
        command += [f"--output={root}/logs/{unit['id']}-%j.out"]
        command += [
            "--wrap",
            f"PYTHONPATH={shlex.quote(deployment)} " + shlex.join(launch),
        ]
        jobs[unit["id"]] = (
            subprocess.check_output(command, text=True).strip().split(";")[0]
        )
    launch = [PYTHON, "-m", "benchmarks.gather_engine.cohort", "collect", str(root)]
    command = ["sbatch", "--parsable", "--partition=benchmark", "--ntasks=1"]
    command += ["--cpus-per-task=1", "--mem=3000M", "--time=00:10:00", "--chdir=/tmp"]
    command += [
        "--dependency=afterany:" + ":".join(jobs.values()),
        "--job-name=collect",
    ]
    command += [f"--output={root}/logs/collect-%j.out", f"--nice={len(units)}"]
    command += ["--wrap", f"PYTHONPATH={shlex.quote(deployment)} " + shlex.join(launch)]
    review = subprocess.check_output(command, text=True).strip().split(";")[0]
    write_json_atomic(
        Path(head(root), "submissions.json"), dict(jobs=jobs, collect=review)
    )
    print(
        json.dumps(
            dict(
                units=len(units),
                first=min(jobs.values(), key=int),
                last=max(jobs.values(), key=int),
                collect=review,
            )
        )
    )


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


def timed(deployment, scratch, evidence, name, arm, cell, source):
    """Run one arm as its own eight-CPU Slurm step."""
    directory = evidence / name
    directory.mkdir(parents=True)
    output = scratch / (name + ".bcf")
    stats = scratch / (name + "-stats")
    params = scratch / (name + ".yaml")
    (scratch / "tmp").mkdir(exist_ok=True)
    params.write_text(
        yaml.safe_dump(
            dict(cell["parameters"], threads=8, temp_dir=str(scratch / "tmp"))
        )
    )
    environment = dict(
        os.environ,
        **cell["environment"],
        PYTHONPATH=deployment,
        VCFCACHE_ROOT=deployment,
        LC_ALL="C",
    )
    if arm == "plain":
        command = [
            PYTHON,
            "-m",
            "benchmarks.gather_engine.plain",
            "--cache",
            cell["cache"],
        ]
        command += ["--input", str(source), "--output", str(output)]
        command += ["--params", str(params), "--work", str(stats)]
    else:
        level = 1 if arm == "cached_l1" else cell["recipe_level"]
        command = [PYTHON, "-c", "from vcfcache.cli import main; main()", "annotate"]
        command += ["-a", cell["cache"], "-i", str(source), "-o", str(output)]
        command += [
            "-y",
            str(params),
            "--statistics",
            "light",
            "--stats-dir",
            str(stats),
        ]
        command += ["--output-compression-level", str(level), "-v"]
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
            launch,
            env=environment,
            stdout=log,
            stderr=log,
            check=True,
            timeout=6 * 3600,
        )
    resources = json.loads((directory / "resources.json").read_text())
    assert resources["exit_code"] == 0 and len(resources["affinity"]) == 8
    row = dict(
        seconds=resources["wall_seconds"],
        cpu_seconds=resources["cpu_seconds"],
        peak_gib=resources["cgroup_after"]["memory_peak_bytes"] / 2**30,
        output_bytes=output.stat().st_size,
        output_records=int(
            subprocess.check_output([BCFTOOLS, "index", "-n", str(output)])
        ),
    )
    if arm != "plain":
        found = list(stats.rglob("compare_stats.yaml"))
        assert len(found) == 1
        counters = yaml.safe_load(found[0].read_text())
        row.update(
            missing=counters["variant_counts"].get("missing_variants"),
            shards=counters.get("shards"),
            bypassed_shards=counters.get("bypassed_shards"),
            cache_access=counters.get("cache_access"),
            hit_marker=counters.get("hit_marker"),
            output_level=counters.get("output_level"),
        )
    shutil.rmtree(stats, ignore_errors=True)
    return output, row


def unit(root, identifier):
    """All cells of one unit, on whichever node Slurm chose."""
    data = json.loads((root / "manifest.json").read_text())
    definition = next(u for u in data["units"] if u["id"] == identifier)
    evidence = root / "units" / identifier
    evidence.mkdir(parents=True, exist_ok=False)
    lscpu = subprocess.check_output(["lscpu"], text=True)
    state = dict(
        state="running",
        unit=identifier,
        node=socket.gethostname(),
        job_id=os.environ.get("SLURM_JOB_ID"),
        cpu_model=re.search(r"Model name:\s*(.+)", lscpu).group(1),
        bcftools=subprocess.check_output(
            [BCFTOOLS, "--version"], text=True
        ).splitlines()[0],
        rows=[],
        audits=[],
        failures=[],
        started=time.time(),
    )
    path = evidence / "state.json"
    write_json_atomic(path, state)
    scratch = Path(tempfile.mkdtemp(prefix="cohort-", dir="/mnt/data/tmp"))
    try:
        assert shutil.disk_usage(scratch).free > 40 * 2**30, "scratch nearly full"
        sources = {}
        for cell in definition["cells"]:
            key = (cell["sample"], cell["assembly"], cell["assay"])
            group = "-".join(
                [
                    cell["tool"],
                    cell["level"],
                    cell["assembly"],
                    cell["assay"],
                    cell["sample"],
                ]
            )
            state["phase"] = group
            write_json_atomic(path, state)
            try:
                if key not in sources:
                    # Assay extraction and BCF encoding only; nothing is
                    # split, normalized or removed.
                    source = scratch / ("-".join(key) + ".bcf")
                    command = [BCFTOOLS, "view", cell["input"]["path"], "-Ob1", "-W"]
                    command += ["--threads", "8", "-o", str(source)]
                    if cell["bed"]:
                        command += ["-R", cell["bed"]]
                    subprocess.run(command, check=True)
                    records = int(
                        subprocess.check_output([BCFTOOLS, "index", "-n", str(source)])
                    )
                    assert (
                        records == cell["input"]["records"]
                    ), "extraction differs from frozen assay"
                    source.chmod(0o444)
                    sources[key] = (source, records)
                source, records = sources[key]
                outputs = {}
                for repeat in range(1, cell["repeats"] + 1):
                    for arm in cell["arms"] if repeat % 2 else cell["arms"][::-1]:
                        output, row = timed(
                            data["deployment"],
                            scratch,
                            evidence,
                            f"{group}-{arm}-r{repeat}",
                            arm,
                            cell,
                            source,
                        )
                        assert row["output_records"] == records
                        outputs[(arm, repeat)] = output
                        state["rows"].append(
                            dict(
                                tool=cell["tool"],
                                level=cell["level"],
                                assembly=cell["assembly"],
                                assay=cell["assay"],
                                sample=cell["sample"],
                                arm=arm,
                                repeat=repeat,
                                records=records,
                                **row,
                            )
                        )
                        write_json_atomic(path, state)
                state["phase"] = group + " audit"
                write_json_atomic(path, state)
                comparison = strict_semantic_compare(
                    outputs[("cached", 1)], outputs[("plain", 1)]
                )
                write_json_atomic(evidence / f"{group}-against_plain.json", comparison)
                digest = record_digest(outputs[("cached", 1)])
                state["audits"].append(
                    dict(
                        tool=cell["tool"],
                        level=cell["level"],
                        assembly=cell["assembly"],
                        assay=cell["assay"],
                        sample=cell["sample"],
                        semantic_pass=comparison.get("semantic_pass"),
                        mismatches=comparison.get("record_mismatches"),
                        records_equal_plain=digest
                        == record_digest(outputs[("plain", 1)]),
                        other_runs_identical=all(
                            record_digest(output) == digest
                            for (arm, _), output in outputs.items()
                            if arm != "plain"
                        ),
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
        state.pop("phase", None)
    except Exception as error:
        state.update(state="failed", error=str(error), traceback=traceback.format_exc())
        raise
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
        write_json_atomic(path, state)


def collect(root):
    """Summarise whatever exists; safe to rerun while jobs are still going."""
    data = json.loads((root / "manifest.json").read_text())
    rows, audits, failures, states = [], [], [], {}
    for definition in data["units"]:
        path = root / "units" / definition["id"] / "state.json"
        if not path.exists():
            states[definition["id"]] = "pending"
            continue
        state = json.loads(path.read_text())
        states[definition["id"]] = state["state"]
        extra = dict(
            unit=state["unit"], node=state["node"], cpu_model=state["cpu_model"]
        )
        rows += [dict(extra, **r) for r in state["rows"]]
        audits += [dict(extra, **a) for a in state["audits"]]
        failures += [dict(extra, **f) for f in state["failures"]]
    for name, table in (("conditions.tsv", rows), ("audits.tsv", audits)):
        if table:
            with (root / name).open("w") as stream:
                fields = list(dict.fromkeys(k for r in table for k in r))
                writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
                writer.writeheader()
                writer.writerows(table)
    count = {
        s: sum(v == s for v in states.values()) for s in sorted(set(states.values()))
    }
    report = [
        "# Gather engine cohort: cached against the plain recipe",
        "",
        "Same-node plain control for every sample; eight CPUs; system bcftools.",
        "Speedup is the median over samples of plain / cached (per-sample means",
        "over repeats). `cached` writes the recipe's compression level,",
        "`level 1` is the vcfcache default. Strict comparisons are untimed.",
        "",
        f"Units: {count}. Failures: {len(failures)}.",
        "",
        "| Tool | Recipe | Assembly | Assay | Samples | Plain s | Cached s | Speedup | Range | Faster | Level 1 | Strict pass | Lines equal |",
        "|---|---|---|---|---:|---:|---:|---:|---|---:|---:|---:|---:|",
    ]
    kinds = sorted(
        {
            (r["tool"], r["level"], r["assembly"], ASSAY_ORDER.index(r["assay"]))
            for r in rows
        }
    )
    for tool, level, assembly, position in kinds:
        assay = ASSAY_ORDER[position]
        key = (tool, level, assembly, assay)
        own = [
            r for r in rows if (r["tool"], r["level"], r["assembly"], r["assay"]) == key
        ]
        ratios, fast, plain, cached = [], [], [], []
        for sample in sorted({r["sample"] for r in own}):
            mean = {
                arm: statistics.mean(
                    r["seconds"]
                    for r in own
                    if (r["sample"], r["arm"]) == (sample, arm)
                )
                for arm in {r["arm"] for r in own if r["sample"] == sample}
            }
            if "plain" in mean and "cached" in mean:
                ratios.append(mean["plain"] / mean["cached"])
                plain.append(mean["plain"])
                cached.append(mean["cached"])
            if "plain" in mean and "cached_l1" in mean:
                fast.append(mean["plain"] / mean["cached_l1"])
        checked = [
            a
            for a in audits
            if (a["tool"], a["level"], a["assembly"], a["assay"]) == key
        ]
        if not ratios:
            continue
        report.append(
            f"| {tool} | {level} | {assembly} | {assay} | {len(ratios)} | "
            f"{statistics.median(plain):.1f} | {statistics.median(cached):.1f} | "
            f"{statistics.median(ratios):.2f}x | {min(ratios):.2f}-{max(ratios):.2f}x | "
            f"{sum(r > 1 for r in ratios)}/{len(ratios)} | "
            + (f"{statistics.median(fast):.2f}x" if fast else "-")
            + f" | {sum(bool(a['semantic_pass']) for a in checked)}/{len(checked)} | "
            f"{sum(bool(a['records_equal_plain']) for a in checked)}/{len(checked)} |"
        )
    report += [""] + [
        f"- FAILED {f['unit']} {f['group']} on {f['node']}: "
        + f["error"].splitlines()[-1][:300]
        for f in failures
    ]
    (root / "RESULTS.md").write_text("\n".join(report) + "\n")
    print("\n".join(report))


def main():
    """Plan and submit, run one unit, or summarise."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=["plan", "unit", "collect"])
    parser.add_argument("root", type=Path)
    parser.add_argument("--deployment")
    parser.add_argument("--commit")
    parser.add_argument("--unit")
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    if args.action == "plan":
        plan(args.root, args.deployment, args.commit, args.smoke)
    elif args.action == "unit":
        unit(args.root, args.unit)
    else:
        collect(Path(head(args.root)) if not args.root.exists() else args.root)


if __name__ == "__main__":
    main()
