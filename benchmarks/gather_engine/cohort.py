"""Cohort benchmark of the gather engine against the plain annotation recipe.

Every unit is one exclusive Slurm job that runs on whichever node is free and
times the plain recipe and the cached run of the same sample on that node, one
eight-CPU step at a time. System bcftools, no wrappers, no engine switches.

    plain      the recipe alone (benchmarks/gather_engine/plain.py)
    cached     vcfcache annotate at the recipe's own output compression level
    cached_l1  vcfcache annotate at its default level 1 (fast groups only)
    exome      WES only: the same cache restricted to the exome capture BED

Tool versions: Ensembl VEP 116.3, fastVEP 0.5.0, vep-rs 0.5.0, snpEff as in
the parent manifest.

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
VEP116 = (
    "apptainer exec -B /mnt/data:/mnt/data /mnt/data/apps/ensembl-vep/116/vep.sif vep"
)
FASTVEP_OLD = "/mnt/data/apps/fastvep/0.3.0-publication"
FASTVEP = "/mnt/data/apps/fastvep/0.5.0"
EXOME_CACHES = "/results/exome-caches"
EXOME_BEDS = dict(
    GRCh38="/results/richness/regions/twist_human_core_exome_hg38.chr1-22-X.pad125.merged.bed",
    GRCh37="/results/richness/regions/twist_human_core_exome_hg19.chr1-22-X.pad125.merged_GRCh37.bed",
)


def upgraded(cell):
    """Point a parent-manifest cell at the current tool release and its cache."""
    tool, level, assembly = cell["tool"], cell["level"], cell["assembly"]
    if tool == "vep":
        cell["cache"] = (
            f"/results/vep116/caches/gnomad_v4.1_{assembly}_joint_af001"
            "/cache/vep116.3_everything"
        )
        cell["parameters"] = dict(
            cell["parameters"],
            annotation_tool_cmd=VEP116,
            tool_version_command=VEP116 + " | awk '/ensembl-vep/ {print $NF}'",
            vep_cache="/mnt/data/apps/ensembl-vep/116/cachedir",
        )
    elif tool == "fastvep":
        name = f"fastvep050_{level}_{assembly.lower()}_af001"
        cell["cache"] = f"/results/fastvep05/caches/{name}/cache/{name}"
        cell["parameters"] = {
            key: (
                value.replace(FASTVEP_OLD, FASTVEP) if isinstance(value, str) else value
            )
            for key, value in cell["parameters"].items()
        }
        assert FASTVEP_OLD not in json.dumps(cell["parameters"])
    elif tool == "veprs":
        name = f"veprs05_{level}"
        cell.update(
            cache=f"/results/veprs05/caches/{name}/cache/{name}",
            parameters=VEPRS_PARAMETERS,
            environment={},
            expected=(
                dict(core=60, clinical=800)[level] if cell["assay"] == "wgs" else 60
            ),
        )
    if cell["assay"] == "wes":
        name = Path(cell["cache"]).name + "_twist"
        cell["exome_cache"] = f"{EXOME_CACHES}/{name}/cache/{name}"
        cell["exome_bed"] = cell["bed"] or EXOME_BEDS[assembly]
    return cell


def head(path):
    """Controller spelling of a worker-visible shared path; unchanged on workers."""
    # Workers mount the share at /results; only the controller lacks it.
    if Path("/results").is_dir():
        return str(path)
    return str(path).replace("/results/", "/mnt/data/slurm-results/", 1)


def plan(root, deployment, commit, smoke, after=None):
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
        upgraded(cell)
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
        if "exome_cache" in cell:
            cell["exome_sha256"] = sha256sum(
                Path(head(cell["exome_cache"]), "vcfcache_annotated.bcf.csi")
            )
            cell["arms"].append("exome")
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
        # One sample per assembly: every panel cell (each tool and recipe
        # once) and the short WES cells, which also exercise the exome caches.
        picked = {}
        for candidate in units:
            assembly = candidate["id"].split("-")[1]
            if candidate["tier"] == "A" and assembly not in picked:
                picked[assembly] = [
                    c
                    for c in candidate["cells"]
                    if c["assay"] == "panel"
                    or (c["assay"] == "wes" and c["expected"] < 120)
                ]
        units = [
            dict(id="smoke", tier="A", cells=[c for v in picked.values() for c in v])
        ]
    submit(root, deployment, commit, units, after)


def submit(root, deployment, commit, units, after=None, **extra):
    """Freeze the units in a manifest and submit one job per unit plus a summary."""
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
            **extra,
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
        if after:
            command.append("--dependency=afterok:" + after)
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


def retry(
    old, root, deployment, commit, repeatability, only=None, overrides=None, rerun=None
):
    """Rerun failed cells of a finished cohort, plus plain-repeat controls.

    Each failed cell is rerun unchanged in a new unit. For every tool/assay
    named in ``repeatability`` the plain recipe is also run twice on that many
    samples, to tell annotator non-determinism from cache effects.
    """
    data = json.loads(Path(head(old), "manifest.json").read_text())
    units = []
    for definition in data["units"]:
        path = Path(head(old), "units", definition["id"], "state.json")
        state = json.loads(path.read_text()) if path.exists() else dict(audits=[])
        # Every cell without a completed comparison: failed, cut short or never
        # started.
        done = {group_name(a) for a in state["audits"]}
        cells = [c for c in definition["cells"] if group_name(c) not in done]
        if only and only not in definition["id"]:
            continue
        # "tool,level:parameter=value" replaces one recipe parameter, for both
        # the plain and the cached arm of the retried cells.
        for spec in overrides or []:
            kind, assignment = spec.split(":", 1)
            name, value = assignment.split("=", 1)
            for c in cells:
                if f"{c['tool']},{c['level']}" == kind:
                    assert name in c["parameters"], name
                    c["parameters"] = dict(c["parameters"], **{name: value})
                    c.setdefault("overrides", {})[name] = value
        if cells:
            units.append(
                dict(id="R-" + definition["id"], tier=definition["tier"], cells=cells)
            )
    for kind in rerun or []:
        # "tool,level,assembly,assay": every cell of that kind again, unchanged.
        tool, level, assembly, assay = kind.split(",")
        units += [
            dict(id=f"X-{tool}-{assay}-{c['sample']}", tier="A", cells=[c])
            for u in data["units"]
            for c in u["cells"]
            if (c["tool"], c["level"], c["assembly"], c["assay"])
            == (tool, level, assembly, assay)
        ]
    for spec in repeatability or []:
        kind, count = spec.rsplit(":", 1)
        tool, level, assembly, assay = kind.split(",")
        picked = [
            dict(c, arms=["plain", "plain_repeat"], repeats=1)
            for u in data["units"]
            for c in u["cells"]
            if (c["tool"], c["level"], c["assembly"], c["assay"])
            == (tool, level, assembly, assay)
        ][: int(count)]
        units += [
            dict(id=f"P-{tool}-{assay}-{c['sample']}", tier="A", cells=[c])
            for c in picked
        ]
    submit(root, deployment, commit, units, retry_of=str(old))


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


def loci(path):
    """Yield (CHROM, POS) and the record lines of each locus, in file order."""
    process = subprocess.Popen(
        [BCFTOOLS, "view", "-H", "--threads", "4", str(path)], stdout=subprocess.PIPE
    )
    assert process.stdout is not None
    key, lines = None, []
    with process.stdout:
        for line in process.stdout:
            first = line.find(b"\t")
            current = line[: line.find(b"\t", first + 1)]
            if current != key:
                if lines:
                    yield key, lines
                key, lines = current, []
            lines.append(line)
    if lines:
        yield key, lines
    if process.wait():
        raise RuntimeError(f"Cannot read {path}")


def canonical(line, unordered):
    """One record with INFO tags and annotation entries in sorted order.

    With ``unordered`` the ``&``-joined values inside every ``|`` field are
    sorted too: annotators emit some of these lists (protein domains, existing
    variant identifiers) in an order that varies between identical runs.
    """
    columns = line.rstrip(b"\n").split(b"\t")
    tags = []
    for tag in columns[7].split(b";"):
        name, _, value = tag.partition(b"=")
        entries = value.split(b",")
        if unordered:
            entries = [
                b"|".join(b"&".join(sorted(f.split(b"&"))) for f in e.split(b"|"))
                for e in entries
            ]
        tags.append(name + b"=" + b",".join(sorted(entries)))
    return (*columns[:7], tuple(sorted(tags)), *columns[8:])


def mismatch_report(candidate, reference, destination, limit=500):
    """Classify records that differ and keep them in full for later review.

    Returns the number of loci whose records differ after sorting INFO tags
    and annotation entries, and the number that still differ once ``&``-joined
    value lists are compared as sets.
    """
    import gzip

    ordered = unordered = kept = 0
    with gzip.open(destination, "wt") as stream:
        for (key_a, a), (key_b, b) in zip(
            loci(candidate), loci(reference), strict=True
        ):
            if key_a != key_b:
                return dict(report_error="locus sequence differs")
            if sorted(a) == sorted(b):
                continue
            strict = sorted(canonical(x, False) for x in a) != sorted(
                canonical(x, False) for x in b
            )
            loose = sorted(canonical(x, True) for x in a) != sorted(
                canonical(x, True) for x in b
            )
            ordered += strict
            unordered += loose
            if strict and kept < limit:
                kept += 1
                stream.write(
                    json.dumps(
                        dict(
                            locus=key_a.decode(),
                            differs_as_sets=loose,
                            candidate=[x.decode() for x in a],
                            reference=[x.decode() for x in b],
                        )
                    )
                    + "\n"
                )
    return dict(report_loci=ordered, report_loci_as_sets=unordered)


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
    if arm in PLAIN_ARMS:
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
        cache = cell["exome_cache"] if arm == "exome" else cell["cache"]
        command = [PYTHON, "-c", "from vcfcache.cli import main; main()", "annotate"]
        command += ["-a", cache, "-i", str(source), "-o", str(output)]
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
    if arm not in PLAIN_ARMS:
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


def group_name(cell):
    """Identifier of one cell within its unit."""
    return "-".join(
        [cell["tool"], cell["level"], cell["assembly"], cell["assay"], cell["sample"]]
    )


PLAIN_ARMS = ("plain", "plain_repeat")


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
                    if records != cell["input"]["records"]:
                        # Accept only the frozen file itself: its frozen record
                        # count may come from a since-dropped duplicate repair.
                        assert (
                            sha256sum(Path(cell["input"]["path"]))
                            == cell["input"]["sha256"]
                        ), "extraction differs from frozen assay"
                        state.setdefault("input_notes", []).append(
                            dict(
                                key=list(key),
                                frozen_records=cell["input"]["records"],
                                records=records,
                                input_sha256=cell["input"]["sha256"],
                            )
                        )
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
                candidate = (
                    outputs[("cached", 1)]
                    if ("cached", 1) in outputs
                    else outputs[("plain_repeat", 1)]
                )
                comparison = strict_semantic_compare(candidate, outputs[("plain", 1)])
                write_json_atomic(evidence / f"{group}-against_plain.json", comparison)
                digest = record_digest(candidate)
                report = {}
                if not comparison.get("semantic_pass"):
                    report = mismatch_report(
                        candidate,
                        outputs[("plain", 1)],
                        evidence / f"{group}-mismatches.jsonl.gz",
                    )
                state["audits"].append(
                    dict(
                        tool=cell["tool"],
                        level=cell["level"],
                        assembly=cell["assembly"],
                        assay=cell["assay"],
                        sample=cell["sample"],
                        semantic_pass=comparison.get("semantic_pass"),
                        mismatches=comparison.get("record_mismatches"),
                        **report,
                        records_equal_plain=digest
                        == record_digest(outputs[("plain", 1)]),
                        other_runs_identical=all(
                            record_digest(output) == digest
                            for (arm, _), output in outputs.items()
                            if arm not in PLAIN_ARMS
                        ),
                        index=Path(str(candidate) + ".csi").exists(),
                        candidate=(
                            "cached" if ("cached", 1) in outputs else "plain_repeat"
                        ),
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
        "`level 1` is the vcfcache default. `Exome cache` (WES only) uses the",
        "same cache restricted to the capture BED. Strict comparisons are untimed;",
        "`Lines equal` also requires every cached run, exome included, to match.",
        "",
        f"Units: {count}. Failures: {len(failures)}.",
        "",
        "| Tool | Recipe | Assembly | Assay | Samples | Plain s | Cached s | Speedup | Range | Faster | Level 1 | Exome cache | Strict pass | Lines equal |",
        "|---|---|---|---|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|",
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
        ratios, fast, plain, cached, exome = [], [], [], [], []
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
            if "plain" in mean and "exome" in mean:
                exome.append(mean["plain"] / mean["exome"])
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
            + " | "
            + (f"{statistics.median(exome):.2f}x" if exome else "-")
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


def gate(root):
    """Fail unless every unit finished without failures and every audit passed."""
    data = json.loads((root / "manifest.json").read_text())
    problems = []
    for definition in data["units"]:
        path = root / "units" / definition["id"] / "state.json"
        state = json.loads(path.read_text()) if path.exists() else dict(state="missing")
        if state["state"] != "complete":
            problems.append(f"{definition['id']}: {state['state']}")
            continue
        expected = sum(1 for _ in definition["cells"])
        problems += [f"{f['group']}: {f['error'][:200]}" for f in state["failures"]]
        problems += [
            f"{a['tool']} {a['level']} {a['assembly']} {a['assay']}: comparison failed"
            for a in state["audits"]
            if not (a["semantic_pass"] and a["other_runs_identical"] and a["index"])
        ]
        if len(state["audits"]) != expected:
            problems.append(
                f"{definition['id']}: {len(state['audits'])}/{expected} audits"
            )
    collect(root)
    if problems:
        print("GATE FAILED\n" + "\n".join(problems))
        raise SystemExit(1)
    print("GATE PASSED")


def main():
    """Plan and submit, run one unit, or summarise."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=["plan", "unit", "collect", "gate", "retry"])
    parser.add_argument("root", type=Path)
    parser.add_argument("--deployment")
    parser.add_argument("--commit")
    parser.add_argument("--unit")
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--after")
    parser.add_argument("--retry-of", type=Path)
    parser.add_argument("--repeatability", nargs="*")
    parser.add_argument("--only")
    parser.add_argument("--override", nargs="*")
    parser.add_argument("--rerun", nargs="*")
    args = parser.parse_args()
    if args.action == "plan":
        plan(args.root, args.deployment, args.commit, args.smoke, args.after)
    elif args.action == "unit":
        unit(args.root, args.unit)
    elif args.action == "retry":
        retry(
            args.retry_of,
            args.root,
            args.deployment,
            args.commit,
            args.repeatability,
            args.only,
            args.override,
            args.rerun,
        )
    elif args.action == "gate":
        gate(Path(head(args.root)) if not args.root.exists() else args.root)
    else:
        collect(Path(head(args.root)) if not args.root.exists() else args.root)


if __name__ == "__main__":
    main()
