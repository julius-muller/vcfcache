#!/usr/bin/env python3
"""Collect existing evidence and write a flat, auditable decision workbook.

No benchmark jobs, source changes, or remote writes. uv invocation is in README.
Raw exports are SHA-bound; aliases, fallbacks and technical repeats stay explicit.
"""

import argparse
import csv
import hashlib
import io
import json
import math
import subprocess
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

from openpyxl import Workbook, load_workbook
from openpyxl.comments import Comment
from openpyxl.formatting.rule import ColorScaleRule, FormulaRule
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
PY = "367100833"
RUST = "e4d278ce1"
FUSION = "9ab232208"
HOST = "ubuntu@10.133.241.128"
CONFIGS, SOURCES, DATA = {}, {}, []
ASSAYS = ["Panel", "WES", "WGS"]


def number(value):
    """Missing remains missing, never silently zero."""
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def truth(value):
    """Parse explicit boolean exports without treating nonempty False as True."""
    return str(value).lower() == "true"


def sha(raw):
    """Bind an original export by its complete byte content."""
    return hashlib.sha256(raw).hexdigest()


def read_git(commit, path):
    """Read pinned repository objects, independent of dirty worktrees."""
    raw = subprocess.check_output(["git", "show", commit + ":" + path], cwd=REPO)
    key = commit + ":" + path
    SOURCES[key] = dict(source_id=key, sha256=sha(raw), location="git:" + key)
    return raw.decode()


def table(commit, path):
    """Parse an original tabular export from its pinned Git object."""
    return list(csv.DictReader(io.StringIO(read_git(commit, path)), delimiter="\t"))


def config(
    label,
    study,
    engine,
    cpus,
    policy="original",
    tools="historical frozen",
    notes="",
    group="other",
):
    """Define one explicitly named configuration, never merge evidence sets."""
    key = sha((study + "|" + label).encode())[:12]
    value = dict(
        config_id=key,
        label=label,
        study=study,
        engine=engine,
        cpus=cpus,
        input_policy=policy,
        tools=tools,
        notes=notes,
        plot_group=group,
    )
    if key in CONFIGS:
        assert CONFIGS[key] == value, label
    CONFIGS[key] = value
    return key


def recipe(tool, level):
    """Use explicit recipe names rather than imply annotators do identical work."""
    if tool == "veprs":
        return "vep-rs " + level
    if tool == "fastvep":
        return "fastVEP " + ("basic" if level == "original_external_recipe" else level)
    if tool == "snpeff":
        return "SnpEff " + level
    if level.startswith("everything_delay_"):
        return "VEP + wait " + str(int(level.split("_")[-1].replace("us", ""))) + " us"
    return "VEP everything"


def add(cid, source, assay, rec, sample, cached, direct=None, **more):
    """Append a measured complete-CLI observation with its own denominator."""
    assay = {"panel": "Panel", "wes": "WES", "wgs": "WGS"}.get(assay.lower(), assay)
    assert assay in ASSAYS, assay
    cached, direct = number(cached), number(direct)
    if cached is None or cached <= 0:
        return
    row = dict(
        config_id=cid,
        source_id=source,
        assay=assay,
        recipe=rec,
        sample=sample,
        node="",
        repetition=1,
        cached_s=cached,
        direct_s=direct,
        speedup=direct / cached if direct else None,
        saved_s=direct - cached if direct else None,
        saved_fraction=1 - cached / direct if direct else None,
        input_n=None,
        hit_fraction=None,
        cpu_s=None,
        peak_gib=None,
        validation="unknown",
        validation_scope="not exported",
        mismatches=None,
        condition_id="",
        backend=CONFIGS[cid]["engine"],
        source_location=source,
    )
    row.update(more)
    DATA.append(row)


def load_cohorts():
    """Keep assemblies/cache strategies separate; use aliases only as views."""
    wanted = {"native_clinical_assays", "vep_matched_assays", "veprs_extension"}
    for engine, folder, branch in [
        ("Python-adaptive", "python_adaptive_results", "feature/python-adaptive"),
        (
            "Rust-adaptive",
            "rust_adaptive_results",
            "feature/performance-sprint-20261027",
        ),
    ]:
        path = f"benchmarks/{folder}/final-20261103/tables/plot_views.tsv"
        src = PY + ":" + path
        for r in table(PY, path):
            assembly = r["assembly"]
            family = r["family"]
            strategy = r["strategy"]
            if family in wanted:
                study = (
                    "Cohort GRCh38 (12 inputs/recipe)"
                    if assembly == "GRCh38"
                    else "Cohort GRCh37 (12 PGP inputs/recipe)"
                )
                suffix = "8CPU; cohort " + assembly
                group = "cohort-primary" if assembly == "GRCh38" else "cohort-pgp"
            elif family == "full52_wgs_strategies":
                study = "WGS strategies " + assembly
                suffix = "8CPU; " + assembly + "; " + strategy
                group = "wgs-strategies"
            elif family == "synthetic_pipeline_spectrum":
                study = "Controlled WGS delay (one genome)"
                suffix = "8CPU; delay; " + strategy
                group = "delay"
            else:
                continue
            cid = config(
                branch
                + " ["
                + ("Rust-adaptive; " if engine == "Rust-adaptive" else "")
                + suffix
                + "]",
                study,
                engine,
                8,
                notes="Historical paired denominators reused; comparator scopes retained; not a randomized engine A/B",
                group=group,
            )
            status = (
                "accepted"
                if truth(r["validated"])
                else (
                    "timing-only"
                    if r["state"] == "complete_with_direct_difference"
                    else "failed audit"
                )
            )
            peak = number(r.get("cached_peak_gib"))
            if peak is None and number(r.get("cached_peak_bytes")) is not None:
                peak = float(r["cached_peak_bytes"]) / 1024**3
            add(
                cid,
                src,
                r["assay"],
                recipe(r["tool"], r["level"]),
                r["sample"],
                r["cached_seconds"],
                r["direct_seconds"],
                node=r["node"],
                input_n=number(r["input_records"]),
                hit_fraction=number(r["cache_hit_rate"]),
                cpu_s=number(r.get("cached_cpu_seconds")),
                peak_gib=peak,
                validation=status,
                validation_scope=r["validation_scope"],
                mismatches=number(
                    r.get("raw_record_mismatches", r.get("mismatched_records"))
                ),
                condition_id=r["measurement_id"],
                source_location=r["source_evidence"],
            )


def load_matrix():
    """Retain primary real-assay cells and distinguish executed fallbacks."""
    path = "benchmarks/branch_eval/integration-review-20261023/cells_backend_audit.tsv"
    rows = table(RUST, path)
    direct = {
        (r["recipe"], r["job"], r["scale"]): r
        for r in rows
        if truth(r["primary"]) and r["requested_arm"] == "direct"
    }
    arms = {
        "main-legacy": ("main@136bafa9", "legacy"),
        "legacy": ("feature/native-bcf-routing", "legacy@head"),
        "original-fused": ("feature/fused-routing", "original fused"),
        "auto": ("feature/adaptive-cache-routing", "Rust-adaptive"),
        "native": ("feature/native-bcf-routing", "Resident Rust"),
    }
    for r in rows:
        if (
            not truth(r["primary"])
            or r["scale"] not in {"panel", "wes", "wgs"}
            or r["requested_arm"] == "direct"
        ):
            continue
        branch, engine = arms[r["requested_arm"]]
        fallback = r["requested_arm"] == "native" and r["effective_backend"] != "native"
        engine = "adaptive fallback" if fallback else engine
        cid = config(
            branch + " [16CPU; " + engine + "]",
            "Archived branch matrix (one SGDP input)",
            engine,
            16,
            notes="Original fused is obsolete; native fallbacks are a separate row, not resident performance",
            group="matrix",
        )
        d = direct.get((r["recipe"], r["job"], r["scale"]))
        # VEP WES inherited a direct from its preceding job; resolve only an
        # unambiguous same-recipe/scale reference, never a different input size.
        if d is None:
            candidates = [
                x
                for x in direct.values()
                if x["recipe"] == r["recipe"] and x["scale"] == r["scale"]
            ]
            if len(candidates) == 1:
                d = candidates[0]
        rec = {
            "fastvep": "fastVEP clinical",
            "veprs": "vep-rs core",
            "veprs-clinical": "vep-rs clinical",
            "vep": "VEP everything",
        }[r["recipe"]]
        add(
            cid,
            RUST + ":" + path,
            r["scale"],
            rec,
            "LP6005519-DNA_A06",
            r["wall_seconds"],
            d["wall_seconds"] if d else None,
            input_n=number(r["records"]),
            hit_fraction=number(r["hit_fraction"]),
            validation="strict" if truth(r["semantic_pass"]) else "timing-only",
            mismatches=number(r["mismatched_records"]),
            validation_scope="complete-record audit; known HGNC_ID/ClinVar differences are NOT ignored here",
            backend=r["effective_backend"],
            condition_id=r["job"] + ":" + r["scale"] + ":" + r["requested_arm"],
            source_location=r["comparison"],
        )


def load_wes():
    """Import the finished resident, worker and decoder screens separately."""
    base = "benchmarks/wes_pilot/"
    src = base + "results-20261103/tables/cli_conditions.tsv"
    modes = {
        "adaptive": ("Rust-adaptive", "adaptive"),
        "resident": ("Resident Rust", "resident w4"),
        "coalesced": ("Resident Rust", "resident w4; CSI coalesced"),
    }
    for r in table(RUST, src):
        if r["mode"] == "direct":
            continue
        engine, mode = modes[r["mode"]]
        cid = config(
            "feature/rust-wes-pilot-20261103 [8CPU; " + mode + "]",
            "WES pilot (12 inputs/recipe)",
            engine,
            8,
            group="wes-cohort",
        )
        add(
            cid,
            RUST + ":" + src,
            "WES",
            recipe(r["tool"], "core" if r["tool"] == "veprs" else "clinical"),
            r["sample"],
            r["cached_seconds"],
            r["direct_seconds"],
            input_n=number(r["input_variants"]),
            hit_fraction=number(r["hit_fraction"]),
            cpu_s=number(r["cpu_seconds"]),
            peak_gib=number(r["peak_gib"]),
            validation="strict",
            validation_scope="complete-record/header, no ignored fields; retained direct",
            condition_id=r["id"] + ":" + r["mode"],
            source_location=r["campaign"] + ":" + r["job_id"],
        )
    src = base + "results-20261103/tables/worker_tuning.tsv"
    for r in table(RUST, src):
        cid = config(
            "feature/rust-wes-pilot-20261103 [8CPU; " + r["mode"] + "]",
            "WES worker allocation",
            "Resident Rust",
            8,
            notes="w1/w2 cover only two difficult inputs; w8 covers 12/recipe; not a sample-specific best-of",
            group="wes-workers",
        )
        add(
            cid,
            RUST + ":" + src,
            "WES",
            recipe(r["tool"], "core" if r["tool"] == "veprs" else "clinical"),
            r["sample"],
            r["cached_seconds"],
            r["direct_seconds"],
            input_n=number(r["input_variants"]),
            hit_fraction=number(r["hit_fraction"]),
            cpu_s=number(r["cpu_seconds"]),
            peak_gib=number(r["peak_gib"]),
            validation="strict",
            validation_scope="complete-record/header; retained direct",
            condition_id=r["id"] + ":" + r["mode"],
            source_location=r["campaign"] + ":" + r["job_id"],
        )
    src = base + "decoder-results-20261103/conditions.tsv"
    for r in table(RUST, src):
        cid = config(
            "feature/rust-wes-pilot-20261103 [8CPU; " + r["mode"] + "]",
            "WES fixed-worker decoder screen (3 inputs/recipe)",
            "Resident Rust",
            8,
            notes="Separate decoder screen; not pooled with earlier worker-allocation observations",
            group="decoder",
        )
        add(
            cid,
            RUST + ":" + src,
            "WES",
            recipe(r["tool"], "core" if r["tool"] == "veprs" else "clinical"),
            r["sample"],
            r["cli_seconds"],
            r["direct_seconds"],
            input_n=number(r["input_variants"]),
            hit_fraction=1 - float(r["missing_variants"]) / float(r["input_variants"]),
            cpu_s=number(r["cpu_seconds"]),
            peak_gib=number(r["peak_gib"]),
            validation="strict",
            validation_scope="complete-record/header; retained direct",
            condition_id=r["id"] + ":" + r["mode"],
            source_location=r["job_id"],
        )
    config(
        "feature/rust-wes-pilot-20261103 [whole-cache RAM; 20GiB stop]",
        "Capacity failure; no accepted assay timing",
        "Resident RAM",
        8,
        notes="Eager and lazy whole-cache loads hit the 20-GiB guard. All assay cells deliberately blank: loading is not annotation runtime.",
        group="capacity",
    )


def load_native_projection():
    """Import real WES projection experiments, not cache-loading microbenchmarks."""
    path = "benchmarks/native_bcf/summary-20261023-final/conditions.tsv"
    rows = [
        r
        for r in table(RUST, path)
        if r["role"] == "primary" and r["sample"].endswith(".wes")
    ]
    refs = {(r["experiment"], r["sample"]): r for r in rows if r["arm"] == "direct"}
    for r in rows:
        if r["arm"] == "direct":
            continue
        kind = (
            "projected WES cache"
            if r["experiment"].endswith("projected")
            else "full cache"
        )
        engine = "Resident Rust" if r["arm"] == "native" else "Rust-adaptive"
        cid = config(
            "feature/native-bcf-routing [16CPU; " + engine + "; " + kind + "]",
            "Resident/projection pilot (two SGDP inputs)",
            engine,
            16,
            notes="Projection build cost excluded; full and projected cache payloads differ",
            group="projection",
        )
        rec = (
            "fastVEP clinical"
            if r["experiment"].startswith("fastvep")
            else ("vep-rs clinical" if "clinical" in r["experiment"] else "vep-rs core")
        )
        direct = refs[(r["experiment"], r["sample"])]
        add(
            cid,
            RUST + ":" + path,
            "WES",
            rec,
            r["sample"].replace(".wes", ""),
            r["cli_seconds"],
            direct["cli_seconds"],
            input_n=number(r["input_records"]),
            hit_fraction=number(r["cache_hit_rate"]),
            cpu_s=number(r["cpu_seconds"]),
            peak_gib=float(r["peak_step_memory_bytes"]) / 1024**3,
            validation="strict" if truth(direct["semantic_pass"]) else "timing-only",
            validation_scope="Native-vs-adaptive alone is not independent validation; use the direct audit",
            condition_id=r["experiment"] + ":" + r["sample"] + ":" + r["arm"],
        )


def load_profile():
    """Retain labelled engineering configurations without inventing direct times."""
    path = "benchmarks/engine_profile/consolidated-20261027/canonical_conditions.tsv"
    for r in table(RUST, path):
        if r["role"] == "direct" or r["state"] == "skipped":
            continue
        # Sampling/instrumentation cells diagnose profiling, not optimization.
        if "instrumentation" in r["purpose"] or r["purpose"] == "stack-sampling":
            continue
        suffix = r["requested"] + "; w" + r["configured_workers"] + "; " + r["purpose"]
        if r["profile"] == "True":
            suffix += "; profiled"
        if r["form"] != "bcf":
            suffix += "; " + r["form"]
        if truth(r["skip_split"]):
            suffix += "; skip split"
        if r["purpose"] == "memory":
            suffix += "; mem=" + r["memory_mb"]
        if r["effective"] != r["requested"]:
            suffix += "; effective=" + r["effective"]
        cid = config(
            "feature/engine-overhead-profile [" + r["cpus"] + "CPU; " + suffix + "]",
            "Profiling/tuning (two SGDP inputs)",
            r["effective"],
            int(r["cpus"]),
            notes="Profiled/unprofiled and recovery VEP fork/buffer settings are distinct. CPU-scaling rows lack matched-CPU direct references.",
            group="profile",
        )
        rec = {
            "fastvep": "fastVEP clinical",
            "veprs-core": "vep-rs core",
            "vep": "VEP everything",
        }.get(r["recipe"], r["recipe"])
        direct = r["direct_seconds"] if r["cpus"] == "16" else None
        if r["recipe"] == "vep" and r["assay"] == "wgs":
            # Recovery VEP settings must not be pooled with the original WES.
            rec = "VEP everything (4 forks)"
        add(
            cid,
            RUST + ":" + path,
            r["assay"],
            rec,
            r["sample"],
            r["wall_seconds"],
            direct,
            node=r["group"],
            input_n=number(r["records"]),
            hit_fraction=number(r["hit_fraction"]),
            cpu_s=number(r["cpu_seconds"]),
            peak_gib=number(r["peak_gib"]),
            validation="strict" if truth(r["direct_pass"]) else "timing-only",
            validation_scope="profile canonical export; known VEP differences retained",
            condition_id=r["campaign"] + ":" + r["cell"],
            backend=r["effective"],
        )
    path = "benchmarks/performance_sprint/evidence-final-engineering-20261027/cells.tsv"
    for r in table(RUST, path):
        cid = config(
            "feature/performance-sprint-20261027 [16CPU; " + r["condition"] + "]",
            "Engineering sprint (two SGDP inputs)",
            (
                "Resident Rust"
                if "native" in r["condition"] or "window" in r["condition"]
                else "Rust-adaptive"
            ),
            16,
            notes="Cached timings available; this compact export does not bind direct elapsed seconds, so no direct speedup is invented",
            group="sprint",
        )
        add(
            cid,
            RUST + ":" + path,
            "WES" if r["recipe"] == "fastvep" else "WGS",
            "fastVEP clinical" if r["recipe"] == "fastvep" else "vep-rs core",
            r["sample"],
            r["wall_seconds"],
            validation=(
                "strict"
                if truth(r["direct_pass"]) and truth(r["cached_reference_pass"])
                else "failed audit"
            ),
            validation_scope="strict direct/cached reference; elapsed direct not included in source table",
            condition_id=r["campaign"] + ":" + r["cell"],
        )


def load_fusion():
    """Import accepted paired repeats, including their actual worker identities."""
    path = "benchmarks/pass_fusion_pilot/results-20261104-v2/conditions.tsv"
    rows = table(FUSION, path)
    refs = {(r["case"], r["repeat"]): r for r in rows if r["mode"] == "direct"}
    for r in rows:
        if r["mode"] == "direct":
            continue
        cid = config(
            "feature/python-pass-fusion-20261104 [8CPU; " + r["mode"] + "]",
            "Pass-fusion paired screen (two genomes, six nodes)",
            "Python-adaptive",
            8,
            notes="Two technical repeats per sample/node. Includes the explicitly labelled AMD sl-w5; not pooled with Intel-only pilots",
            group="fusion",
        )
        d = refs[(r["case"], r["repeat"])]
        add(
            cid,
            FUSION + ":" + path,
            "WGS",
            "vep-rs core",
            r["sample"],
            r["wall_seconds"],
            d["wall_seconds"],
            node=r["node"],
            repetition=int(r["repeat"]),
            input_n=number(r["input_variants"]),
            hit_fraction=number(r["hit_fraction"]),
            cpu_s=number(r["cpu_seconds"]),
            peak_gib=number(r["cgroup_peak_gib"]),
            validation="strict",
            validation_scope="64 strict audits accepted, including documented audit-only recovery",
            condition_id=r["case"] + ":" + r["repeat"] + ":" + r["mode"],
        )


def recent_snapshot(refresh=False):
    """Download small read-only checkpoints, not BCFs or remote modifications."""
    file = ROOT / "source_recent_pilots.json"
    if file.exists() and not refresh:
        return json.loads(file.read_text())
    script = r"""
import pathlib,json,hashlib
base=pathlib.Path('/mnt/data/slurm-results')
out={}
for name in ['input-preserving-20261105','adaptive-hts125-20261107','rust-hts125-20261107']:
 r=base/name; m=json.loads((r/'manifest.json').read_text()); cases=[]
 for c in m['cases']:
  if c['id']=='reader': continue
  p=r/'cases'/c['id']; state=json.loads((p/'state.json').read_text()); task=next(t for t in m['tasks'] if t['id']==c['task_id']); cells=[]
  for row in state['rows']:
   cell=p/row['mode']; payload=dict(row=row)
   for kind,fn in [('audit','against_direct.json'),('resources','resources.json'),('statistics','statistics.json')]:
    f=cell/fn
    if f.exists():
     raw=f.read_bytes(); x=json.loads(raw); payload[kind]=({k:v for k,v in x.items() if k!='examples'} if kind=='audit' else x); payload[kind+'_sha256']=hashlib.sha256(raw).hexdigest()
   cells.append(payload)
  cases.append(dict(definition=c,state={k:v for k,v in state.items() if k!='rows'},task={k:task[k] for k in ['id','tool','level','sample','assay','assembly','cohort','parameters','input','cache']},cells=cells))
 out[name]=dict(manifest_sha256=hashlib.sha256((r/'manifest.json').read_bytes()).hexdigest(),commit=m['commit'],cases=cases,source=str(r))
print(json.dumps(out,separators=(',',':')))
"""
    raw = subprocess.check_output(
        [
            "ssh",
            "-o",
            "BatchMode=yes",
            "-o",
            "ConnectTimeout=15",
            HOST,
            "cd /tmp && sudo -u appuser python3 -c "
            + __import__("shlex").quote(script),
        ]
    )
    data = json.loads(raw)
    file.write_text(json.dumps(data, indent=2) + "\n")
    return data


def load_recent(refresh=False):
    """Import verified main/Python/Rust toolchain pilots without cherry-picking."""
    recent = recent_snapshot(refresh)
    raw = (ROOT / "source_recent_pilots.json").read_bytes()
    SOURCES["recent-pilots"] = dict(
        source_id="recent-pilots",
        sha256=sha(raw),
        location=HOST + ":/mnt/data/slurm-results/{three immutable pilot roots}",
    )
    for name, pilot in recent.items():
        for case in pilot["cases"]:
            task = case["task"]
            refs = {
                x["row"]["mode"]: x
                for x in case["cells"]
                if "direct" in x["row"]["mode"]
            }
            for c in case["cells"]:
                mode = c["row"]["mode"]
                if "direct" in mode:
                    continue
                if name.startswith("input-"):
                    label = "main [8CPU; " + mode + "]"
                    engine = "legacy"
                    tools = (
                        "bcftools 1.21" if mode.startswith("old") else "bcftools 1.25"
                    )
                    policy = (
                        "original"
                        if mode.startswith("old")
                        else "input-preserving; timed identity/duplicate checks"
                    )
                    study, group = (
                        "Main input-preserving pilot (one KPGP input)",
                        "main125",
                    )
                    d = refs["old-direct" if mode == "old-cached" else "new-direct"]
                elif name.startswith("adaptive-"):
                    label = (
                        "feature/python-adaptive-hts125-pilot-20261107 [8CPU; "
                        + mode
                        + "]"
                    )
                    engine = "Python-adaptive"
                    tools = (
                        "bcftools 1.21" if mode == "baseline-121" else "bcftools 1.25"
                    )
                    policy = (
                        "original"
                        if mode.startswith("baseline")
                        else "prepared assertion; untimed preflight"
                    )
                    study, group = "Python HTSlib pilot (one KPGP input)", "python125"
                    d = refs["direct-125"]
                else:
                    label = "feature/rust-hts125-pilot-20261107 [8CPU; " + mode + "]"
                    engine = (
                        "Resident Rust"
                        if mode.startswith("resident")
                        else "Rust-adaptive"
                    )
                    tools = (
                        "external 1.21 / linked 1.19.1"
                        if "control121" in mode
                        else (
                            "external 1.25 / linked 1.19.1"
                            if "external125" in mode
                            else "external 1.25 / linked 1.25"
                        )
                    )
                    policy = (
                        "prepared assertion; untimed preflight"
                        if "prepared" in mode
                        else "original"
                    )
                    study, group = (
                        "Rust HTSlib pilot (two WES, one WGS input)",
                        "rust125",
                    )
                    d = refs["direct-prepared125"]
                cid = config(
                    label,
                    study,
                    engine,
                    8,
                    policy,
                    tools,
                    notes="Intel-only, one measurement/condition; do not attribute cross-pilot differences solely to engine",
                    group=group,
                )
                resource = c["resources"]
                assert resource["exit_code"] == 0 and len(resource["affinity"]) == 8
                audit = c["audit"]
                add(
                    cid,
                    "recent-pilots",
                    task["assay"],
                    recipe(task["tool"], task["level"]),
                    task["sample"],
                    c["row"]["seconds"],
                    d["row"]["seconds"],
                    node=case["state"]["node"],
                    input_n=number(c["row"]["records"]),
                    hit_fraction=(
                        1 - float(c["row"]["missing"]) / float(c["row"]["records"])
                        if c["row"].get("missing") is not None
                        else None
                    ),
                    cpu_s=number(resource["cpu_seconds"]),
                    peak_gib=(
                        number(
                            resource.get("cgroup_after", {}).get("memory_peak_bytes")
                        )
                        / 1024**3
                        if number(
                            resource.get("cgroup_after", {}).get("memory_peak_bytes")
                        )
                        is not None
                        else None
                    ),
                    validation="strict" if audit["semantic_pass"] else "timing-only",
                    validation_scope="complete-record/header, zero ignored fields",
                    mismatches=audit["record_mismatches"],
                    condition_id=case["definition"]["id"] + ":" + mode,
                    backend=c["row"].get("effective_engine", engine),
                    source_location=pilot["source"]
                    + "/cases/"
                    + case["definition"]["id"]
                    + "/"
                    + mode,
                )


def export(path, rows):
    """Export one flat, typed-on-import TSV with true empty missing values."""
    columns = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=columns, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def quantile(values, p):
    """Use linear interpolation, matching R's default empirical quantiles."""
    values = sorted(x for x in values if x is not None)
    if not values:
        return None
    index = (len(values) - 1) * p
    lo = math.floor(index)
    hi = math.ceil(index)
    return values[lo] + (values[hi] - values[lo]) * (index - lo)


def summarize():
    """Reduce technical repeats within sample/node before reporting distributions."""
    groups = defaultdict(list)
    for r in DATA:
        groups[(r["config_id"], r["assay"], r["recipe"])].append(r)
    result = []
    for (cid, assay, rec), rows in groups.items():
        units = defaultdict(list)
        for r in rows:
            units[(r["sample"], r["node"])].append(r)
        per = []
        for values in units.values():
            per.append(
                {
                    key: quantile([r.get(key) for r in values], 0.5)
                    for key in [
                        "cached_s",
                        "direct_s",
                        "speedup",
                        "saved_s",
                        "saved_fraction",
                        "input_n",
                        "hit_fraction",
                        "cpu_s",
                        "peak_gib",
                    ]
                }
            )

        def nums(field, values=per):
            return [r[field] for r in values if r[field] is not None]

        statuses = {r["validation"] for r in rows}
        validity = (
            "strict"
            if statuses == {"strict"}
            else (
                "accepted (recorded scope)"
                if statuses <= {"strict", "accepted"}
                else (
                    "timing-only / audit differences"
                    if "timing-only" in statuses
                    else "failed/unknown audit"
                )
            )
        )
        d = dict(
            config_id=cid,
            assay=assay,
            recipe=rec,
            n_samples=len({r["sample"] for r in rows}),
            n_units=len(per),
            n_timings=len(rows),
            runtime_median_s=quantile(nums("cached_s"), 0.5),
            speedup_median=quantile(nums("speedup"), 0.5),
            validation=validity,
            runtime_min_s=min(nums("cached_s")),
            runtime_max_s=max(nums("cached_s")),
            runtime_p95_s=quantile(nums("cached_s"), 0.95) if len(per) >= 5 else None,
            runtime_sum_s=sum(r["cached_s"] for r in rows),
            direct_median_s=quantile(nums("direct_s"), 0.5),
            saved_median_s=quantile(nums("saved_s"), 0.5),
            saved_fraction_median=quantile(nums("saved_fraction"), 0.5),
            speedup_min=min(nums("speedup")) if nums("speedup") else None,
            speedup_max=max(nums("speedup")) if nums("speedup") else None,
            speedup_p95=(
                quantile(nums("speedup"), 0.95) if len(nums("speedup")) >= 5 else None
            ),
            hit_fraction_median=quantile(nums("hit_fraction"), 0.5),
            cpu_median_s=quantile(nums("cpu_s"), 0.5),
            peak_gib_max=max(nums("peak_gib")) if nums("peak_gib") else None,
            input_n_median=quantile(nums("input_n"), 0.5),
            n_audit_differences=sum(
                r["validation"] in {"timing-only", "failed audit"} for r in rows
            ),
            effective_backends="; ".join(sorted({r["backend"] for r in rows})),
        )
        costs = [
            r["direct_s"] * 1e6 / r["input_n"]
            for r in per
            if r["direct_s"] and r["input_n"]
        ]
        d["direct_us_per_variant_median"] = quantile(costs, 0.5)
        result.append(d)
    return result


def workbook(summary):
    """Exactly three assay sheets; identical configuration and recipe ordering."""
    configs = list(CONFIGS.values())
    preferred = [
        "vep-rs core",
        "fastVEP basic",
        "fastVEP clinical",
        "SnpEff core",
        "vep-rs clinical",
        "SnpEff clinical",
        "VEP everything",
        "VEP everything (4 forks)",
    ]
    recipes = sorted(
        {r["recipe"] for r in summary},
        key=lambda x: (preferred.index(x) if x in preferred else 100, x),
    )
    visible = [
        ("runtime_median_s", "Runtime median (s)"),
        ("speedup_median", "Speedup median"),
        ("n_samples", "Samples"),
        ("validation", "Validation"),
    ]
    detail = [
        ("runtime_min_s", "Runtime min (s)"),
        ("runtime_max_s", "Runtime max (s)"),
        ("runtime_p95_s", "Runtime p95 (s)"),
        ("direct_median_s", "Direct median (s)"),
        ("saved_median_s", "Saved median (s)"),
        ("speedup_min", "Speedup min"),
        ("speedup_max", "Speedup max"),
        ("speedup_p95", "Speedup p95"),
        ("n_units", "Sample-node units"),
        ("n_timings", "Timed runs"),
        ("runtime_sum_s", "Sum measured wall (s)"),
        ("hit_fraction_median", "Hit fraction median"),
        ("cpu_median_s", "CPU median (s)"),
        ("peak_gib_max", "Peak cgroup memory (GiB)"),
        ("input_n_median", "Variants median"),
        ("n_audit_differences", "Audit-different timings"),
        ("effective_backends", "Actual backends"),
    ]
    metrics = visible + detail
    lookup = {(r["config_id"], r["assay"], r["recipe"]): r for r in summary}
    wb = Workbook()
    wb.remove(wb.active)
    navy = "17354A"
    blue = "DCEAF3"
    grey = "F0F0F0"
    amber = "FFF0C2"
    green = "D8EBD5"
    red = "F6D4D4"
    for assay in ASSAYS:
        ws = wb.create_sheet(assay)
        headers = [
            "Branch / exact configuration",
            "Evidence set",
            "Configured CPUs",
            "Input policy",
            "Toolchain",
            "Notes",
        ]
        last = len(headers) + len(metrics) * len(recipes)
        ws.merge_cells(start_row=1, start_column=1, end_row=1, end_column=min(last, 14))
        ws.cell(1, 1, "VCFcache decision overview — " + assay)
        ws.cell(1, 1).font = Font(size=16, bold=True, color=navy)
        ws.merge_cells(start_row=2, start_column=1, end_row=2, end_column=min(last, 20))
        ws.cell(
            2,
            1,
            "Runtime = complete CLI elapsed seconds, build excluded. Speedup = own recorded direct / cached; green >1, red <1. Expand + for details.",
        )
        ws.merge_cells(start_row=3, start_column=1, end_row=3, end_column=min(last, 20))
        ws.cell(
            3,
            1,
            "Blank = not measured. p95 shown only for >=5 sample-node units. Technical repeats reduced within sample/node. Historical/pilot/cohort strata NEVER pooled.",
        )
        for j, name in enumerate(headers, 1):
            ws.cell(5, j, name)
        for k, rec in enumerate(recipes):
            start = 7 + k * len(metrics)
            end = start + len(metrics) - 1
            ws.merge_cells(start_row=4, start_column=start, end_row=4, end_column=end)
            ws.cell(4, start, rec)
            ws.cell(4, start).fill = PatternFill("solid", fgColor=navy)
            ws.cell(4, start).font = Font(color="FFFFFF", bold=True)
            for j, (_, label) in enumerate(metrics, start):
                ws.cell(5, j, rec + " | " + label)
            ws.column_dimensions.group(
                get_column_letter(start + len(visible)),
                get_column_letter(end),
                outline_level=1,
                hidden=True,
            )
            ws.column_dimensions[get_column_letter(end)].collapsed = True
        for i, c in enumerate(configs, 6):
            values = [
                c["label"],
                c["study"],
                c["cpus"],
                c["input_policy"],
                c["tools"],
                c["notes"],
            ]
            for j, value in enumerate(values, 1):
                ws.cell(i, j, value)
            ws.row_dimensions[i].height = 44
            for k, rec in enumerate(recipes):
                row = lookup.get((c["config_id"], assay, rec))
                for j, (key, _) in enumerate(metrics, 7 + k * len(metrics)):
                    cell = ws.cell(i, j, row.get(key) if row else None)
                    if not row:
                        cell.fill = PatternFill("solid", fgColor=grey)
                    if isinstance(cell.value, (int, float)):
                        cell.number_format = (
                            "0.00"
                            if key
                            not in {
                                "n_samples",
                                "n_units",
                                "n_timings",
                                "input_n_median",
                                "n_audit_differences",
                            }
                            else "0"
                        )
                    if key == "validation" and row:
                        cell.fill = PatternFill(
                            "solid",
                            fgColor=(
                                green
                                if row[key] == "strict"
                                else (
                                    blue
                                    if row[key] == "accepted (recorded scope)"
                                    else amber
                                )
                            ),
                        )
                    if key == "runtime_median_s" and row:
                        cell.comment = Comment(
                            "Full end-to-end elapsed seconds. Quantiles describe exported observations, not confidence intervals. Actual backend: "
                            + row["effective_backends"]
                            + ". Evidence: "
                            + c["study"],
                            "VCFcache",
                        )
        for col in range(1, last + 1):
            ws.cell(5, col).fill = PatternFill("solid", fgColor=blue)
            ws.cell(5, col).font = Font(bold=True, color=navy)
            ws.cell(5, col).alignment = Alignment(wrap_text=True, vertical="center")
            ws.column_dimensions[get_column_letter(col)].width = (
                18 if col >= 7 else [64, 39, 13, 34, 33, 63][col - 1]
            )
        # Individual width definitions must not undo the hidden group ranges.
        for k in range(len(recipes)):
            start = 7 + k * len(metrics)
            for col in range(start + len(visible), start + len(metrics)):
                ws.column_dimensions[get_column_letter(col)].hidden = True
                ws.column_dimensions[get_column_letter(col)].outlineLevel = 1
        # Keep provenance/settings available, but do not consume the entire
        # initial screen before the first runtime column appears.
        for letter in ("D", "E", "F"):
            ws.column_dimensions[letter].hidden = True
            ws.column_dimensions[letter].outlineLevel = 1
        ws.column_dimensions["F"].collapsed = True
        for row in ws.iter_rows(min_row=6):
            for cell in row:
                cell.alignment = Alignment(wrap_text=True, vertical="center")
        for k in range(len(recipes)):
            start = 7 + k * len(metrics)
            for offset in [
                1,
                visible.__len__() + 5,
                visible.__len__() + 6,
                visible.__len__() + 7,
            ]:
                col = get_column_letter(start + offset)
                area = f"{col}6:{col}{5+len(configs)}"
                ws.conditional_formatting.add(
                    area,
                    FormulaRule(
                        formula=[f"AND(ISNUMBER({col}6),{col}6<1)"],
                        fill=PatternFill("solid", fgColor=red),
                    ),
                )
                ws.conditional_formatting.add(
                    area,
                    FormulaRule(
                        formula=[f"AND(ISNUMBER({col}6),{col}6>1.05)"],
                        fill=PatternFill("solid", fgColor=green),
                    ),
                )
            # Runtime colouring is within the SAME evidence set and recipe only.
            bystudy = defaultdict(list)
            for i, c in enumerate(configs, 6):
                if (c["config_id"], assay, recipes[k]) in lookup:
                    bystudy[c["study"]].append(i)
            for indices in bystudy.values():
                if len(indices) > 1:
                    area = " ".join(f"{get_column_letter(start)}{i}" for i in indices)
                    ws.conditional_formatting.add(
                        area,
                        ColorScaleRule(
                            start_type="min",
                            start_color="D8EBD5",
                            end_type="max",
                            end_color="F6D4D4",
                        ),
                    )
        ws.freeze_panes = "C6"
        ws.auto_filter.ref = f"A5:{get_column_letter(last)}{5+len(configs)}"
        ws.row_dimensions[5].height = 62
        ws.sheet_properties.outlinePr.summaryRight = True
        ws.sheet_view.zoomScale = 75
    path = ROOT / "decision_runtime_comparison.xlsx"
    wb.save(path)
    check = load_workbook(path, data_only=True)
    assert check.sheetnames == ASSAYS
    orders = [[ws.cell(i, 1).value for i in range(6, 6 + len(configs))] for ws in check]
    assert orders[0] == orders[1] == orders[2] and len(orders[0]) == len(set(orders[0]))
    return len(configs), len(recipes)


def main():
    """Collect read-only evidence, verify identities and create the three sheets."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--refresh",
        action="store_true",
        help="Read current remote immutable pilot metadata again",
    )
    args = parser.parse_args()
    load_cohorts()
    load_matrix()
    load_wes()
    load_native_projection()
    load_profile()
    load_fusion()
    load_recent(args.refresh)
    # Current decision candidates first, older exhaustive tuning screens later.
    group_order = [
        "main125",
        "cohort-primary",
        "wes-cohort",
        "wes-workers",
        "rust125",
        "python125",
        "matrix",
        "fusion",
        "decoder",
        "projection",
        "cohort-pgp",
        "wgs-strategies",
        "delay",
        "profile",
        "sprint",
        "capacity",
    ]
    ordered = sorted(
        CONFIGS.items(), key=lambda item: group_order.index(item[1]["plot_group"])
    )
    CONFIGS.clear()
    CONFIGS.update(ordered)
    # Each condition is unique within its explicitly separated analysis stratum.
    keys = [
        (
            r["config_id"],
            r["assay"],
            r["recipe"],
            r["condition_id"],
            r["node"],
            r["repetition"],
        )
        for r in DATA
    ]
    assert len(keys) == len(set(keys)), "Duplicate observations would distort summaries"
    assert not any("1000g" in r["source_location"].lower() for r in DATA)
    summary = summarize()
    export(ROOT / "measurements.tsv", DATA)
    export(ROOT / "summary.tsv", summary)
    export(ROOT / "configurations.tsv", list(CONFIGS.values()))
    export(ROOT / "sources.tsv", list(SOURCES.values()))
    n, rec = workbook(summary)
    result = dict(
        created_utc=datetime.now(timezone.utc).isoformat(),
        configurations=n,
        recipes=rec,
        exported_timing_rows=len(DATA),
        summary_cells=len(summary),
        workbook_sha256=sha((ROOT / "decision_runtime_comparison.xlsx").read_bytes()),
        notes="Exported rows/strata overlap; not independent biological sample counts. No new annotation or remote writes.",
    )
    (ROOT / "BUILD.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
