#!/usr/bin/env python3
"""Reanalyse existing measurements, excluding sl-w5 without new benchmarks."""

import argparse
import importlib.util
import json
import math
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

sys.dont_write_bytecode = True
from node_provenance import digest, key, load, read_tsv  # noqa: E402 (keep output flat)

ROOT = Path(__file__).resolve().parent
V1 = ROOT.parent / "decision_20261109"
INTEL = {"sl-w1", "sl-w2", "sl-w3", "sl-w4", "sl-w6"}
NUMERIC = {
    "repetition",
    "cached_s",
    "direct_s",
    "speedup",
    "saved_s",
    "saved_fraction",
    "input_n",
    "hit_fraction",
    "cpu_s",
    "peak_gib",
    "mismatches",
}

METRICS = {
    "Runtime median (s)": "runtime_median_s",
    "Speedup median": "speedup_median",
    "Samples": "n_samples",
    "Validation": "validation",
    "Runtime min (s)": "runtime_min_s",
    "Runtime max (s)": "runtime_max_s",
    "Runtime p95 (s)": "runtime_p95_s",
    "Direct median (s)": "direct_median_s",
    "Saved median (s)": "saved_median_s",
    "Speedup min": "speedup_min",
    "Speedup max": "speedup_max",
    "Speedup p95": "speedup_p95",
    "Sample-node units": "n_units",
    "Timed runs": "n_timings",
    "Sum measured wall (s)": "runtime_sum_s",
    "Hit fraction median": "hit_fraction_median",
    "CPU median (s)": "cpu_median_s",
    "Peak cgroup memory (GiB)": "peak_gib_max",
    "Variants median": "input_n_median",
    "Audit-different timings": "n_audit_differences",
    "Actual backends": "effective_backends",
}


def workbook(base, summary):
    """Keep v1's exact 98 rows and 14 recipe blocks, replacing EVERY metric.

    Excluded-only recipes/configurations remain genuinely blank. Styles, numeric
    formats and expandable columns are inherited, never historical cell values.
    """
    wb = base.load_workbook(V1 / "decision_runtime_comparison.xlsx")
    lookup = {(r["config_id"], r["assay"], r["recipe"]): r for r in summary}
    recipes = set()
    for ws in wb:
        ws.cell(1, 1, "VCFcache decision v2 — Intel-only — " + ws.title)
        ws.cell(
            3,
            1,
            "sl-w5/unknown cached or direct worker excluded; common cohort samples across engines. Blank = absent or excluded; see retained_counts.tsv. CPU placement NOT controlled by this filter.",
        )
        for i, cfg in enumerate(base.CONFIGS.values(), 6):
            for j, field in enumerate(
                ["label", "study", "cpus", "input_policy", "tools", "notes"], 1
            ):
                ws.cell(i, j, cfg[field])
            for j in range(7, ws.max_column + 1):
                rec, metric = ws.cell(5, j).value.split(" | ")
                recipes.add(rec)
                data = lookup.get((cfg["config_id"], ws.title, rec))
                cell = ws.cell(i, j)
                cell.value = data[METRICS[metric]] if data else None
                cell.comment = None
                if data and metric == "Runtime median (s)":
                    cell.comment = base.Comment(
                        "Intel-only reanalysis; "
                        + cfg["config_id"]
                        + "; retained samples="
                        + str(data["n_samples"])
                        + ". See measurements.tsv and excluded_measurements.tsv for node/reference proofs.",
                        "VCFcache",
                    )
    wb.save(ROOT / "decision_runtime_comparison.xlsx")
    return len(base.CONFIGS), len(recipes)


def main():
    """Filter frozen observations, preserve v1, and refresh derived exports."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--refresh",
        action="store_true",
        help="Read remote metadata only; never run annotation",
    )
    args = parser.parse_args()
    before = {p.name: digest(p.read_bytes()) for p in V1.iterdir() if p.is_file()}
    provenance = load(args.refresh)
    spec = importlib.util.spec_from_file_location(
        "decision_v1", V1 / "build_decision.py"
    )
    base = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(base)
    base.ROOT = ROOT
    configs = read_tsv(V1 / "configurations.tsv")
    for c in configs:
        c["cpus"] = int(c["cpus"])
        c["study"] = (
            c["study"]
            .replace("12 inputs/recipe", "filtered inputs/recipe")
            .replace("12 PGP inputs/recipe", "filtered PGP inputs/recipe")
        )
        if c["plot_group"] == "fusion":
            c["notes"] = (
                "Intel workers only; two technical repeats/sample/node, reduced within sample/node"
            )
        c[
            "notes"
        ] += "; v2 excludes sl-w5 in cached/direct arms and unknown placement; see exclusion audit"
        base.CONFIGS[c["config_id"]] = c
    rows = read_tsv(V1 / "measurements.tsv")
    reasons = {}
    for r in rows:
        p = provenance["rows"][key(r)]
        for field in NUMERIC:
            r[field] = float(r[field]) if r[field] else None
        r["node"] = p["node"]
        r["direct_node"] = p["direct_node"]
        r["cached_node_proof"] = p["cached_proof"]
        r["direct_node_proof"] = p["direct_proof"]
        bad = []
        if r["node"] not in INTEL:
            bad.append(
                "cached_sl-w5" if r["node"] == "sl-w5" else "cached_node_unknown"
            )
        if r["direct_s"] is not None:
            if r["direct_node"] not in INTEL:
                bad.append(
                    "direct_sl-w5"
                    if r["direct_node"] == "sl-w5"
                    else "direct_node_unknown"
                )
            if p["direct_seconds"] is None or not math.isclose(
                r["direct_s"], p["direct_seconds"], rel_tol=1e-9, abs_tol=1e-6
            ):
                bad.append("direct_duration_not_bound_to_node_proof")
        if bad:
            reasons[key(r)] = "; ".join(bad)

    # Use the common surviving sample set across engines, not separately
    # filtered cohorts that could bias an engine comparison by composition.
    pairs = defaultdict(list)
    for r in rows:
        group = base.CONFIGS[r["config_id"]]["plot_group"]
        if group in {"cohort-primary", "cohort-pgp", "wgs-strategies", "delay"}:
            cfg = base.CONFIGS[r["config_id"]]
            strategy = cfg["label"].rsplit("; ", 1)[-1]
            pairs[(group, r["condition_id"], r["assay"], r["recipe"], strategy)].append(
                r
            )
    for pair in pairs.values():
        assert len(pair) == 2, "Expected one Python and one Rust analysis view"
        if any(key(r) in reasons for r in pair):
            for r in pair:
                reasons.setdefault(
                    key(r), "paired_engine_excluded_to_keep_common_sample_set"
                )

    excluded = []
    for r in rows:
        if key(r) in reasons:
            excluded.append(dict(**r, exclusion_reason=reasons[key(r)]))
        else:
            base.DATA.append(r)
    for r in read_tsv(V1 / "sources.tsv"):
        base.SOURCES[r["source_id"]] = r
    for name in (
        "measurements.tsv",
        "configurations.tsv",
        "build_decision.py",
        "decision_runtime_comparison.xlsx",
    ):
        base.SOURCES["v1:" + name] = dict(
            source_id="v1:" + name, sha256=before[name], location=str(V1 / name)
        )
    base.SOURCES["worker-provenance"] = dict(
        source_id="worker-provenance",
        sha256=digest((ROOT / "node_provenance.json").read_bytes()),
        location="node_provenance.json",
    )
    audit = ROOT / "SLURM_CPU_PLACEMENT_AUDIT_20261109.md"
    if audit.exists():
        base.SOURCES["cpu-placement-audit"] = dict(
            source_id="cpu-placement-audit",
            sha256=digest(audit.read_bytes()),
            location="ubuntu@10.133.242.60:/home/ubuntu/portal-services/" + audit.name,
        )
    summary = base.summarize()
    base.export(ROOT / "measurements.tsv", base.DATA)
    base.export(ROOT / "excluded_measurements.tsv", excluded)
    base.export(ROOT / "summary.tsv", summary)
    base.export(ROOT / "configurations.tsv", list(base.CONFIGS.values()))
    base.export(ROOT / "sources.tsv", list(base.SOURCES.values()))
    # Keep source facts separate from the filtering rationale and actual counts.
    counts = []
    for cid, cfg in base.CONFIGS.items():
        original = [r for r in rows if r["config_id"] == cid]
        for assay, recipe in sorted({(r["assay"], r["recipe"]) for r in original}):
            old = [r for r in original if r["assay"] == assay and r["recipe"] == recipe]
            new = [
                r
                for r in base.DATA
                if r["config_id"] == cid
                and r["assay"] == assay
                and r["recipe"] == recipe
            ]
            counts.append(
                dict(
                    config_id=cid,
                    label=cfg["label"],
                    assay=assay,
                    recipe=recipe,
                    original_samples=len({r["sample"] for r in old}),
                    retained_samples=len({r["sample"] for r in new}),
                    original_timings=len(old),
                    retained_timings=len(new),
                )
            )
    base.export(ROOT / "retained_counts.tsv", counts)
    n, rec = workbook(base, summary)
    after = {p.name: digest(p.read_bytes()) for p in V1.iterdir() if p.is_file()}
    assert before == after, "v1 must not change"
    report = dict(
        created_utc=datetime.now(timezone.utc).isoformat(),
        configurations=n,
        recipes=rec,
        original_timing_rows=len(rows),
        exported_timing_rows=len(base.DATA),
        excluded_timing_rows=len(excluded),
        summary_cells=len(summary),
        excluded_reasons=dict(Counter(r["exclusion_reason"] for r in excluded)),
        workbook_sha256=digest(
            (ROOT / "decision_runtime_comparison.xlsx").read_bytes()
        ),
        v1_file_sha256=before,
        new_benchmark_jobs=0,
        policy=provenance["policy"],
    )
    (ROOT / "BUILD.json").write_text(json.dumps(report, indent=2) + "\n")
    print(
        json.dumps({k: v for k, v in report.items() if k != "v1_file_sha256"}, indent=2)
    )


if __name__ == "__main__":
    main()
