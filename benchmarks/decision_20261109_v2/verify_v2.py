"""Verify the flat deliverables and their numerical/assay identities."""

import csv
import hashlib
import json
import math
import statistics
import struct
import sys
from collections import defaultdict
from pathlib import Path

from openpyxl import load_workbook

sys.dont_write_bytecode = True
from build_v2 import METRICS  # noqa: E402 (keep output flat)
from node_provenance import key as observation_key  # noqa: E402 (keep output flat)

ROOT = Path(__file__).resolve().parent


def read(name):
    """Read the exact flat exported data."""
    return list(csv.DictReader((ROOT / name).open(), delimiter="\t"))


def read_tsv_v1(name):
    """Read original observations solely for preservation and pairing checks."""
    return list(
        csv.DictReader(
            (ROOT.parent / "decision_20261109" / name).open(), delimiter="\t"
        )
    )


def main():
    """Check sheet identity, blank semantics, plots, sources and key references."""
    configs = read("configurations.tsv")
    summary = read("summary.tsv")
    rows = read("measurements.tsv")
    sources = read("sources.tsv")
    wb = load_workbook(ROOT / "decision_runtime_comparison.xlsx", data_only=True)
    assert wb.sheetnames == ["Panel", "WES", "WGS"]
    expected = [r["label"] for r in configs]
    assert len(expected) == len(set(expected)) == 98
    lookup = {(r["config_id"], r["assay"], r["recipe"]): r for r in summary}
    # Independent recomputation from observations, not just workbook-vs-summary.
    grouped = defaultdict(list)
    for r in rows:
        grouped[(r["config_id"], r["assay"], r["recipe"])].append(r)
    for key, observed in grouped.items():
        cell = lookup[key]
        units = defaultdict(list)
        for r in observed:
            units[(r["sample"], r["node"])].append(r)
        runtimes = [
            statistics.median(float(r["cached_s"]) for r in unit)
            for unit in units.values()
        ]
        ratios = [
            statistics.median(
                float(r["direct_s"]) / float(r["cached_s"])
                for r in unit
                if r["direct_s"]
            )
            for unit in units.values()
            if any(r["direct_s"] for r in unit)
        ]
        assert int(cell["n_samples"]) == len({r["sample"] for r in observed})
        assert int(cell["n_timings"]) == len(observed)
        assert int(cell["n_units"]) == len(units)
        assert math.isclose(
            float(cell["runtime_median_s"]), statistics.median(runtimes), rel_tol=1e-12
        )
        assert math.isclose(
            float(cell["runtime_sum_s"]),
            sum(float(r["cached_s"]) for r in observed),
            rel_tol=1e-12,
        )
        assert math.isclose(float(cell["runtime_min_s"]), min(runtimes), rel_tol=1e-12)
        assert math.isclose(float(cell["runtime_max_s"]), max(runtimes), rel_tol=1e-12)
        if ratios:
            assert math.isclose(
                float(cell["speedup_median"]), statistics.median(ratios), rel_tol=1e-12
            )
    for ws in wb:
        assert [ws.cell(i, 1).value for i in range(6, 6 + len(configs))] == expected
        assert ws.max_row == len(configs) + 5
        assert ws.freeze_panes == "C6" and ws.auto_filter.ref
        assert all(
            ws.column_dimensions[x].hidden for x in ["D", "E", "F", "K", "L", "M"]
        )
        assert len(ws.conditional_formatting) > 0
        for col in range(7, ws.max_column + 1):
            rec, metric = ws.cell(5, col).value.split(" | ")
            for i, config in enumerate(configs, 6):
                entry = lookup.get((config["config_id"], ws.title, rec))
                observed = ws.cell(i, col).value
                expected_value = entry[METRICS[metric]] if entry else ""
                if expected_value == "":
                    assert observed is None
                elif METRICS[metric] in {"validation", "effective_backends"}:
                    assert observed == expected_value
                else:
                    assert math.isclose(
                        float(observed), float(expected_value), rel_tol=1e-12
                    )
        for col in range(7, ws.max_column + 1):
            header = ws.cell(5, col).value
            if not header.endswith(" | Runtime median (s)"):
                continue
            rec = header.split(" | ")[0]
            for i, cfg in enumerate(configs, 6):
                observed = ws.cell(i, col).value
                entry = lookup.get((cfg["config_id"], ws.title, rec))
                if not entry:
                    assert observed is None
                else:
                    assert isinstance(observed, (int, float)) and observed > 0
                    assert math.isclose(
                        observed, float(entry["runtime_median_s"]), rel_tol=1e-12
                    )
                    assert ws.cell(i, col + 2).value == int(entry["n_samples"])
    assert not any(p.is_dir() for p in ROOT.iterdir()), "Folder must stay flat"
    assert not list(ROOT.glob("*.pdf")) and not list(ROOT.glob("*.svg"))
    pngs = [x for x in (ROOT / "PNG_FILES.txt").read_text().splitlines() if x]
    assert len(pngs) == 15 and len(set(pngs)) == 15
    assert {p.name for p in ROOT.glob("*.png")} == set(pngs)
    for name in pngs:
        data = (ROOT / name).read_bytes()
        assert data[:8] == b"\x89PNG\r\n\x1a\n"
        width, height = struct.unpack(">II", data[16:24])
        assert width >= 2000 and height >= 1000
    for prefix in [
        "01_cohort_recipe_speedup",
        "02_historical_branch_speedup",
        "03_observed_cost_speedup",
        "07_latency_vs_direct",
    ]:
        hashes = {
            hashlib.sha256(
                (ROOT / (prefix + "_" + assay + ".png")).read_bytes()
            ).hexdigest()
            for assay in wb.sheetnames
        }
        assert (
            len(hashes) == 3
        ), "Assays must not accidentally share pooled/identical graphics"
    cfg = {r["config_id"]: r for r in configs}
    assert all(
        r["backend"] != "native"
        for r in rows
        if r["assay"] == "WGS" and cfg[r["config_id"]]["plot_group"] == "matrix"
    )
    assert all(
        1 <= int(r["n_samples"]) <= 12
        for r in summary
        if cfg[r["config_id"]]["plot_group"] == "cohort-primary"
    )
    assert all(not r["runtime_p95_s"] for r in summary if int(r["n_units"]) < 5)
    assert all(not r["speedup"] for r in rows if not r["direct_s"])
    assert all("1000g" not in r["source_location"].lower() for r in rows)
    assert all(float(r["cached_s"]) > 0 for r in rows)
    recent = next(r for r in sources if r["source_id"] == "recent-pilots")
    assert (
        recent["sha256"]
        == hashlib.sha256(
            (
                ROOT.parent / "decision_20261109" / "source_recent_pilots.json"
            ).read_bytes()
        ).hexdigest()
    )
    build = json.loads((ROOT / "BUILD.json").read_text())
    assert (
        build["workbook_sha256"]
        == hashlib.sha256(
            (ROOT / "decision_runtime_comparison.xlsx").read_bytes()
        ).hexdigest()
    )
    expected_wgs = {
        "Python-adaptive": 0.8611374471781069,
        "Rust-adaptive": 1.319682105167947,
    }
    for engine, value in expected_wgs.items():
        cell = next(
            r
            for r in summary
            if cfg[r["config_id"]]["plot_group"] == "cohort-primary"
            and cfg[r["config_id"]]["engine"] == engine
            and r["assay"] == "WGS"
            and r["recipe"] == "vep-rs core"
        )
        assert math.isclose(float(cell["speedup_median"]), value, rel_tol=1e-12)

    intel = {"sl-w1", "sl-w2", "sl-w3", "sl-w4", "sl-w6"}
    original = read_tsv_v1("measurements.tsv")
    omitted = read("excluded_measurements.tsv")
    originals = {observation_key(r): r for r in original}
    assert len(rows) + len(omitted) == len(original) == 1664
    assert {observation_key(r) for r in rows}.isdisjoint(
        {observation_key(r) for r in omitted}
    )
    assert {observation_key(r) for r in rows + omitted} == set(originals)
    metadata = json.loads((ROOT / "node_provenance.json").read_text())["rows"]
    for r in rows:
        assert r["node"] in intel
        assert not r["direct_s"] or r["direct_node"] in intel
        assert math.isclose(
            float(r["cached_s"]),
            float(originals[observation_key(r)]["cached_s"]),
            rel_tol=1e-12,
        )
        assert r["direct_s"] == originals[observation_key(r)]["direct_s"]
        proof = metadata[observation_key(r)]
        assert proof["node"] == r["node"]
        if r["direct_s"]:
            assert proof["direct_node"] == r["direct_node"]
            assert math.isclose(
                float(r["direct_s"]),
                proof["direct_seconds"],
                rel_tol=1e-9,
                abs_tol=1e-6,
            )
    paired = defaultdict(dict)
    for r in rows:
        if cfg[r["config_id"]]["plot_group"] == "cohort-primary":
            paired[(r["condition_id"], r["assay"], r["recipe"])][
                cfg[r["config_id"]]["engine"]
            ] = r
    assert len(paired) == 88
    assert all(
        set(pair) == {"Python-adaptive", "Rust-adaptive"} for pair in paired.values()
    )
    snpeff = [
        p
        for (_, assay, recipe), p in paired.items()
        if assay == "WGS" and recipe == "SnpEff clinical"
    ]
    assert len(snpeff) == 4
    gains = [
        float(p["Python-adaptive"]["cached_s"]) / float(p["Rust-adaptive"]["cached_s"])
        for p in snpeff
    ]
    assert 1.02 < min(gains) < max(gains) < 1.10
    for name, expected_sha in build["v1_file_sha256"].items():
        assert (
            hashlib.sha256(
                (ROOT.parent / "decision_20261109" / name).read_bytes()
            ).hexdigest()
            == expected_sha
        ), name
    assert build["new_benchmark_jobs"] == 0
    assert (
        len(
            {
                wb["WGS"].cell(5, j).value.split(" | ")[0]
                for j in range(7, wb["WGS"].max_column + 1)
            }
        )
        == 14
    )
    result = dict(
        passed=True,
        sheets=wb.sheetnames,
        configurations=len(configs),
        exported_timing_rows=len(rows),
        summary_cells=len(summary),
        png_figures=len(pngs),
        checks=[
            "sl-w5 excluded from cached runs AND direct denominators",
            "unknown worker proofs excluded",
            "88 common Intel-only cohort pairs",
            "SnpEff WGS four pairs with no hardware-class bimodality",
            "v1 files unchanged",
            "all original timings preserved numerically",
            "identical row order",
            "numeric runtime cells",
            "true missing cells",
            "expandable detail columns",
            "15 distinct assay-specific/overall PNGs",
            "no subfolders",
            "fallback labels",
            "p95 support",
            "source/workbook hashes",
            "filtered cohort reference medians",
            "independent per-sample/node runtime and speedup recomputation",
        ],
    )
    (ROOT / "VALIDATION.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
