"""Verify the flat deliverables and their numerical/assay identities."""

import csv
import hashlib
import json
import math
import statistics
import struct
from collections import defaultdict
from pathlib import Path

from openpyxl import load_workbook

ROOT = Path(__file__).resolve().parent


def read(name):
    """Read the exact flat exported data."""
    return list(csv.DictReader((ROOT / name).open(), delimiter="\t"))


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
        int(r["n_samples"]) == 12
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
        == hashlib.sha256((ROOT / "source_recent_pilots.json").read_bytes()).hexdigest()
    )
    build = json.loads((ROOT / "BUILD.json").read_text())
    assert (
        build["workbook_sha256"]
        == hashlib.sha256(
            (ROOT / "decision_runtime_comparison.xlsx").read_bytes()
        ).hexdigest()
    )
    expected_wgs = {
        "Python-adaptive": 0.8880397698862494,
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
    result = dict(
        passed=True,
        sheets=wb.sheetnames,
        configurations=len(configs),
        exported_timing_rows=len(rows),
        summary_cells=len(summary),
        png_figures=len(pngs),
        checks=[
            "identical row order",
            "numeric runtime cells",
            "true missing cells",
            "expandable detail columns",
            "15 distinct assay-specific/overall PNGs",
            "no subfolders",
            "fallback labels",
            "p95 support",
            "source/workbook hashes",
            "known cohort reference medians",
            "independent per-sample/node runtime and speedup recomputation",
        ],
    )
    (ROOT / "VALIDATION.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
