"""Final analyses preserve denominators, failed attempts and validation scope."""

import copy
import csv
import json
import tarfile
from collections import Counter

import pytest

from benchmarks.performance_sprint.audit_final import build_threshold, seed_audit
from benchmarks.performance_sprint.snapshot_final import audit_cell, collect, cpu, peak
from benchmarks.run_pilot import REPO_ROOT, sha256sum

ROOT = REPO_ROOT / "benchmarks/performance_sprint/final-analysis-20261029"


def rows(name):
    with (ROOT / "tables" / name).open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def test_frozen_inventory_has_no_alias_inflation():
    cells, views, attempts = (
        rows("cells.tsv"),
        rows("plot_views.tsv"),
        rows("attempts.tsv"),
    )
    assert len(cells) == len({r["condition_id"] for r in cells}) == 600
    assert len(views) == 612
    assert len({r["measurement_id"] for r in views}) == 600
    assert sum(r["is_alias"] == "True" for r in views) == 12
    assert len(attempts) == 604
    assert sum(r["latest"] == "True" for r in attempts) == 600
    assert len({r["sample"] for r in cells}) == 54
    full = [r for r in views if r["family"] == "full52_wgs_strategies"]
    assert set(Counter((r["tool"], r["strategy"]) for r in full).values()) == {52}


def test_exclusions_never_become_validated_speedups():
    cells = rows("cells.tsv")
    assert Counter(r["evidence_class"] for r in cells) == {
        "validated": 572,
        "veprs_csq_timing_only": 24,
        "snpeff_warning_failed_strict": 4,
    }
    excluded = [r for r in cells if r["validated"] == "False"]
    assert len(excluded) == 28
    assert all(int(r["mismatched_records"]) > 0 for r in excluded)
    assert all(r["headers_equal"] == "True" for r in cells)
    assert all(r["engine_effective"] == "adaptive" for r in cells)


def test_measured_phase_wall_times_add_to_elapsed():
    fields = [
        "cli_other_seconds",
        "preprocess_seconds",
        "route_seconds",
        "assemble_seconds",
        "engine_other_seconds",
        "annotate_misses_seconds",
    ]
    for r in rows("cells.tsv"):
        assert sum(float(r[k]) for k in fields) == pytest.approx(
            float(r["cached_seconds"])
        )
        assert float(r["speedup"]) == pytest.approx(
            float(r["direct_seconds"]) / float(r["cached_seconds"])
        )
        assert float(r["cache_hit_rate"]) == pytest.approx(
            1 - int(r["missing_variants"]) / int(r["input_records"])
        )


def test_snapshot_archive_digest_is_frozen():
    snapshot = json.loads((ROOT / "source/SNAPSHOT.json").read_text())
    assert (
        sha256sum(ROOT / "source/terminal-evidence.tar.gz")
        == snapshot["archive_sha256"]
    )


def test_seed_audit_rejects_leakage_and_does_not_claim_gnomad_disjointness():
    strategies = dict(
        cohort_strategies=dict(kpgp=dict(training_samples=["seed"], assembly="GRCh38"))
    )
    cells = [dict(family="full52_wgs_strategies", cohort="kpgp", sample="test")]
    result = seed_audit(strategies, cells)
    assert result[0]["n_overlapping_ids"] == 0
    assert "not a gnomAD" in result[0]["definition"]
    cells[0]["sample"] = "seed"
    with pytest.raises(ValueError, match="leakage"):
        seed_audit(strategies, cells)


@pytest.mark.parametrize(
    "build,direct,cached,expected",
    [
        (100, 30, 20, 11),
        (99, 30, 20, 10),
        (0, 30, 20, 1),
        (100, 20, 20, None),
        (100, 19, 20, None),
    ],
)
def test_recorded_build_repayment_is_strictly_positive(build, direct, cached, expected):
    assert build_threshold(build, direct, cached) == expected


def test_resource_counters_use_step_delta_and_aggregate_memory():
    metrics = dict(
        cpu_seconds=99,
        cgroup_before=dict(cpu_stat=dict(usage_usec=100)),
        cgroup_after=dict(cpu_stat=dict(usage_usec=1_500_100), memory_peak_bytes=123),
    )
    assert cpu(metrics) == 1.5
    assert peak(metrics) == 123
    assert cpu(dict(cpu_seconds=4)) == 4


def test_complete_snapshot_cannot_be_overwritten_or_resumed(tmp_path):
    (tmp_path / "source").mkdir()
    (tmp_path / "source/SNAPSHOT.json").write_text("{}")
    with pytest.raises(FileExistsError, match="immutable"):
        collect(tmp_path, resume_preparation=True)


@pytest.fixture
def accepted_cell(tmp_path):
    r = rows("cells.tsv")[0]
    directory = tmp_path / "a001" / r["condition_id"]
    directory.mkdir(parents=True)
    with tarfile.open(ROOT / "source/terminal-evidence.tar.gz") as archive:
        manifest = json.load(archive.extractfile("main/manifest.json"))
        for name in [
            "state.json",
            "command.json",
            "statistics.json",
            "resources.json",
            "comparison.json",
        ]:
            (directory / name).write_bytes(
                archive.extractfile(r["source_archive_member"] + "/" + name).read()
            )
    state = json.loads((directory / "state.json").read_text())
    task = json.loads((directory / "command.json").read_text())["task"]
    return task, directory, state, manifest


def test_cell_audit_fails_closed_on_binding_counts_and_equivalence(accepted_cell):
    task, directory, state, manifest = accepted_cell
    audit_cell(task, directory, state, manifest, manifest, {})
    wrong = copy.deepcopy(task)
    wrong["input"]["records"] += 1
    with pytest.raises(ValueError, match="binding"):
        audit_cell(wrong, directory, state, manifest, manifest, {})
    comparison = json.loads((directory / "comparison.json").read_text())
    comparison["semantic_pass"] = False
    (directory / "comparison.json").write_text(json.dumps(comparison))
    with pytest.raises(ValueError, match="validation"):
        audit_cell(task, directory, state, manifest, manifest, {})


def test_payload_bound_builds_override_stale_strategy_costs():
    builds = rows("recorded_build_amortization.tsv")
    fast = [r for r in builds if r["tool"] == "fastvep"]
    assert len(fast) == 9
    assert all(r["provenance_archive_member"] for r in fast)
    assert all(float(r["payload_bound_historical_build_seconds"]) > 100 for r in fast)
    assert any(abs(float(r["recorded_cost_difference_seconds"])) > 40 for r in fast)
    assert all(r["kind"].startswith("locally_built_fastvep_") for r in fast)
