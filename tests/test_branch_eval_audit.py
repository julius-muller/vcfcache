"""Backend-aware reporting never silently replaces superseded measurements."""

import json

import pytest
import yaml

from benchmarks.branch_eval.audit_engines import collect, effective_backend


@pytest.mark.parametrize(
    "arm,stats,expected",
    [
        ("native", {"engine": "fused"}, "adaptive"),
        ("native", {"engine": "native"}, "native"),
        ("original-fused", {"engine": "fused"}, "original_fused"),
        ("main-legacy", {}, "legacy"),
        ("auto", {}, "unrecorded"),
        ("direct", {}, "uncached"),
    ],
)
def test_backend_is_observed_not_requested(arm, stats, expected):
    assert effective_backend(arm, stats) == expected


def test_collect_preserves_supersession_failure_and_heavy_records(tmp_path):
    jobs = [
        {"id": 1, "status": "ABORTED: superseded by 2"},
        {"id": 2, "status": "complete"},
    ]
    (tmp_path / "RUN_MANIFEST.json").write_text(json.dumps({"jobs": jobs}))
    for job in jobs:
        condition = tmp_path / f"heavy-{job['id']}"
        cell_dir = condition / "benchmark/wes-native-1"
        stats_dir = cell_dir / "wes_vcstats"
        stats_dir.mkdir(parents=True)
        (condition / "CONDITIONS.txt").write_text(
            f"job={job['id']}\npipeline=synthetic_heavy\ncpus=8\n"
        )
        (cell_dir.parent / "measurements.json").write_text(
            json.dumps(
                [
                    dict(
                        input="/inputs/wes.bcf",
                        arm="native",
                        rep=1,
                        wall_seconds=12.0,
                        semantic_pass=False,
                    )
                ]
            )
        )
        (stats_dir / "compare_stats.yaml").write_text(
            yaml.safe_dump(
                dict(
                    engine="fused",
                    requested_engine="native",
                    fallback_reason="budget",
                    variant_counts=dict(input_variants=100, missing_variants=10),
                    stage_timings=[dict(stage="annotate-misses", wall_seconds=8.0)],
                )
            )
        )
        (cell_dir / "comparison.json").write_text(
            json.dumps(dict(record_mismatches=22))
        )
    rows = collect(tmp_path)
    assert len(rows) == 2
    assert [r["primary"] for r in rows] == [False, True]
    assert all(r["effective_backend"] == "adaptive" for r in rows)
    assert all(r["mismatched_records"] == 22 for r in rows)
    assert all(not r["semantic_pass"] for r in rows)
    assert all(r["family"] == "heavy" for r in rows)
    assert all(r["miss_annotation_seconds"] == 8.0 for r in rows)
    assert all(r["hit_fraction"] == 0.9 for r in rows)
