from __future__ import annotations

import json

import pytest

from benchmarks.analyze_assay_benchmark import bootstrap_interval, summarize
from benchmarks.collect_paired_benchmark import assay_label, collect_rows


def test_assay_labels_are_path_specific():
    assert assay_label("/samples/GRCh38/panel_acmg_sf_v3.3/S.vcf.gz") == "panel"
    assert assay_label("/samples/GRCh38/wes_twist_core/S.vcf.gz") == "wes"
    assert (
        assay_label("/samples/GRCh38/wes_twist_core_targets/S.vcf.gz")
        == "wes_target_control"
    )
    assert assay_label("/samples/GRCh38/1000g/EUR/S.vcf.gz") == "wgs"


def test_assay_summary_uses_samples_and_deterministic_bootstrap():
    rows = []
    for assay in ("panel", "wes", "wgs"):
        for index, relative in enumerate((0.1, 0.2, 0.3)):
            rows.append(
                {
                    "assay": assay,
                    "sample": f"{assay}-{index}",
                    "relative_runtime": str(relative),
                    "speedup": str(1 / relative),
                    "cache_hit_rate": str(1 - relative),
                    "input_records": "100",
                    "wall_seconds_saved": "80",
                    "semantic_pass": "true",
                    "annotation_mismatches": "0",
                    "ignored_annotation_mismatches": "1",
                    "annotation_order_only": "2",
                }
            )
    values = summarize(rows)
    assert [row["samples"] for row in values] == [3, 3, 3]
    assert values[0]["median_relative_runtime"] == pytest.approx(0.2)
    assert values[0]["unexpected_annotation_mismatches"] == 0
    assert bootstrap_interval([1, 2, 3], seed=1, replicates=100) == bootstrap_interval(
        [1, 2, 3], seed=1, replicates=100
    )


def test_paired_collector_requires_complete_semantically_valid_tasks(tmp_path):
    campaign = tmp_path / "campaign"
    (campaign / "measured/tasks/task-0/pilot/S/commit").mkdir(parents=True)
    (campaign / "campaign.json").write_text(
        json.dumps(
            {
                "campaign_id": "run",
                "phases": {"measured": {"tasks": 1}},
            }
        )
    )
    task = campaign / "measured/tasks/task-0"
    (task / "slurm-task.json").write_text(
        json.dumps(
            {
                "task_id": 0,
                "input_vcf": "/samples/GRCh38/wes_twist_core/S.vcf.gz",
                "input_sha256": "sha",
                "population": "POP",
                "superpopulation": "EUR",
                "first_mode": "cached",
                "slurm_job_id": "1",
                "slurm_node": "worker",
            }
        )
    )
    summary = {
        "sample": "S",
        "commit": "abc",
        "replicate": 1,
        "input_records": 100,
        "cache_hit_rate": 0.9,
        "cached_wall_seconds": 10,
        "uncached_wall_seconds": 100,
        "semantic_comparison": {
            "semantic_pass": True,
            "records_compared": 100,
            "key_mismatches": 0,
            "annotation_mismatches": 0,
        },
    }
    summary_path = task / "pilot/S/commit/summary_r01.json"
    summary_path.write_text(json.dumps(summary))
    rows = collect_rows(campaign, "measured")
    assert len(rows) == 1
    assert rows[0]["assay"] == "wes"
    assert rows[0]["relative_runtime"] == pytest.approx(0.1)
