import csv
from pathlib import Path

import pytest

from benchmarks.summarize_final_campaign import read_and_pair, summarize


FIELDS = [
    "sample",
    "cohort",
    "tool",
    "level",
    "assembly",
    "assay",
    "condition",
    "wall_seconds",
    "input_records",
    "output_records",
    "missing_variants",
    "digest",
    "concordant",
]


def _write(path: Path, rows: list[list[str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(FIELDS)
        writer.writerows(rows)


def test_summary_uses_paired_speedups_and_excludes_failed_cells(tmp_path: Path) -> None:
    metrics = tmp_path / "metrics.tsv"
    common = ["cohort", "fastvep", "clinical", "GRCh38", "wgs"]
    _write(
        metrics,
        [
            ["s1", *common, "direct", "100", "100", "100", "", "a", ""],
            ["s1", *common, "cached", "20", "100", "100", "10", "a", "PASS"],
            ["s2", *common, "direct", "90", "200", "200", "", "b", ""],
            ["s2", *common, "cached", "30", "200", "200", "40", "c", "FAIL_CONTENT"],
        ],
    )
    summary, failures = summarize(read_and_pair(metrics))
    assert len(summary) == 1
    assert summary[0]["n_valid"] == 1
    assert summary[0]["n_failed"] == 1
    assert summary[0]["median_speedup"] == 5.0
    assert summary[0]["median_hit_rate"] == pytest.approx(0.9)
    assert failures[0]["sample"] == "s2"
    assert failures[0]["reason"] == "FAIL_CONTENT"


def test_incomplete_pair_is_rejected(tmp_path: Path) -> None:
    metrics = tmp_path / "metrics.tsv"
    _write(
        metrics,
        [["s1", "c", "t", "l", "a", "wgs", "direct", "1", "1", "1", "", "x", ""]],
    )
    with pytest.raises(ValueError, match="incomplete"):
        read_and_pair(metrics)
