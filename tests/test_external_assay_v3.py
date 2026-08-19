from __future__ import annotations

import csv
from pathlib import Path

import pytest

from benchmarks.collect_external_assay_v3 import DEFAULT_WGS, matched_wgs_rows
from benchmarks.prepare_external_assay_v3 import (
    DEFAULT_METRICS,
    selected_samples,
    write_tsv,
)
from benchmarks.run_external_assay_task import (
    config_for,
    task_row,
    warmup_config_for,
)


def test_selected_samples_are_independent_grch38_and_deterministic() -> None:
    first = selected_samples(DEFAULT_METRICS, 6)
    second = selected_samples(DEFAULT_METRICS, 6)

    assert first == second
    assert len(first) == 12
    assert {row["cohort"] for row in first} == {"kpgp", "sgdp"}
    assert {row["assembly"] for row in first} == {"GRCh38"}
    assert len({(row["cohort"], row["sample"]) for row in first}) == 12


def test_matched_wgs_rows_form_two_tool_matrix() -> None:
    selected = selected_samples(DEFAULT_METRICS, 6)
    sample_keys = {(row["cohort"], row["sample"]) for row in selected}

    rows = matched_wgs_rows(DEFAULT_WGS, sample_keys)

    assert len(rows) == 24
    assert {row["tool"] for row in rows} == {"vep", "fastvep"}
    assert {row["assay"] for row in rows} == {"wgs"}
    assert {(row["cohort"], row["sample"]) for row in rows} == sample_keys


def test_manifest_writer_uses_unix_line_endings(tmp_path: Path) -> None:
    path = tmp_path / "manifest.tsv"
    write_tsv(path, [{"task_id": 0, "sample": "S1"}])

    assert b"\r" not in path.read_bytes()
    assert path.read_text() == "task_id\tsample\n0\tS1\n"


def test_task_row_is_position_stable_and_builds_config(tmp_path: Path) -> None:
    manifest = tmp_path / "tasks.tsv"
    fields = (
        "task_id", "tool", "assay", "cohort", "sample", "assembly",
        "input_vcf", "cache_strategy", "cache_kind", "cache_alias",
        "warmup_input_vcf", "cache_dir", "params_file", "condition_order",
        "replicate",
    )
    row = dict.fromkeys(fields, "")
    row.update(
        {
            "task_id": "0",
            "tool": "vep",
            "input_vcf": "/tmp/sample.vcf.gz",
            "warmup_input_vcf": "/tmp/warmup.vcf.gz",
            "cache_dir": "/tmp/cache",
            "params_file": "/tmp/params.yaml",
            "condition_order": "uncached,cached",
            "replicate": "1",
        }
    )
    with manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerow(row)

    observed = task_row(manifest, 0)
    config = config_for(observed, tmp_path / "run")
    warmup = warmup_config_for(observed, tmp_path / "run")

    assert config.input_vcf == Path("/tmp/sample.vcf.gz")
    assert config.cache_dir == Path("/tmp/cache")
    assert config.params_file == Path("/tmp/params.yaml")
    assert config.replicate == 1
    assert warmup.input_vcf == Path("/tmp/warmup.vcf.gz")
    assert warmup.data_root == tmp_path / "run" / "untimed_warmup"
    with pytest.raises(RuntimeError, match="absent or unstable"):
        task_row(manifest, 1)
