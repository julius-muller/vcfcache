from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from benchmarks.analyze_controlled_runtime import fit_models
from benchmarks.prepare_controlled_runtime import (
    HIT_RATES,
    PIPELINES,
    _placeholder_cache,
    annotation_yaml,
)
from benchmarks.run_controlled_cohort import prepare

ANNOTATION = """genome_build: GRCh38
annotation_cmd: |
  vep \\
    --offline \\
    --hgvsg \\
    --everything \\
    --stats_file ${AUXILIARY_DIR}/stats.html \\
    -o ${OUTPUT_BCF}
must_contain_info_tag: CSQ
"""


def test_controlled_recipes_remove_everything_and_add_only_requested_delay():
    vanilla = annotation_yaml(ANNOTATION)
    delayed = annotation_yaml(ANNOTATION, 5_000)
    assert "--everything" not in vanilla
    assert "--hgvsg" not in vanilla
    assert "--plugin" not in vanilla
    assert "--plugin SyntheticDelay,delay_us=5000" in delayed
    assert delayed.count("--dir_plugins") == 1
    assert "--stats_file" in delayed


def test_placeholder_cache_has_a_valid_database_layout(tmp_path):
    bundled_root = tmp_path / "bundled"
    bundled = bundled_root / "cache/everything"
    bundled.mkdir(parents=True)
    (bundled / "vcfcache_annotated.bcf").write_bytes(b"bcf")
    (bundled / "vcfcache_annotated.bcf.csi").write_bytes(b"csi")
    (bundled_root / "workflow").mkdir()
    (bundled_root / "workflow/init.yaml").write_text("genome_build: GRCh38\n")
    params = tmp_path / "params.yaml"
    params.write_text("genome_build: GRCh38\n")

    database = tmp_path / "database"
    cache = _placeholder_cache(database, ANNOTATION, params, bundled)

    assert cache == database / "cache/placeholder"
    assert (database / "blueprint").is_dir()
    assert (database / "workflow/init.yaml").is_file()
    assert (cache / "vcfcache_annotated.bcf").is_symlink()


def _controlled_root(root: Path) -> Path:
    controlled = root / "controlled"
    input_vcf = root / "HG02374.vcf.gz"
    input_vcf.write_bytes(b"input")
    (controlled / "READY.json").parent.mkdir(parents=True)
    (controlled / "READY.json").write_text(
        json.dumps(
            {
                "input_vcf": str(input_vcf),
                "input_sha256": "input-sha",
            }
        )
    )
    for pipeline in PIPELINES:
        for hit_rate in HIT_RATES:
            cache = controlled / "caches" / pipeline / f"hit-{hit_rate:03d}"
            cache.mkdir(parents=True)
            for name in (
                "annotation.yaml",
                "params.snapshot.yaml",
                "vcfcache_annotated.bcf",
                "vcfcache_annotated.bcf.csi",
            ):
                (cache / name).write_text(name)
            (cache / "controlled_cache.json").write_text(
                json.dumps({"pipeline": pipeline, "target_hit_rate": hit_rate / 100})
            )
    return controlled


def test_prepare_freezes_four_baselines_and_twenty_cached_tasks(
    tmp_path, monkeypatch, capsys
):
    controlled = _controlled_root(tmp_path)
    controller = tmp_path / "results"
    monkeypatch.setattr(
        "benchmarks.run_controlled_cohort.git_output", lambda *_args: "commit"
    )
    monkeypatch.setattr(
        "benchmarks.run_controlled_cohort.sha256sum", lambda _path: "ready-sha"
    )
    args = type(
        "Args",
        (),
        {
            "campaign_id": "controlled-test",
            "controlled_root": controlled,
            "controller_results": controller,
            "worker_results": Path("/results"),
        },
    )()
    prepare(args)
    capsys.readouterr()
    root = controller / "campaigns/controlled-test"
    with (root / "manifests/baseline.tsv").open(newline="") as handle:
        baselines = list(csv.DictReader(handle, delimiter="\t"))
    with (root / "manifests/cached.tsv").open(newline="") as handle:
        cached = list(csv.DictReader(handle, delimiter="\t"))
    assert len(baselines) == len(PIPELINES) == 4
    assert len(cached) == len(PIPELINES) * len(HIT_RATES) == 20
    assert (root / "logs").is_dir()
    assert {float(row["target_hit_rate"]) for row in cached} == {
        value / 100 for value in HIT_RATES
    }
    assert all(row["baseline_result"].startswith("/results/") for row in cached)
    campaign = json.loads((root / "campaign.json").read_text())
    assert campaign["phases"] == {"baseline": 4, "cached": 20}


def test_controlled_model_uses_prespecified_unit_slope_and_robust_overhead():
    rows = []
    for pipeline in PIPELINES:
        for hit_rate in HIT_RATES:
            observed = hit_rate / 100
            rows.append(
                {
                    "pipeline": pipeline,
                    "delay_us": str(PIPELINES[pipeline] or 0),
                    "target_hit_rate": str(observed),
                    "observed_hit_rate": str(observed),
                    "cached_wall_seconds": str(5 + (1 - observed) * 100),
                    "uncached_wall_seconds": "100",
                    "relative_runtime": str(0.05 + (1 - observed)),
                    "speedup": str(100 / (5 + (1 - observed) * 100)),
                    "semantic_pass": "True",
                }
            )
    modeled, model = fit_models(rows)
    assert len(modeled) == 20
    assert model["pipelines"]["everything"][
        "lookup_preprocessing_overhead_seconds"
    ] == pytest.approx(5)
    assert all(row["residual_seconds"] == pytest.approx(0) for row in modeled)
