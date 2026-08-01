from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from benchmarks.run_strategy_comparison import (
    BUNDLED_CACHE_SPECS,
    VEP_CACHE_NAME,
    make_svg,
    public_strategies,
    select_samples,
)


def test_strategy_selection_is_deterministic_and_disjoint(tmp_path):
    selection = tmp_path / "selected.tsv"
    selection.write_text(
        "sample\tpopulation\tsuperpopulation\tsex\tseed\n"
        + "".join(
            f"{pop}{index}\tP\t{pop}\tfemale\tseed\n"
            for pop in ("AFR", "AMR", "EAS", "EUR", "SAS")
            for index in range(3)
        )
    )
    first = select_samples(selection)
    second = select_samples(selection)
    assert first == second
    assert len(first["training"]) == 3
    assert len(first["evaluation"]) == 3
    assert {row["sample"] for row in first["training"]}.isdisjoint(
        row["sample"] for row in first["evaluation"]
    )
    assert [row["superpopulation"] for row in first["training"]] == first[
        "superpopulations"
    ]
    assert [row["superpopulation"] for row in first["evaluation"]] == first[
        "superpopulations"
    ]


def test_strategy_svg_contains_all_strategies_and_no_external_assets(tmp_path):
    rows = []
    strategies = (
        "gnomad_af_0.1",
        "gnomad_af_0.01",
        "cohort_3_genomes",
    )
    for strategy_index, strategy in enumerate(strategies):
        for sample_index in range(3):
            rows.append(
                {
                    "strategy": strategy,
                    "cache_hit_rate": 0.5 + 0.05 * strategy_index,
                    "speedup": 2 + strategy_index + sample_index / 10,
                }
            )
    output = tmp_path / "figure.svg"
    make_svg(output, rows)
    svg = output.read_text()
    assert svg.startswith('<svg xmlns="http://www.w3.org/2000/svg"')
    assert "gnomAD" in svg
    assert "3-genome" in svg
    assert "held-out genomes" in svg
    assert "href=" not in svg


def test_public_strategies_require_verified_zenodo_bundles(tmp_path):
    with pytest.raises(FileNotFoundError):
        public_strategies(tmp_path)

    for spec in BUNDLED_CACHE_SPECS:
        root = tmp_path / spec.root_name
        (root / "blueprint").mkdir(parents=True)
        (root / "blueprint/vcfcache.bcf").write_text("blueprint")
        cache = root / "cache" / VEP_CACHE_NAME
        cache.mkdir(parents=True)
        (cache / "vcfcache_annotated.bcf").write_text("cache")
        (cache / "params.snapshot.yaml").write_text("genome_build: GRCh38\n")
        (root / "zenodo_provenance.json").write_text(
            json.dumps(
                {
                    "alias": spec.alias,
                    "doi": spec.doi,
                    "source": "zenodo_production",
                    "archive_name": spec.archive_name,
                    "archive_md5": spec.archive_md5,
                    "archive_md5_verified": True,
                }
            )
        )

    strategies = public_strategies(tmp_path)
    assert [strategy.name for strategy in strategies] == [
        "gnomad_af_0.1",
        "gnomad_af_0.01",
    ]
    assert all(strategy.kind == "bundled_zenodo" for strategy in strategies)
    assert all(strategy.doi for strategy in strategies)


def test_tracked_bundle_manifest_matches_runner_allow_list():
    manifest = (
        Path(__file__).resolve().parents[1] / "benchmarks/manifests/bundled_caches.tsv"
    )
    with manifest.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["alias"] for row in rows] == [
        spec.alias for spec in BUNDLED_CACHE_SPECS
    ]
    assert [row["doi"] for row in rows] == [spec.doi for spec in BUNDLED_CACHE_SPECS]
    assert [row["archive_md5"] for row in rows] == [
        spec.archive_md5 for spec in BUNDLED_CACHE_SPECS
    ]
