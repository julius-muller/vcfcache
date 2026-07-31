from __future__ import annotations

import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

from benchmarks.prepare_assay_data import (
    ENSEMBL_115_GTF_URL,
    HPRC_QUOTAS,
    HPRC_R2_URL,
    HPRC_SELECTION_SEED,
    THOUSAND_GENOMES_X_NAME,
    TWIST_HG38_BED_URL,
)
from benchmarks.prepare_data import GIAB_SOURCES, SELECTION_SEED

MANIFESTS = Path(__file__).parents[1] / "benchmarks" / "manifests"


def _rows(name: str) -> list[dict[str, str]]:
    with (MANIFESTS / name).open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_tracked_1000g_selection_is_complete_balanced_and_unique():
    rows = _rows("selected_1000g_samples.tsv")

    assert len(rows) == 50
    assert len({row["sample"] for row in rows}) == 50
    assert {row["seed"] for row in rows} == {SELECTION_SEED}
    strata = Counter((row["superpopulation"], row["sex"]) for row in rows)
    assert strata == Counter(
        {
            (superpopulation, sex): 5
            for superpopulation in ("AFR", "AMR", "EAS", "EUR", "SAS")
            for sex in ("female", "male")
        }
    )


def test_tracked_hprc_selection_matches_quotas_and_does_not_overlap():
    primary = {row["sample"] for row in _rows("selected_1000g_samples.tsv")}
    rows = _rows("selected_hprc_r2_samples.tsv")
    hprc = {row["sample"] for row in rows}

    assert len(rows) == 20
    assert len(hprc) == 20
    assert not (primary & hprc)
    assert {row["seed"] for row in rows} == {HPRC_SELECTION_SEED}
    assert Counter(row["superpopulation"] for row in rows) == HPRC_QUOTAS
    sexes = defaultdict(set)
    for row in rows:
        sexes[row["superpopulation"]].add(row["sex"])
    assert all(values == {"female", "male"} for values in sexes.values())


def test_tracked_source_file_manifest_covers_every_primary_download():
    rows = _rows("source_files.tsv")
    assert len(rows) == 58  # 22 autosomal VCF/index pairs plus 7 GIAB pairs
    assert Counter(row["cohort"] for row in rows) == {"1000g": 44, "giab": 14}
    assert {row["sample_or_chromosome"] for row in rows if row["cohort"] == "giab"} == {
        source.sample for source in GIAB_SOURCES
    }
    for row in rows:
        assert row["url"].startswith("https://")
        if row["upstream_md5"]:
            assert re.fullmatch(r"[0-9a-f]{32}", row["upstream_md5"])


def test_tracked_assay_sources_match_pipeline_constants_and_are_portable():
    rows = _rows("assay_sources.tsv")
    urls = {row["url"] for row in rows}

    assert HPRC_R2_URL in urls
    assert TWIST_HG38_BED_URL in urls
    assert ENSEMBL_115_GTF_URL in urls
    assert any(THOUSAND_GENOMES_X_NAME in url for url in urls)
    for row in rows:
        assert not row["path"].startswith("/")
        assert int(row["bytes"]) > 0
        assert re.fullmatch(r"[0-9a-f]{64}", row["sha256"])


def test_tracked_region_and_resource_catalog_values_are_frozen():
    regions = {row["assay"]: row for row in _rows("assay_regions.tsv")}
    assert {
        key: (int(row["regions"]), int(row["target_bases"]), row["sha256"])
        for key, row in regions.items()
    } == {
        "wes_twist_core_targets": (
            191_723,
            33_074_111,
            "b1c4f29837061526d8524035524cdd745886c767c8f6819a47e4a75ea63c2221",
        ),
        "wes_twist_core": (
            165_014,
            78_026_576,
            "f007730dc6165a7ca74ea62f4b889d7b4e90be9e52583da56417ffc37f8833e3",
        ),
        "panel_acmg_sf_v3.3": (
            1_823,
            423_617,
            "11f6395868425252215de90d0ca485e97b3cb2c462d2ec468be009370b571450",
        ),
    }

    catalog_rows = _rows("source_catalog.tsv")
    catalog = {row["resource_id"]: row for row in catalog_rows}
    assert len(catalog) == len(catalog_rows)
    assert catalog["hprc_r2_wave"]["local_sha256"] == (
        "6316681cb75a418d518c657160e0e700b3252af6e01965d5304d1dccd11b365b"
    )
    assert catalog["gnomad_joint_ht"]["source_uri"].endswith(
        "/release/4.1/ht/joint/gnomad.joint.v4.1.sites.ht"
    )
    assert catalog["vcfcache_vep_cache"]["local_sha256"] == (
        "e9ae4a550c42a7cb1e68ccc1ef358c5e3558dab998c7bf833faa4f0571e1da7e"
    )
