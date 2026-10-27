"""Candidate experiments must preserve complete records and remain disabled by default."""

import pytest

from benchmarks.run_pilot import strict_semantic_compare
from tests.test_fused_engine import dataset, run_case
from tests.test_fused_engine import native as native_fixture

native = native_fixture


@pytest.mark.parametrize("hits", [0, 1, 2, 10])
@pytest.mark.parametrize("kind", ["normal", "multiallelic", "spanning", "preannotated"])
def test_indexed_miss_spool(native, tmp_path, monkeypatch, hits, kind):
    files = dataset(tmp_path / "data", hit_every=hits, kind=kind)
    before, old = run_case(tmp_path / "before", files, "fused", preserve=True)
    assert old.last_run_stats["experimental_indexed_misses"] is False
    monkeypatch.setenv("VCFCACHE_EXPERIMENT_INDEXED_MISSES", "1")
    after, new = run_case(tmp_path / "after", files, "fused", preserve=True)
    assert strict_semantic_compare(after, before)["semantic_pass"]
    assert new.last_run_stats["experimental_indexed_misses"] is True
    assert list(new.work_dir.rglob("misses.raw.bcf.csi"))
    assert not list((new.work_dir / "fused").glob("shard-*/misses.bcf"))


@pytest.mark.parametrize("window", [1024, 4096])
@pytest.mark.parametrize("workers", [1, 3])
def test_native_windows(native, tmp_path, monkeypatch, window, workers):
    files = dataset(tmp_path / "data", records=200)
    before, old = run_case(tmp_path / "before", files, "native", workers=workers)
    assert old.last_run_stats["experimental_window_bp"] == 16384
    monkeypatch.setenv("VCFCACHE_EXPERIMENT_WINDOW_BP", str(window))
    after, new = run_case(tmp_path / "after", files, "native", workers=workers)
    assert strict_semantic_compare(after, before)["semantic_pass"]
    assert new.last_run_stats["experimental_window_bp"] == window


def test_screening_is_cached_only_and_bounded():
    from benchmarks.engine_profile.matrix import design as original
    from benchmarks.performance_sprint.pilot import design

    groups = design(dict(groups=original()))
    cells = [c for g in groups for c in g["cells"]]
    assert len(cells) == 14
    assert all(c["role"] == "measure" for c in cells)
    assert {c["purpose"] for c in cells} == {
        "control-adaptive",
        "control-native",
        "indexed-misses",
        "window-1024",
        "window-4096",
    }
    assert all(c["cpus"] == 16 for c in cells)
    assert {g["sample"] for g in groups} == {"LP6005519-DNA_A06", "LP6005592-DNA_C05"}
    second = design(dict(groups=original()), "trusted-hits")
    assert sum(len(g["cells"]) for g in second) == 8
    third = design(dict(groups=original()), "memory-preflight")
    assert sum(len(g["cells"]) for g in third) == 4
    assert all(c["memory_mb"] == 8192 for g in third for c in g["cells"])


@pytest.mark.parametrize("hits", [0, 1, 2, 10])
@pytest.mark.parametrize("kind", ["normal", "multiallelic", "spanning", "preannotated"])
def test_trusted_routing_avoids_only_redundant_hit_inspection(
    native, tmp_path, monkeypatch, hits, kind
):
    import json

    files = dataset(tmp_path / "data", hit_every=hits, kind=kind)
    before, old = run_case(tmp_path / "before", files, "fused", preserve=True)
    monkeypatch.setenv("VCFCACHE_EXPERIMENT_TRUST_ROUTED_HITS", "1")
    after, new = run_case(tmp_path / "after", files, "fused", preserve=True)
    assert strict_semantic_compare(after, before)["semantic_pass"]
    assert new.last_run_stats["experimental_trust_routed_hits"] is True
    for key in [
        "input_variants",
        "missing_variants",
        "output_variants",
        "missing_annotated",
    ]:
        assert old.last_run_stats[key] == new.last_run_stats[key]
    trusted = sum(
        json.loads(p.read_text())["trusted_hit_records"]
        for p in new.work_dir.glob("fused/shard-*/assemble.json")
    )
    assert (
        trusted
        == new.last_run_stats["input_variants"] - new.last_run_stats["missing_variants"]
    )


def test_payload_prediction_is_bounded_and_abstains(tmp_path):
    from vcfcache.database.memory_preflight import predict_buffer

    sample, cache, *_ = dataset(tmp_path / "data", records=100, payload=2000)
    prediction = predict_buffer(sample, cache, 5_000_000, max_seconds=5)
    assert prediction["sampled_records"] == 16
    assert prediction["exact_cache_matches"] >= 8
    assert prediction["estimated_buffer_bytes"] > 5_000_000 * 1024
    assert predict_buffer(sample, cache, 5_000_000, max_seconds=0) is None
    assert predict_buffer(sample, tmp_path / "absent", 5_000_000) is None


@pytest.mark.parametrize("hits", [0, 1, 2, 10])
@pytest.mark.parametrize("kind", ["normal", "multiallelic", "spanning", "preannotated"])
def test_bounded_native_route(native, tmp_path, monkeypatch, hits, kind):
    files = dataset(tmp_path / "data", hit_every=hits, kind=kind)
    before, _ = run_case(tmp_path / "before", files, "fused", preserve=True)
    monkeypatch.setenv("VCFCACHE_EXPERIMENT_NATIVE_ROUTE", "1")
    after, new = run_case(tmp_path / "after", files, "fused", preserve=True)
    assert strict_semantic_compare(after, before)["semantic_pass"]
    assert new.last_run_stats.get("experimental_native_route") is True


def test_bounded_native_route_chunk_boundaries(native, tmp_path, monkeypatch):
    files = dataset(tmp_path / "data", records=40000, payload=100)
    before, _ = run_case(tmp_path / "before", files, "fused", preserve=True)
    monkeypatch.setenv("VCFCACHE_EXPERIMENT_NATIVE_ROUTE", "1")
    after, new = run_case(tmp_path / "after", files, "fused", preserve=True)
    assert strict_semantic_compare(after, before)["semantic_pass"]
    assert new.last_run_stats["input_variants"] == 40000
