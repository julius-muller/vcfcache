"""Experimental CSI-window coalescing must not change any record fields."""

import pytest
import yaml

from benchmarks.run_pilot import strict_semantic_compare
from tests.test_fused_engine import dataset, run_case
from tests.test_fused_engine import native as native_fixture

native = native_fixture


@pytest.mark.parametrize("hits", [0, 1, 2, 10])
@pytest.mark.parametrize("workers", [1, 4, 8])
def test_coalesced_complete_records(native, tmp_path, monkeypatch, hits, workers):
    files = dataset(tmp_path / "data", records=80, hit_every=hits)
    params = yaml.safe_load(files[2].read_text())
    params["threads"] = 8
    files[2].write_text(yaml.safe_dump(params))
    legacy, _ = run_case(tmp_path, files, "legacy")
    monkeypatch.setenv("VCFCACHE_EXPERIMENT_COALESCE_CHUNKS", "1")
    candidate, manager = run_case(tmp_path, files, "native", workers=workers)
    assert manager.last_run_stats["coalesced_chunks"] is True
    assert manager.last_run_stats["effective_engine"] == "native"
    assert strict_semantic_compare(candidate, legacy)["semantic_pass"]
