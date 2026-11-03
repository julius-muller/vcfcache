"""Decoder-pool overrides retain complete native output semantics and CPU bounds."""

import pytest
import yaml

from benchmarks.run_pilot import strict_semantic_compare
from tests.test_fused_engine import dataset, run_case
from tests.test_fused_engine import native as native_fixture

native = native_fixture


@pytest.mark.parametrize("hits", [0, 1, 2])
@pytest.mark.parametrize("codec_threads", [0, 4])
def test_decoder_pool_equivalence(native, tmp_path, monkeypatch, hits, codec_threads):
    files = dataset(tmp_path / "data", records=80, hit_every=hits)
    params = yaml.safe_load(files[2].read_text())
    params["threads"] = 8
    files[2].write_text(yaml.safe_dump(params))
    direct, _ = run_case(tmp_path, files, "legacy")
    monkeypatch.setenv("VCFCACHE_EXPERIMENT_CODEC_THREADS", str(codec_threads))
    candidate, manager = run_case(tmp_path, files, "native", workers=4)
    assert manager.last_run_stats["lookup_codec_threads"] == codec_threads
    assert manager.last_run_stats["effective_engine"] == "native"
    assert strict_semantic_compare(candidate, direct)["semantic_pass"]


@pytest.mark.parametrize("codec_threads", ["-1", "9", "invalid"])
def test_invalid_decoder_pool(native, tmp_path, monkeypatch, codec_threads):
    files = dataset(tmp_path / "data", records=10)
    monkeypatch.setenv("VCFCACHE_EXPERIMENT_CODEC_THREADS", codec_threads)
    with pytest.raises(ValueError):
        run_case(tmp_path, files, "native", workers=4)
