"""Shared optional-engine selection and output publication contracts."""

import os
import subprocess
from pathlib import Path

import pytest

from benchmarks.run_pilot import strict_semantic_compare
from tests.test_fused_engine import dataset, run_case
from tests.test_fused_engine import native as native_fixture
from vcfcache.database import engine_runtime as runtime

native = native_fixture


@pytest.mark.parametrize("engine", ["legacy", "fused", "adaptive", "native"])
def test_explicit_selection_never_changes_request(engine):
    assert runtime.select_engine(engine) == (engine, None)


def test_auto_missing_helper_is_legacy(monkeypatch):
    monkeypatch.delenv("VCFCACHE_ENGINE", raising=False)

    def absent():
        raise FileNotFoundError("not installed")

    monkeypatch.setattr(runtime, "native_binary", absent)
    engine, reason = runtime.select_engine("auto")
    assert engine == "legacy"
    assert "not installed" in reason


def test_auto_respects_invalid_explicit_helper(monkeypatch, tmp_path):
    monkeypatch.setenv("VCFCACHE_ENGINE", str(tmp_path / "absent"))
    with pytest.raises(FileNotFoundError, match="VCFCACHE_ENGINE"):
        runtime.select_engine("auto")


def test_auto_invalid_override_is_not_hidden_by_platform(monkeypatch, tmp_path):
    monkeypatch.delattr(os, "sched_getaffinity")
    monkeypatch.setenv("VCFCACHE_ENGINE", str(tmp_path / "absent"))
    with pytest.raises(FileNotFoundError, match="VCFCACHE_ENGINE"):
        runtime.select_engine("auto")


def test_auto_checks_platform_before_optional_helper(monkeypatch):
    monkeypatch.delattr(os, "sched_getaffinity")
    assert runtime.select_engine("auto")[0] == "legacy"
    with pytest.raises(RuntimeError, match="Linux"):
        with runtime.cpu_affinity(1):
            pytest.fail("unsupported platform entered affinity context")


@pytest.mark.parametrize("version", ["vcfcache-engine 0.1.0", "unexpected 99.0.0"])
def test_auto_broken_or_stale_helper_does_not_silently_fallback(monkeypatch, version):
    monkeypatch.setattr(runtime, "native_binary", lambda: "helper")
    monkeypatch.setattr(subprocess, "check_output", lambda *a, **kw: version)
    with pytest.raises(RuntimeError, match="rebuild"):
        runtime.select_engine("auto")


def test_auto_selects_adaptive_not_native(monkeypatch):
    monkeypatch.setattr(runtime, "native_binary", lambda: "helper")
    monkeypatch.setattr(
        subprocess, "check_output", lambda *a, **kw: "vcfcache-engine 0.4.0"
    )
    assert runtime.select_engine("auto")[0] == "adaptive"


@pytest.mark.parametrize("outcome", ["success", "abort", "error", "missing-index"])
def test_staging_publishes_only_completed_outputs(tmp_path, outcome):
    target = tmp_path / "result.bcf"
    index = Path(str(target) + ".csi")
    target.write_bytes(b"original")
    index.write_bytes(b"original-index")

    def exercise():
        with runtime.staged_output(target, ".test-") as stage:
            stage.path.write_bytes(b"new")
            if outcome != "missing-index":
                Path(str(stage.path) + ".csi").write_bytes(b"new-index")
            if outcome == "error":
                raise ValueError("cancelled")
            if outcome == "abort":
                stage.aborted = True

    if outcome in {"error", "missing-index"}:
        with pytest.raises((ValueError, RuntimeError)):
            exercise()
    else:
        exercise()
    assert target.read_bytes() == (b"new" if outcome == "success" else b"original")
    assert index.read_bytes() == (
        b"new-index" if outcome == "success" else b"original-index"
    )
    assert not list(tmp_path.glob(".test-*"))


@pytest.mark.parametrize("engine", ["adaptive", "auto"])
def test_public_alias_matches_legacy_complete_records(native, tmp_path, engine):
    files = dataset(tmp_path / "data")
    reference, _ = run_case(tmp_path, files, "legacy")
    result, manager = run_case(tmp_path, files, engine)
    assert strict_semantic_compare(result, reference)["semantic_pass"]
    assert manager.last_run_stats["requested_engine"] == engine
    assert manager.last_run_stats["effective_engine"] == "adaptive"


def test_auto_legacy_fallback_is_reported(native, tmp_path, monkeypatch):
    files = dataset(tmp_path / "data")
    monkeypatch.setattr(
        runtime, "select_engine", lambda requested: ("legacy", "test: unavailable")
    )
    _, manager = run_case(tmp_path, files, "auto")
    assert manager.last_run_stats["requested_engine"] == "auto"
    assert manager.last_run_stats["effective_engine"] == "legacy"
    assert manager.last_run_stats["selection_reason"] == "test: unavailable"
