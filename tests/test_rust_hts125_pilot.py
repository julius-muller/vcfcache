"""Prepared-input variants must preserve complete records on their verified domain."""

import json

import pytest

from benchmarks.run_pilot import sha256sum, strict_semantic_compare, write_json_atomic
from benchmarks.rust_hts125_pilot.common import preflight_prepared
from benchmarks.rust_hts125_pilot.run import CONDITIONS, INTEL, collect, modes_for
from tests.test_fused_engine import checked, dataset, run_case
from tests.test_fused_engine import native as native_fixture

native = native_fixture


@pytest.mark.parametrize("engine", ["adaptive", "native"])
@pytest.mark.parametrize("hits", [0, 1, 2])
def test_prepared_equivalence(native, tmp_path, monkeypatch, engine, hits):
    """No normalization subprocess is needed for already-biallelic input."""
    files = dataset(tmp_path / "data", hit_every=hits)
    before = sha256sum(files[0])
    assert preflight_prepared("bcftools", files[0], "CSQ")["records"] == 20
    expected, _ = run_case(tmp_path / "control", files, "legacy", preserve=True)
    monkeypatch.setenv("VCFCACHE_PILOT_PREPARED", "1")
    output, manager = run_case(tmp_path / "candidate", files, engine, preserve=True)
    assert strict_semantic_compare(output, expected)["semantic_pass"]
    assert manager.last_run_stats["pilot_prepared"]
    assert sha256sum(files[0]) == before
    assert not list(manager.work_dir.rglob("input.bcf.pipe"))
    if engine == "native":
        spec = json.loads(next(manager.work_dir.rglob("job.json")).read_text())
        assert spec["input"] == str(files[0])


def test_linked_build_info(native):
    """Versions come from the helper's linked library, not external tools."""
    info = json.loads(checked(native, "build-info").stdout)
    assert info["htslib_version"] in {"1.19.1", "1.25"}
    assert info["libdeflate"]


@pytest.mark.parametrize("kind", ["multiallelic", "spanning", "preannotated"])
def test_prepared_preflight_rejects_unsupported(tmp_path, kind):
    """Unsupported input is rejected before timing, not silently rewritten."""
    files = dataset(tmp_path / "data", kind=kind)
    with pytest.raises(ValueError):
        preflight_prepared("bcftools", files[0], "CSQ")


def test_pilot_is_bounded():
    """No hidden resident WGS fallbacks, full cohort or AMD worker."""
    assert len(CONDITIONS) == 9
    assert 4 * len(modes_for("wes")) + 2 * len(modes_for("wgs")) == 46
    assert not any(m.startswith("resident-") for m in modes_for("wgs"))
    assert "sl-w5" not in INTEL


@pytest.mark.parametrize("terminal", [False, True])
def test_terminal_review_preserves_incomplete_evidence(tmp_path, terminal):
    """A scheduler timeout must not become an eternal running success claim."""
    write_json_atomic(
        tmp_path / "manifest.json", dict(cases=[dict(id="case", node="sl-w1")])
    )
    (tmp_path / "READY").write_text(sha256sum(tmp_path / "manifest.json") + "\n")
    collect(tmp_path, terminal=terminal)
    review = json.loads((tmp_path / "REVIEW.json").read_text())
    assert review["complete"] == 0
    assert review["cases"][0]["state"] == ("interrupted" if terminal else "not_started")
    assert (tmp_path / "terminal-evidence.tar.gz").exists() == terminal
