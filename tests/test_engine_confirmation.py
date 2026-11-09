"""Regression checks for the bounded confirmation pilot, not annotation engines."""

import json
from collections import Counter
from pathlib import Path

import pytest

from benchmarks.engine_confirmation import run
from benchmarks.run_pilot import sha256sum, write_json_atomic


def freeze_minimal(root, cases):
    root.mkdir()
    write_json_atomic(root / "manifest.json", dict(cases=cases))
    (root / "READY").write_text(sha256sum(root / "manifest.json") + "\n")


def test_balanced_positions_and_distinct_cells():
    for offset in range(3):
        positions = Counter(
            (mode, position)
            for rep in range(1, 4)
            for position, mode in enumerate(run.order_for(rep, offset))
        )
        assert len(positions) == 9
        assert set(positions.values()) == {1}
    assert len(run.SELECTION) == 6
    assert {c[-1] for c in run.SELECTION} == {f"sl-w{i}" for i in range(1, 7)}


@pytest.mark.parametrize("mode", run.MODES)
def test_light_statistics_and_correct_direct_engine(mode):
    command = run.command_for(
        mode, *map(Path, ["cache", "input", "output", "params", "stats"])
    )
    assert command[command.index("--statistics") + 1] == "light"
    assert ("--uncached" in command) == (mode == "direct")
    assert command[command.index("--engine") + 1] == (
        "legacy" if mode == "direct" else "adaptive"
    )


def test_backend_environment_isolated(monkeypatch):
    monkeypatch.setenv("VCFCACHE_EXPERIMENT_NATIVE_ROUTE", "1")
    monkeypatch.setenv("VCFCACHE_PILOT_PREPARED", "1")
    monkeypatch.setenv("VCFCACHE_PROFILE_DIR", "/unwanted")
    data = dict(
        backends={m: dict(path="/frozen/" + m) for m in ("python", "rust")},
        tools=dict(helper=dict(path="/frozen/helper")),
    )
    for mode in run.MODES:
        env = run.environment_for(data, mode, Path("/tools/bcftools"))
        assert not any(
            k.startswith(
                ("VCFCACHE_EXPERIMENT_", "VCFCACHE_PILOT_", "VCFCACHE_PROFILE")
            )
            for k in env
        )
        assert env["PYTHONPATH"] == "/frozen/" + (
            "rust" if mode == "rust" else "python"
        )
        assert ("VCFCACHE_ENGINE" in env) == (mode == "rust")


def test_definition_rejects_modified_manifest(tmp_path):
    root = tmp_path / "root"
    freeze_minimal(root, [])
    assert run.definition(root) == dict(cases=[])
    (root / "manifest.json").write_text('{"changed":true}')
    with pytest.raises(RuntimeError, match="changed"):
        run.definition(root)


def test_terminal_review_records_missing_cases_without_acceptance(tmp_path):
    root = tmp_path / "root"
    freeze_minimal(root, [dict(id="unstarted", node="sl-w1")])
    run.collect(root, terminal=True)
    review = json.loads((root / "REVIEW.json").read_text())
    assert review["accepted_cases"] == 0
    assert review["cases"][0]["state"]["state"] == "interrupted"
    assert (root / "terminal-evidence.tar.gz").exists()


def test_submission_refuses_duplicates_without_sbatch(tmp_path, monkeypatch):
    root = tmp_path / "root"
    freeze_minimal(root, [])
    (root / "submission_intent.json").write_text("{}")
    monkeypatch.setattr(
        run.subprocess, "check_output", lambda *a, **k: pytest.fail("must not resubmit")
    )
    with pytest.raises(RuntimeError, match="already attempted"):
        run.submit(root)
