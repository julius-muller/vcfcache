"""Regression checks for the bounded confirmation pilot, not annotation engines."""

import io
import json
from collections import Counter
from pathlib import Path
from types import SimpleNamespace

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


def test_candidate_matrix_is_bounded_and_has_no_native_wgs():
    total = 0
    for _, _, _, assay, _ in run.CANDIDATE_SELECTION:
        modes = run.candidate_modes(assay)
        total += 3 * len(modes)
        assert modes[0] == "direct"
        assert all(run.CANDIDATES[m]["prepared"] for m in modes)
        for rep in range(1, 4):
            assert sorted(run.order_for(rep, 0, modes)) == sorted(modes)
        if assay == "wgs":
            assert all(run.CANDIDATES[m]["engine"] != "native" for m in modes)
    assert total == 96


@pytest.mark.parametrize("mode", list(run.CANDIDATES))
def test_explicit_candidate_configuration(mode):
    spec = run.CANDIDATES[mode]
    data = dict(
        configurations=run.CANDIDATES,
        backends={b: dict(path="/frozen/" + b) for b in ("python", "rust")},
        tools=dict(helper=dict(path="/frozen/helper")),
    )
    command = run.command_for(
        mode, *map(Path, ["cache", "in", "out", "params", "stats"]), spec
    )
    assert command[command.index("--engine") + 1] == spec["engine"]
    env = run.environment_for(data, mode, Path("/tools/bcftools"))
    assert env["VCFCACHE_PILOT_PREPARED"] == "1"
    assert env["PYTHONPATH"] == "/frozen/" + spec["backend"]
    if spec["engine"] == "native":
        assert "--engine-memory-mb" in command
        assert int(env["VCFCACHE_EXPERIMENT_CODEC_THREADS"]) + spec["workers"] <= 8
    if spec.get("coalesced"):
        assert env["VCFCACHE_EXPERIMENT_COALESCE_CHUNKS"] == "1"


@pytest.mark.parametrize(
    "records,raises",
    [
        ("chr1\t1\tA\tG\t.\n", False),
        ("chr1\t1\tA\tG,C\t.\n", True),
        ("chr1\t1\tA\t*\t.\n", True),
        ("chr1\t1\tA\tG\tannotation\n", True),
        ("chr1\t1\tA\tG\t.\nchr1\t1\tA\tG\t.\n", True),
    ],
)
def test_prepared_domain_is_checked_not_transformed(monkeypatch, records, raises):
    monkeypatch.setattr(
        run.subprocess,
        "check_output",
        lambda *a, **k: "##INFO=<ID=CSQ,Number=.,Type=String>\n",
    )
    monkeypatch.setattr(
        run.subprocess,
        "Popen",
        lambda *a, **k: SimpleNamespace(
            stdout=io.StringIO(records), wait=lambda: 0, poll=lambda: 0
        ),
    )
    if raises:
        with pytest.raises(ValueError):
            run.preflight_prepared("bcftools", Path("in.bcf"), "CSQ")
    else:
        assert run.preflight_prepared("bcftools", Path("in.bcf"), "CSQ")["records"] == 1


def test_candidate_review_keeps_all_configurations(tmp_path):
    root = tmp_path / "root"
    case = dict(id="veprs-wes-demo", node="sl-w1", modes=run.candidate_modes("wes"))
    freeze_minimal(root, [case])
    doc = json.loads((root / "manifest.json").read_text())
    doc["profile"] = "candidates"
    write_json_atomic(root / "manifest.json", doc)
    (root / "READY").write_text(sha256sum(root / "manifest.json") + "\n")
    folder = root / "cases" / case["id"]
    folder.mkdir(parents=True)
    rows = [
        dict(cell="r01-" + m, mode=m, repetition=1, seconds=10.0 + i)
        for i, m in enumerate(case["modes"])
    ]
    write_json_atomic(
        folder / "state.json", dict(state="running", rows=rows, audits={})
    )
    run.collect(root)
    text = (root / "RESULTS.md").read_text()
    assert all(m in text for m in case["modes"])
    assert json.loads((root / "REVIEW.json").read_text())["accepted_cases"] == 0
