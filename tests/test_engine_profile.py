"""Profiling must preserve records and leave durable, honest partial reviews."""

import argparse
import json
from collections import Counter

import pytest

from benchmarks.engine_profile import review, run
from benchmarks.engine_profile.matrix import design
from benchmarks.run_pilot import sha256sum, strict_semantic_compare, write_json_atomic
from tests.test_fused_engine import dataset, run_case
from tests.test_fused_engine import native as native_fixture
from vcfcache.database.profiling import phase

native = native_fixture


def test_design_bound():
    groups = design()
    cells = [c for g in groups for c in g["cells"]]
    assert len(cells) == len({c["id"] for c in cells}) == 86
    assert Counter(c["purpose"] for c in cells)["primary"] == 24
    assert Counter(c["purpose"] for c in cells)["reference"] == 24
    assert max(c["memory_mb"] for c in cells) == 8192
    assert set(c["cpus"] for c in cells) == {1, 8, 16}
    assert all(c["recipe"] != "veprs-clinical" for c in cells)


def test_union_overlapping_nested_and_disjoint():
    assert review.union_seconds([(0, 10), (2, 3), (5, 12), (14, 20)]) == 18e-9
    with pytest.raises(ValueError):
        review.union_seconds([(5, 4)])


def test_python_trace_disabled_and_exception(tmp_path, monkeypatch):
    monkeypatch.delenv("VCFCACHE_PROFILE_DIR", raising=False)
    with phase("off"):
        pass
    assert not list(tmp_path.iterdir())
    monkeypatch.setenv("VCFCACHE_PROFILE_DIR", str(tmp_path))
    with pytest.raises(ValueError), phase("failed"):
        raise ValueError("expected")
    record = json.loads(next(tmp_path.glob("*.jsonl")).read_text())
    assert record["outcome"] == "ValueError"
    assert record["duration_ns"] >= 0


@pytest.mark.parametrize("engine", ["native", "adaptive"])
def test_profile_record_equivalence(native, tmp_path, monkeypatch, engine):
    files = dataset(tmp_path / "data")
    off, _ = run_case(tmp_path / "off", files, engine)
    traces = tmp_path / "trace"
    traces.mkdir()
    monkeypatch.setenv("VCFCACHE_PROFILE_DIR", str(traces))
    on, _ = run_case(tmp_path / "on", files, engine)
    assert strict_semantic_compare(on, off)["semantic_pass"]
    events = [
        json.loads(line)
        for path in traces.glob("*.jsonl")
        for line in path.read_text().splitlines()
    ]
    assert any(e["name"] == f"{engine}.attempt" for e in events)
    assert any(p.name.startswith("rust-") for p in traces.iterdir())
    assert all(e["duration_ns"] >= 0 for e in events)


def campaign(tmp_path):
    groups = design()[:1]
    groups[0]["cells"] = groups[0]["cells"][:1]
    manifest = dict(
        groups=groups,
        root=str(tmp_path),
        builds={"profile": {"root": "/results/profile"}},
    )
    write_json_atomic(tmp_path / "manifest.json", manifest)
    (tmp_path / "READY").write_text(sha256sum(tmp_path / "manifest.json"))
    return groups[0]["cells"][0]


def test_native_abandoned_attempt_is_traced(native, tmp_path, monkeypatch):
    from vcfcache.database.workflow_manager import WorkflowManager

    files = dataset(tmp_path / "data", records=200, payload=20000)
    traces = tmp_path / "trace"
    traces.mkdir()
    monkeypatch.setenv("VCFCACHE_PROFILE_DIR", str(traces))
    manager = WorkflowManager(
        input_file=files[0],
        output_dir=tmp_path / "run",
        name="fallback",
        params_file=files[2],
        anno_config_file=files[3],
        output_file=tmp_path / "output.bcf",
    )
    manager.run("annotate", db_bcf=files[1], engine="native", engine_memory_mb=1)
    names = {
        json.loads(line)["name"]
        for path in traces.glob("*.jsonl")
        for line in path.read_text().splitlines()
    }
    assert {"native.attempt", "native.execute", "adaptive.attempt"} <= names
    assert "memory" in manager.last_run_stats["fallback_reason"]


def test_partial_collection_and_idempotence(tmp_path, monkeypatch):
    cell = campaign(tmp_path)
    monkeypatch.setattr(
        review, "job_states", lambda _: {"12": dict(state="TIMEOUT", group="g00")}
    )
    monkeypatch.setattr(review.shutil, "which", lambda _: None)
    write_json_atomic(
        tmp_path / "submissions.json",
        [dict(jobs=[dict(job_id="12", group="g00", cells=[cell["id"]])])],
    )
    assert review.collect(tmp_path, write=False)["counts"] == {"incomplete_timeout": 1}
    assert not (tmp_path / "reports").exists()
    review.collect(tmp_path)
    review.collect(tmp_path)
    assert len(list((tmp_path / "reports").iterdir())) == 1
    assert (tmp_path / "scheduler-status.json").exists()


def test_partial_submission_has_afterany_collector(tmp_path, monkeypatch):
    campaign(tmp_path)
    manifest = run.read(tmp_path / "manifest.json")
    manifest["groups"] += design()[1:2]
    write_json_atomic(tmp_path / "manifest.json", manifest)
    (tmp_path / "READY").write_text(sha256sum(tmp_path / "manifest.json"))
    calls = []

    def submit(command, **kwargs):
        calls.append(command)
        if len(calls) == 2:
            raise RuntimeError("scheduler failure")
        return "123\n" if len(calls) == 1 else "125\n"

    monkeypatch.setattr(run.subprocess, "check_output", submit)
    with pytest.raises(RuntimeError, match="scheduler failure"):
        run.submit(argparse.Namespace(campaign=tmp_path, command="submit", cells=None))
    batch = run.read(tmp_path / "submissions.json")[0]
    assert batch["jobs"][0]["job_id"] == "123"
    assert batch["collector_job_id"] == "125"
    assert "--dependency=afterany:123" in calls[-1]
    assert batch["submission_error"] == "scheduler failure"


def test_manifest_drift_rejected(tmp_path):
    campaign(tmp_path)
    (tmp_path / "manifest.json").write_text("{}")
    with pytest.raises(RuntimeError, match="modified"):
        run.verify(tmp_path)
