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


@pytest.mark.parametrize(
    "marker,expected",
    [(None, "unavailable"), ("", "GROUP_FINISHED"), ("-failure", "GROUP_FAILED")],
)
def test_expired_running_snapshot(tmp_path, marker, expected):
    write_json_atomic(
        tmp_path / "submissions.json",
        [dict(attempt=1, jobs=[dict(job_id="7", group="g02", cells=[])])],
    )
    if marker is not None:
        write_json_atomic(
            tmp_path / f"group-g02-a1{marker}.json", dict(finished_unix=123)
        )
    result = review.resolve_scheduler(
        tmp_path,
        {"7": dict(group="g02", state="unavailable")},
        {"7": dict(group="g02", state="RUNNING")},
    )
    assert result["7"]["state"] == expected
    assert "exit_code" not in result["7"]
    assert (
        review.resolve_scheduler(
            tmp_path, {"7": dict(group="g02", state="RUNNING")}, {}
        )["7"]["state"]
        == "RUNNING"
    )


def test_terminal_snapshot_retained(tmp_path):
    write_json_atomic(
        tmp_path / "submissions.json",
        [dict(jobs=[dict(job_id="7", group="g02", cells=[])])],
    )
    assert (
        review.resolve_scheduler(tmp_path, {}, {"7": dict(state="COMPLETED")})["7"][
            "state"
        ]
        == "COMPLETED"
    )


def test_vep_recovery_controls():
    manifest = dict(
        recipes={"vep": dict(params={"vep_forks": 16, "vep_buffer": 500000})}
    )
    cell = dict(
        recipe="vep",
        cpus=16,
        annotation_overrides={"vep_forks": 4, "vep_buffer": 50000},
    )
    assert run.cell_params(manifest, cell) == dict(
        threads=16, vep_forks=4, vep_buffer=50000
    )
    cell["annotation_overrides"]["vep_forks"] = 17
    with pytest.raises(ValueError):
        run.cell_params(manifest, cell)


def test_recovery_only_unsuccessful_wgs(tmp_path):
    from benchmarks.engine_profile.recover_vep import create

    parent = tmp_path / "parent"
    parent.mkdir()
    manifest = dict(groups=design(), root=str(parent))
    write_json_atomic(parent / "manifest.json", manifest)
    (parent / "READY").write_text(sha256sum(parent / "manifest.json"))
    original = (parent / "manifest.json").read_bytes()
    for group in manifest["groups"]:
        for cell in group["cells"]:
            if cell["recipe"] == "vep" and cell["assay"] == "wgs":
                file = parent / "attempts/a001" / cell["id"] / "state.json"
                file.parent.mkdir(parents=True)
                write_json_atomic(
                    file,
                    dict(state="failed" if cell["role"] == "direct" else "skipped"),
                )
    harness = tmp_path / "harness"
    script = harness / "benchmarks/engine_profile/job.sh"
    script.parent.mkdir(parents=True)
    script.write_text("test")
    result = create(parent, tmp_path / "retry", str(harness), "test-commit")
    assert [g["id"] for g in result["groups"]] == ["g02", "g05"]
    assert sum(len(g["cells"]) for g in result["groups"]) == 8
    assert (parent / "manifest.json").read_bytes() == original
    file = parent / "attempts/a001/g02-c004/state.json"
    write_json_atomic(file, dict(state="complete"))
    with pytest.raises(ValueError, match="unsuccessful"):
        create(parent, tmp_path / "invalid", str(harness), "test-commit")


def test_consolidation_keeps_configuration_and_failed_evidence(tmp_path):
    from benchmarks.engine_profile.consolidate import consolidate

    group = design()[2]
    group["cells"] = group["cells"][4:8]
    for campaign_name in ("original", "recovery"):
        root = tmp_path / "source" / campaign_name
        root.mkdir(parents=True)
        manifest = dict(
            groups=[group], recipes={"vep": dict(params={"vep_buffer": 500000})}
        )
        if campaign_name == "recovery":
            manifest["recovery"] = dict(
                parent_manifest_sha256=sha256sum(
                    tmp_path / "source/original/manifest.json"
                )
            )
        write_json_atomic(root / "manifest.json", manifest)
        (root / "READY").write_text(sha256sum(root / "manifest.json"))
        for cell in group["cells"]:
            path = root / "attempts/a001" / cell["id"]
            path.mkdir(parents=True)
            state = dict(state="failed")
            if campaign_name == "recovery":
                state.update(
                    state="complete",
                    effective_engine="adaptive",
                    direct_pass=True,
                    metrics=dict(
                        wall_seconds=10 if cell["role"] == "direct" else 5,
                        cpu_seconds=8,
                    ),
                )
                write_json_atomic(
                    path / "statistics.json",
                    dict(
                        variant_counts=dict(input_variants=10, missing_variants=2),
                        stage_timings=[dict(stage="annotate-misses", wall_seconds=2)],
                    ),
                )
            write_json_atomic(path / "state.json", state)
    rows = consolidate(tmp_path)
    assert len(rows) == 8
    assert sum(r["state"] == "failed" for r in rows) == 4
    assert all(
        r.get("speedup") == 2
        for r in rows
        if r["campaign"] == "recovery" and r["role"] != "direct"
    )
    manifest_path = tmp_path / "source/recovery/manifest.json"
    manifest = run.read(manifest_path)
    manifest["recovery"]["parent_manifest_sha256"] = "wrong"
    write_json_atomic(manifest_path, manifest)
    (manifest_path.parent / "READY").write_text(sha256sum(manifest_path))
    with pytest.raises(ValueError, match="does not belong"):
        consolidate(tmp_path)
