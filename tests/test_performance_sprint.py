"""Frozen sprint bookkeeping must not silently expand the engineering budget."""

import json
from pathlib import Path

import pytest

from benchmarks.performance_sprint.deploy import check_budget, reserved_cells
from benchmarks.performance_sprint.references import audit, request
from benchmarks.run_pilot import REPO_ROOT, sha256sum


def test_budget_counts_uncertain_submissions(tmp_path):
    for index, screen in enumerate(["initial", "trusted-hits", "memory-preflight"]):
        (tmp_path / f"launch-{index}.json").write_text(
            json.dumps(dict(campaign=f"/results/{index}", request=dict(screen=screen)))
        )
    assert reserved_cells(tmp_path) == 26
    assert check_budget(tmp_path, "native-route") == 34
    (tmp_path / "launch-last.json").write_text(
        json.dumps(dict(campaign="/results/last", request=dict(screen="native-route")))
    )
    with pytest.raises(ValueError, match="budget exceeded"):
        check_budget(tmp_path, "memory-preflight")
    assert check_budget(tmp_path, "memory-confirm") == 36


def test_first_receipt_without_screen_uses_original_design(tmp_path):
    (tmp_path / "launch-old.json").write_text(
        json.dumps(dict(campaign="/results/first", request=dict(profile={})))
    )
    assert reserved_cells(tmp_path) == 14


def test_duplicate_receipts_are_not_silently_counted(tmp_path):
    value = json.dumps(dict(campaign="/results/one", request=dict(screen="initial")))
    for name in ["launch-one.json", "launch-two.json"]:
        (tmp_path / name).write_text(value)
    with pytest.raises(ValueError, match="duplicate campaign"):
        reserved_cells(tmp_path)


def test_all_manuscript_families_have_explicit_reference_resolution():
    document = request(REPO_ROOT)
    assert len(document["rows"]) == 540
    limited = [
        r
        for r in document["rows"]
        if r["reference_kind"] == "historical_annotation_key_digest_only"
    ]
    assert len(limited) == 180
    assert all(r["full_direct_output_retained"] is False for r in limited)
    assert all(r["cpu_limit"] == 8 for r in document["rows"])
    reused = [
        r
        for r in document["rows"]
        if r["family"] == "synthetic_pipeline_spectrum"
        and r["strategy"] == "AF ≥ 1% cache"
    ]
    assert len(reused) == 6
    assert all("wgs-pipeline-spectrum-4f854e5" in r["direct_root"] for r in reused)
    assert all(
        r["reference_kind"] == "historical_semantic_metadata_only" for r in reused
    )
    snapshot_path = "benchmarks/figures/source_data/final/2026-08-09T223718Z-light-matched-final/SNAPSHOT.json"
    snapshot = json.loads((REPO_ROOT / snapshot_path).read_text())
    fast = [
        r
        for r in document["rows"]
        if r["family"] == "full52_wgs_strategies" and r["tool"] == "fastvep"
    ]
    assert len(fast) == 156
    assert all(snapshot["fastvep"]["campaign_id"] in r["direct_root"] for r in fast)
    assert snapshot_path in document["sources"]


def test_historical_metadata_is_not_full_output_validation(tmp_path):
    from benchmarks.performance_sprint.references import historical_metadata
    from benchmarks.run_pilot import write_json_atomic

    baseline = tmp_path / "baseline.json"
    comparison = tmp_path / "comparison.json"
    write_json_atomic(
        baseline,
        dict(
            semantic_pass=True,
            uncached_wall_seconds=100,
            cached_wall_seconds=20,
            input_records=1000,
        ),
    )
    write_json_atomic(comparison, dict(semantic_pass=True))
    row = dict(
        baseline_summary=dict(path=str(baseline), sha256=sha256sum(baseline)),
        cached_summary=str(baseline),
        baseline_comparison=str(comparison),
        historical_direct_seconds=100,
        historical_cached_seconds=20,
        input_records=1000,
    )
    result = historical_metadata(row)
    assert result["state"] == "historical_semantic_metadata_only"
    assert result["full_output_validation"] is False
    row["historical_direct_seconds"] = 101
    with pytest.raises(RuntimeError, match="timing does not match"):
        historical_metadata(row)


def test_supplemental_audit_does_not_repeat_resolved_rows(tmp_path):
    from benchmarks.run_pilot import write_json_atomic

    request_file = tmp_path / "request.json"
    write_json_atomic(
        request_file,
        dict(
            rows=[
                dict(
                    id=name,
                    family="fixture",
                    reference_kind="historical_annotation_key_digest_only",
                    expected_direct_digest="abc",
                )
                for name in ["done", "retry"]
            ]
        ),
    )
    old = tmp_path / "old.jsonl"
    old.write_text(
        json.dumps(dict(id="done", state="digest_only"))
        + "\n"
        + json.dumps(dict(id="retry", state="unresolved"))
        + "\n"
    )
    result = audit(request_file, tmp_path / "supplement", old)
    assert result["rows"] == 1 and result["full_request_rows"] == 2
    assert result["unresolved_checkpoint_sha256"] == sha256sum(old)
    assert result["digest_only"] == 1


def test_final_confirmation_uses_only_two_remaining_budget_cells():
    from benchmarks.engine_profile.matrix import design as original
    from benchmarks.performance_sprint.pilot import design

    groups = design(dict(groups=original()), "memory-confirm")
    assert len(groups) == 1 and groups[0]["id"] == "g01"
    assert [c["purpose"] for c in groups[0]["cells"]] == [
        "memory-preflight",
        "control-high-memory",
    ]
    assert all(c["memory_mb"] == 8192 for c in groups[0]["cells"])


def test_reference_audit_hashes_retained_output_and_detects_drift(tmp_path):
    from benchmarks.run_pilot import write_json_atomic

    root = tmp_path / "references"
    for role in ["direct", "cached"]:
        directory = root / role
        directory.mkdir(parents=True)
        for name in ["output.bcf", "output.bcf.csi"]:
            (directory / name).write_bytes(b"fixture")
        write_json_atomic(
            directory / "metrics.json",
            dict(
                output_sha256=sha256sum(directory / "output.bcf"),
                output_index_sha256=sha256sum(directory / "output.bcf.csi"),
                command=[
                    "vcfcache",
                    "annotate",
                    *(["--uncached"] if role == "direct" else []),
                ],
                input_vcf="same.input.bcf",
                output_records=100,
                wall_seconds=1,
            ),
        )
    document = tmp_path / "request.json"
    write_json_atomic(
        document,
        dict(
            rows=[
                dict(
                    id="one",
                    family="fixture",
                    reference_kind="retained_complete_record_output",
                    direct_root=str(root),
                    cached_root=str(root),
                    direct_metrics_glob="direct/metrics.json",
                    cached_metrics_glob="cached/metrics.json",
                    historical_direct_seconds=1,
                    historical_cached_seconds=1,
                )
            ]
        ),
    )
    assert audit(document, tmp_path / "valid")["retained_outputs_verified"] == 1
    (root / "cached/output.bcf").write_bytes(b"drift")
    assert audit(document, tmp_path / "drift")["unresolved"] == 1


def manuscript_task():
    """A complete structural task; file hashes are checked by the worker."""
    return dict(
        id="cell-one",
        family="fixture",
        cpus=8,
        engine="adaptive",
        tool="fastvep",
        parameters=dict(threads=8, rayon_threads=8),
        historical_cached_seconds=300,
        input=dict(path="input.bcf", sha256="fixture", records=100),
        cache=dict(
            directory="cache",
            files={
                name: "fixture"
                for name in [
                    "annotation.yaml",
                    "vcfcache_annotated.bcf",
                    "vcfcache_annotated.bcf.csi",
                ]
            },
        ),
        annotation_assets=[dict(path="tool", sha256="fixture")],
        baseline=dict(
            cpus=8,
            reuse_decision="approved",
            direct_seconds=1000,
            approval_evidence=[dict(path="audit.json", sha256="fixture")],
            kind="complete_records",
            direct=dict(path="direct.bcf", sha256="fixture"),
            cached=dict(path="cached.bcf", sha256="fixture"),
        ),
    )


def test_manuscript_gate_checks_cpu_and_reuse_approval():
    from benchmarks.performance_sprint.manuscript import validate_task

    task = manuscript_task()
    validate_task(task)
    task["baseline"]["reuse_decision"] = "pending_provenance_audit"
    with pytest.raises(ValueError, match="not been approved"):
        validate_task(task)
    task["baseline"]["reuse_decision"] = "approved"
    task["parameters"]["rayon_threads"] = 16
    with pytest.raises(ValueError, match="threads exceed"):
        validate_task(task)


def test_freeze_requires_exact_task_document_and_complete_ids(tmp_path):
    from benchmarks.performance_sprint.manuscript import prepare
    from benchmarks.run_pilot import write_json_atomic

    tasks = tmp_path / "tasks.json"
    write_json_atomic(tasks, [manuscript_task()])
    freeze = dict(
        approved=False,
        tests_passed=True,
        tasks_sha256=sha256sum(tasks),
        required_task_ids=["cell-one"],
    )
    with pytest.raises(ValueError, match="not approved"):
        prepare(tmp_path / "blocked", tasks, freeze)
    freeze["approved"] = True
    freeze["required_task_ids"].append("missing")
    with pytest.raises(ValueError, match="were dropped"):
        prepare(tmp_path / "blocked", tasks, freeze)
    freeze["required_task_ids"] = ["cell-one"]
    prepare(tmp_path / "allowed", tasks, freeze)
    manifest = json.loads((tmp_path / "allowed/manifest.json").read_text())
    assert manifest["mode"] == "cached_only_no_implicit_direct_reruns"


def test_historical_digest_matches_original_shell_protocol(tmp_path):
    import subprocess

    from benchmarks.performance_sprint.manuscript import annotation_digest
    from tests.test_fused_engine import dataset

    sample, *_ = dataset(tmp_path / "fixture", records=30)
    blocksort = (
        REPO_ROOT
        / "benchmarks/reproduction/final_clinical_campaign/blocksort.frozen.awk"
    )
    query_format = "%CHROM\t%POS\t%REF\t%ALT\n"
    reference = dict(
        query_format=query_format,
        blocksort=dict(path=str(blocksort), sha256=sha256sum(blocksort)),
    )
    query = subprocess.check_output(
        ["bcftools", "query", "-f", query_format, str(sample)]
    )
    sorted_rows = subprocess.check_output(["awk", "-f", str(blocksort)], input=query)
    for line in sorted_rows.splitlines(keepends=True):
        with (tmp_path / line.split(b"\t", 1)[0].decode()).open("ab") as stream:
            stream.write(line)
    command = (
        'for c in 1 2; do printf \'%s %s\\n\' "$c" "$(md5sum < "'
        + str(tmp_path)
        + "/$c\" | cut -d' ' -f1)\"; done | md5sum | cut -d' ' -f1"
    )
    expected = subprocess.check_output(
        ["bash", "-o", "pipefail", "-c", command], text=True
    ).strip()
    assert annotation_digest(sample, reference) == expected


def test_manuscript_collector_keeps_missing_jobs_visible(tmp_path):
    from benchmarks.performance_sprint.manuscript import collect, prepare
    from benchmarks.run_pilot import write_json_atomic

    tasks = tmp_path / "tasks.json"
    write_json_atomic(tasks, [manuscript_task()])
    prepare(
        tmp_path / "campaign",
        tasks,
        dict(
            approved=True,
            tests_passed=True,
            tasks_sha256=sha256sum(tasks),
            required_task_ids=["cell-one"],
        ),
    )
    result = collect(tmp_path / "campaign")
    assert result["complete"] == 0
    assert len(result["rows"]) == 1
    assert len(list((tmp_path / "campaign/reviews").glob("*.json"))) == 1


@pytest.mark.parametrize("new_pair", [False, True])
def test_cached_runner_executes_real_cli_with_mock_slurm_only(
    tmp_path, monkeypatch, new_pair
):
    import shutil
    import subprocess
    import sys

    import yaml

    from benchmarks.performance_sprint import manuscript
    from benchmarks.run_pilot import write_json_atomic
    from tests.test_fused_engine import dataset, run_case

    files = dataset(tmp_path / "fixture", records=20)
    reference, _ = run_case(tmp_path / "reference", files, "legacy")
    base = tmp_path / "cache-root"
    cache = base / "cache/test"
    cache.mkdir(parents=True)
    (base / "blueprint").mkdir()
    (base / "workflow").mkdir()
    (base / "workflow/init.yaml").write_text("{}\n")
    for source, name in [
        (files[1], "vcfcache_annotated.bcf"),
        (Path(str(files[1]) + ".csi"), "vcfcache_annotated.bcf.csi"),
        (files[3], "annotation.yaml"),
        (files[2], "params.snapshot.yaml"),
    ]:
        shutil.copy2(source, cache / name)
    task = manuscript_task()
    task.update(
        engine="legacy", cpus=6, parameters=yaml.safe_load(files[2].read_text())
    )
    task["baseline"]["cpus"] = 6
    task["input"] = dict(path=str(files[0]), sha256=sha256sum(files[0]), records=20)
    task["cache"] = dict(
        directory=str(cache), files={p.name: sha256sum(p) for p in cache.iterdir()}
    )
    artifact = dict(path=str(files[3]), sha256=sha256sum(files[3]))
    task["annotation_assets"] = [artifact]
    task["baseline"]["approval_evidence"] = [artifact]
    task["baseline"]["direct"] = task["baseline"]["cached"] = dict(
        path=str(reference), sha256=sha256sum(reference)
    )
    if new_pair:
        task.update(family="veprs_extension", tool="veprs", level="core")
        task["baseline"].update(
            kind="new_paired_reference",
            reuse_decision="new_direct_required",
            direct_seconds=None,
        )
    tasks = tmp_path / "tasks.json"
    write_json_atomic(tasks, [task])
    build = dict(
        root=str(REPO_ROOT),
        files={"vcfcache/cli.py": sha256sum(REPO_ROOT / "vcfcache/cli.py")},
    )
    campaign = tmp_path / "campaign"
    manuscript.prepare(
        campaign,
        tasks,
        dict(
            approved=True,
            tests_passed=True,
            tasks_sha256=sha256sum(tasks),
            required_task_ids=[task["id"]],
            build=build,
        ),
    )
    monkeypatch.setattr(manuscript, "PYTHON", sys.executable)
    original_temp = manuscript.tempfile.TemporaryDirectory
    monkeypatch.setattr(
        manuscript.tempfile,
        "TemporaryDirectory",
        lambda **kwargs: original_temp(prefix=kwargs["prefix"], dir=tmp_path),
    )
    original_run = subprocess.run

    def local_step(command, **kwargs):
        if command[0] != "srun":
            return original_run(command, **kwargs)
        directory = Path(command[command.index("--directory") + 1])
        result = original_run(command[command.index("--") + 1 :], **kwargs)
        write_json_atomic(
            directory / "resources.json",
            dict(wall_seconds=1, cgroup_after=dict(test_mock=True)),
        )
        return result

    monkeypatch.setattr(subprocess, "run", local_step)
    result = manuscript.run(campaign, task["id"], 1)
    assert result["state"] == "complete", result.get("traceback")
    assert result["input_records"] == 20
    assert result["cache_hit_rate"] == 0.5
    if new_pair:
        assert (campaign / "attempts/a001/cell-one/direct.bcf").exists()
        assert result["direct_seconds"] == 1
    assert (campaign / "attempts/a001/cell-one/comparison.json").exists()
    # Explicit replay is not allowed to overwrite its first attempt.
    with pytest.raises(FileExistsError):
        manuscript.run(campaign, task["id"], 1)


def test_partial_slurm_submission_still_registers_collector(tmp_path, monkeypatch):
    from benchmarks.performance_sprint import manuscript
    from benchmarks.run_pilot import write_json_atomic

    tasks = [manuscript_task(), dict(manuscript_task(), id="cell-two")]
    document = tmp_path / "tasks.json"
    write_json_atomic(document, tasks)
    campaign = tmp_path / "campaign"
    manuscript.prepare(
        campaign,
        document,
        dict(
            approved=True,
            tests_passed=True,
            tasks_sha256=sha256sum(document),
            required_task_ids=[t["id"] for t in tasks],
            build=dict(root="/results/test"),
        ),
    )
    submitted = []

    def sbatch(command, **kwargs):
        submitted.append(command)
        if len(submitted) == 2:
            raise RuntimeError("scheduler unavailable")
        return "123\n" if len(submitted) == 1 else "125\n"

    monkeypatch.setattr(manuscript.subprocess, "check_output", sbatch)
    with pytest.raises(RuntimeError, match="scheduler unavailable"):
        manuscript.submit(campaign)
    history = json.loads((campaign / "submissions.json").read_text())
    assert history[0]["jobs"][0]["job_id"] == "123"
    assert history[0]["collector"] == "125"
    assert "--dependency=afterany:123" in submitted[-1]
    assert all("--exclude=sl-w5" not in command for command in submitted)
    with pytest.raises(ValueError, match="already submitted"):
        manuscript.submit(campaign)


def test_worker_exclusions_must_be_explicit_and_justified():
    from benchmarks.performance_sprint.manuscript import exclusion_args

    assert exclusion_args({}) == []
    assert exclusion_args(
        dict(node_exclusions={"sl-w5": "failed current preflight"})
    ) == ["--exclude=sl-w5"]
    with pytest.raises(ValueError, match="explicit reason"):
        exclusion_args(dict(node_exclusions={"sl-w5": ""}))


def test_verified_metadata_never_ignores_a_changed_file(tmp_path):
    from benchmarks.performance_sprint.assemble import Files
    from benchmarks.performance_sprint.manuscript import verified_file

    path = tmp_path / "asset"
    path.write_text("original")
    files = Files()
    record = files.get(path)
    assert verified_file(record, files.known) == path
    path.write_text("changed contents")
    with pytest.raises(RuntimeError, match="frozen file changed"):
        verified_file(record, files.known)


def test_audit_merge_preserves_coverage_and_rejects_unexpected_rows(tmp_path):
    from benchmarks.performance_sprint.assemble import merge_audits

    original, supplement = tmp_path / "original.jsonl", tmp_path / "supplement.jsonl"
    original.write_text(
        json.dumps(dict(id="a", state="retained_outputs_verified"))
        + "\n"
        + json.dumps(dict(id="b", state="unresolved"))
        + "\n"
    )
    supplement.write_text(
        json.dumps(dict(id="b", state="historical_semantic_metadata_only")) + "\n"
    )
    assert set(merge_audits(original, supplement)) == {"a", "b"}
    supplement.write_text(
        json.dumps(dict(id="a", state="retained_outputs_verified")) + "\n"
    )
    with pytest.raises(ValueError, match="original unresolved"):
        merge_audits(original, supplement)


def test_dispatch_stops_on_failed_first_cells(tmp_path, monkeypatch):
    from benchmarks.performance_sprint import manuscript

    monkeypatch.setattr(
        manuscript, "verify", lambda root: dict(freeze=dict(first_cells=["a", "b"]))
    )
    monkeypatch.setattr(
        manuscript,
        "latest",
        lambda root, cell: dict(state="complete" if cell == "a" else "failed"),
    )
    with pytest.raises(RuntimeError, match="blocked automatic dispatch"):
        manuscript.dispatch(tmp_path)
    assert (tmp_path / "DISPATCH_BLOCKED.json").exists()


def test_dispatch_schedules_only_remaining_primary_cells(tmp_path, monkeypatch):
    from benchmarks.performance_sprint import manuscript

    monkeypatch.setattr(
        manuscript,
        "verify",
        lambda root: dict(
            freeze=dict(first_cells=["a"]), tasks=[dict(id="a"), dict(id="b")]
        ),
    )
    monkeypatch.setattr(manuscript, "latest", lambda root, cell: dict(state="complete"))
    calls = []
    monkeypatch.setattr(
        manuscript,
        "submit",
        lambda root, cells: calls.append(cells) or dict(jobs=["123"]),
    )
    manuscript.dispatch(tmp_path)
    assert calls == [["b"]]


def test_execution_deduplication_retains_source_aliases():
    import copy

    from benchmarks.performance_sprint.assemble import deduplicate

    original = manuscript_task()
    duplicate = copy.deepcopy(original)
    duplicate["id"] = "plot-copy"
    duplicate["family"] = "another_plot"
    tasks, aliases = deduplicate([original, duplicate])
    assert len(tasks) == 1
    assert aliases == {"plot-copy": "cell-one"}
    assert tasks[0]["source_aliases"] == ["plot-copy"]


def test_c2_pgp_does_not_require_membership_in_full52():
    from benchmarks.performance_sprint.assemble import c2_input

    assert c2_input(
        dict(assembly="GRCh37", cohort="pgp", sample="hu24385B"), {}
    ).endswith("/GRCh37/pgp/hu24385B.GRCh37.small_variants.vcf.gz")
    assert (
        c2_input(
            dict(assembly="GRCh38", cohort="sgdp", sample="sample"),
            {"sample": dict(input_vcf="original.vcf.gz")},
        )
        == "original.vcf.gz"
    )
