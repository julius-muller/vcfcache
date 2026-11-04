"""The adaptive engine must reproduce the legacy path record for record."""

import subprocess
from pathlib import Path

import pytest
import yaml

from vcfcache.database.adaptive import plan_shards
from vcfcache.database.workflow_manager import WorkflowManager

DATA = Path(__file__).resolve().parent / "data" / "nodata"
CONFIG = Path(__file__).resolve().parent / "config"


@pytest.fixture(scope="module")
def cache(tmp_path_factory):
    """Build a cache from the committed fixtures, as cache-build would."""
    root = tmp_path_factory.mktemp("adaptive-cache")
    manager = WorkflowManager(
        input_file=DATA / "sample4.bcf",
        output_dir=root,
        name="adaptive_cache",
        anno_config_file=CONFIG / "test_annotation.yaml",
        params_file=CONFIG / "test_params.yaml",
        verbosity=0,
    )
    manager.run("cache-build", db_bcf=DATA / "gnomad_test.bcf")
    built = list(root.rglob("vcfcache_annotated.bcf"))
    assert built, "cache-build produced no cache BCF"
    return built[0]


def annotate(
    root, cache, engine, *, threads=8, input_file=None, recipe_cmd=None, **options
):
    """Annotate the shared sample with one engine and return (output, stats)."""
    root.mkdir(parents=True, exist_ok=True)
    output = root / "result.bcf"
    parameters = yaml.safe_load((CONFIG / "test_params.yaml").read_text())
    parameters["threads"] = threads
    params_file = root / "params.yaml"
    params_file.write_text(yaml.safe_dump(parameters))
    manager = WorkflowManager(
        input_file=input_file or DATA / "sample4.bcf",
        output_dir=root,
        name="adaptive_run",
        anno_config_file=CONFIG / "test_annotation.yaml",
        params_file=params_file,
        output_file=output,
        verbosity=0,
    )
    if recipe_cmd is not None:
        manager.nfa_config_content["annotation_cmd"] = recipe_cmd
    manager.run("annotate", db_bcf=cache, engine=engine, **options)
    return output, manager.last_run_stats


def records(path):
    """Return the record block, which must match the legacy output exactly."""
    return subprocess.run(
        ["bcftools", "view", "-H", str(path)],
        capture_output=True,
        text=True,
        check=True,
    ).stdout


@pytest.mark.parametrize("workers", [1, 2, 4])
@pytest.mark.parametrize("preserve", [False, True])
@pytest.mark.parametrize("skip_split", [False, True])
@pytest.mark.parametrize("fused", [False, True])
def test_matches_legacy(
    cache, tmp_path, monkeypatch, workers, preserve, skip_split, fused
):
    """Multiallelics, spanning deletions and cache-absent contigs must agree."""
    options = dict(preserve_unannotated=preserve, skip_split_multiallelic=skip_split)
    monkeypatch.setenv("VCFCACHE_EXPERIMENT_FUSED_PUBLISH", str(int(fused)))
    expected, baseline = annotate(tmp_path / "legacy", cache, "legacy", **options)
    actual, stats = annotate(
        tmp_path / "adaptive", cache, "adaptive", engine_workers=workers, **options
    )
    assert records(actual) == records(expected)
    for key in ("input_variants", "output_variants", "missing_variants"):
        assert stats[key] == baseline[key]
    assert stats["engine"] == "adaptive"
    assert stats["fused_publish"] is fused
    if fused:
        assert stats["shards"] * (2 + stats["publish_threads_per_shard"]) <= 8
        assert not list((tmp_path / "adaptive").rglob("step4.bcf"))
        assert not list((tmp_path / "adaptive").rglob("*_merged.bcf"))


def test_fused_single_cpu_retains_original_path(cache, tmp_path, monkeypatch):
    """A two-process merge pipe must not be launched on the one-CPU path."""
    monkeypatch.setenv("VCFCACHE_EXPERIMENT_FUSED_PUBLISH", "1")
    output, stats = annotate(tmp_path, cache, "adaptive", threads=1)
    assert output.exists() and stats["fused_publish"] is False


@pytest.mark.parametrize("preserve", [False, True])
def test_fused_filters_after_overlay(cache, tmp_path, monkeypatch, preserve):
    """A recipe may keep records without adding the required annotation value."""
    command = "bcftools view ${INPUT_BCF} -Ob -o ${OUTPUT_BCF} -W"
    expected, _ = annotate(
        tmp_path / "legacy",
        cache,
        "legacy",
        recipe_cmd=command,
        preserve_unannotated=preserve,
    )
    monkeypatch.setenv("VCFCACHE_EXPERIMENT_FUSED_PUBLISH", "1")
    actual, stats = annotate(
        tmp_path / "fused",
        cache,
        "adaptive",
        recipe_cmd=command,
        preserve_unannotated=preserve,
        engine_workers=2,
    )
    assert stats["missing_variants"] > 0 and stats["fused_publish"] is True
    assert records(actual) == records(expected)


@pytest.mark.parametrize("kind", ["all-hits", "empty-cache", "empty-input"])
def test_fused_degenerate_inputs(cache, tmp_path, monkeypatch, kind):
    """Zero misses and empty streams still publish valid, indexed, sorted BCF."""
    fully_annotated, _ = annotate(tmp_path / "reference", cache, "legacy")
    empty = tmp_path / "empty.bcf"
    subprocess.run(
        [
            "bcftools",
            "view",
            "-i",
            "POS<0",
            str(fully_annotated),
            "-Ob",
            "-W",
            "-o",
            str(empty),
        ],
        check=True,
    )
    source = empty if kind == "empty-input" else DATA / "sample4.bcf"
    database = empty if kind == "empty-cache" else fully_annotated
    expected, _ = annotate(tmp_path / "legacy", database, "legacy", input_file=source)
    monkeypatch.setenv("VCFCACHE_EXPERIMENT_FUSED_PUBLISH", "1")
    actual, stats = annotate(
        tmp_path / "fused", database, "adaptive", input_file=source
    )
    assert stats["fused_publish"] is True
    assert records(actual) == records(expected)
    if kind != "empty-cache":
        assert stats["missing_variants"] == 0


def test_incompatible_naive_headers_preserve_existing_output(tmp_path):
    """Do not force concatenation across incompatible BCF tag dictionaries."""
    import logging

    import pysam

    from vcfcache.database.adaptive import publish_shards

    parts = []
    for index, order in enumerate([("FIRST", "SECOND"), ("SECOND", "FIRST")]):
        header = pysam.VariantHeader()
        header.contigs.add("1", length=1000)
        for tag in order:
            header.info.add(tag, 1, "Integer", "test")
        path = tmp_path / f"part-{index}.bcf"
        with pysam.VariantFile(path, "wb", header=header) as writer:
            record = writer.new_record(contig="1", start=index, alleles=("A", "C"))
            record.info["FIRST"] = 1
            record.info["SECOND"] = 2
            writer.write(record)
        subprocess.run(["bcftools", "index", str(path)], check=True)
        parts.append(path)
    output = tmp_path / "result.bcf"
    output.write_bytes(b"existing output")
    manager = WorkflowManager.__new__(WorkflowManager)
    manager.logger = logging.getLogger(__name__)
    with pytest.raises(subprocess.CalledProcessError):
        publish_shards(manager, parts, output, tmp_path, "bcftools")
    assert output.read_bytes() == b"existing output"


def test_fused_binary_stdout(cache, tmp_path):
    """BCF stdout must not be decoded as text or mixed with log output."""
    import os
    import sys

    expected, _ = annotate(tmp_path / "legacy", cache, "legacy")
    directory = tmp_path / "stdout"
    directory.mkdir()
    program = """
import sys
from pathlib import Path
from vcfcache.database.workflow_manager import WorkflowManager
m=WorkflowManager(input_file=Path(sys.argv[1]),output_dir=Path(sys.argv[2]),name="stdout",
    anno_config_file=Path(sys.argv[3]),params_file=Path(sys.argv[4]),output_file="-",verbosity=0)
m.run("annotate",db_bcf=Path(sys.argv[5]),engine="adaptive",engine_workers=2)
"""
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            program,
            str(DATA / "sample4.bcf"),
            str(directory),
            str(CONFIG / "test_annotation.yaml"),
            str(tmp_path / "legacy/params.yaml"),
            str(cache),
        ],
        env=dict(os.environ, VCFCACHE_EXPERIMENT_FUSED_PUBLISH="1"),
        capture_output=True,
        check=True,
    )
    decoded = subprocess.run(
        ["bcftools", "view", "-H", "-"],
        input=result.stdout,
        capture_output=True,
        check=True,
    )
    assert decoded.stdout.decode() == records(expected)


def test_fused_pipe_failure_preserves_previous_output(cache, tmp_path, monkeypatch):
    """An upstream failure must not replace an existing requested deliverable."""
    from vcfcache.database.workflow_manager import BcftoolsCommand

    output = tmp_path / "result.bcf"
    output.write_bytes(b"previous output")
    original = BcftoolsCommand.run

    def failing_pipeline(command, *args, **kwargs):
        if " -c INFO -Ou | " in command.cmd:
            command.cmd = "false | true"
        return original(command, *args, **kwargs)

    monkeypatch.setattr(BcftoolsCommand, "run", failing_pipeline)
    monkeypatch.setenv("VCFCACHE_EXPERIMENT_FUSED_PUBLISH", "1")
    with pytest.raises(RuntimeError, match="Workflow execution failed"):
        annotate(tmp_path, cache, "adaptive", engine_workers=2)
    assert output.read_bytes() == b"previous output"


def test_gathered_shards_stay_sorted():
    """Shards must be consecutive in index order, or a naive gather misorders."""
    occupied = [(f"chr{index}", 1000, 10 * index) for index in range(1, 9)]
    for workers in (0, 1, 3, 8, 64):
        shards = plan_shards(occupied, workers)
        assert [contig for shard in shards for contig, _ in shard] == [
            contig for contig, _, _ in occupied
        ]
        assert all(shard for shard in shards)


@pytest.mark.parametrize("engine", ["legacy", "adaptive"])
@pytest.mark.parametrize("uncached", [False, True])
def test_zero_exit_truncated_annotation_is_rejected(cache, tmp_path, engine, uncached):
    """An indexable, zero-exit partial recipe must not become a successful run."""
    manager = WorkflowManager(
        input_file=DATA / "sample4.bcf",
        output_dir=tmp_path,
        name="truncated",
        anno_config_file=CONFIG / "test_annotation.yaml",
        params_file=CONFIG / "test_params.yaml",
        output_file=tmp_path / "result.bcf",
        verbosity=0,
    )
    original = manager.nfa_config_content["annotation_cmd"]
    manager.nfa_config_content["annotation_cmd"] = (
        original
        + '\nbcftools view -i "POS<0" ${OUTPUT_BCF} -Ob -o ${OUTPUT_BCF}.partial -W'
        + "\nmv ${OUTPUT_BCF}.partial ${OUTPUT_BCF}"
        + "\nmv ${OUTPUT_BCF}.partial.csi ${OUTPUT_BCF}.csi"
    )
    with pytest.raises(RuntimeError, match="Annotation record count changed"):
        manager.run(
            "annotate-nocache" if uncached else "annotate", db_bcf=cache, engine=engine
        )
