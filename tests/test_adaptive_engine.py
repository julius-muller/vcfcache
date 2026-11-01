"""The adaptive engine must reproduce the legacy path record for record."""

import subprocess
from pathlib import Path

import pytest

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


def annotate(root, cache, engine, **options):
    """Annotate the shared sample with one engine and return (output, stats)."""
    root.mkdir(parents=True, exist_ok=True)
    output = root / "result.bcf"
    manager = WorkflowManager(
        input_file=DATA / "sample4.bcf",
        output_dir=root,
        name="adaptive_run",
        anno_config_file=CONFIG / "test_annotation.yaml",
        params_file=CONFIG / "test_params.yaml",
        output_file=output,
        verbosity=0,
    )
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
def test_matches_legacy(cache, tmp_path, workers, preserve, skip_split):
    """Multiallelics, spanning deletions and cache-absent contigs must agree."""
    options = dict(preserve_unannotated=preserve, skip_split_multiallelic=skip_split)
    expected, baseline = annotate(tmp_path / "legacy", cache, "legacy", **options)
    actual, stats = annotate(
        tmp_path / "adaptive", cache, "adaptive", engine_workers=workers, **options
    )
    assert records(actual) == records(expected)
    for key in ("input_variants", "output_variants", "missing_variants"):
        assert stats[key] == baseline[key]
    assert stats["engine"] == "adaptive"


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
