"""Prepared-input experiments must remain equivalent on their verified domain."""

import hashlib
import subprocess
from pathlib import Path

import pytest
import yaml

from benchmarks.adaptive_hts125_pilot.run import (
    INTEL,
    MODES,
    preflight_prepared,
    wrapper,
)
from vcfcache.database.adaptive import write_regions
from vcfcache.database.workflow_manager import WorkflowManager

DATA = Path(__file__).resolve().parent / "data/nodata"
CONFIG = Path(__file__).resolve().parent / "config"


def annotate(root, cache, engine, input_file, **options):
    """Use eight configured CPUs and the selected immutable fixture input."""
    root.mkdir()
    params = yaml.safe_load((CONFIG / "test_params.yaml").read_text())
    params["threads"] = 8
    parameters = root / "params.yaml"
    parameters.write_text(yaml.safe_dump(params))
    output = root / "result.bcf"
    manager = WorkflowManager(
        input_file=input_file,
        output_dir=root,
        name="pilot",
        anno_config_file=CONFIG / "test_annotation.yaml",
        params_file=parameters,
        output_file=output,
    )
    manager.run("annotate", db_bcf=cache, engine=engine, **options)
    return output, manager.last_run_stats


def records(path):
    """Headers may include timestamps; complete record values must match."""
    return subprocess.check_output(["bcftools", "view", "-H", str(path)])


@pytest.fixture(scope="module")
def cache(tmp_path_factory):
    """Build a small immutable annotation cache for equivalence checks."""
    root = tmp_path_factory.mktemp("adaptive125-cache")
    manager = WorkflowManager(
        input_file=DATA / "sample4.bcf",
        output_dir=root,
        name="pilot_cache",
        anno_config_file=CONFIG / "test_annotation.yaml",
        params_file=CONFIG / "test_params.yaml",
    )
    manager.run("cache-build", db_bcf=DATA / "gnomad_test.bcf")
    return next(root.rglob("vcfcache_annotated.bcf"))


@pytest.fixture
def prepared(tmp_path):
    """Select a domain-valid test subset; never preprocess real pilot inputs."""
    path = tmp_path / "prepared.bcf"
    subprocess.run(
        [
            "bcftools",
            "view",
            "-m2",
            "-M2",
            "-e",
            'ALT="*"',
            str(DATA / "sample4.bcf"),
            "-Ob",
            "-o",
            str(path),
            "-W",
        ],
        check=True,
    )
    return path


@pytest.mark.parametrize("workers", [2, 4])
@pytest.mark.parametrize("lookup", ["contigs", "windows"])
def test_prepared_routing_is_equivalent(
    cache, prepared, tmp_path, monkeypatch, workers, lookup
):
    """No new full-file validation pass belongs inside these engineering cells."""
    before = hashlib.sha256(prepared.read_bytes()).hexdigest()
    expected, _ = annotate(
        tmp_path / "reference",
        cache,
        "adaptive",
        input_file=prepared,
        preserve_unannotated=True,
    )
    monkeypatch.setenv("VCFCACHE_PILOT_PREPARED", "1")
    monkeypatch.setenv("VCFCACHE_PILOT_LOOKUP", lookup)
    monkeypatch.setenv("VCFCACHE_PILOT_THREAD_BUDGET", "1")
    monkeypatch.setenv("VCFCACHE_PILOT_MERGE_WINDOWS", "1")
    actual, stats = annotate(
        tmp_path / "pilot",
        cache,
        "adaptive",
        input_file=prepared,
        preserve_unannotated=True,
        engine_workers=workers,
    )
    assert records(actual) == records(expected)
    assert stats["pilot_prepared"] and stats["cache_access"] == lookup
    assert stats["shards"] * (1 + stats["pilot_codec_threads_per_shard"]) <= 8
    assert not list((tmp_path / "pilot").rglob("*normalized.bcf"))
    assert hashlib.sha256(prepared.read_bytes()).hexdigest() == before


def test_prepared_direct_bypasses_preprocessing(prepared, tmp_path, monkeypatch):
    """Direct and cached candidates must use the same prepared-input contract."""
    monkeypatch.setenv("VCFCACHE_PILOT_PREPARED", "1")
    manager = WorkflowManager(
        input_file=prepared,
        output_dir=tmp_path,
        name="direct",
        anno_config_file=CONFIG / "test_annotation.yaml",
        params_file=CONFIG / "test_params.yaml",
        output_file=tmp_path / "direct.bcf",
    )
    manager.run("annotate-nocache")
    assert (
        manager.last_run_stats["input_variants"]
        == manager.last_run_stats["output_variants"]
    )
    assert not list(tmp_path.rglob("*normalized.bcf"))


@pytest.mark.parametrize(
    "rows",
    [
        ["chr1\t1\t.\tA\tC,G\t.\tPASS\t."],
        ["chr1\t1\t.\tA\t*\t.\tPASS\t."],
        ["chr1\t1\t.\tA\tC\t.\tPASS\t."] * 2,
        ["chr1\t1\t.\tA\tC\t.\tPASS\tCSQ=old"],
    ],
)
def test_preflight_rejects_unsupported_input(tmp_path, rows):
    """The trusted fast path is not silently enabled on arbitrary VCFs."""
    source = tmp_path / "input.vcf"
    source.write_text(
        '##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n##INFO=<ID=CSQ,Number=.,Type=String,Description="Annotation">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        + "\n".join(rows)
        + "\n"
    )
    path = tmp_path / "input.bcf"
    subprocess.run(
        ["bcftools", "view", str(source), "-Ob", "-o", str(path), "-W"], check=True
    )
    with pytest.raises(ValueError):
        preflight_prepared("bcftools", path, "CSQ")


def test_thread_wrapper_does_not_override_shard_budgets(tmp_path):
    """A global eight-thread recipe wrapper must respect per-shard codec pools."""
    binary = tmp_path / "echo-args"
    binary.write_text('#!/bin/bash\nprintf "%s\\n" "$@"\n')
    binary.chmod(0o755)
    command = wrapper(tmp_path / "wrapped", binary, True)
    args = subprocess.check_output(
        [str(command), "annotate", "--threads", "3", "input"], text=True
    ).splitlines()
    assert args == ["annotate", "--threads", "3", "input"]
    assert "sl-w5" not in INTEL and len(MODES) == 7


def test_windows_are_coalesced_without_extra_records(tmp_path, monkeypatch):
    """Adjacent 16-kb blocks should form one region before calling HTSlib."""
    monkeypatch.setenv("VCFCACHE_PILOT_MERGE_WINDOWS", "1")
    monkeypatch.setattr(
        subprocess,
        "run",
        lambda *_a, **_k: type(
            "R", (), {"stdout": "chr1\t1\nchr1\t17000\nchr1\t50000\n"}
        )(),
    )
    paths = write_regions(
        "bcftools", tmp_path / "unused", [[("chr1", 100000)]], tmp_path, True
    )
    assert paths[0].read_text() == "chr1\t0\t32768\nchr1\t49152\t65536\n"
