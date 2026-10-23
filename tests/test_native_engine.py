"""Differential tests for resident BCF/CSI lookup, no external annotator needed."""

import os
import shutil
import subprocess
import sys

import pysam
import pytest
import yaml

from benchmarks.run_pilot import strict_semantic_compare
from tests import test_fused_engine as differential
from tests.test_fused_engine import checked, dataset, run_case
from tests.test_fused_engine import native as native_fixture

native = native_fixture


@pytest.mark.parametrize(
    "check,options",
    [
        ("test_reordered_info_dictionary", {}),
        ("test_missing_cache_contig", {"lookup": "regions"}),
        ("test_missing_cache_contig", {"lookup": "contigs"}),
        ("test_lookup_modes_and_boundaries", {"lookup": "auto", "workers": 1}),
        ("test_lookup_modes_and_boundaries", {"lookup": "regions", "workers": 2}),
        ("test_lookup_modes_and_boundaries", {"lookup": "contigs", "workers": 2}),
        ("test_missing_required_values", {"preserve": False}),
        ("test_missing_required_values", {"preserve": True}),
    ],
)
def test_shared_native_differential_cases(
    native, tmp_path, monkeypatch, check, options
):
    original = differential.run_case

    def resident_case(root, files, engine, **kwargs):
        return original(
            root, files, "native" if engine == "fused" else engine, **kwargs
        )

    # Reuse identical real-bcftools fixtures, substituting only the engine.
    monkeypatch.setattr(differential, "run_case", resident_case)
    getattr(differential, check)(native, tmp_path, **options)


@pytest.mark.parametrize("hits", [0, 1, 2, 10])
@pytest.mark.parametrize("workers", [1, 3])
def test_native_complete_records(native, tmp_path, hits, workers):
    files = dataset(tmp_path / "data", hit_every=hits)
    legacy, old = run_case(tmp_path, files, "legacy")
    output, new = run_case(tmp_path, files, "native", workers=workers)
    assert strict_semantic_compare(output, legacy)["semantic_pass"]
    assert new.last_run_stats["engine"] == "native"
    for key in (
        "input_variants",
        "missing_variants",
        "output_variants",
        "missing_annotated",
    ):
        assert new.last_run_stats[key] == old.last_run_stats[key]
    assert not list(new.work_dir.rglob("hits.bcf"))


@pytest.mark.parametrize("kind", ["normal", "multiallelic", "spanning", "preannotated"])
@pytest.mark.parametrize("preserve", [False, True])
def test_native_representation(native, tmp_path, kind, preserve):
    files = dataset(tmp_path / "data", kind=kind)
    legacy, _ = run_case(tmp_path, files, "legacy", preserve=preserve)
    output, _ = run_case(tmp_path, files, "native", preserve=preserve)
    assert strict_semantic_compare(output, legacy)["semantic_pass"]


def test_native_empty(native, tmp_path):
    files = dataset(tmp_path / "data", records=0)
    legacy, _ = run_case(tmp_path, files, "legacy")
    output, manager = run_case(tmp_path, files, "native")
    assert strict_semantic_compare(output, legacy)["semantic_pass"]
    assert manager.last_run_stats["candidate_records"] == 0


def test_native_memory_fallback(native, tmp_path):
    files = dataset(tmp_path / "data", records=2000, payload=100)
    legacy, _ = run_case(tmp_path, files, "legacy")
    from vcfcache.database.workflow_manager import WorkflowManager

    output = tmp_path / "fallback.bcf"
    manager = WorkflowManager(
        input_file=files[0],
        output_dir=tmp_path / "run",
        name="memory",
        params_file=files[2],
        anno_config_file=files[3],
        output_file=output,
    )
    manager.run("annotate", db_bcf=files[1], engine="native", engine_memory_mb=1)
    assert strict_semantic_compare(output, legacy)["semantic_pass"]
    assert manager.last_run_stats["requested_engine"] == "native"
    assert "memory" in manager.last_run_stats["fallback_reason"]


@pytest.mark.parametrize("behavior", ["drop", "format", "reverse", "header"])
def test_native_recipe_contract(native, tmp_path, behavior):
    files = dataset(tmp_path / "data")
    recipe = yaml.safe_load(files[3].read_text())
    command = recipe["annotation_cmd"].replace("-Ob -W -o ${OUTPUT_BCF}", "-Ou")
    if behavior == "drop":
        command += " | bcftools view -e 'POS%4==0' -Ob -W -o ${OUTPUT_BCF}"
    elif behavior == "format":
        command += " | bcftools annotate -x FORMAT/DP -Ob -W -o ${OUTPUT_BCF}"
    elif behavior == "reverse":
        command += " | bcftools view -Ov | awk '/^#/ {print; next} {line[++n]=$0} END {for(i=n;i>0;i--) print line[i]}' | bcftools view -Ob -o ${OUTPUT_BCF}"
    else:
        header = tmp_path / "header.txt"
        header.write_text('##INFO=<ID=NEW,Number=1,Type=String,Description="new">\n')
        command += f" | bcftools annotate -h {header} -Ob -W -o ${{OUTPUT_BCF}}"
    recipe["annotation_cmd"] = command
    files[3].write_text(yaml.safe_dump(recipe))
    legacy, _ = run_case(
        tmp_path, files, "fused" if behavior == "reverse" else "legacy"
    )
    output, _ = run_case(tmp_path, files, "native")
    assert strict_semantic_compare(output, legacy)["semantic_pass"]


def test_native_failure_no_public_output(native, tmp_path):
    files = dataset(tmp_path / "data", hit_every=0)
    recipe = yaml.safe_load(files[3].read_text())
    recipe["annotation_cmd"] = "exit 19"
    files[3].write_text(yaml.safe_dump(recipe))
    with pytest.raises(RuntimeError, match="19"):
        run_case(tmp_path, files, "native")
    assert not (tmp_path / "native/output.bcf").exists()


def test_native_cli_budget(tmp_path):
    files = dataset(tmp_path / "data")
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            "from vcfcache.cli import main; main()",
            "annotate",
            "-i",
            str(files[0]),
            "-a",
            str(files[1]),
            "-o",
            str(tmp_path / "output.bcf"),
            "--engine",
            "native",
            "--engine-memory-mb",
            "0",
        ],
        text=True,
        capture_output=True,
    )
    assert result.returncode == 2
    assert "--engine-memory-mb must be positive" in result.stderr


@pytest.mark.parametrize(
    "field,value",
    [("AN", None), ("AF", (None,)), ("SUP", "."), ("FLAG", False), ("SUP", "")],
)
def test_native_missing_overlay_values(native, tmp_path, field, value):
    files = dataset(tmp_path / "data", hit_every=1)
    for path, is_cache in ((files[0], False), (files[1], True)):
        replacement = path.with_suffix(".new.bcf")
        with (
            pysam.VariantFile(str(path)) as source,
            pysam.VariantFile(str(replacement), "wb", header=source.header) as out,
        ):
            for record in source:
                if is_cache:
                    record.info[field] = value
                elif field == "FLAG":
                    record.info[field] = True
                elif field == "SUP":
                    record.info[field] = "original"
                out.write(record)
        replacement.replace(path)
        checked("bcftools", "index", "-f", path)
    legacy, _ = run_case(tmp_path, files, "legacy")
    output, _ = run_case(tmp_path, files, "native")
    assert strict_semantic_compare(output, legacy)["semantic_pass"]


def test_native_large_header_index(native, tmp_path):
    files = dataset(tmp_path / "data", records=400, payload=400)
    replacement = tmp_path / "large-header.bcf"
    with pysam.VariantFile(str(files[0])) as source:
        header = source.header.copy()
        for i in range(3000):
            header.contigs.add(f"unused_scaffold_{i}", length=100000)
        with pysam.VariantFile(str(replacement), "wb", header=header) as out:
            for record in source:
                record.translate(header)
                out.write(record)
    checked("bcftools", "index", replacement)
    files = replacement, *files[1:]
    legacy, _ = run_case(tmp_path, files, "legacy")
    output, _ = run_case(tmp_path, files, "native")
    assert strict_semantic_compare(output, legacy)["semantic_pass"]


def test_native_recipe_padded_pairing_reuses_recipe(native, tmp_path):
    files = dataset(tmp_path / "data", hit_every=0)
    truth = tmp_path / "data/truth.bcf"
    replacement = tmp_path / "padded.bcf"
    with (
        pysam.VariantFile(str(truth)) as source,
        pysam.VariantFile(str(replacement), "wb", header=source.header) as out,
    ):
        for i, record in enumerate(source):
            if i == 0:
                record.alleles = ("AC", "CC")
            out.write(record)
    checked("bcftools", "index", replacement)
    count = tmp_path / "recipe-invocations.txt"
    recipe = yaml.safe_load(files[3].read_text())
    recipe["annotation_cmd"] = (
        f"printf x >> {count}; bcftools view {replacement} -Ob -W -o ${{OUTPUT_BCF}}"
    )
    files[3].write_text(yaml.safe_dump(recipe))
    legacy, _ = run_case(tmp_path, files, "legacy")
    output, manager = run_case(tmp_path, files, "native")
    assert strict_semantic_compare(output, legacy)["semantic_pass"]
    assert count.read_text() == "xx"
    assert "compatibility_merge_reason" in manager.last_run_stats


def test_native_padded_cache_fallback(native, tmp_path):
    files = dataset(tmp_path / "data", hit_every=1)
    replacement = tmp_path / "padded-cache.bcf"
    with (
        pysam.VariantFile(str(files[1])) as source,
        pysam.VariantFile(str(replacement), "wb", header=source.header) as out,
    ):
        for i, record in enumerate(source):
            if i == 0:
                record.alleles = ("AC", "CC")
            out.write(record)
    checked("bcftools", "index", replacement)
    files = files[0], replacement, *files[2:]
    legacy, _ = run_case(tmp_path, files, "legacy")
    output, manager = run_case(tmp_path, files, "native")
    assert strict_semantic_compare(output, legacy)["semantic_pass"]
    assert "non-exact" in manager.last_run_stats["fallback_reason"]


@pytest.mark.parametrize("statistics", ["light", "full", "none"])
def test_resident_public_stdout(native, tmp_path, statistics):
    files = dataset(tmp_path / "data")
    base = tmp_path / "cache-root"
    cache = base / "cache/test"
    cache.mkdir(parents=True)
    (base / "blueprint").mkdir()
    (base / "workflow").mkdir()
    (base / "workflow/init.yaml").write_text("{}\n")
    for source, name in (
        (files[1], "vcfcache_annotated.bcf"),
        (files[1].with_suffix(".bcf.csi"), "vcfcache_annotated.bcf.csi"),
        (files[2], "params.snapshot.yaml"),
        (files[3], "annotation.yaml"),
    ):
        shutil.copyfile(source, cache / name)
    command = [
        sys.executable,
        "-m",
        "vcfcache",
        "annotate",
        "-a",
        str(cache),
        "-i",
        str(files[0]),
        "-o",
        "-",
        "--engine",
        "native",
        "--statistics",
        statistics,
    ]
    if statistics != "none":
        command += ["--stats-dir", str(tmp_path / "stats")]
    result = subprocess.run(command, capture_output=True)
    assert result.returncode == 0, result.stderr.decode(errors="replace")
    output = tmp_path / "stdout.bcf"
    output.write_bytes(result.stdout)
    checked("bcftools", "index", output)
    assert int(checked("bcftools", "index", "-n", output).stdout) == 20
    if statistics != "none":
        stats = yaml.safe_load(
            (tmp_path / "stats/sample_vcstats/compare_stats.yaml").read_text()
        )
        assert stats["engine"] == "native"
        assert stats["buffer_bytes"] > 0
        assert stats["variant_counts"]["total_output"] == 20


@pytest.mark.parametrize("fails", [False, True])
def test_resident_restores_affinity(native, tmp_path, fails):
    files = dataset(tmp_path / "data", hit_every=0)
    if fails:
        recipe = yaml.safe_load(files[3].read_text())
        recipe["annotation_cmd"] = "exit 21"
        files[3].write_text(yaml.safe_dump(recipe))
    before = os.sched_getaffinity(0)
    try:
        if fails:
            with pytest.raises(RuntimeError):
                run_case(tmp_path, files, "native")
        else:
            run_case(tmp_path, files, "native")
    finally:
        assert os.sched_getaffinity(0) == before
