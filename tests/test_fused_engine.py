"""Differential tests of the optional native engine, using real bcftools."""

import json
import subprocess
from pathlib import Path

import pysam
import pytest
import yaml

from benchmarks.run_pilot import strict_semantic_compare
from vcfcache.database.fused import Pipelines, native_binary, validate_info_headers
from vcfcache.database.workflow_manager import WorkflowManager


@pytest.fixture
def native():
    try:
        return native_binary()
    except FileNotFoundError:
        pytest.skip("optional native executable not installed")


def checked(*args, **kwargs):
    return subprocess.run(
        [str(arg) for arg in args], check=True, capture_output=True, **kwargs
    )


def dataset(directory, *, hit_every=2, records=20, payload=20, kind="normal"):
    """Generate deterministic annotations and genotypes, independent of timing."""
    directory.mkdir(parents=True, exist_ok=True)
    header = pysam.VariantHeader()
    for chromosome in ("2", "1"):
        header.contigs.add(chromosome, length=100000000)
    for tag, number, dtype in (
        ("CSQ", ".", "String"),
        ("SUP", 1, "String"),
        ("AF", "A", "Float"),
        ("AC", "A", "Integer"),
        ("AN", 1, "Integer"),
        ("FLAG", 0, "Flag"),
    ):
        header.info.add(tag, number, dtype, f"test {tag}")
    for tag, number, dtype in (
        ("GT", 1, "String"),
        ("DP", 1, "Integer"),
        ("AD", "R", "Integer"),
    ):
        header.formats.add(tag, number, dtype, f"test {tag}")
    for sample in ("sample-A", "sample-B"):
        header.add_sample(sample)
    sample = directory / "sample.bcf"
    cache = directory / "cache.bcf"
    truth = directory / "truth.bcf"
    with (
        pysam.VariantFile(str(sample), "wb", header=header) as inp,
        pysam.VariantFile(str(cache), "wb", header=header) as cached,
        pysam.VariantFile(str(truth), "wb", header=header) as complete,
    ):
        for i in range(records):
            chromosome = "2" if i < records // 2 else "1"
            pos = (i % max(1, records // 2)) * 3 + 1
            alleles = ("A", "C")
            if kind == "multiallelic" and i == 1:
                alleles = ("A", "C", "G")
            elif kind == "spanning" and i == 1:
                alleles = ("A", "*")
            record = inp.new_record(
                contig=chromosome,
                start=pos,
                alleles=alleles,
                id=f"id-{i}",
                qual=40,
                filter="PASS",
            )
            record.info["AF"] = tuple(0.25 for _ in alleles[1:])
            record.info["AC"] = tuple(1 for _ in alleles[1:])
            record.info["AN"] = 4
            for sample_name in header.samples:
                record.samples[sample_name]["GT"] = (0, 1)
                record.samples[sample_name]["DP"] = 30 + i
                record.samples[sample_name]["AD"] = tuple(10 for _ in alleles)
            if kind == "preannotated" and i == 1:
                record.info["CSQ"] = "preexisting"
            inp.write(record)
            annotation = record.copy()
            annotation.info["CSQ"] = f"C|consequence_{i}|" + "x" * payload
            annotation.info["SUP"] = f"supplement_{i}"
            annotation.info["FLAG"] = True
            annotation.info["AF"] = tuple(0.5 for _ in alleles[1:])
            complete.write(annotation)
            if hit_every and i % hit_every == 0:
                cached.write(annotation)
    for path in (sample, cache, truth):
        checked("bcftools", "index", path)
    params = directory / "params.yaml"
    params.write_text(
        yaml.safe_dump(
            {
                "bcftools_cmd": "bcftools",
                "annotation_tool_cmd": "bcftools",
                "temp_dir": str(directory),
                "threads": 6,
                "genome_build": "GRCh38",
            }
        )
    )
    recipe = directory / "annotation.yaml"
    recipe.write_text(
        yaml.safe_dump(
            {
                "annotation_cmd": f"bcftools annotate -a {truth} ${{INPUT_BCF}} -c INFO -Ob -W -o ${{OUTPUT_BCF}}",
                "must_contain_info_tag": "CSQ",
                "required_tool_version": "test",
                "genome_build": "GRCh38",
            }
        )
    )
    return sample, cache, params, recipe


def run_case(
    root, files, engine, *, workers=1, preserve=False, skip_split=False, lookup="auto"
):
    sample, cache, params, recipe = files
    output = root / engine / "output.bcf"
    output.parent.mkdir(parents=True, exist_ok=True)
    manager = WorkflowManager(
        input_file=sample,
        output_dir=output.parent,
        name="test",
        params_file=params,
        anno_config_file=recipe,
        output_file=output,
    )
    manager.run(
        "annotate",
        db_bcf=cache,
        engine=engine,
        engine_workers=workers if engine in {"fused", "native"} else 0,
        engine_lookup=lookup,
        preserve_unannotated=preserve,
        skip_split_multiallelic=skip_split,
    )
    return output, manager


@pytest.mark.parametrize("hits", [0, 1, 2, 10])
@pytest.mark.parametrize("workers", [1, 2])
def test_complete_record_parity(native, tmp_path, hits, workers):
    files = dataset(tmp_path / "data", hit_every=hits)
    legacy, old = run_case(tmp_path, files, "legacy")
    fused, new = run_case(tmp_path, files, "fused", workers=workers)
    assert strict_semantic_compare(fused, legacy)["semantic_pass"]
    assert (
        checked("bcftools", "view", "-H", fused).stdout
        == checked("bcftools", "view", "-H", legacy).stdout
    )
    for key in (
        "input_variants",
        "output_variants",
        "missing_variants",
        "missing_annotated",
        "dropped_variants",
    ):
        assert new.last_run_stats[key] == old.last_run_stats[key]
    assert new.last_run_stats["engine"] == "fused"


@pytest.mark.parametrize("kind", ["normal", "multiallelic", "spanning", "preannotated"])
@pytest.mark.parametrize("preserve", [False, True])
def test_variant_representation_and_fields(native, tmp_path, kind, preserve):
    files = dataset(tmp_path / "data", kind=kind)
    legacy, _ = run_case(tmp_path, files, "legacy", preserve=preserve)
    fused, _ = run_case(tmp_path, files, "fused", preserve=preserve)
    assert strict_semantic_compare(fused, legacy)["semantic_pass"]


def test_empty(native, tmp_path):
    files = dataset(tmp_path / "data", records=0)
    legacy, _ = run_case(tmp_path, files, "legacy")
    fused, manager = run_case(tmp_path, files, "fused")
    assert strict_semantic_compare(fused, legacy)["semantic_pass"]
    assert manager.last_run_stats["output_variants"] == 0


@pytest.mark.parametrize("lookup", ["auto", "regions", "contigs"])
@pytest.mark.parametrize("workers", [1, 2])
def test_lookup_modes_and_boundaries(native, tmp_path, lookup, workers):
    files = dataset(tmp_path / "data")
    # Same-coordinate groups and long deletions crossing window/task boundaries
    # must not be duplicated by record-overlap selection.
    positions = [100, 100, 16383, 16384, 16384, 32767, 32768, 40000, 49999, 50000]
    for path in (files[0], files[1], tmp_path / "data/truth.bcf"):
        replacement = path.with_name("replacement.bcf")
        with (
            pysam.VariantFile(str(path)) as inp,
            pysam.VariantFile(str(replacement), "wb", header=inp.header) as out,
        ):
            for record in inp:
                i = int(record.id.removeprefix("id-")) % 10
                record.start = positions[i]
                record.alleles = (
                    ("A" * 20000, "A") if i == 2 else ("A", "G" if i in {1, 4} else "C")
                )
                out.write(record)
        replacement.replace(path)
        checked("bcftools", "index", "-f", path)
    legacy, _ = run_case(tmp_path, files, "legacy", preserve=True)
    fused, manager = run_case(
        tmp_path, files, "fused", workers=workers, preserve=True, lookup=lookup
    )
    assert strict_semantic_compare(fused, legacy)["semantic_pass"]
    # bcftools' overlays may reorder different alleles within one coordinate.
    # The strict comparator checks complete records, duplicates and all fields;
    # byte-identical order within a coordinate is not the sorted-output contract.
    assert manager.last_run_stats["engine_lookup"] == lookup
    # Every overlay, even with a single worker, retains explicit selection.
    for command in (manager.work_dir / "fused").glob("shard-*/route.sh"):
        assert "-R regions.bed --regions-overlap 0" in command.read_text()
        assert "--threads " in command.read_text()


@pytest.mark.parametrize("lookup", ["regions", "contigs"])
def test_missing_cache_contig(native, tmp_path, lookup):
    files = dataset(tmp_path / "data")
    replacement = tmp_path / "short-cache.bcf"
    with (pysam.VariantFile(str(files[1])) as inp,):
        header = pysam.VariantHeader()
        for declaration in inp.header.records:
            if declaration.type == "CONTIG" and declaration["ID"] == "1":
                continue
            header.add_line(str(declaration))
        for name in inp.header.samples:
            header.add_sample(name)
        with pysam.VariantFile(str(replacement), "wb", header=header) as out:
            for record in inp:
                if record.contig == "2":
                    record.translate(header)
                    out.write(record)
    replacement.replace(files[1])
    checked("bcftools", "index", "-f", files[1])
    legacy, _ = run_case(tmp_path, files, "legacy")
    fused, _ = run_case(tmp_path, files, "fused", lookup=lookup)
    assert strict_semantic_compare(fused, legacy)["semantic_pass"]


def test_cache_work_planner_splits_single_contig():
    from vcfcache.database.adaptive import plan_tasks

    windows = [
        dict(
            contig="chr1",
            start=i * 16384,
            end=(i + 1) * 16384,
            records=10,
            estimated_cache_bytes=2**20,
        )
        for i in range(32)
    ]
    metadata = {
        "windows": windows,
        "contigs": [
            dict(name="chr1", estimated_cache_bytes=2**25, requested_cache_bytes=2**25)
        ],
    }
    tasks, strategies = plan_tasks(
        [(Path("normalized.bcf"), metadata)], {"chr1": 32 * 16384}, 4, "auto"
    )
    assert len(tasks) >= 4
    assert strategies[0]["strategy"] == "contigs"
    regions = [region for task in tasks for region in task["regions"]]
    assert regions[0][1] == 0 and regions[-1][2] == 32 * 16384
    assert all(
        left[2] == right[1] for left, right in zip(regions, regions[1:], strict=False)
    )


def test_sparse_planner_uses_cache_bytes_not_sample_counts():
    from vcfcache.database.adaptive import plan_tasks

    metadata = {
        "windows": [
            dict(
                contig="chr1",
                start=0,
                end=16384,
                records=250,
                estimated_cache_bytes=2**20,
            )
        ],
        "contigs": [
            dict(name="chr1", estimated_cache_bytes=2**30, requested_cache_bytes=2**20)
        ],
    }
    tasks, strategies = plan_tasks(
        [(Path("normalized.bcf"), metadata)], {"chr1": 10**8}, 1, "auto"
    )
    assert strategies[0]["strategy"] == "regions"
    assert tasks[0]["regions"] == [("chr1", 0, 16384)]


@pytest.mark.parametrize(
    "input_records,cache_records,expected",
    [(10000, 9000, 8), (10000, 100000, 4), (250, 9000, 1)],
)
def test_worker_policy_reserves_codecs_only_when_useful(
    input_records, cache_records, expected
):
    from vcfcache.database.adaptive import select_workers

    metadata = {
        "windows": [{"records": input_records}],
        "contigs": [
            dict(
                estimated_cache_bytes=100,
                requested_cache_bytes=100,
                cache_records=cache_records,
            )
        ],
    }
    assert select_workers([(None, metadata)], 16, "auto") == expected


@pytest.mark.parametrize("fails", [False, True])
def test_cpu_budget_restores_callers_affinity(monkeypatch, fails):
    import os
    from types import SimpleNamespace

    import vcfcache.database.fused as module

    original = os.sched_getaffinity(0)

    def workflow(*args):
        assert os.sched_getaffinity(0) == {min(original)}
        if fails:
            raise RuntimeError("test failure")
        return "ok"

    monkeypatch.setattr(module, "_run_fused", workflow)
    manager = SimpleNamespace(params_file_content={"threads": 1})
    if fails:
        with pytest.raises(RuntimeError, match="test failure"):
            module.run_fused(manager, None)
    else:
        assert module.run_fused(manager, None) == "ok"
    assert os.sched_getaffinity(0) == original


def test_benchmark_claim_preserves_prior_output(tmp_path):
    from benchmarks.adaptive_routing import claim_directory

    destination = tmp_path / "pre-created"
    destination.mkdir()
    claim_directory(destination)
    marker = (destination / ".benchmark_owner").read_text()
    with pytest.raises(FileExistsError):
        claim_directory(destination)
    assert (destination / ".benchmark_owner").read_text() == marker


def test_single_contig_multiple_tasks_annotates_once(native, tmp_path, monkeypatch):
    import vcfcache.database.fused as module

    files = dataset(tmp_path / "data")
    for path in (files[0], files[1], tmp_path / "data/truth.bcf"):
        replacement = path.with_name("replacement.bcf")
        with (
            pysam.VariantFile(str(path)) as inp,
            pysam.VariantFile(str(replacement), "wb", header=inp.header) as out,
        ):
            for record in inp:
                record.contig = "2"
                record.start = int(record.id.removeprefix("id-")) * 4_000_000
                out.write(record)
        replacement.replace(path)
        checked("bcftools", "index", "-f", path)
    planner = module.plan_tasks

    def expensive_cache(prepared, lengths, workers, lookup):
        # Inflate ONLY work estimates to exercise real task splitting on a tiny
        # fixture. The records, cache lookup and annotation recipe remain real.
        for _, metadata in prepared:
            for contig in metadata["contigs"]:
                contig["estimated_cache_bytes"] = 2**30
                contig["requested_cache_bytes"] = 2**30
        return planner(prepared, lengths, workers, lookup)

    monkeypatch.setattr(module, "plan_tasks", expensive_cache)
    legacy, _ = run_case(tmp_path, files, "legacy")
    fused, manager = run_case(tmp_path, files, "fused", workers=2)
    assert strict_semantic_compare(fused, legacy)["semantic_pass"]
    assert manager.last_run_stats["lookup_tasks"] >= 4
    assert (
        sum(
            row["stage"] == "annotate-misses"
            for row in manager.last_run_stats["stage_timings"]
        )
        == 1
    )


@pytest.mark.parametrize("hits", [0, 1, 2])
def test_experimental_store_candidates(native, tmp_path, hits):
    if "store-build" not in subprocess.check_output([native, "--help"], text=True):
        pytest.skip("helper built without experimental-store feature")
    files = dataset(tmp_path / "data", hit_every=hits)
    store = tmp_path / "store"
    checked(native, "store-build", "--cache", files[1], "--output", store)
    selected = tmp_path / "selected.bcf"
    with selected.open("wb") as stream:
        subprocess.run(
            [native, "store-query", "--store", str(store), "--input", str(files[0])],
            stdout=stream,
            check=True,
        )
    indexed = tmp_path / "selected-indexed.bcf"
    checked("bcftools", "view", selected, "-Ob", "-W", "-o", indexed)
    assert strict_semantic_compare(indexed, files[1])["semantic_pass"]
    # Candidate extraction never substitutes target GT or other sample fields.
    assert json.loads((store / "manifest.json").read_text())["experimental"]
    again = subprocess.run(
        [native, "store-build", "--cache", str(files[1]), "--output", str(store)],
        capture_output=True,
    )
    assert again.returncode != 0  # never overwrite an existing store


def test_cli_rejects_lookup_for_legacy(tmp_path):
    import sys

    response = subprocess.run(
        [
            sys.executable,
            "-m",
            "vcfcache.cli",
            "annotate",
            "-a",
            "cache",
            "-i",
            "x.bcf",
            "-o",
            "y.bcf",
            "--engine-lookup",
            "regions",
        ],
        capture_output=True,
        text=True,
    )
    assert response.returncode == 2
    assert "requires --engine fused" in response.stderr


def test_multiallelic_replay(native, tmp_path):
    files = dataset(tmp_path / "data", kind="multiallelic")
    legacy, _ = run_case(tmp_path, files, "legacy", skip_split=True)
    fused, manager = run_case(tmp_path, files, "fused", skip_split=True)
    assert strict_semantic_compare(fused, legacy)["semantic_pass"]
    assert manager.last_run_stats["engine"] == "legacy"
    assert "multiallelic" in manager.last_run_stats["fallback_reason"]


def test_failed_recipe_does_not_publish(native, tmp_path):
    files = dataset(tmp_path / "data")
    recipe = yaml.safe_load(files[3].read_text())
    recipe["annotation_cmd"] = "printf 'deliberate failure\\n' >&2; exit 17"
    files[3].write_text(yaml.safe_dump(recipe))
    with pytest.raises(RuntimeError, match="deliberate failure"):
        run_case(tmp_path, files, "fused", workers=2)
    assert not (tmp_path / "fused/output.bcf").exists()
    assert not (tmp_path / "fused/output.bcf.csi").exists()


def test_bad_executable_override(monkeypatch, tmp_path):
    monkeypatch.setenv("VCFCACHE_ENGINE", str(tmp_path / "missing"))
    with pytest.raises(FileNotFoundError, match="not executable"):
        native_binary()


def test_incompatible_info_header_rejected(tmp_path):
    files = dataset(tmp_path / "data")
    bad = tmp_path / "bad.vcf"
    bad.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=2>\n"
        '##INFO=<ID=CSQ,Number=1,Type=Integer,Description="wrong type">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )
    with pytest.raises(ValueError, match="Incompatible INFO header for CSQ"):
        validate_info_headers(files[0], bad)


def test_bad_routing_manifest(native, tmp_path):
    files = dataset(tmp_path / "data")
    order = tmp_path / "bad.bin"
    order.write_bytes(b"\x02")
    result = subprocess.run(
        [
            native,
            "assemble",
            "--hits",
            str(files[0]),
            "--misses",
            str(files[0]),
            "--order",
            str(order),
            "--tag",
            "CSQ",
            "--stats",
            str(tmp_path / "stats.json"),
        ],
        capture_output=True,
    )
    assert result.returncode != 0
    assert b"invalid routing manifest" in result.stderr
    assert not (tmp_path / "stats.json").exists()


@pytest.mark.parametrize("preserve", [False, True])
@pytest.mark.parametrize(
    "recipe_kind", ["drop", "strip", "reorder", "extra-header", "change-format"]
)
def test_recipe_behavior(native, tmp_path, preserve, recipe_kind):
    files = dataset(tmp_path / "data")
    recipe = yaml.safe_load(files[3].read_text())
    command = recipe["annotation_cmd"].replace("-Ob -W -o ${OUTPUT_BCF}", "-Ou")
    if recipe_kind == "drop":
        command += " | bcftools view -e 'POS%4==0' -Ob -W -o ${OUTPUT_BCF}"
    elif recipe_kind == "strip":
        command += " | bcftools annotate -x INFO/CSQ -Ob -W -o ${OUTPUT_BCF}"
        # Keep the declaration so the legacy required-tag expression remains valid.
        command = recipe["annotation_cmd"].replace("-Ob -W -o ${OUTPUT_BCF}", "-Ou")
        command += " | bcftools +setGT -Ou -- -t a -n . | bcftools view -Ob -W -o ${OUTPUT_BCF}"
    elif recipe_kind == "reorder":
        command += " | bcftools view -Ov | awk '/^#/ {print; next} {line[++n]=$0} END {for(i=n;i>0;i--) print line[i]}' | bcftools sort -Ob -o ${OUTPUT_BCF}; bcftools index ${OUTPUT_BCF}"
    elif recipe_kind == "extra-header":
        new_header = tmp_path / "header.txt"
        new_header.write_text(
            '##INFO=<ID=NEW,Number=1,Type=String,Description="new field">\n'
        )
        command += f" | bcftools annotate -h {new_header} -Ob -W -o ${{OUTPUT_BCF}}"
    else:
        command += " | bcftools annotate -x FORMAT/DP -Ob -W -o ${OUTPUT_BCF}"
    recipe["annotation_cmd"] = command
    files[3].write_text(yaml.safe_dump(recipe))
    legacy, old = run_case(tmp_path, files, "legacy", preserve=preserve)
    fused, new = run_case(tmp_path, files, "fused", workers=2, preserve=preserve)
    assert strict_semantic_compare(fused, legacy)["semantic_pass"]
    for key in (
        "input_variants",
        "output_variants",
        "missing_variants",
        "missing_annotated",
        "dropped_variants",
    ):
        assert new.last_run_stats[key] == old.last_run_stats[key]


def test_vector_missing_predicate(native, tmp_path):
    files = dataset(tmp_path / "data", hit_every=0)
    truth = tmp_path / "data/vector.bcf"
    values = [None, ".", "", ".,a", "a,.", "a,", ",a", ",", "a,,b", "a"]
    with (
        pysam.VariantFile(str(files[0])) as source,
        pysam.VariantFile(str(truth), "wb", header=source.header) as out,
    ):
        for index, record in enumerate(source):
            if values[index % len(values)] is not None:
                record.info["CSQ"] = values[index % len(values)].split(",")
            out.write(record)
    checked("bcftools", "index", truth)
    command = [
        native,
        "route",
        "--hits",
        str(tmp_path / "hits.bcf"),
        "--misses",
        str(tmp_path / "miss.bcf"),
        "--order",
        str(tmp_path / "order.bin"),
        "--tag",
        "CSQ",
        "--stats",
        str(tmp_path / "stats.json"),
    ]
    stream = checked("bcftools", "view", truth, "-Ou").stdout
    subprocess.run(command, input=stream, check=True, capture_output=True)
    expected = checked("bcftools", "view", "-H", "-i", 'INFO/CSQ!=""', truth).stdout
    assert checked("bcftools", "view", "-H", tmp_path / "hits.bcf").stdout == expected
    assert json.loads((tmp_path / "stats.json").read_text())["input_variants"] == 20


@pytest.mark.parametrize("preserve", [False, True])
def test_missing_required_values(native, tmp_path, preserve):
    files = dataset(tmp_path / "data")
    replacement = tmp_path / "data/incomplete-truth.bcf"
    with (
        pysam.VariantFile(str(tmp_path / "data/truth.bcf")) as source,
        pysam.VariantFile(str(replacement), "wb", header=source.header) as out,
    ):
        for record in source:
            if record.pos % 4:
                del record.info["CSQ"]
            out.write(record)
    checked("bcftools", "index", replacement)
    recipe = yaml.safe_load(files[3].read_text())
    recipe["annotation_cmd"] = recipe["annotation_cmd"].replace(
        str(tmp_path / "data/truth.bcf"), str(replacement)
    )
    files[3].write_text(yaml.safe_dump(recipe))
    legacy, old = run_case(tmp_path, files, "legacy", preserve=preserve)
    fused, new = run_case(tmp_path, files, "fused", workers=2, preserve=preserve)
    assert strict_semantic_compare(fused, legacy)["semantic_pass"]
    assert (
        new.last_run_stats["missing_annotated"]
        == old.last_run_stats["missing_annotated"]
    )


def test_actual_unsorted_recipe_is_sorted(native, tmp_path):
    files = dataset(tmp_path / "data")
    recipe = yaml.safe_load(files[3].read_text())
    command = recipe["annotation_cmd"].replace("-Ob -W -o ${OUTPUT_BCF}", "-Ov")
    recipe["annotation_cmd"] = (
        command
        + " | awk '/^#/ {print; next} {line[++n]=$0} END {for(i=n;i>0;i--) print line[i]}' | bcftools view -Ob -o ${OUTPUT_BCF}"
    )
    files[3].write_text(yaml.safe_dump(recipe))
    fused, manager = run_case(tmp_path, files, "fused", workers=2)
    assert int(checked("bcftools", "index", "-n", fused).stdout) == 20
    assert any(
        row["stage"] == "sort-annotator-output"
        for row in manager.last_run_stats["stage_timings"]
    )


def test_duplicate_keys_replay(native, tmp_path):
    files = dataset(tmp_path / "data", hit_every=0)
    replacement = tmp_path / "data/duplicates.bcf"
    with (
        pysam.VariantFile(str(files[0])) as source,
        pysam.VariantFile(str(replacement), "wb", header=source.header) as out,
    ):
        for index, record in enumerate(source):
            out.write(record)
            if index == 0:
                duplicate = record.copy()
                duplicate.info["CSQ"] = "preexisting"
                out.write(duplicate)
    checked("bcftools", "index", replacement)
    files = replacement, *files[1:]
    legacy, _ = run_case(tmp_path, files, "legacy")
    fused, manager = run_case(tmp_path, files, "fused", workers=2)
    assert strict_semantic_compare(fused, legacy)["semantic_pass"]
    assert manager.last_run_stats["engine"] == "legacy"
    assert "duplicate" in manager.last_run_stats["fallback_reason"]


def test_reordered_info_dictionary(native, tmp_path):
    files = dataset(tmp_path / "data")
    replacement = tmp_path / "data/reordered.bcf"
    with pysam.VariantFile(str(files[1])) as source:
        header = pysam.VariantHeader()
        for name, contig in source.header.contigs.items():
            header.contigs.add(name, length=contig.length)
        for kind in ("INFO", "FORMAT"):
            for declaration in reversed(list(source.header.records)):
                if declaration.type == kind:
                    header.add_line(str(declaration))
        for sample in source.header.samples:
            header.add_sample(sample)
        with pysam.VariantFile(str(replacement), "wb", header=header) as out:
            for record in source:
                record.translate(header)
                out.write(record)
    checked("bcftools", "index", replacement)
    files = files[0], replacement, *files[2:]
    legacy, _ = run_case(tmp_path, files, "legacy")
    fused, _ = run_case(tmp_path, files, "fused", workers=2)
    assert strict_semantic_compare(fused, legacy)["semantic_pass"]


def test_numeric_required_tag_preserves_legacy_error(native, tmp_path):
    files = dataset(tmp_path / "data")
    recipe = yaml.safe_load(files[3].read_text())
    recipe["must_contain_info_tag"] = "AF"
    files[3].write_text(yaml.safe_dump(recipe))
    # The legacy empty-string predicate does not accept numeric tags either.
    for engine in ("legacy", "fused"):
        with pytest.raises(RuntimeError, match="failed"):
            run_case(tmp_path, files, engine)
        assert not (tmp_path / engine / "output.bcf").exists()


@pytest.mark.parametrize("statistics", ["light", "full", "none"])
def test_native_stdout_is_valid_bcf(native, tmp_path, statistics):
    files = dataset(tmp_path / "data")
    output = tmp_path / "stdout.bcf"
    # Exercise the public CLI, including construction/status/logging messages.
    cache_root = tmp_path / "cache-root"
    cache_dir = cache_root / "cache/test"
    cache_dir.mkdir(parents=True)
    (cache_root / "blueprint").mkdir()
    (cache_root / "workflow").mkdir()
    (cache_root / "workflow/init.yaml").write_text("{}\n")
    import shutil

    for source, name in (
        (files[1], "vcfcache_annotated.bcf"),
        (files[1].with_suffix(".bcf.csi"), "vcfcache_annotated.bcf.csi"),
        (files[2], "params.snapshot.yaml"),
        (files[3], "annotation.yaml"),
    ):
        shutil.copyfile(source, cache_dir / name)
    import sys

    command = [
        sys.executable,
        "-m",
        "vcfcache",
        "annotate",
        "--debug",
        "-a",
        str(cache_dir),
        "-i",
        str(files[0]),
        "-o",
        "-",
        "--engine",
        "fused",
        "--statistics",
        statistics,
    ]
    if statistics != "none":
        command.extend(["--stats-dir", str(tmp_path / "statistics")])
    result = subprocess.run(command, capture_output=True)
    assert result.returncode == 0, result.stderr.decode(errors="replace")
    output.write_bytes(result.stdout)
    checked("bcftools", "index", output)
    assert int(checked("bcftools", "index", "-n", output).stdout) == 20
    if statistics != "none":
        stats = yaml.safe_load(
            (tmp_path / "statistics/sample_vcstats/compare_stats.yaml").read_text()
        )
        assert stats["engine"] == "fused"
        assert stats["variant_counts"]["total_output"] == 20
        assert stats["engine_lookup"] == "auto"
        assert stats["lookup_regions"] >= 1
        assert sum(item["cache_records"] for item in stats["lookup_strategies"]) == 10
        assert len(stats["phase_timings"]) == 3
        assert stats["engine_wall_seconds"] > 0


def test_pipeline_failure_does_not_deadlock(tmp_path):
    import logging

    pipeline = Pipelines(logging.getLogger("test"))
    with pytest.raises(RuntimeError, match="failed"):
        pipeline.run("yes x | (exit 9)", tmp_path, "broken-pipe")
    assert not pipeline.active


def test_pipeline_cancel_kills_group(tmp_path):
    import logging
    import threading

    pipeline = Pipelines(logging.getLogger("test"))
    errors = []
    started = threading.Event()

    def execute():
        started.set()
        try:
            pipeline.run("sleep 60", tmp_path, "cancel")
        except RuntimeError as error:
            errors.append(error)

    thread = threading.Thread(target=execute)
    thread.start()
    assert started.wait(2)
    # Wait until creation is visible under the same lock used by cancellation.
    import time

    deadline = time.monotonic() + 2
    while not pipeline.active and time.monotonic() < deadline:
        time.sleep(0.01)
    pipeline.stop()
    thread.join(3)
    assert not thread.is_alive()
    assert errors
    assert not pipeline.active
