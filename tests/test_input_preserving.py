"""Real bcftools checks for the immutable-record annotation contract."""

import hashlib
import subprocess

import pytest
import yaml

from vcfcache.database import preserving
from vcfcache.database.preserving import HIT, identity
from vcfcache.database.workflow_manager import BcftoolsCommand, WorkflowManager

HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=MOCK_ANNO,Number=1,Type=String,Description="Test annotation">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
"""
ROWS = [
    "chr1\t10\ta\tA\tC\t50\tPASS\t.\tGT:DP\t0/1:12",
    "chr1\t20\tb\tA\tC,G\t51\tPASS\t.\tGT:DP\t1/2:20",
    "chr1\t30\tc\tA\t*\t52\tPASS\t.\tGT:DP\t0/1:30",
    "chr1\t40\td\tT\tC\t53\tPASS\tMOCK_ANNO=old\tGT:DP\t0/1:40",
    "chr1\t50\te\tA\tC,G\t54\tPASS\t.\tGT:DP\t1|2:50",
    "chr1\t60\tf\tA\tC\t55\tPASS\t.\tGT:DP\t0/1:60",
    "chr1\t60\tg\tA\tC\t56\tPASS\t.\tGT:DP\t1/1:61",
    "chr1\t70\th\tCAA\tCA\t57\tPASS\t.\tGT:DP\t0/1:70",
]
RECIPE = """echo '##INFO=<ID=MOCK_ANNO,Number=1,Type=String,Description="Test annotation">' > extra.hdr
${params.bcftools_cmd} query -i 'ALT!="*"' -f '%CHROM\\t%POS\\ttest\\n' ${INPUT_BCF} > extra.tsv
bgzip -c extra.tsv > extra.tsv.gz
tabix -s1 -b2 -e2 extra.tsv.gz
${params.bcftools_cmd} annotate -a extra.tsv.gz -h extra.hdr -c CHROM,POS,INFO/MOCK_ANNO ${INPUT_BCF} -Ob -o ${OUTPUT_BCF} -W
"""


def bcf(folder, name, rows, header=HEADER):
    """Write test-only VCF then index its BCF representation."""
    vcf = folder / f"{name}.vcf"
    vcf.write_text(header + "\n".join(rows) + ("\n" if rows else ""))
    output = folder / f"{name}.bcf"
    subprocess.run(
        ["bcftools", "view", str(vcf), "-Ob", "-o", str(output), "-W"], check=True
    )
    return output


def manager(folder, source, recipe=RECIPE):
    """Use the normal manager, not a stub of bcftools behaviour."""
    folder.mkdir()
    params = folder / "params.yaml"
    params.write_text(
        yaml.safe_dump(
            dict(
                annotation_tool_cmd="bcftools",
                bcftools_cmd="bcftools",
                threads=2,
                genome_build="GRCh38",
                temp_dir=str(folder),
            )
        )
    )
    config = folder / "annotation.yaml"
    config.write_text(
        yaml.safe_dump(
            dict(
                annotation_cmd=recipe,
                must_contain_info_tag="MOCK_ANNO",
                genome_build="GRCh38",
                required_tool_version="test",
            )
        )
    )
    return WorkflowManager(
        input_file=source,
        output_dir=folder,
        name="test",
        anno_config_file=config,
        params_file=params,
        output_file=folder / "result.bcf",
        verbosity=0,
    )


def records(path):
    """Record lines only; headers legitimately differ between routes."""
    return subprocess.check_output(["bcftools", "view", "-H", str(path)])


@pytest.fixture
def inputs(tmp_path):
    source = bcf(tmp_path, "input", [r for i, r in enumerate(ROWS) if i != 6])
    rows = [ROWS[i].replace("\t.\tGT:", "\tMOCK_ANNO=test\tGT:") for i in (0, 5)]
    rows.insert(
        1, ROWS[1].replace("C,G", "C").replace("\t.\tGT:", "\tMOCK_ANNO=WRONG\tGT:")
    )
    rows.insert(
        2, ROWS[4].replace("C,G", "G,C").replace("\t.\tGT:", "\tMOCK_ANNO=WRONG\tGT:")
    )
    return source, bcf(tmp_path, "cache", rows)


@pytest.mark.parametrize("ratio", [0, 10**9])
def test_exact_preserves_all_records_and_annotates_partial_hits(
    tmp_path, inputs, monkeypatch, ratio
):
    source, cache = inputs
    before = hashlib.sha256(source.read_bytes()).hexdigest()
    direct = manager(tmp_path / "direct", source)
    direct.run("annotate-nocache")
    monkeypatch.setattr(preserving, "WINDOW_RATIO", ratio)
    cached = manager(tmp_path / "cached", source)
    cached.run("annotate", db_bcf=cache)
    assert identity("bcftools", source) == identity("bcftools", cached.output_file)

    assert records(cached.output_file) == records(direct.output_file)
    assert cached.last_run_stats["input_variants"] == 7
    assert cached.last_run_stats["missing_variants"] == 5
    assert cached.last_run_stats["cache_access"] == ("contigs" if ratio else "windows")
    assert HIT.encode() not in subprocess.check_output(
        ["bcftools", "view", "-h", str(cached.output_file)]
    )
    assert hashlib.sha256(source.read_bytes()).hexdigest() == before


@pytest.mark.parametrize("empty", ["input", "cache"])
def test_empty_inputs_and_empty_cache(tmp_path, inputs, empty):
    source, cache = inputs
    if empty == "input":
        source = bcf(tmp_path, "empty_input", [])
    else:
        cache = bcf(tmp_path, "empty_cache", [])
    run = manager(tmp_path / "run", source)
    run.run("annotate", db_bcf=cache)
    assert run.last_run_stats["output_variants"] == (0 if empty == "input" else 7)


def test_duplicate_keys_bypass_the_cache_and_equal_direct(tmp_path, inputs):
    _, cache = inputs
    source = bcf(tmp_path, "duplicates", ROWS)
    direct = manager(tmp_path / "direct", source)
    direct.run("annotate-nocache")
    run = manager(tmp_path / "run", source)
    run.run("annotate", db_bcf=cache)
    assert identity("bcftools", source) == identity("bcftools", run.output_file)
    assert records(run.output_file) == records(direct.output_file)
    assert run.last_run_stats["output_variants"] == 8
    assert run.last_run_stats["bypassed_shards"] == 1
    assert run.last_run_stats["missing_variants"] == 8


def test_all_hits_skip_recipe(tmp_path, inputs):
    _, cache = inputs
    source = bcf(tmp_path, "all_hits", [ROWS[0], ROWS[5]])
    run = manager(tmp_path / "run", source, recipe="exit 17")
    run.run("annotate", db_bcf=cache)
    assert run.last_run_stats["missing_variants"] == 0


@pytest.mark.parametrize(
    "recipe",
    [
        "${params.bcftools_cmd} view -i 'POS=10' ${INPUT_BCF} -Ob -o ${OUTPUT_BCF} -W",
        "${params.bcftools_cmd} annotate --set-id CHANGED ${INPUT_BCF} -Ob -o ${OUTPUT_BCF} -W",
        "${params.bcftools_cmd} +setGT ${INPUT_BCF} -Ob -o ${OUTPUT_BCF} -W -- -t a -n 0",
    ],
)
def test_rejects_changed_recipe_output_and_preserves_existing_target(
    tmp_path, inputs, recipe
):
    source, _ = inputs
    run = manager(tmp_path / "run", source, recipe=recipe)
    run.output_file.write_bytes(b"previous deliverable")
    with pytest.raises(RuntimeError, match="(count changed|changed variant identity)"):
        run.run("annotate-nocache")
    assert run.output_file.read_bytes() == b"previous deliverable"


def test_pipeline_checks_producer_failure(tmp_path, inputs):
    run = manager(tmp_path / "run", inputs[0])
    with pytest.raises(subprocess.CalledProcessError):
        BcftoolsCommand("false | true", run.logger, tmp_path).run()


def test_reserved_hit_marker_is_rejected(tmp_path, inputs):
    source, cache = inputs
    run = manager(tmp_path / "run", source)
    marked = tmp_path / "marker.bcf"
    subprocess.run(
        [
            "bcftools",
            "annotate",
            "-a",
            str(cache),
            "-m",
            HIT,
            str(source),
            "-Ob",
            "-o",
            str(marked),
            "-W",
        ],
        check=True,
    )
    run.input_file = marked
    with pytest.raises(ValueError, match="Reserved INFO"):
        run.run("annotate", db_bcf=cache)


def test_binary_stdout_is_a_valid_bcf(tmp_path, inputs, capfdbinary):
    """No completion/logging text may corrupt the binary output stream."""
    source, cache = inputs
    run = manager(tmp_path / "stdout", source)
    run.output_file = "-"
    run.run("annotate", db_bcf=cache)
    captured = capfdbinary.readouterr()
    decoded = subprocess.run(
        ["bcftools", "view", "-H", "-"],
        input=captured.out,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )
    assert len(decoded.stdout.splitlines()) == 7


def test_output_cannot_overwrite_input(tmp_path, inputs):
    source, cache = inputs
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    run = manager(tmp_path / "same_path", source)
    run.output_file = source
    with pytest.raises(ValueError, match="overwrite the input"):
        run.run("annotate", db_bcf=cache)
    assert hashlib.sha256(source.read_bytes()).hexdigest() == digest


WIDE = HEADER.replace(
    "##contig=<ID=chr1,length=1000>",
    "".join(f"##contig=<ID={c},length=1000>\n" for c in ("chr1", "chr2", "chr3"))[:-1],
)
SITES = WIDE.replace("\tFORMAT\tSAMPLE", "")
UNUSUAL = [
    ("chr1", 10, "A", "C", "."),  # plain hit
    ("chr1", 20, "A", "C,G", "."),  # multiallelic; cache holds A>C
    ("chr1", 30, "A", "*", "."),  # spanning deletion
    ("chr1", 40, "T", "C", "MOCK_ANNO=old"),  # annotated already, not cached
    ("chr1", 45, "T", "C", "MOCK_ANNO=old"),  # annotated already, cached
    ("chr1", 70, "CAA", "CA", "."),  # indel absent from the cache
    ("chr2", 10, "A", "C", "."),  # plain hit
    ("chr2", 20, "a", "g", "."),  # lower case alleles, cached in upper case
    ("chr2", 30, "A", "<DEL>", "."),  # symbolic allele
    ("chr2", 40, "A", ".", "."),  # no alternate allele
    ("chr2", 50, "A", "C", "."),  # cached without annotation
    ("chr3", 60, "A", "C", "."),  # duplicate key: would be hit ...
    ("chr3", 60, "A", "C", "."),  # ... then miss
    ("chr3", 80, "G", "T", "."),  # duplicate key absent from the cache
    ("chr3", 80, "G", "T", "."),
    ("chr3", 90, "A", "C", "."),  # cached, but its shard bypasses the cache
    ("chrUn", 5, "A", "C", "."),  # contig unknown to the cache
]
CACHED = [
    ("chr1", 10, "A", "C", "MOCK_ANNO=test"),
    ("chr1", 20, "A", "C", "MOCK_ANNO=WRONG"),
    ("chr1", 45, "T", "C", "MOCK_ANNO=test"),
    ("chr2", 10, "A", "C", "MOCK_ANNO=test"),
    ("chr2", 20, "A", "G", "MOCK_ANNO=test"),
    ("chr2", 50, "A", "C", "."),
    ("chr3", 60, "A", "C", "MOCK_ANNO=test"),
    ("chr3", 90, "A", "C", "MOCK_ANNO=test"),
]


@pytest.fixture
def unusual(tmp_path):
    """Every input irregularity the cache must hand to the recipe unchanged."""
    header = WIDE.replace("##INFO", "##contig=<ID=chrUn,length=1000>\n##INFO", 1)
    source = bcf(
        tmp_path,
        "unusual",
        [
            f"{c}\t{p}\tid{n}\t{r}\t{a}\t50\tPASS\t{i}\tGT:DP\t0/1:{n}"
            for n, (c, p, r, a, i) in enumerate(UNUSUAL)
        ],
        header,
    )
    cache = bcf(
        tmp_path,
        "wide_cache",
        [f"{c}\t{p}\t.\t{r}\t{a}\t.\t.\t{i}" for c, p, r, a, i in CACHED],
        SITES,
    )
    return source, cache


@pytest.mark.parametrize("ratio", [0, 10**9])
@pytest.mark.parametrize("per_thread", [4, 1])
def test_unusual_input_equals_direct_annotation(
    tmp_path, unusual, monkeypatch, ratio, per_thread
):
    source, cache = unusual
    monkeypatch.setattr(preserving, "WINDOW_RATIO", ratio)
    monkeypatch.setattr(preserving, "SHARDS_PER_THREAD", per_thread)
    direct = manager(tmp_path / "direct", source)
    direct.run("annotate-nocache")
    cached = manager(tmp_path / "cached", source)
    cached.run("annotate", db_bcf=cache)
    assert records(cached.output_file) == records(direct.output_file)
    assert identity("bcftools", source) == identity("bcftools", cached.output_file)
    stats = cached.last_run_stats
    assert stats["shards"] == (4 if per_thread == 4 else 2)
    assert stats["bypassed_shards"] == 1
    if per_thread == 4:
        # Hits: chr1:10, chr1:45, chr2:10 and the lower-case chr2:20 only.
        assert stats["missing_variants"] == len(UNUSUAL) - 4
    assert (cached.output_file.parent / "result.bcf.csi").exists()


def test_duplicates_are_found_across_read_chunks(tmp_path, unusual, monkeypatch):
    source, _ = unusual
    bed = tmp_path / "all.bed"
    bed.write_text("".join(f"{c}\t0\t1000\n" for c in ("chr1", "chr2", "chr3")))
    clean = tmp_path / "clean.bed"
    clean.write_text("chr1\t0\t1000\nchr2\t0\t1000\n")
    for size in (1 << 24, 16, 1):
        monkeypatch.setattr(preserving, "CHUNK", size)
        assert preserving.inspect_keys("bcftools", source, bed, False) == (True, None)
        assert preserving.inspect_keys("bcftools", source, clean, True) == (
            False,
            ["chr1\t0\t16384", "chr2\t0\t16384"],
        )


@pytest.mark.parametrize("level", [0, 1, 9])
def test_output_level_changes_bytes_not_records(tmp_path, inputs, level):
    source, cache = inputs
    reference = manager(tmp_path / "reference", source)
    reference.run("annotate", db_bcf=cache)
    run = manager(tmp_path / "run", source)
    run.output_level = level
    run.run("annotate", db_bcf=cache)
    assert records(run.output_file) == records(reference.output_file)
    assert run.last_run_stats["output_level"] == level
    assert reference.last_run_stats["output_level"] == 1


def test_invalid_output_level_is_rejected(tmp_path, inputs):
    source, cache = inputs
    run = manager(tmp_path / "run", source)
    run.output_level = 10
    with pytest.raises(ValueError, match="compression level"):
        run.run("annotate", db_bcf=cache)


def test_shards_keep_index_order_and_respect_the_limit():
    occupied = [("chr1", 100, 50), ("chr2", 100, 5), ("chr3", 100, 40), ("chrX", 9, 5)]
    for limit in (1, 2, 3, 4, 99):
        shards = preserving.plan_shards(occupied, limit)
        assert len(shards) <= min(limit, len(occupied))
        assert [c for shard in shards for c, _ in shard] == [c for c, _, _ in occupied]
    # Small contigs share a shard, so shards carry comparable record counts.
    assert preserving.plan_shards(occupied, 99)[1] == [("chr2", 100), ("chr3", 100)]
    assert preserving.plan_shards([], 8) == []


def test_cached_run_rejects_changed_recipe_output(tmp_path, inputs):
    """The recipe check runs beside the merge, yet nothing may be published."""
    source, cache = inputs
    recipe = "${params.bcftools_cmd} annotate --set-id CHANGED ${INPUT_BCF} -Ob -o ${OUTPUT_BCF} -W"
    run = manager(tmp_path / "run", source, recipe=recipe)
    run.output_file.write_bytes(b"previous deliverable")
    with pytest.raises(RuntimeError, match="changed variant identity"):
        run.run("annotate", db_bcf=cache)
    assert run.output_file.read_bytes() == b"previous deliverable"


@pytest.mark.parametrize("ratio", [0, 10**9])
@pytest.mark.parametrize("subset", ["all", "hits"])
def test_untagged_input_needs_no_hit_marker(
    tmp_path, unusual, monkeypatch, ratio, subset
):
    """An input that never declared the tag cannot hide annotation of its own."""
    _, cache = unusual
    monkeypatch.setattr(preserving, "WINDOW_RATIO", ratio)
    header = "\n".join(
        line for line in WIDE.splitlines() if "MOCK_ANNO" not in line
    ).replace("##FORMAT", "##contig=<ID=chrUn,length=1000>\n##FORMAT", 1)
    keep = [
        (n, row)
        for n, row in enumerate(UNUSUAL)
        if row[4] == "."
        and (subset == "all" or row[:2] in {("chr1", 10), ("chr2", 10)})
    ]
    source = bcf(
        tmp_path,
        "untagged",
        [
            f"{c}\t{p}\tid{n}\t{r}\t{a}\t50\tPASS\t.\tGT:DP\t0/1:{n}"
            for n, (c, p, r, a, _) in keep
        ],
        header + "\n",
    )
    cached = manager(tmp_path / "cached", source)
    cached.run("annotate", db_bcf=cache)
    direct = manager(tmp_path / "direct", source)
    direct.run("annotate-nocache")
    assert records(cached.output_file) == records(direct.output_file)
    stats = cached.last_run_stats
    assert stats["hit_marker"] is False
    assert stats["missing_variants"] == (len(keep) - 3 if subset == "all" else 0)
    assert HIT.encode() not in subprocess.check_output(
        ["bcftools", "view", "-h", str(cached.output_file)]
    )
