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


def test_duplicate_keys_go_to_the_recipe_and_equal_direct(tmp_path, inputs):
    _, cache = inputs
    source = bcf(tmp_path, "duplicates", ROWS)
    direct = manager(tmp_path / "direct", source)
    direct.run("annotate-nocache")
    run = manager(tmp_path / "run", source)
    run.run("annotate", db_bcf=cache)
    assert identity("bcftools", source) == identity("bcftools", run.output_file)
    assert records(run.output_file) == records(direct.output_file)
    assert run.last_run_stats["output_variants"] == 8
    assert run.last_run_stats["duplicate_sites"] == 1
    # The hit at chr1:10 still uses the cache; both copies at chr1:60 do not.
    assert run.last_run_stats["missing_variants"] == 7


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
@pytest.mark.parametrize("per_thread", [None, 1, 8])
def test_unusual_input_equals_direct_annotation(
    tmp_path, unusual, monkeypatch, ratio, per_thread
):
    """Identical records whether the input is one shard or cut everywhere."""
    source, cache = unusual
    monkeypatch.setattr(preserving, "WINDOW_RATIO", ratio)
    if per_thread:
        # Cut at every sampled record, inside contigs and next to duplicates.
        monkeypatch.setattr(preserving, "MIN_SHARD", 1)
        monkeypatch.setattr(preserving, "SAMPLE", 1)
        monkeypatch.setattr(preserving, "SHARDS_PER_THREAD", per_thread)
    direct = manager(tmp_path / "direct", source)
    direct.run("annotate-nocache")
    cached = manager(tmp_path / "cached", source)
    cached.run("annotate", db_bcf=cache)
    assert records(cached.output_file) == records(direct.output_file)
    assert identity("bcftools", source) == identity("bcftools", cached.output_file)
    stats = cached.last_run_stats
    # 17 records at 15 distinct positions: with a cut tried at every record,
    # the two duplicated positions each stay in one shard.
    if ratio:
        # Not windowed: planned from the index, whose whole test file sits in
        # one compressed block, so no cut can be placed inside it.
        assert stats["planning"] == "index" if per_thread else stats["shards"] == 1
    else:
        assert stats["shards"] == {None: 1, 1: 2, 8: 15}[per_thread]
    assert stats["duplicate_sites"] == 2
    # Hits: chr1:10, chr1:45, chr2:10, the lower-case chr2:20 and chr3:90,
    # which shares its contig with duplicates but not their position.
    assert stats["missing_variants"] == len(UNUSUAL) - 5
    assert (cached.output_file.parent / "result.bcf.csi").exists()


def test_key_scan_finds_duplicates_samples_and_windows(unusual, monkeypatch):
    source, _ = unusual
    monkeypatch.setattr(preserving, "SAMPLE", 2)
    chr3 = preserving.scan_contig("bcftools", source, "chr3", False)
    assert chr3 == dict(
        records=5,
        samples=[(0, 60), (2, 80), (4, 90)],
        duplicates=[60, 80],
        blocks=None,
    )
    chr1 = preserving.scan_contig("bcftools", source, "chr1", True)
    assert (chr1["records"], chr1["duplicates"], chr1["blocks"]) == (6, [], {0})
    # chr2:20 holds a lower-case key; no other key there, so no duplicate.
    assert preserving.scan_contig("bcftools", source, "chr2", False)["duplicates"] == []


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


def scans_of(**contigs):
    """Synthetic scans: one sample per record at the given positions."""
    return {
        c: dict(
            records=len(p),
            samples=list(enumerate(p)),
            duplicates=[],
            blocks=None,
        )
        for c, p in contigs.items()
    }


def test_shards_have_equal_records_and_cut_only_between_positions():
    U = preserving.UNBOUNDED
    scans = scans_of(chr1=[10, 20, 20, 20, 30, 40], chr2=[5, 6], chr3=[7, 8, 9, 10])
    order = ["chr1", "chr2", "chr3"]
    assert preserving.plan_ranges(order, scans, 1) == [
        [("chr1", 1, U), ("chr2", 1, U), ("chr3", 1, U)]
    ]
    # Twelve records in three parts; the middle shard spans chr1 and chr2.
    assert preserving.plan_ranges(order, scans, 3) == [
        [("chr1", 1, 29)],
        [("chr1", 30, U), ("chr2", 1, U)],
        [("chr3", 1, U)],
    ]
    # A file whose contigs are not in header order: every contig on its own.
    assert preserving.plan_ranges(order, scans, 3, span=False) == [
        [("chr1", 1, 29)],
        [("chr1", 30, U)],
        [("chr2", 1, U)],
        [("chr3", 1, U)],
    ]
    # A cut that would fall inside the group at chr1:20 moves before it.
    assert preserving.plan_ranges(["chr1"], scans, 2) == [
        [("chr1", 1, 19)],
        [("chr1", 20, U)],
    ]
    shards = preserving.plan_ranges(order, scans, 12)
    covered = [(c, a, b) for shard in shards for c, a, b in shard]
    assert [c for c, _, _ in covered] == sorted(
        [c for c, _, _ in covered], key=order.index
    )
    # Every position, the group at chr1:20 included, falls in exactly one range.
    for contig, positions in (("chr1", scans["chr1"]["samples"]),):
        for _, position in positions:
            assert sum(c == contig and a <= position <= b for c, a, b in covered) == 1
    assert preserving.plan_ranges(["chr1"], scans_of(chr1=[]), 4) == []


def test_windows_are_clipped_to_shard_ranges():
    scans = {
        "chr1": dict(records=3, samples=[], duplicates=[16380, 16390], blocks={0, 1})
    }
    first = [("chr1", 1, 16384)]
    second = [("chr1", 16385, preserving.UNBOUNDED)]
    assert preserving.shard_rows(first, scans) == ["chr1\t0\t16384"]
    assert preserving.shard_rows(second, scans) == ["chr1\t16384\t32768"]
    assert preserving.shard_duplicates(first, scans) == ["chr1\t16380"]
    assert preserving.shard_duplicates(second, scans) == ["chr1\t16390"]


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
    assert stats["missing_variants"] == (len(keep) - 4 if subset == "all" else 0)
    assert HIT.encode() not in subprocess.check_output(
        ["bcftools", "view", "-h", str(cached.output_file)]
    )


# Annotation that differs per record (from its ID), so a value paired with the
# wrong copy of a duplicated key, or a cached value left in place, is visible.
PER_RECORD = r"""${params.bcftools_cmd} view ${INPUT_BCF} -Ov | awk -F '\t' -v OFS='\t' '
  BEGIN { while ((getline id < "skip.txt") > 0) skip[id] = 1 }
  /^##INFO=<ID=MOCK_ANNO,/ { declared = 1 }
  /^#CHROM/ && !declared {
    print "##INFO=<ID=MOCK_ANNO,Number=1,Type=String,Description=\"Test annotation\">" }
  /^#/ { print; next }
  { if (!($3 in skip)) $8 = ($8 == "." ? "" : $8 ";") "MOCK_ANNO=v_" $3; print }
' | ${params.bcftools_cmd} view -Ob -o ${OUTPUT_BCF} -W
"""


@pytest.mark.parametrize("cuts", [False, True])
@pytest.mark.parametrize("declared", [True, False])
@pytest.mark.parametrize("skip", ["none", "f", "h"])
def test_duplicate_positions_never_keep_cached_values(
    tmp_path, monkeypatch, declared, skip, cuts
):
    """Every record at a duplicate position gets its own recipe result.

    The cache holds WRONG for each duplicated key; the recipe writes a value
    derived from each record's ID and, with ``skip``, writes nothing for one
    record that the cache would otherwise have answered. With ``cuts`` a shard
    boundary is attempted at every record.
    """
    if cuts:
        monkeypatch.setattr(preserving, "MIN_SHARD", 1)
        monkeypatch.setattr(preserving, "SAMPLE", 1)
        monkeypatch.setattr(preserving, "SHARDS_PER_THREAD", 8)
    header = (
        HEADER
        if declared
        else HEADER.replace(
            '##INFO=<ID=MOCK_ANNO,Number=1,Type=String,Description="Test annotation">\n',
            "",
        )
    )
    rows = [
        ("chr1", 10, "a", "A", "C"),
        ("chr1", 60, "f", "A", "C"),
        ("chr1", 60, "g", "A", "C"),
        ("chr1", 60, "h", "A", "G"),
        ("chr1", 70, "x", "A", "T"),
        ("chr1", 70, "y", "A", "T"),
        ("chr1", 80, "z", "G", "T"),
    ]
    source = bcf(
        tmp_path,
        "dups",
        [
            f"{c}\t{p}\t{i}\t{r}\t{a}\t50\tPASS\t.\tGT:DP\t0/1:9"
            for c, p, i, r, a in rows
        ],
        header,
    )
    cache = bcf(
        tmp_path,
        "dup_cache",
        [
            "chr1\t60\t.\tA\tC\t.\t.\tMOCK_ANNO=WRONG",
            "chr1\t60\t.\tA\tG\t.\t.\tMOCK_ANNO=WRONG",
            "chr1\t70\t.\tA\tT\t.\t.\tMOCK_ANNO=WRONG",
        ],
        SITES.replace("chr2", "chrA").replace("chr3", "chrB"),
    )
    runs = {}
    for name, step in (("direct", "annotate-nocache"), ("cached", "annotate")):
        run = manager(tmp_path / name, source, recipe=PER_RECORD)
        (run.output_dir / "work" / "annotate").mkdir(parents=True)
        (run.output_dir / "work" / "annotate-nocache").mkdir(parents=True)
        for folder in ("annotate", "annotate-nocache"):
            (run.output_dir / "work" / folder / "skip.txt").write_text(
                ("" if skip == "none" else skip) + "\n"
            )
        run.run(step, db_bcf=cache)
        runs[name] = run
    cached, direct = (records(runs[n].output_file) for n in ("cached", "direct"))
    # bcftools annotate may reorder records that share a position (a known,
    # documented property); content and the order of positions must match.
    assert sorted(cached.splitlines()) == sorted(direct.splitlines())
    assert [line.split(b"\t")[:2] for line in cached.splitlines()] == [
        line.split(b"\t")[:2] for line in direct.splitlines()
    ]
    assert b"WRONG" not in cached
    stats = runs["cached"].last_run_stats
    assert stats["duplicate_sites"] == 2
    # f and h sit at a duplicate position and would hold a cached WRONG.
    assert stats["duplicate_clean_base"] == (skip in ("f", "h"))


@pytest.mark.parametrize("index", ["csi", "tbi", "bcf"])
def test_index_planning_keeps_file_order_and_cuts_large_inputs(
    tmp_path, monkeypatch, index
):
    """A file whose contigs are not in header order, large enough for cuts.

    Records come out in the file's own order, cuts fall inside contigs, and
    duplicates on either side of a cut still go to the recipe.
    """
    monkeypatch.setattr(preserving, "WINDOW_RATIO", 10**9)
    monkeypatch.setattr(preserving, "MIN_SHARD", 1)
    monkeypatch.setattr(preserving, "SHARDS_PER_THREAD", 4)
    header = HEADER.replace(
        "##contig=<ID=chr1,length=1000>",
        "##contig=<ID=chr1,length=50000000>\n##contig=<ID=chr2,length=50000000>",
    )
    rows = []
    for contig in ("chr2", "chr1"):  # file order differs from header order
        for n, position in enumerate(range(100, 40_000_000, 2_000)):
            rows.append(
                f"{contig}\t{position}\tr{contig}{n}\tA\tC\t50\tPASS\t.\tGT:DP\t0/1:{n % 90}"
            )
            if n % 997 == 0:  # an exact duplicate now and then
                rows.append(rows[-1].replace(f"r{contig}{n}", f"d{contig}{n}"))
    vcf = tmp_path / "big.vcf"
    vcf.write_text(header + "\n".join(rows) + "\n")
    if index == "bcf":
        source = tmp_path / "big.bcf"
        subprocess.run(
            ["bcftools", "view", str(vcf), "-Ob", "-o", str(source), "-W"], check=True
        )
    else:
        source = tmp_path / "big.vcf.gz"
        subprocess.run(
            ["bcftools", "view", str(vcf), "-Oz", "-o", str(source)], check=True
        )
        subprocess.run(
            ["bcftools", "index", "-t" if index == "tbi" else "-c", str(source)],
            check=True,
        )
    cache = bcf(
        tmp_path,
        "big_cache",
        [
            line.split("\t", 5)[0]
            + "\t"
            + line.split("\t", 5)[1]
            + "\t.\tA\tC\t.\t.\tMOCK_ANNO=test"
            for line in rows[::3]
        ],
        header.replace("\tFORMAT\tSAMPLE", ""),
    )
    direct = manager(tmp_path / "direct", source)
    direct.run("annotate-nocache")
    cached = manager(tmp_path / "cached", source)
    cached.run("annotate", db_bcf=cache)
    stats = cached.last_run_stats
    assert stats["planning"] == "index"
    assert stats["shards"] > 2
    # The joined shard indexes must answer exactly as bcftools' own index.
    assert stats["joined_index"]
    own = tmp_path / "reindexed.bcf"
    own.write_bytes(cached.output_file.read_bytes())
    subprocess.run(["bcftools", "index", str(own)], check=True)
    for region in ("chr2", "chr1", "chr2:1-15000000", "chr1:20000000-39000000"):
        query = ["bcftools", "view", "-H", "-r", region]
        assert subprocess.check_output(query + [str(cached.output_file)]) == (
            subprocess.check_output(query + [str(own)])
        )
    assert stats["duplicate_sites"] > 0
    out = records(cached.output_file)
    assert out.split(b"\n", 1)[0].startswith(b"chr2\t")
    assert sorted(out.splitlines()) == sorted(records(direct.output_file).splitlines())
    assert [l.split(b"\t")[:2] for l in out.splitlines()] == [
        l.split(b"\t")[:2] for l in records(direct.output_file).splitlines()
    ]
