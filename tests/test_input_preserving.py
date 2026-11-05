"""Real bcftools checks for the immutable-record annotation contract."""

import hashlib
import subprocess

import pytest
import yaml

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


def bcf(folder, name, rows):
    """Write test-only VCF then index its BCF representation."""
    vcf = folder / f"{name}.vcf"
    vcf.write_text(HEADER + "\n".join(rows) + ("\n" if rows else ""))
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


@pytest.mark.parametrize("regions", ["off", "windows"])
def test_exact_preserves_all_records_and_annotates_partial_hits(
    tmp_path, inputs, monkeypatch, regions
):
    source, cache = inputs
    before = hashlib.sha256(source.read_bytes()).hexdigest()
    direct = manager(tmp_path / "direct", source)
    direct.run("annotate-nocache")
    monkeypatch.setenv("VCFCACHE_PILOT_REGION_MODE", regions)
    cached = manager(tmp_path / "cached", source)
    cached.run("annotate", db_bcf=cache)
    assert identity("bcftools", source) == identity("bcftools", cached.output_file)

    def records(path):
        return subprocess.check_output(["bcftools", "view", "-H", str(path)])

    assert records(cached.output_file) == records(direct.output_file)
    assert cached.last_run_stats["input_variants"] == 7
    assert cached.last_run_stats["missing_variants"] == 5
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


def test_duplicate_keys_are_preserved_using_direct_recipe(tmp_path, inputs):
    _, cache = inputs
    source = bcf(tmp_path, "duplicates", ROWS)
    run = manager(tmp_path / "run", source)
    run.run("annotate", db_bcf=cache)
    assert identity("bcftools", source) == identity("bcftools", run.output_file)
    assert run.last_run_stats["output_variants"] == 8
    assert run.last_run_stats["cache_matching"] == "bypassed-duplicate-keys"


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
