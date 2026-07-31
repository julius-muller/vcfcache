from __future__ import annotations

import csv
import subprocess
from pathlib import Path

import pytest

from benchmarks.prepare_assay_data import (
    ACMG_GENE_FILE,
    HPRC_QUOTAS,
    AssaySample,
    build_mane_panel_intervals,
    deterministic_hprc_selection,
    load_acmg_genes,
    merge_intervals,
    pad_intervals,
    parse_gtf_attributes,
    prepare_one_interval_sample,
    validate_interval_vcf,
)


def _candidate(
    sample: str, superpopulation: str, sex: str, population: str = "POP"
) -> AssaySample:
    return AssaySample(sample, population, superpopulation, sex)


def test_acmg_v33_gene_list_is_exact_and_frozen():
    genes = load_acmg_genes(ACMG_GENE_FILE)
    assert len(genes) == 84
    assert len(set(genes)) == 84
    assert {"ABCD1", "CYP27A1", "PLN"} <= set(genes)
    assert not ({"ADA2", "GCK", "RUNX1", "SLC4A1"} & set(genes))


def test_hprc_selection_is_deterministic_stratified_and_excludes_existing():
    candidates = []
    for superpopulation, quota in HPRC_QUOTAS.items():
        for number in range(quota + 3):
            sex = "female" if number % 2 == 0 else "male"
            candidates.append(
                _candidate(
                    f"{superpopulation}{number:02d}",
                    superpopulation,
                    sex,
                    f"P{number % 3}",
                )
            )
    excluded = {"AFR00", "EUR00"}

    first = deterministic_hprc_selection(candidates, exclusions=excluded)
    second = deterministic_hprc_selection(candidates, exclusions=excluded)

    assert first == second
    assert len(first) == 20
    assert not ({sample.sample for sample in first} & excluded)
    for superpopulation, quota in HPRC_QUOTAS.items():
        selected = [s for s in first if s.superpopulation == superpopulation]
        assert len(selected) == quota
        assert {s.sex for s in selected} == {"female", "male"}


def test_hprc_selection_reports_insufficient_stratum():
    with pytest.raises(ValueError, match="AFR"):
        deterministic_hprc_selection(
            [_candidate("AMR1", "AMR", "female")], exclusions=set()
        )


def test_parse_gtf_attributes():
    parsed = parse_gtf_attributes(
        'gene_id "ENSG1"; transcript_id "ENST1"; gene_name "APC"; '
        'tag "basic"; tag "MANE_Select";'
    )
    assert parsed["gene_name"] == ["APC"]
    assert parsed["tag"] == ["basic", "MANE_Select"]


def test_merge_intervals_merges_overlap_and_bookended_regions():
    assert merge_intervals([("chr2", 20, 30), ("chr1", 20, 30), ("chr1", 10, 20)]) == [
        ("chr1", 10, 30),
        ("chr2", 20, 30),
    ]


def test_pad_intervals_clips_to_contigs_then_merges():
    assert pad_intervals(
        [("chr1", 5, 10), ("chr1", 20, 25), ("chrX", 95, 100)],
        padding=10,
        contig_lengths={"chr1": 100, "chrX": 100},
    ) == [("chr1", 0, 35), ("chrX", 85, 100)]


def test_build_mane_panel_intervals_filters_pads_merges_and_maps_contigs(tmp_path):
    gtf = tmp_path / "mini.gtf"
    gtf.write_text(
        '1\tsrc\tCDS\t101\t110\t.\t+\t0\tgene_name "APC"; '
        'transcript_id "T1"; tag "MANE_Select";\n'
        '1\tsrc\tCDS\t131\t140\t.\t+\t0\tgene_name "APC"; '
        'transcript_id "T1"; tag "MANE_Select";\n'
        'X\tsrc\tCDS\t5\t10\t.\t-\t0\tgene_name "ABCD1"; '
        'transcript_id "T2"; tag "MANE_Select";\n'
        '2\tsrc\texon\t1\t50\t.\t+\t.\tgene_name "APC"; '
        'transcript_id "T1"; tag "MANE_Select";\n'
        '3\tsrc\tCDS\t1\t50\t.\t+\t0\tgene_name "APC"; '
        'transcript_id "T3"; tag "basic";\n'
        '4\tsrc\tCDS\t1\t50\t.\t+\t0\tgene_name "NOT_PANEL"; '
        'transcript_id "T4"; tag "MANE_Select";\n'
    )

    intervals, observed = build_mane_panel_intervals(gtf, {"APC", "ABCD1"}, padding=20)

    assert observed == {"APC", "ABCD1"}
    assert intervals == [("chr1", 80, 160), ("chrX", 0, 30)]


def _write_vcf(path: Path, sample: str, records: list[str]) -> Path:
    plain = path.with_suffix("")
    plain.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=1000>\n"
        "##contig=<ID=chrX,length=1000>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"
        + "".join(records)
    )
    with path.open("wb") as output:
        subprocess.run(["bgzip", "-c", str(plain)], check=True, stdout=output)
    subprocess.run(["tabix", "-p", "vcf", str(path)], check=True)
    return path


def test_prepare_interval_sample_combines_autosomes_and_x_as_vcfgz(tmp_path):
    sample = "S1"
    autosomal = _write_vcf(
        tmp_path / "autosomal.vcf.gz",
        sample,
        [
            "chr1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n",
            "chr1\t900\t.\tC\tT\t.\tPASS\t.\tGT\t0/1\n",
        ],
    )
    chromosome_x = _write_vcf(
        tmp_path / "x.vcf.gz",
        sample,
        ["chrX\t20\t.\tG\tA\t.\tPASS\t.\tGT\t1\n"],
    )
    bed = tmp_path / "targets.bed"
    bed.write_text("chr1\t0\t100\nchrX\t0\t100\n")
    destination = tmp_path / "output/sample.vcf.gz"

    result = prepare_one_interval_sample(
        sample=sample,
        autosomal_vcf=autosomal,
        chromosome_x_vcf=chromosome_x,
        bed=bed,
        destination=destination,
        work_dir=tmp_path / "work",
    )

    assert result == destination
    assert destination.exists()
    assert Path(f"{destination}.tbi").exists()
    validation = validate_interval_vcf(destination, allowed_contigs={"chr1", "chrX"})
    assert validation["status"] == "PASS"
    assert validation["records"] == 2
    assert validation["sample"] == sample
    queried = subprocess.check_output(
        ["bcftools", "query", "-f", "%CHROM:%POS\\n", str(destination)], text=True
    ).splitlines()
    assert queried == ["chr1:10", "chrX:20"]


def test_readable_manifest_schema(tmp_path):
    manifest = tmp_path / "samples.tsv"
    manifest.write_text(
        "sample\tpopulation\tsuperpopulation\tsex\tseed\n" "S1\tP1\tAFR\tfemale\tseed\n"
    )
    rows = list(csv.DictReader(manifest.open(), delimiter="\t"))
    assert rows[0]["sample"] == "S1"
