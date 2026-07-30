from __future__ import annotations

import hashlib
import os
import subprocess
from pathlib import Path

import pytest

from benchmarks.prepare_data import (
    GIAB_EXCLUSIONS,
    SelectedSample,
    SourceFile,
    build_flat_sample_view,
    deterministic_sample_selection,
    download_file,
    finalize_1000g_sample,
    parse_1000g_manifest,
    parse_population_panel,
    prepare_1000g_chromosome,
    smoke_test,
    validate_prepared_vcf,
    write_selected_samples,
)


def test_parse_1000g_manifest_pairs():
    parsed = parse_1000g_manifest(
        "chr1.vcf.gz 0123456789abcdef0123456789abcdef "
        "chr1.vcf.gz.tbi fedcba9876543210fedcba9876543210"
    )
    assert parsed == {
        "chr1.vcf.gz": "0123456789abcdef0123456789abcdef",
        "chr1.vcf.gz.tbi": "fedcba9876543210fedcba9876543210",
    }


def test_parse_1000g_manifest_rejects_malformed():
    with pytest.raises(ValueError, match="odd number"):
        parse_1000g_manifest("file checksum extra")
    with pytest.raises(ValueError, match="Invalid MD5"):
        parse_1000g_manifest("file not-an-md5")


def test_population_selection_is_deterministic_sex_balanced_and_excludes_giab():
    rows = ["sample\tpop\tsuper_pop\tgender"]
    for superpopulation in ("AFR", "AMR", "EAS", "EUR", "SAS"):
        for sex in ("female", "male"):
            for number in range(8):
                rows.append(
                    f"{superpopulation}_{sex}_{number}\tP{number}\t"
                    f"{superpopulation}\t{sex}"
                )
    rows.append("NA12878\tCEU\tEUR\tfemale")
    candidates = parse_population_panel("\n".join(rows))
    first = deterministic_sample_selection(candidates)
    second = deterministic_sample_selection(candidates)
    assert first == second
    assert len(first) == 50
    assert not ({sample.sample for sample in first} & GIAB_EXCLUSIONS)
    for superpopulation in ("AFR", "AMR", "EAS", "EUR", "SAS"):
        subset = [s for s in first if s.superpopulation == superpopulation]
        assert len(subset) == 10
        assert sum(s.sex == "female" for s in subset) == 5
        assert sum(s.sex == "male" for s in subset) == 5


def test_selection_requires_even_group_size():
    with pytest.raises(ValueError, match="must be even"):
        deterministic_sample_selection([], per_superpopulation=9)


def test_download_file_resumes_partial_verifies_and_renames_atomically(
    tmp_path, monkeypatch
):
    destination = tmp_path / "source.vcf.gz"
    partial = tmp_path / "source.vcf.gz.partial"
    partial.write_bytes(b"partial")
    expected = b"complete source"
    source = SourceFile(
        "test",
        "sample",
        "https://example.invalid/source.vcf.gz",
        hashlib.md5(expected).hexdigest(),
    )

    def fake_run(args, **_kwargs):
        assert "--continue-at" in args
        assert args[args.index("--continue-at") + 1] == "-"
        partial.write_bytes(expected)
        return subprocess.CompletedProcess(args, 0)

    monkeypatch.setattr("benchmarks.prepare_data.run_command", fake_run)

    assert download_file(source, destination) == destination
    assert destination.read_bytes() == expected
    assert not partial.exists()


def _write_compact_vcf(tmp_path: Path) -> Path:
    plain = tmp_path / "sample.vcf"
    plain.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=1000>\n"
        '##INFO=<ID=AF,Number=A,Type=Float,Description="AF">\n'
        '##INFO=<ID=AC,Number=A,Type=Integer,Description="AC">\n'
        '##INFO=<ID=AN,Number=1,Type=Integer,Description="AN">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr1\t10\t.\tA\tG\t.\tPASS\tAF=0.1;AC=2;AN=20\tGT\t0/1\n"
        "chr1\t20\t.\tAT\tA\t.\tPASS\tAF=0.2;AC=4;AN=20\tGT\t1/1\n"
    )
    compressed = tmp_path / "sample.vcf.gz"
    with compressed.open("wb") as output:
        subprocess.run(["bgzip", "--stdout", str(plain)], check=True, stdout=output)
    subprocess.run(["tabix", "--preset", "vcf", str(compressed)], check=True)
    return compressed


def test_validate_prepared_vcf_compact_1000g(tmp_path):
    vcf = _write_compact_vcf(tmp_path)
    result = validate_prepared_vcf(vcf, cohort="1000g")
    assert result["status"] == "PASS"
    assert result["sample"] == "S1"
    assert result["records"] == 2
    assert result["snps"] == 1
    assert result["indels"] == 1
    assert result["sha256"]


def test_prepare_1000g_chromosome_splits_filters_and_compacts(tmp_path):
    root = tmp_path / "benchmark"
    source_dir = root / "sources/1000g/vcf"
    source_dir.mkdir(parents=True)
    (root / "manifests").mkdir()
    selected = [
        SelectedSample("S1", "P1", "AFR", "female"),
        SelectedSample("S2", "P2", "EUR", "male"),
    ]
    write_selected_samples(root, selected)
    plain = tmp_path / "cohort.vcf"
    plain.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=1000>\n"
        '##INFO=<ID=AF,Number=A,Type=Float,Description="AF">\n'
        '##INFO=<ID=AC,Number=A,Type=Integer,Description="AC">\n'
        '##INFO=<ID=AN,Number=1,Type=Integer,Description="AN">\n'
        '##INFO=<ID=EXTRA,Number=1,Type=String,Description="remove me">\n'
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="remove me">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
        "chr1\t10\t.\tA\tG\t.\tPASS\tAF=0.1;AC=2;AN=4;EXTRA=x\t"
        "GT:DP\t0/1:20\t0/0:18\n"
        "chr1\t20\t.\tAT\tA\t.\tPASS\tAF=0.2;AC=2;AN=4;EXTRA=x\t"
        "GT:DP\t0/0:20\t1/1:18\n"
        "chr1\t30\t.\tN\t<DEL>\t.\tPASS\tAF=0.3;AC=2;AN=4;SVTYPE=DEL\t"
        "GT:DP\t0/1:20\t0/1:18\n"
        "chr1\t35\t.\tA\tG,<DEL>\t.\tPASS\t"
        "AF=0.2,0.1;AC=1,1;AN=4;SVTYPE=DEL\tGT:DP\t1/2:20\t0/0:18\n"
        "chr1\t40\t.\tC\tT\t.\tPASS\tAF=0.4;AC=0;AN=4;EXTRA=x\t"
        "GT:DP\t0/0:20\t0/0:18\n"
    )
    source = source_dir / (
        "1kGP_high_coverage_Illumina.chr1." "filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    )
    with source.open("wb") as output:
        subprocess.run(["bgzip", "--stdout", str(plain)], check=True, stdout=output)
    subprocess.run(["tabix", "--preset", "vcf", str(source)], check=True)

    output_dir = prepare_1000g_chromosome(root, "1")

    s1 = output_dir / "S1.vcf.gz"
    s2 = output_dir / "S2.vcf.gz"
    assert validate_prepared_vcf(s1, cohort="1000g")["status"] == "PASS"
    assert validate_prepared_vcf(s2, cohort="1000g")["status"] == "PASS"
    assert validate_prepared_vcf(s1, cohort="1000g")["records"] == 1
    assert validate_prepared_vcf(s2, cohort="1000g")["records"] == 1


def test_finalize_1000g_sample_concatenates_sorts_compresses_and_indexes(tmp_path):
    root = tmp_path / "benchmark"
    (root / "manifests").mkdir(parents=True)
    sample = SelectedSample("S1", "P1", "AFR", "female")
    write_selected_samples(root, [sample])
    contigs = "".join(
        f"##contig=<ID=chr{chromosome},length=1000>\n" for chromosome in range(1, 23)
    )
    header = (
        "##fileformat=VCFv4.2\n"
        f"{contigs}"
        '##INFO=<ID=AF,Number=A,Type=Float,Description="AF">\n'
        '##INFO=<ID=AC,Number=A,Type=Integer,Description="AC">\n'
        '##INFO=<ID=AN,Number=1,Type=Integer,Description="AN">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
    )
    for chromosome in range(1, 23):
        shard_dir = root / f"work/chromosome_shards/chr{chromosome}"
        shard_dir.mkdir(parents=True)
        shard = shard_dir / "S1.vcf.gz"
        record = (
            f"chr{chromosome}\t10\t.\tA\tG\t.\tPASS\t" "AF=0.1;AC=1;AN=2\tGT\t0/1\n"
        )
        with shard.open("wb") as output:
            subprocess.run(
                ["bgzip", "--stdout"],
                input=(header + record).encode(),
                check=True,
                stdout=output,
            )
        subprocess.run(["tabix", "--preset", "vcf", str(shard)], check=True)

    output = finalize_1000g_sample(root, sample)

    result = validate_prepared_vcf(output, cohort="1000g")
    assert result["status"] == "PASS"
    assert result["records"] == 22
    assert Path(f"{output}.tbi").exists()


def test_build_flat_sample_view_links_all_vcfs_and_indexes(tmp_path):
    root = tmp_path / "benchmark"
    for number in range(50):
        vcf = root / f"samples/GRCh38/1000g/POP/S{number}.vcf.gz"
        vcf.parent.mkdir(parents=True, exist_ok=True)
        vcf.touch()
        Path(f"{vcf}.tbi").touch()
    for number in range(7):
        vcf = root / f"samples/GRCh38/giab/HG{number}.vcf.gz"
        vcf.parent.mkdir(parents=True, exist_ok=True)
        vcf.touch()
        Path(f"{vcf}.tbi").touch()

    output = build_flat_sample_view(root)

    links = list(output.iterdir())
    assert len(links) == 114
    assert all(link.is_symlink() and link.resolve().exists() for link in links)


@pytest.mark.benchmark_network
@pytest.mark.skipif(
    os.environ.get("VCFCACHE_BENCHMARK_NETWORK") != "1",
    reason="set VCFCACHE_BENCHMARK_NETWORK=1 for the official-source smoke test",
)
def test_official_source_smoke(tmp_path):
    smoke_test(tmp_path / "benchmark_data")
