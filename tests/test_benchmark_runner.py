from __future__ import annotations

import subprocess
from pathlib import Path

from benchmarks.run_pilot import (
    PilotConfig,
    annotation_command,
    calculate_cache_hit_rate,
    parse_elapsed,
    parse_gnu_time,
    semantic_compare,
)


def test_annotation_commands_differ_only_by_uncached_flag(tmp_path):
    config = PilotConfig(
        data_root=tmp_path,
        input_vcf=tmp_path / "sample.vcf.gz",
        cache_dir=tmp_path / "cache",
        params_file=tmp_path / "params.yaml",
        replicate=1,
    )
    run_dir = tmp_path / "run"
    cached = annotation_command(config, "cached", run_dir)
    uncached = annotation_command(config, "uncached", run_dir)
    assert uncached[:-1] == cached
    assert uncached[-1] == "--uncached"


def test_parse_elapsed_formats():
    assert parse_elapsed("7.5") == 7.5
    assert parse_elapsed("2:03.5") == 123.5
    assert parse_elapsed("1:02:03") == 3723


def test_cache_hit_rate_falls_back_to_normalized_output_total():
    counts = {
        "input_variants": None,
        "total_output": 1_772,
        "tool_annotated": 71,
    }
    assert calculate_cache_hit_rate("uncached", counts, 1_772) is None
    assert calculate_cache_hit_rate("cached", counts, 1_772) == 1 - (71 / 1_772)


def test_parse_gnu_time(tmp_path):
    output = tmp_path / "time.txt"
    output.write_text(
        "User time (seconds): 10.5\n"
        "System time (seconds): 1.5\n"
        "Percent of CPU this job got: 240%\n"
        "Elapsed (wall clock) time (h:mm:ss or m:ss): 0:05.25\n"
        "Maximum resident set size (kbytes): 12345\n"
        "File system inputs: 10\n"
        "File system outputs: 20\n"
        "Voluntary context switches: 30\n"
        "Involuntary context switches: 40\n"
    )
    parsed = parse_gnu_time(output)
    assert parsed["wall_seconds_gnu"] == 5.25
    assert parsed["user_seconds"] == 10.5
    assert parsed["cpu_percent"] == 240
    assert parsed["max_rss_kib"] == 12345


def _write_annotated_vcf(
    path: Path, *, csq: str, position: int = 10, af: str = "0.1"
) -> Path:
    plain = path.with_suffix("")
    plain.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=1000>\n"
        '##INFO=<ID=AF,Number=A,Type=Float,Description="AF">\n'
        '##INFO=<ID=AC,Number=A,Type=Integer,Description="AC">\n'
        '##INFO=<ID=AN,Number=1,Type=Integer,Description="AN">\n'
        '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        f"chr1\t{position}\t.\tA\tG\t.\tPASS\t"
        f"AF={af};AC=1;AN=2;CSQ={csq}\tGT\t0/1\n"
    )
    with path.open("wb") as output:
        subprocess.run(["bgzip", "--stdout", str(plain)], check=True, stdout=output)
    subprocess.run(["tabix", "--preset", "vcf", str(path)], check=True)
    return path


def _write_same_locus_vcf(path: Path, alts: list[str]) -> Path:
    plain = path.with_suffix("")
    records = []
    values = {"C": ("0.1", "1"), "G": ("0.2", "2")}
    for alt in alts:
        af, ac = values[alt]
        records.append(
            f"chr1\t10\t.\tA\t{alt}\t.\tPASS\t"
            f"AF={af};AC={ac};AN=4;CSQ={alt}|effect\tGT\t0/1\n"
        )
    plain.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=1000>\n"
        '##INFO=<ID=AF,Number=A,Type=Float,Description="AF">\n'
        '##INFO=<ID=AC,Number=A,Type=Integer,Description="AC">\n'
        '##INFO=<ID=AN,Number=1,Type=Integer,Description="AN">\n'
        '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n" + "".join(records)
    )
    with path.open("wb") as output:
        subprocess.run(["bgzip", "--stdout", str(plain)], check=True, stdout=output)
    subprocess.run(["tabix", "--preset", "vcf", str(path)], check=True)
    return path


def test_semantic_compare_ignores_csq_item_order(tmp_path):
    cached = _write_annotated_vcf(tmp_path / "cached.vcf.gz", csq="a|1,b|2")
    uncached = _write_annotated_vcf(tmp_path / "uncached.vcf.gz", csq="b|2,a|1")
    report = semantic_compare(cached, uncached)
    assert report["semantic_pass"] is True
    assert report["annotation_order_only"] == 1
    assert report["annotation_mismatches"] == 0


def test_semantic_compare_ignores_split_allele_order_within_locus(tmp_path):
    cached = _write_same_locus_vcf(tmp_path / "cached.vcf.gz", ["C", "G"])
    uncached = _write_same_locus_vcf(tmp_path / "uncached.vcf.gz", ["G", "C"])
    report = semantic_compare(cached, uncached)
    assert report["semantic_pass"] is True
    assert report["record_order_only_loci"] == 1
    assert report["key_mismatches"] == 0
    assert report["annotation_mismatches"] == 0


def test_semantic_compare_detects_annotation_difference(tmp_path):
    cached = _write_annotated_vcf(tmp_path / "cached.vcf.gz", csq="a|1")
    uncached = _write_annotated_vcf(tmp_path / "uncached.vcf.gz", csq="b|2")
    report = semantic_compare(cached, uncached)
    assert report["semantic_pass"] is False
    assert report["annotation_mismatches"] == 1
