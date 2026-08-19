from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest

from benchmarks.run_cohort import (
    build_tasks,
    mode_order,
    prepare_campaign,
    submit_phase,
    worker_path,
    write_tasks,
)
from benchmarks.run_pilot import (
    DEFAULT_CACHE_PROVENANCE_EXPECTED,
    PilotConfig,
    annotation_command,
    calculate_cache_hit_rate,
    cgroup_v2_snapshot,
    parse_elapsed,
    parse_gnu_time,
    semantic_compare,
    strict_semantic_compare,
    validate_default_cache_provenance,
    validate_requirements_output,
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
    assert cached[cached.index("--statistics") + 1] == "light"


def test_default_cache_provenance_is_fail_closed(tmp_path):
    provenance = tmp_path / "zenodo_provenance.json"
    provenance.write_text(json.dumps(DEFAULT_CACHE_PROVENANCE_EXPECTED))
    assert validate_default_cache_provenance(provenance)["source"] == (
        "zenodo_production"
    )
    provenance.write_text(
        json.dumps({**DEFAULT_CACHE_PROVENANCE_EXPECTED, "doi": "wrong"})
    )
    with pytest.raises(RuntimeError, match="Invalid bundled-cache provenance"):
        validate_default_cache_provenance(provenance)


@pytest.mark.parametrize(
    ("requirements", "expected_version"),
    [
        (
            "Tool checks\n  bcftools: 1.21+htslib-1.21  ✓\n"
            "  annotation tool: 115.2  ✓\n",
            "115.2",
        ),
        (
            "Tool checks\n  bcftools: 1.21+htslib-1.21  ✓\n"
            "  annotation tool: fastvep 0.3.0  ✓\n",
            "0.3.0",
        ),
    ],
)
def test_requirements_validation_accepts_cache_declared_tool_version(
    requirements, expected_version
):
    validate_requirements_output(requirements, expected_tool_version=expected_version)


def test_requirements_validation_rejects_wrong_annotation_tool_version():
    requirements = (
        "Tool checks\n  bcftools: 1.21+htslib-1.21  ✓\n"
        "  annotation tool: fastvep 0.3.0  ✓\n"
    )
    with pytest.raises(RuntimeError, match="expected version '115.2'"):
        validate_requirements_output(requirements, expected_tool_version="115.2")


def test_requirements_validation_rejects_failed_bcftools_check():
    requirements = (
        "Tool checks\n  bcftools: ERROR (not found)  ✗\n"
        "  annotation tool: fastvep 0.3.0  ✓\n"
    )
    with pytest.raises(RuntimeError, match="bcftools requirement check"):
        validate_requirements_output(requirements, expected_tool_version="0.3.0")


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


def test_replicates_use_distinct_report_paths(tmp_path):
    first = PilotConfig(tmp_path, tmp_path / "sample.vcf.gz", tmp_path, tmp_path, 1)
    second = PilotConfig(tmp_path, tmp_path / "sample.vcf.gz", tmp_path, tmp_path, 2)
    assert first.comparison_path.name == "semantic_comparison_r01.json"
    assert second.comparison_path.name == "semantic_comparison_r02.json"
    assert first.summary_path.name == "summary_r01.json"
    assert second.summary_path.name == "summary_r02.json"


def test_cohort_tasks_are_deterministic_and_single_pass(tmp_path):
    qc = tmp_path / "sample_qc.tsv"
    qc.write_text(
        "cohort\tsample\tpopulation\tsuperpopulation\tsex\tpath\trecords\tsha256\tstatus\n"
        "1000g\tS1\tPOP\tEUR\tfemale\t/mnt/data/S1.vcf.gz\t10\tabc\tPASS\n"
    )
    first = build_tasks(qc, seed="paper", selected_sample="S1")
    second = build_tasks(qc, seed="paper", selected_sample="S1")
    assert first == second
    assert [task.task_id for task in first] == [0]
    assert [task.replicate for task in first] == [1]
    assert all(
        {task.first_mode, task.second_mode} == {"cached", "uncached"} for task in first
    )
    output = tmp_path / "tasks.tsv"
    write_tasks(output, first)
    assert (
        output.read_text()
        .splitlines()[0]
        .startswith("task_id\tphase\tmeasured\tsample")
    )
    with pytest.raises(ValueError, match="exactly one run"):
        build_tasks(qc, replicates=2, seed="paper", selected_sample="S1")


def test_mode_order_changes_with_auditable_key():
    first, second, key = mode_order("HG02079", 1, "vcfcache-paper-v1")
    assert {first, second} == {"cached", "uncached"}
    assert len(key) == 64


def test_measured_full_cohort_has_exact_balanced_order(tmp_path):
    qc = tmp_path / "sample_qc.tsv"
    header = (
        "cohort\tsample\tpopulation\tsuperpopulation\tsex\tpath\t"
        "records\tsha256\tstatus\n"
    )
    rows = [
        f"1000g\tS{index:02d}\tPOP\tEUR\tfemale\t/mnt/data/S{index:02d}.vcf.gz\t"
        f"10\tsha{index:02d}\tPASS\n"
        for index in range(50)
    ]
    qc.write_text(header + "".join(rows))
    tasks = build_tasks(qc, phase="measured", seed="paper")
    assert len(tasks) == 50
    assert sum(task.first_mode == "cached" for task in tasks) == 25
    assert len({task.sample for task in tasks}) == 50


def test_worker_path_translates_export_mount():
    assert worker_path(
        Path("/mnt/data/slurm-results/campaigns/run/manifests/smoke.tsv"),
        Path("/mnt/data/slurm-results"),
        Path("/results"),
    ) == Path("/results/campaigns/run/manifests/smoke.tsv")


def test_submit_phase_sets_worker_paths_and_working_directory(tmp_path, monkeypatch):
    controller = tmp_path / "controller"
    manifest = controller / "campaigns/run/manifests/smoke.tsv"
    manifest.parent.mkdir(parents=True)
    manifest.write_text("header\nrow\n")
    captured = {}

    def _run(command, **_kwargs):
        captured["command"] = command
        return subprocess.CompletedProcess(command, 0, stdout="42\n", stderr="")

    monkeypatch.setattr(
        "benchmarks.run_cohort.shutil.which", lambda _name: "/bin/sbatch"
    )
    monkeypatch.setattr("benchmarks.run_cohort.subprocess.run", _run)
    job_id, command = submit_phase(
        campaign_id="run",
        phase="smoke",
        controller_results=controller,
        worker_results=Path("/results"),
        concurrency=1,
    )
    assert job_id == "42"
    assert command == captured["command"]
    assert any(value.startswith("--chdir=") for value in command)
    export = next(value for value in command if value.startswith("--export="))
    assert "VCFCACHE_TASK_MANIFEST=/results/campaigns/run/manifests/smoke.tsv" in export


def test_prepare_campaign_precreates_shared_phase_directories(
    tmp_path, monkeypatch, capsys
):
    qc = tmp_path / "sample_qc.tsv"
    header = (
        "cohort\tsample\tpopulation\tsuperpopulation\tsex\tpath\t"
        "records\tsha256\tstatus\n"
    )
    rows = [
        f"1000g\tS{index:02d}\tPOP\tEUR\tfemale\t/mnt/data/S{index:02d}.vcf.gz\t"
        f"10\tsha{index:02d}\tPASS\n"
        for index in range(49)
    ]
    rows.append(
        "1000g\tHG02079\tPOP\tEAS\tmale\t/mnt/data/HG02079.vcf.gz\t10\tsmoke\tPASS\n"
    )
    qc.write_text(header + "".join(rows))
    controller = tmp_path / "controller"
    monkeypatch.setattr(
        "benchmarks.run_cohort.git_output",
        lambda *args: "abcdef123456" if "--short=12" in args else "abcdef1234567890",
    )
    monkeypatch.setattr(
        "benchmarks.run_cohort.subprocess.run",
        lambda *_args, **_kwargs: subprocess.CompletedProcess([], 0),
    )
    args = type(
        "Args",
        (),
        {
            "campaign_id": "run",
            "controller_results": controller,
            "worker_results": Path("/results"),
            "qc": qc,
            "seed": "paper",
        },
    )()
    prepare_campaign(args)
    capsys.readouterr()
    for phase in ("measured",):
        assert (controller / "campaigns/run" / phase / "tasks").is_dir()
        assert (controller / "campaigns/run" / phase / "attempts").is_dir()
    assert not (controller / "campaigns/run/smoke").exists()
    assert not (controller / "campaigns/run/warmup").exists()
    campaign = json.loads((controller / "campaigns/run/campaign.json").read_text())
    assert campaign["phases"]["measured"]["replicates"] == 1
    assert campaign["phases"]["measured"]["tasks"] == 50


def test_cgroup_v2_snapshot_parses_peak_and_counters(tmp_path):
    proc = tmp_path / "proc-cgroup"
    root = tmp_path / "cgroup"
    group = root / "slurm/job/step"
    group.mkdir(parents=True)
    proc.write_text("0::/slurm/job/step\n")
    (group / "memory.peak").write_text("123456\n")
    (group / "memory.current").write_text("123\n")
    (group / "memory.events").write_text("oom 0\noom_kill 0\n")
    (group / "cpu.stat").write_text("usage_usec 999\n")
    (group / "io.stat").write_text("8:0 rbytes=10 wbytes=20\n")
    snapshot = cgroup_v2_snapshot(proc_cgroup=proc, cgroup_root=root)
    assert snapshot is not None
    assert snapshot["memory_peak_bytes"] == 123456
    assert snapshot["cpu_stat"] == {"usage_usec": 999}


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
    path: Path,
    *,
    csq: str,
    position: int = 10,
    af: str = "0.1",
    csq_format: str | None = None,
) -> Path:
    plain = path.with_suffix("")
    csq_description = "Consequence"
    if csq_format:
        csq_description += f". Format: {csq_format}"
    plain.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=1000>\n"
        '##INFO=<ID=AF,Number=A,Type=Float,Description="AF">\n'
        '##INFO=<ID=AC,Number=A,Type=Integer,Description="AC">\n'
        '##INFO=<ID=AN,Number=1,Type=Integer,Description="AN">\n'
        f'##INFO=<ID=CSQ,Number=.,Type=String,Description="{csq_description}">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        f"chr1\t{position}\t.\tA\tG\t.\tPASS\t"
        f"AF={af};AC=1;AN=2;CSQ={csq}\tGT\t0/1\n"
    )
    with path.open("wb") as output:
        subprocess.run(["bgzip", "--stdout", str(plain)], check=True, stdout=output)
    subprocess.run(["tabix", "--preset", "vcf", str(path)], check=True)
    return path


def _write_same_locus_vcf(
    path: Path, alts: list[str], *, csq_suffix: dict[str, str] | None = None
) -> Path:
    plain = path.with_suffix("")
    records = []
    values = {"C": ("0.1", "1"), "G": ("0.2", "2")}
    csq_suffix = csq_suffix or {}
    for alt in alts:
        af, ac = values[alt]
        records.append(
            f"chr1\t10\t.\tA\t{alt}\t.\tPASS\t"
            f"AF={af};AC={ac};AN=4;CSQ={alt}|{csq_suffix.get(alt, 'effect')}"
            "\tGT\t0/1\n"
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


def _write_multicontig_vcf(path: Path, contigs: list[str]) -> Path:
    plain = path.with_suffix("")
    records = {
        "chr1": "chr1\t10\t.\tA\tG\t.\tPASS\tCSQ=G|effect\tGT\t0/1\n",
        "chr2": "chr2\t20\t.\tC\tT\t.\tPASS\tCSQ=T|effect\tGT\t0/1\n",
    }
    plain.write_text(
        "##fileformat=VCFv4.2\n"
        + "".join(f"##contig=<ID={contig},length=1000>\n" for contig in contigs)
        + '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence">\n'
        + '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        + "".join(records[contig] for contig in contigs)
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


def test_semantic_compare_ignores_domains_subfield_order(tmp_path):
    cached = _write_annotated_vcf(
        tmp_path / "cached.vcf.gz",
        csq="G|missense_variant|Pfam:1&CDD:2&SFLD:3",
        csq_format="Allele|Consequence|DOMAINS",
    )
    uncached = _write_annotated_vcf(
        tmp_path / "uncached.vcf.gz",
        csq="G|missense_variant|SFLD:3&Pfam:1&CDD:2",
        csq_format="Allele|Consequence|DOMAINS",
    )
    report = semantic_compare(cached, uncached)
    assert report["semantic_pass"] is True
    assert report["annotation_order_only"] == 1
    assert report["unordered_csq_fields"] == ["DOMAINS"]


def test_semantic_compare_detects_domains_member_difference(tmp_path):
    cached = _write_annotated_vcf(
        tmp_path / "cached.vcf.gz",
        csq="G|missense_variant|Pfam:1&CDD:2",
        csq_format="Allele|Consequence|DOMAINS",
    )
    uncached = _write_annotated_vcf(
        tmp_path / "uncached.vcf.gz",
        csq="G|missense_variant|Pfam:1&CDD:9",
        csq_format="Allele|Consequence|DOMAINS",
    )
    report = semantic_compare(cached, uncached)
    assert report["semantic_pass"] is False
    assert report["annotation_mismatches"] == 1


def test_semantic_compare_ignores_split_allele_order_within_locus(tmp_path):
    cached = _write_same_locus_vcf(tmp_path / "cached.vcf.gz", ["C", "G"])
    uncached = _write_same_locus_vcf(tmp_path / "uncached.vcf.gz", ["G", "C"])
    report = semantic_compare(cached, uncached)
    assert report["semantic_pass"] is True
    assert report["record_order_only_loci"] == 1
    assert report["key_mismatches"] == 0
    assert report["annotation_mismatches"] == 0


def test_semantic_compare_canonicalizes_contig_order_and_allows_missing_info(tmp_path):
    cached = _write_multicontig_vcf(tmp_path / "cached.vcf.gz", ["chr1", "chr2"])
    uncached = _write_multicontig_vcf(tmp_path / "uncached.vcf.gz", ["chr2", "chr1"])
    report = semantic_compare(cached, uncached)
    assert report["semantic_pass"] is True
    assert report["records_compared"] == 2
    assert report["key_mismatches"] == 0


def test_semantic_compare_detects_annotation_difference(tmp_path):
    cached = _write_annotated_vcf(tmp_path / "cached.vcf.gz", csq="a|1")
    uncached = _write_annotated_vcf(tmp_path / "uncached.vcf.gz", csq="b|2")
    report = semantic_compare(cached, uncached)
    assert report["semantic_pass"] is False
    assert report["annotation_mismatches"] == 1


def test_semantic_compare_ignores_only_hgnc_id_for_known_vep_bug(tmp_path):
    cached = _write_annotated_vcf(
        tmp_path / "cached.vcf.gz",
        csq="G|missense_variant|HGNC:1",
        csq_format="Allele|Consequence|HGNC_ID",
    )
    uncached = _write_annotated_vcf(
        tmp_path / "uncached.vcf.gz",
        csq="G|missense_variant|",
        csq_format="Allele|Consequence|HGNC_ID",
    )
    report = semantic_compare(cached, uncached)
    assert report["semantic_pass"] is True
    assert report["annotation_mismatches"] == 0
    assert report["raw_annotation_mismatches"] == 1
    assert report["ignored_annotation_mismatches"] == 1
    assert report["ignored_csq_fields"] == ["HGNC_ID"]
    assert report["ignored_difference_reference"].endswith("/issues/1959")
    assert report["examples"][0]["kind"] == "CSQ_ignored_known_vep_issue"


def test_semantic_compare_does_not_hide_other_csq_difference_with_hgnc_rule(
    tmp_path,
):
    cached = _write_annotated_vcf(
        tmp_path / "cached.vcf.gz",
        csq="G|missense_variant|HGNC:1",
        csq_format="Allele|Consequence|HGNC_ID",
    )
    uncached = _write_annotated_vcf(
        tmp_path / "uncached.vcf.gz",
        csq="G|synonymous_variant|",
        csq_format="Allele|Consequence|HGNC_ID",
    )
    report = semantic_compare(cached, uncached)
    assert report["semantic_pass"] is False
    assert report["annotation_mismatches"] == 1
    assert report["ignored_annotation_mismatches"] == 0


def _write_strict_vcf(path: Path, *, info: str, hgnc_description: str = "HGNC") -> Path:
    plain = path.with_suffix("")
    plain.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=1000>\n"
        '##FILTER=<ID=PASS,Description="All filters passed">\n'
        '##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: Allele|HGNC_ID">\n'
        f'##INFO=<ID=FV_CLINVAR,Number=1,Type=String,Description="{hgnc_description}">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        f"chr1\t10\t.\tA\tG\t.\tPASS\t{info}\tGT\t0/1\n"
    )
    with path.open("wb") as output:
        subprocess.run(["bgzip", "--stdout", str(plain)], check=True, stdout=output)
    subprocess.run(["tabix", "--preset", "vcf", str(path)], check=True)
    return path


def test_strict_fastvep_compare_only_canonicalizes_info_and_csq_order(tmp_path):
    cached = _write_strict_vcf(
        tmp_path / "cached.vcf.gz",
        info="FV_CLINVAR=pathogenic;CSQ=G|HGNC:2,G|HGNC:1",
    )
    uncached = _write_strict_vcf(
        tmp_path / "uncached.vcf.gz",
        info="CSQ=G|HGNC:1,G|HGNC:2;FV_CLINVAR=pathogenic",
    )
    report = strict_semantic_compare(cached, uncached)
    assert report["semantic_pass"] is True
    assert report["ignored_fields"] == []


def test_strict_fastvep_compare_allows_complete_record_order_within_locus(tmp_path):
    cached = _write_same_locus_vcf(tmp_path / "cached.vcf.gz", ["C", "G"])
    uncached = _write_same_locus_vcf(tmp_path / "uncached.vcf.gz", ["G", "C"])
    report = strict_semantic_compare(cached, uncached)
    assert report["semantic_pass"] is True
    assert report["record_mismatches"] == 0
    assert report["record_order_only_loci"] == 1
    assert report["ignored_fields"] == []


def test_strict_fastvep_compare_canonicalizes_indexed_contig_order(tmp_path):
    cached = _write_multicontig_vcf(tmp_path / "cached.vcf.gz", ["chr1", "chr2"])
    uncached = _write_multicontig_vcf(tmp_path / "uncached.vcf.gz", ["chr2", "chr1"])
    report = strict_semantic_compare(cached, uncached)
    assert report["semantic_pass"] is True
    assert report["record_mismatches"] == 0
    assert report["locus_mismatches"] == 0


def test_strict_fastvep_compare_does_not_hide_same_locus_value_change(tmp_path):
    cached = _write_same_locus_vcf(tmp_path / "cached.vcf.gz", ["C", "G"])
    uncached = _write_same_locus_vcf(
        tmp_path / "uncached.vcf.gz",
        ["G", "C"],
        csq_suffix={"G": "changed"},
    )
    report = strict_semantic_compare(cached, uncached)
    assert report["semantic_pass"] is False
    assert report["record_mismatches"] == 1
    assert report["record_order_only_loci"] == 0


def test_strict_fastvep_compare_detects_supplementary_and_hgnc_changes(tmp_path):
    cached = _write_strict_vcf(
        tmp_path / "cached.vcf.gz",
        info="CSQ=G|HGNC:2;FV_CLINVAR=pathogenic",
    )
    uncached = _write_strict_vcf(
        tmp_path / "uncached.vcf.gz",
        info="CSQ=G|HGNC:1;FV_CLINVAR=benign",
    )
    report = strict_semantic_compare(cached, uncached)
    assert report["semantic_pass"] is False
    assert report["record_mismatches"] == 1
    assert report["ignored_annotation_mismatches"] == 0


def test_strict_fastvep_compare_detects_relevant_header_change(tmp_path):
    cached = _write_strict_vcf(
        tmp_path / "cached.vcf.gz",
        info="CSQ=G|HGNC:1;FV_CLINVAR=pathogenic",
        hgnc_description="ClinVar A",
    )
    uncached = _write_strict_vcf(
        tmp_path / "uncached.vcf.gz",
        info="CSQ=G|HGNC:1;FV_CLINVAR=pathogenic",
        hgnc_description="ClinVar B",
    )
    report = strict_semantic_compare(cached, uncached)
    assert report["semantic_pass"] is False
    assert report["record_mismatches"] == 0
    assert report["relevant_headers_equal"] is False
