from __future__ import annotations

import subprocess
from pathlib import Path

from benchmarks.run_cohort import (
    build_tasks,
    mode_order,
    submit_phase,
    worker_path,
    write_tasks,
)
from benchmarks.run_pilot import (
    PilotConfig,
    annotation_command,
    calculate_cache_hit_rate,
    cgroup_v2_snapshot,
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


def test_replicates_use_distinct_report_paths(tmp_path):
    first = PilotConfig(tmp_path, tmp_path / "sample.vcf.gz", tmp_path, tmp_path, 1)
    second = PilotConfig(tmp_path, tmp_path / "sample.vcf.gz", tmp_path, tmp_path, 2)
    assert first.comparison_path.name == "semantic_comparison_r01.json"
    assert second.comparison_path.name == "semantic_comparison_r02.json"
    assert first.summary_path.name == "summary_r01.json"
    assert second.summary_path.name == "summary_r02.json"


def test_cohort_tasks_are_deterministic_and_replicate_specific(tmp_path):
    qc = tmp_path / "sample_qc.tsv"
    qc.write_text(
        "cohort\tsample\tpopulation\tsuperpopulation\tsex\tpath\trecords\tsha256\tstatus\n"
        "1000g\tS1\tPOP\tEUR\tfemale\t/mnt/data/S1.vcf.gz\t10\tabc\tPASS\n"
    )
    first = build_tasks(qc, replicates=3, seed="paper", selected_sample="S1")
    second = build_tasks(qc, replicates=3, seed="paper", selected_sample="S1")
    assert first == second
    assert [task.task_id for task in first] == [0, 1, 2]
    assert [task.replicate for task in first] == [1, 2, 3]
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
    tasks = build_tasks(qc, phase="measured", replicates=3, seed="paper")
    assert len(tasks) == 150
    assert sum(task.first_mode == "cached" for task in tasks) == 75
    for sample in {task.sample for task in tasks}:
        orders = {task.first_mode for task in tasks if task.sample == sample}
        assert orders == {"cached", "uncached"}


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
