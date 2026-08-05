from __future__ import annotations

import csv
import gzip
import json
import os
from pathlib import Path

import pytest

from benchmarks.analyze_external_benchmark import scaling_rows, strategy_summary
from benchmarks.prepare_external_wgs import (
    COHORT_RECORD_LIMITS,
    PGP_PROVIDER_CAP,
    Candidate,
    Selected,
    _choose_relatedness_replacement,
    _download_pgp_landing,
    _normalization_command,
    _select_pgp,
    _select_sgdp,
    _source_with_reference_contigs,
    _unique_variant_records,
    qc,
    stable_key,
)
from benchmarks.recover_external_attempts import completed_attempt, promote_attempt
from benchmarks.repair_external_duplicate_inputs import variant_key_summary
from benchmarks.run_external_cohort import (
    STRATEGIES,
    _runtime_params,
    build_tasks,
    condition_order,
    read_evaluation_samples,
)


def candidate(
    cohort: str,
    sample: str,
    *,
    provider: str = "DDBJ_NIG",
    region: str = "",
    assembly: str | None = None,
) -> Candidate:
    return Candidate(
        cohort=cohort,
        sample=sample,
        provider=provider,
        population="population",
        region=region,
        sex="XX",
        assembly=assembly or ("GRCh37" if cohort == "pgp" else "GRCh38"),
        source_kind="gVCF",
        url=f"https://example.invalid/{sample}.vcf.gz",
        index_url=f"https://example.invalid/{sample}.vcf.gz.tbi",
        source_name=f"{sample}.vcf.gz",
        source_bytes=100,
    )


def test_stable_key_is_repeatable_and_namespaced():
    assert stable_key("a", "sample") == stable_key("a", "sample")
    assert stable_key("a", "sample") != stable_key("b", "sample")


def test_relatedness_replacement_preserves_evaluation_set_and_cohort():
    def selected(sample: str, role: str, key: str) -> Selected:
        source = candidate("kpgp", sample)
        values = {
            name: value
            for name, value in source.__dict__.items()
            if name not in {"eligibility", "exclusion_reason"}
        }
        return Selected(
            **values,
            role=role,
            selection_seed="seed",
            selection_key=key,
        )

    training = selected("training", "training", "b")
    evaluation = selected("evaluation", "evaluation", "a")
    reserve = candidate("kpgp", "reserve")
    victim, replacement = _choose_relatedness_replacement(
        [training, evaluation], [reserve], "training", "evaluation"
    )
    assert victim == training
    assert replacement == reserve


def test_pgp_landing_fetch_is_single_request(monkeypatch, tmp_path):
    class Response:
        content = b"<!DOCTYPE HTML><html></html>"

        @staticmethod
        def raise_for_status():
            return None

    calls = []

    def get(url, *, timeout):
        calls.append((url, timeout))
        return Response()

    monkeypatch.setattr("benchmarks.prepare_external_wgs.requests.get", get)
    destination = tmp_path / "landing.html"
    assert (
        _download_pgp_landing("https://example.invalid/file", destination)
        == destination
    )
    assert destination.read_bytes() == Response.content
    assert calls == [("https://example.invalid/file", 90)]


def test_sgdp_selection_has_disjoint_balanced_evaluation():
    rows = []
    for region in ("Africa", "WestEurasia", "EastAsia"):
        rows.extend(
            candidate("sgdp", f"{region}-{index}", region=region) for index in range(6)
        )
    rows.extend(
        candidate("sgdp", f"South-{index}", region="SouthAsia") for index in range(3)
    )
    rows.extend(
        candidate("sgdp", f"Central-{index}", region="CentralAsiaSiberia")
        for index in range(3)
    )
    selected = _select_sgdp(rows)
    training = {row.sample for row, role in selected if role == "training"}
    evaluation = [row for row, role in selected if role == "evaluation"]
    assert len(training) == 3
    assert len(evaluation) == 20
    assert not training & {row.sample for row in evaluation}
    assert sum(row.region == "Africa" for row in evaluation) == 5
    assert sum(row.region == "WestEurasia" for row in evaluation) == 5
    assert sum(row.region == "EastAsia" for row in evaluation) == 5
    assert (
        sum(row.region in {"SouthAsia", "CentralAsiaSiberia"} for row in evaluation)
        == 5
    )


def test_pgp_selection_uses_distinct_training_providers_and_evaluation_cap():
    rows = []
    for provider in ("Nebula", "Dante", "Veritas", "Sequencing.com"):
        rows.extend(
            candidate("pgp", f"{provider}-{index}", provider=provider)
            for index in range(7)
        )
    selected = _select_pgp(rows)
    training = [row for row, role in selected if role == "training"]
    evaluation = [row for row, role in selected if role == "evaluation"]
    assert len({row.provider for row in training}) == 3
    assert {row.assembly for row in [*training, *evaluation]} == {"GRCh37"}
    assert len(evaluation) == 12
    assert all(
        sum(row.provider == provider for row in evaluation) <= PGP_PROVIDER_CAP
        for provider in {row.provider for row in evaluation}
    )


def test_pgp_selection_excludes_grch38_candidates():
    rows = [
        candidate("pgp", f"hg19-{provider}-{index}", provider=provider)
        for provider in ("Nebula", "Dante", "Veritas")
        for index in range(7)
    ]
    rows.extend(
        candidate(
            "pgp",
            f"hg38-{provider}",
            provider=provider,
            assembly="GRCh38",
        )
        for provider in ("Nebula", "Dante", "Gencove")
    )
    selected = _select_pgp(rows)
    assert len(selected) == 15
    assert all(row.assembly == "GRCh37" for row, _role in selected)


def test_pgp_reference_header_repair_preserves_records(tmp_path):
    source = tmp_path / "source.vcf.gz"
    with gzip.open(source, "wt") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write("##contig=<ID=chr1,length=10>\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n")
        handle.write(
            "chr10\t2\t.\tA\tG\t.\tPASS\tcustomer_score=7;flagged\tGT:XY\t0/1:value\n"
        )
    reference = tmp_path / "hg19.fa.gz"
    Path(f"{reference}.fai").write_text("chr1\t20\t0\t0\t0\nchr10\t30\t0\t0\t0\n")
    selected = Selected(
        cohort="pgp",
        sample="sample",
        role="evaluation",
        provider="provider",
        population="population",
        region="region",
        sex="XX",
        assembly="GRCh37",
        source_kind="VCF",
        url="https://example.invalid/sample.vcf.gz",
        index_url="",
        source_name="sample.vcf.gz",
        source_bytes=100,
        upstream_md5="",
        documented_overlap="none",
        landing_url="",
        selection_seed="seed",
        selection_key="selection-key",
    )
    repaired = _source_with_reference_contigs(tmp_path, selected, source, reference)
    with gzip.open(repaired, "rt") as handle:
        value = handle.read()
    assert value.count("##contig=<ID=chr1,length=20>") == 1
    assert value.count("##contig=<ID=1,length=20>") == 1
    assert value.count("##contig=<ID=chr10,length=30>") == 1
    assert value.count("##contig=<ID=10,length=30>") == 1
    assert "##INFO=<ID=customer_score,Number=.,Type=String" in value
    assert "##INFO=<ID=flagged,Number=0,Type=Flag" in value
    assert "##FORMAT=<ID=XY,Number=.,Type=String" in value
    assert "chr10\t2\t.\tA\tG" in value
    assert (
        _source_with_reference_contigs(tmp_path, selected, source, reference)
        == repaired
    )


def test_only_pgp_normalization_forces_malformed_auxiliary_info():
    reference = Path("/reference.fa.gz")
    assert "--force" in _normalization_command(reference, "pgp")
    assert "--force" not in _normalization_command(reference, "kpgp")
    assert "--force" not in _normalization_command(reference, "sgdp")


def test_unique_variant_records_preserves_distinct_alleles_and_first_duplicate():
    lines = [
        "##fileformat=VCFv4.2\n",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        "chr1\t1\tfirst\tA\tC\t.\tPASS\t.\n",
        "chr1\t1\tdistinct\tA\tG\t.\tPASS\t.\n",
        "chr1\t1\tduplicate\tA\tC\t.\tPASS\t.\n",
        "chr1\t2\tnext\tA\tC\t.\tPASS\t.\n",
    ]
    observed = list(_unique_variant_records(lines))
    assert observed == [*lines[:4], lines[5]]


def test_variant_key_summary_counts_adjacent_duplicate_keys(monkeypatch, tmp_path):
    class Stream:
        stdout = iter(
            [
                "chr1\t1\tA\tC\tSNP\n",
                "chr1\t1\tA\tC\tSNP\n",
                "chr1\t2\tA\tAT\tINDEL\n",
            ]
        )

        class Stderr:
            @staticmethod
            def read():
                return ""

        stderr = Stderr()

        @staticmethod
        def wait():
            return 0

    monkeypatch.setattr("subprocess.Popen", lambda *args, **kwargs: Stream())
    result = variant_key_summary(tmp_path / "input.vcf.gz")
    assert result["records"] == 3
    assert result["duplicate_keys"] == 1
    assert result["snps"] == 2
    assert result["indels"] == 1
    assert result["contigs"] == ["chr1"]


def test_record_count_guard_accounts_for_heterogeneous_pgp_callers():
    assert COHORT_RECORD_LIMITS["pgp"] == (2_400_000, 6_500_000)
    assert COHORT_RECORD_LIMITS["kpgp"] == (2_500_000, 6_500_000)
    assert COHORT_RECORD_LIMITS["sgdp"] == (2_500_000, 6_500_000)


def test_external_stager_installs_grch37_vep_cache_on_workers():
    script = (
        Path(__file__).parents[1] / "benchmarks/stage_external_controller.sh"
    ).read_text()
    assert "115_GRCh37/" in script
    assert "115_GRCh37/info.txt" in script
    assert "115_GRCh37/1/all_vars.gz" in script


def test_external_stager_can_wait_without_reserving_a_slurm_worker():
    script = (
        Path(__file__).parents[1] / "benchmarks/stage_external_controller.sh"
    ).read_text()
    assert "gate_job_id=${VCFCACHE_GATE_JOB_ID:-}" in script
    assert 'if [[ -n "$gate_job_id" ]]; then' in script


def test_external_stager_creates_worker_receiver_parents():
    script = (
        Path(__file__).parents[1] / "benchmarks/stage_external_controller.sh"
    ).read_text()
    assert "sudo mkdir -p '$external_root' '$external_root/deployment'" in script


def _complete_attempt_run(root: Path, sample: str, replicate: int = 1) -> Path:
    for condition in (
        "uncached",
        "gnomad_af_0.1",
        "gnomad_af_0.01",
        "cohort_3_genomes",
    ):
        mode = "uncached" if condition == "uncached" else "cached"
        run = (
            root
            / "runs"
            / condition
            / "pilot"
            / sample
            / "commit"
            / f"{mode}_r{replicate:02d}"
        )
        run.mkdir(parents=True)
        (run / "metrics.json").write_text(json.dumps({"wall_seconds": 1}))
        (run / "output.bcf").write_bytes(b"bcf")
        (run / "output.bcf.csi").write_bytes(b"index")
    return root


def test_completed_attempt_selects_newest_four_condition_attempt(tmp_path):
    task = tmp_path / "warmup/attempts/task-3"
    old = _complete_attempt_run(task / "job-10/run", "sample")
    (task / "job-10/failure.json").write_text("{}")
    newest = _complete_attempt_run(task / "job-12/run", "sample")
    (task / "job-12/failure.json").write_text("{}")
    found, runs = completed_attempt(tmp_path, "warmup", 3, "sample", 1)
    assert found == newest
    assert set(runs) == {
        "uncached",
        "gnomad_af_0.1",
        "gnomad_af_0.01",
        "cohort_3_genomes",
    }
    assert found != old


def test_promote_attempt_uses_hardlinks_and_records_recovery(tmp_path, monkeypatch):
    attempt = tmp_path / "warmup/attempts/task-4/job-20/run"
    attempt.mkdir(parents=True)
    source = attempt / "large-output.bcf"
    source.write_bytes(b"content")
    failure = attempt.parent / "failure.json"
    failure.write_text('{"exit_code": 1}\n')
    result = tmp_path / "warmup/tasks/task-4"
    monkeypatch.setenv("SLURM_JOB_ID", "99")
    monkeypatch.setenv("SLURM_ARRAY_JOB_ID", "98")
    monkeypatch.setenv("SLURMD_NODENAME", "sl-w5")
    promote_attempt(
        attempt_run=attempt,
        result_dir=result,
        task_id=4,
        phase="warmup",
        source_failure=failure,
    )
    assert (result / "large-output.bcf").read_bytes() == b"content"
    assert os.stat(source).st_ino == os.stat(result / "large-output.bcf").st_ino
    metadata = json.loads((result / "slurm-task.json").read_text())
    assert metadata["recovered"] is True
    assert metadata["slurm_node"] == "sl-w5"
    assert json.loads((result / "recovery.json").read_text())["storage"] == "hardlinks"


def test_external_campaign_cleanliness_check_is_cwd_independent():
    source = (
        Path(__file__).parents[1] / "benchmarks/run_external_cohort.py"
    ).read_text()
    assert '["git", "-C", str(REPO_ROOT), "diff", "--quiet"]' in source


def test_relatedness_screen_acknowledges_small_assembly_strata():
    source = (
        Path(__file__).parents[1] / "benchmarks/prepare_external_wgs.py"
    ).read_text()
    assert '"--bad-ld",\n                "--indep-pairwise"' in source


def test_relatedness_screen_assigns_stable_allele_specific_variant_ids():
    source = (
        Path(__file__).parents[1] / "benchmarks/prepare_external_wgs.py"
    ).read_text()
    assert source.count('"--set-all-var-ids"') == 2
    assert source.count('"@:#:$r:$a"') == 2
    assert source.count('"--rm-dup"') == 2
    assert "_deduplicate_variant_keys(raw_partial, partial)" in source
    assert source.count('"exclude-all"') == 2
    assert source.count('"--snps-only"') == 2
    assert source.count('"--max-alleles"') == 2


def test_qc_preserves_source_id_separately_from_vcf_sample(
    tmp_path, monkeypatch, capsys
):
    row = Selected(
        cohort="pgp",
        sample="huDEFDD1",
        role="evaluation",
        provider="Veritas",
        population="",
        region="",
        sex="",
        assembly="GRCh37",
        source_kind="VCF",
        url="https://example.invalid/file",
        index_url="",
        source_name="source.vcf.gz",
        source_bytes=100,
        upstream_md5="",
        documented_overlap="none",
        landing_url="",
        selection_seed="seed",
        selection_key="key",
    )
    monkeypatch.setattr(
        "benchmarks.prepare_external_wgs._selected", lambda _path: [row]
    )
    monkeypatch.setattr(
        "benchmarks.prepare_external_wgs._prepared_path",
        lambda _root, _row: tmp_path / "huDEFDD1.vcf.gz",
    )
    monkeypatch.setattr(
        "benchmarks.prepare_external_wgs.validate_prepared_vcf",
        lambda _path, *, cohort: {
            "status": "PASS",
            "errors": "",
            "sample": "Sample1",
            "records": 2_474_642,
            "snps": 2_077_905,
            "indels": 396_737,
            "multiallelic": 0,
            "contigs": "chr1,chr2",
            "bytes": 100,
            "sha256": "sha",
        },
    )
    monkeypatch.setattr(
        "benchmarks.prepare_external_wgs._duplicate_variant_keys", lambda _path: 0
    )
    args = type("Args", (), {"root": tmp_path})()
    output = qc(args)
    capsys.readouterr()
    with output.open(newline="") as handle:
        value = next(csv.DictReader(handle, delimiter="\t"))
    assert value["sample"] == "huDEFDD1"
    assert value["vcf_sample"] == "Sample1"
    assert value["status"] == "PASS"


def test_four_condition_order_is_balanced_across_52_warmups():
    first = [condition_order(index, 1)[0][0] for index in range(52)]
    assert set(first) == set(STRATEGIES)
    assert {name: first.count(name) for name in STRATEGIES} == {
        name: 13 for name in STRATEGIES
    }


def test_runtime_params_lower_buffer_without_mutating_cache_snapshot(tmp_path):
    source = tmp_path / "params.snapshot.yaml"
    source.write_text("genome_build: GRCh38\nvep_buffer: 500000\nvep_forks: 8\n")
    manifest = _runtime_params(tmp_path / "campaign", "GRCh38", source, 100_000)
    runtime = Path(str(manifest["path"]))
    assert "vep_buffer: 100000" in runtime.read_text()
    assert source.read_text().splitlines()[1] == "vep_buffer: 500000"
    assert manifest["vep_buffer"] == 100_000


def _qc(path: Path, *, relatedness: str = "PASS") -> None:
    fields = (
        "cohort",
        "sample",
        "role",
        "assembly",
        "population",
        "superpopulation",
        "sex",
        "provider",
        "path",
        "records",
        "sha256",
        "status",
        "relatedness_status",
        "documented_overlap",
    )
    rows = []
    for cohort, count in (("kpgp", 20), ("sgdp", 20), ("pgp", 12)):
        rows.extend(
            {
                "cohort": cohort,
                "sample": f"{cohort}-{index}",
                "role": "evaluation",
                "assembly": "GRCh37" if cohort == "pgp" else "GRCh38",
                "population": "pop",
                "superpopulation": "region",
                "sex": "XX",
                "provider": "provider",
                "path": f"/mnt/data/{cohort}-{index}.vcf.gz",
                "records": "4000000",
                "sha256": f"sha-{cohort}-{index}",
                "status": "PASS",
                "relatedness_status": relatedness,
                "documented_overlap": "none",
            }
            for index in range(count)
        )
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def test_campaign_tasks_require_relatedness_and_have_one_common_baseline(tmp_path):
    qc = tmp_path / "qc.tsv"
    _qc(qc)
    samples = read_evaluation_samples(qc)
    tasks = build_tasks(samples, phase="measured", replicates=3)
    assert len(tasks) == 156
    assert {task.assembly for task in tasks} == {"GRCh37", "GRCh38"}
    assert all(set(task.strategy_order.split(",")) == set(STRATEGIES) for task in tasks)
    assert all(task.strategy_order.split(",").count("uncached") == 1 for task in tasks)

    _qc(qc, relatedness="PENDING")
    with pytest.raises(RuntimeError, match="incomplete"):
        read_evaluation_samples(qc)


def test_scaling_model_amortizes_only_custom_build_cost():
    rows = []
    for strategy, cached, hit in (
        ("gnomad_af_0.1", 70.0, 0.4),
        ("gnomad_af_0.01", 50.0, 0.6),
        ("cohort_3_genomes", 60.0, 0.5),
    ):
        rows.append(
            {
                "cohort": "kpgp",
                "strategy": strategy,
                "sample": "S1",
                "replicate": "1",
                "cache_hit_rate": str(hit),
                "cached_wall_seconds": str(cached),
                "uncached_wall_seconds": "100",
                "speedup": str(100 / cached),
                "relative_runtime": str(cached / 100),
                "semantic_pass": "true",
            }
        )
    summaries = strategy_summary(rows)
    modeled = scaling_rows(summaries, {"kpgp": 120.0})
    custom = [row for row in modeled if row["strategy"] == "cohort_3_genomes"]
    assert custom[0]["effective_seconds_per_sample"] == 180.0
    assert custom[99]["effective_seconds_per_sample"] == 61.2
    public = [row for row in modeled if row["strategy"] == "gnomad_af_0.01"]
    assert public[0]["effective_seconds_per_sample"] == 50.0
