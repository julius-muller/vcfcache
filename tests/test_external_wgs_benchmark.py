from __future__ import annotations

import csv
from pathlib import Path

import pytest

from benchmarks.analyze_external_benchmark import scaling_rows, strategy_summary
from benchmarks.prepare_external_wgs import (
    PGP_PROVIDER_CAP,
    Candidate,
    _download_pgp_landing,
    _select_pgp,
    _select_sgdp,
    stable_key,
)
from benchmarks.run_external_cohort import (
    STRATEGIES,
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


def test_four_condition_order_is_balanced_across_52_warmups():
    first = [condition_order(index, 1)[0][0] for index in range(52)]
    assert set(first) == set(STRATEGIES)
    assert {name: first.count(name) for name in STRATEGIES} == {
        name: 13 for name in STRATEGIES
    }


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
